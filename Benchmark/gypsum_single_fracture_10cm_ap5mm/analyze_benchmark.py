#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import re
import shutil
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
BENCH_DIR = REPO_ROOT / "Benchmark" / "gypsum_single_fracture_10cm_ap5mm"
RESULTS_DIR = REPO_ROOT / "Output" / "benchmark_gypsum_single_fracture_10cm_ap5mm"
CASE_INDEX = RESULTS_DIR / "case_index.csv"
FIG_DIR = RESULTS_DIR / "figures"
LOCAL_RESULTS_DIR = BENCH_DIR / "results_snapshot"
LOCAL_FIG_DIR = LOCAL_RESULTS_DIR / "figures"

LENGTH = 0.1
APERTURE0 = 5.0e-3
THICKNESS = 0.1
RHO = 1.0e3
MU = 1.0e-3
GLOBAL_INJECTION_TIME = 1.0e5
PROFILE_TIMES = [0.0, 1.0e4, 5.0e4, 1.0e5]
COMBINED_TIMES = [1.0e2, 1.0e3, 1.0e4]
D_EFFECTIVE = 1.0e-9

A2 = 2.50122e-6
K2 = 2.72730e-5
VREF = 1.0e-3


def chemistry_coefficients(travel_time: float) -> dict[str, float]:
    return {"A2": A2, "K2": K2, "L": 0.0, "VREF": VREF}


def effective_diffusion_height_factor(travel_time: float, aperture_height: float = APERTURE0) -> float:
    if aperture_height <= 0.0:
        return 0.0
    return math.sqrt(D_EFFECTIVE * travel_time) / aperture_height


def load_aperture_delta(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def load_snapshot(path: Path) -> tuple[float, np.ndarray]:
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().strip()
    match = re.search(r"time=([0-9eE+.\-]+)", header)
    if match is None:
        raise ValueError(f"Missing time header in {path}")
    time_value = float(match.group(1))
    data = np.loadtxt(path, comments="#", skiprows=2)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return time_value, data


def load_snapshot_nearest(raw_dir: Path, target_time: float) -> tuple[float, np.ndarray]:
    best_match: tuple[float, Path] | None = None
    for path in sorted(raw_dir.glob("DFN_step*.txt")):
        time_value, _ = load_snapshot(path)
        distance = abs(time_value - target_time)
        if best_match is None or distance < best_match[0]:
            best_match = (distance, path)
    if best_match is None:
        raise FileNotFoundError(f"No DFN_step snapshots found in {raw_dir}")
    return load_snapshot(best_match[1])


def compute_case_metrics(row: dict[str, str]) -> dict[str, float | str]:
    case_dir = Path(row["case_output_dir"])
    raw_dir = case_dir / "raw_output"
    final = load_aperture_delta(raw_dir / "DFN_aperture_delta.txt")
    lengths = np.sqrt((final[:, 2] - final[:, 0]) ** 2 + (final[:, 3] - final[:, 1]) ** 2)
    delta_b = final[:, 6]
    dissolved_volume = float(np.sum(2.0 * delta_b * lengths * THICKNESS))
    travel_time = float(row["travel_time_s"])
    velocity = float(row["velocity_m_per_s"])
    chemistry = chemistry_coefficients(travel_time)
    reference_volume = effective_diffusion_height_factor(travel_time) * (velocity * APERTURE0 * THICKNESS / chemistry["VREF"]) * chemistry["A2"] * (
        GLOBAL_INJECTION_TIME
        - (1.0 - math.exp(-chemistry["K2"] * GLOBAL_INJECTION_TIME)) / chemistry["K2"]
    )
    rate = dissolved_volume / travel_time
    relative_error = (dissolved_volume - reference_volume) / reference_volume if reference_volume else math.nan
    final_velocity = float(np.mean(np.abs(final[:, 5])))
    return {
        "case": row["case"],
        "travel_time_s": travel_time,
        "segments": int(row["segments"]),
        "velocity_m_per_s": velocity,
        "final_velocity_m_per_s": final_velocity,
        "reynolds": float(row["reynolds"]),
        "damkohler": float(row["damkohler"]),
        "delta_v_a2_m3": chemistry["A2"],
        "delta_v_k2_s_inv": chemistry["K2"],
        "vref_m3": chemistry["VREF"],
        "effective_diffusion_height_factor": effective_diffusion_height_factor(travel_time),
        "reference_dissolved_volume_m3": reference_volume,
        "dissolved_volume_m3": dissolved_volume,
        "avg_dissolution_rate_m3_per_s": rate,
        "relative_error": relative_error,
    }


def save_summary(rows: list[dict[str, float | str]]) -> None:
    summary_path = RESULTS_DIR / "summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_reaction_law() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    t = np.logspace(1, 5, 400)
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    chemistry = chemistry_coefficients(1.0e2)
    y = chemistry["A2"] * (1.0 - np.exp(-chemistry["K2"] * t))
    ax.plot(t, y * 1e6, color="#b13a1b", lw=2.2, label=r"$\Delta V_{ref}(t)=A_2(1-e^{-k_2 t})$")
    for tt in [1e2, 1e3, 1e4, 1e5]:
        ax.axvline(tt, color="#4a5568", lw=0.9, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("Residence Time [s]")
    ax.set_ylabel(r"Reference $\Delta V$ [$10^{-6}$ m$^3$]")
    ax.set_title("Pure Gypsum Fast2-Only Mineral Volume Law")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_reaction_law_and_targets.png", dpi=220)
    plt.close(fig)


def plot_re_da_map(rows: list[dict[str, float | str]]) -> None:
    unique = {}
    for row in rows:
        unique[row["travel_time_s"]] = (row["reynolds"], row["damkohler"])
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    for tt in sorted(unique):
        re_value, da_value = unique[tt]
        ax.scatter(re_value, da_value, s=70, color="#1f4e79")
        ax.text(re_value * 1.08, da_value * 1.06, f"t={tt:.0e} s", fontsize=9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Re")
    ax.set_ylabel("Da = k2 * t_cross")
    ax.set_title("Benchmark Parameter Points")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_re_da_map.png", dpi=220)
    plt.close(fig)


def plot_vs_segments(rows: list[dict[str, float | str]]) -> None:
    grouped: dict[float, list[dict[str, float | str]]] = defaultdict(list)
    for row in rows:
        grouped[float(row["travel_time_s"])].append(row)

    colors = {1e2: "#9b2226", 1e3: "#ca6702", 1e4: "#005f73", 1e5: "#0a9396"}

    fig, ax = plt.subplots(figsize=(7.8, 5.0))
    for tt in sorted(grouped):
        data = sorted(grouped[tt], key=lambda item: int(item["segments"]))
        x = [int(item["segments"]) for item in data]
        y = [float(item["dissolved_volume_m3"]) * 1e6 for item in data]
        ref_value = float(data[0]["reference_dissolved_volume_m3"]) * 1e6
        ax.plot(x, y, marker="o", lw=1.8, color=colors[tt], label=f"t={tt:.0e} s")
        ax.hlines(ref_value, x[0], x[-1], colors=colors[tt], linestyles="dashed", lw=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("Number Of Fracture Segments")
    ax.set_ylabel(r"Dissolved Gypsum Volume [$10^{-6}$ m$^3$]")
    ax.set_title("Model Dissolved Volume vs Segment Discretization")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_dissolved_volume_vs_segments.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.8, 5.0))
    for tt in sorted(grouped):
        data = sorted(grouped[tt], key=lambda item: int(item["segments"]))
        x = [int(item["segments"]) for item in data]
        y = [100.0 * float(item["relative_error"]) for item in data]
        ax.plot(x, y, marker="o", lw=1.8, color=colors[tt], label=f"t={tt:.0e} s")
    ax.axhline(0.0, color="#222222", lw=0.9)
    ax.set_xscale("log")
    ax.set_xlabel("Number Of Fracture Segments")
    ax.set_ylabel("Relative Error [%]")
    ax.set_title("Error Relative To Analytical Plug-Flow Reference")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_relative_error_vs_segments.png", dpi=220)
    plt.close(fig)

    travel_times = sorted(grouped)
    segments = sorted({int(row["segments"]) for row in rows})
    z = np.zeros((len(travel_times), len(segments)))
    for i, tt in enumerate(travel_times):
        rows_at_time = grouped[tt]
        ref_row = next(item for item in rows_at_time if int(item["segments"]) == 100)
        ref_value = float(ref_row["dissolved_volume_m3"])
        for j, seg in enumerate(segments):
            value = float(next(item for item in rows_at_time if int(item["segments"]) == seg)["dissolved_volume_m3"])
            z[i, j] = abs(100.0 * (value - ref_value) / ref_value) if ref_value else 0.0
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    im = ax.imshow(z, cmap="magma_r", aspect="auto")
    ax.set_xticks(range(len(segments)), [str(seg) for seg in segments])
    ax.set_yticks(range(len(travel_times)), [f"{tt:.0e}" for tt in travel_times])
    ax.set_xlabel("Segments")
    ax.set_ylabel("Crossing Time [s]")
    ax.set_title("Absolute Error Heatmap Relative To n=100 [%]")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Absolute Difference Relative To n=100 [%]")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_error_heatmap.png", dpi=220)
    plt.close(fig)


def plot_custom_t1e5(rows: list[dict[str, float | str]]) -> None:
    rows_t1e5 = sorted([row for row in rows if float(row["travel_time_s"]) == 1.0e5], key=lambda item: int(item["segments"]))
    labels = [str(int(row["segments"])) for row in rows_t1e5]
    volumes = np.array([float(row["dissolved_volume_m3"]) for row in rows_t1e5])
    ref_value = volumes[-1]
    abs_diff = np.abs(volumes - ref_value)
    x = np.arange(len(labels), dtype=float)

    fig, ax1 = plt.subplots(figsize=(7.8, 4.8))
    ax2 = ax1.twinx()
    ax1.plot(x, volumes, color="#1f4e79", marker="o", lw=2.0)
    ax2.bar(x, abs_diff, width=0.55, color="#d9c6a5", alpha=0.75)
    ax1.set_xticks(x, labels)
    ax1.set_xlabel("n")
    ax1.set_ylabel("Dissolved volume [m$^3$]")
    ax2.set_ylabel(r"Absolute difference from $n=100$ [m$^3$]")
    ax1.set_title(r"Dissolved volume at $t=10^5$ s")
    ax1.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_dissolution_volume_t1e5_vs_n_absdiff_bar.png", dpi=220)
    plt.close(fig)

    fig, ax1 = plt.subplots(figsize=(7.8, 4.8))
    ax2 = ax1.twinx()
    diff_percent = np.zeros_like(volumes) if ref_value == 0.0 else 100.0 * (volumes - ref_value) / ref_value
    ax1.plot(x, volumes, color="#1f4e79", marker="o", lw=2.0)
    ax2.plot(x, diff_percent, color="#b5651d", marker="s", lw=1.8)
    ax1.set_xticks(x, labels)
    ax1.set_xlabel("n")
    ax1.set_ylabel("Dissolved volume [m$^3$]")
    ax2.set_ylabel(r"Difference from $n=100$ [%]")
    ax1.set_title(r"Dissolved volume at $t=10^5$ s")
    ax1.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_dissolution_volume_t1e5_vs_n_with_diff.png", dpi=220)
    plt.close(fig)


def plot_aperture_profile_t1e5() -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    colors = ["#9b2226", "#ca6702", "#005f73", "#0a9396", "#3d405b"]
    for seg, color in zip([1, 5, 10, 50, 100], colors):
        raw = RESULTS_DIR / "cases" / f"t1e05_n{seg:03d}" / "raw_output" / "DFN_aperture_delta.txt"
        data = load_aperture_delta(raw)
        x_left = data[:, 0] + 0.05
        x_right = data[:, 2] + 0.05
        aperture_mm = data[:, 4] * 1e3
        order = np.argsort(x_left)
        x_left = x_left[order] * 100.0
        x_right = x_right[order] * 100.0
        aperture_mm = aperture_mm[order]
        x_step = np.empty(2 * len(aperture_mm))
        y_step = np.empty(2 * len(aperture_mm))
        x_step[0::2] = x_left
        x_step[1::2] = x_right
        y_step[0::2] = aperture_mm
        y_step[1::2] = aperture_mm
        ax.plot(x_step, y_step, lw=2.0, color=color, label=f"n={seg}")
    ax.axhline(APERTURE0 * 1e3, color="#666666", lw=1.0, ls="--")
    ax.set_xlim(0.0, 10.0)
    ax.set_xlabel("Position along fracture [cm]")
    ax.set_ylabel("Final aperture [mm]")
    ax.set_title(r"Segmented aperture profile at $t=10^5$ s")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_aperture_profile_t1e5_by_segments.png", dpi=220)
    plt.close(fig)


def plot_aperture_evolution_t1e04_n100() -> None:
    raw_dir = RESULTS_DIR / "cases" / "t1e04_n100" / "raw_output"
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    colors = ["#6c757d", "#9b2226", "#ca6702", "#0a9396"]
    for t_value, color in zip(PROFILE_TIMES, colors):
        _, data = load_snapshot_nearest(raw_dir, t_value)
        x_left = (data[:, 0] + 0.05) * 100.0
        x_right = (data[:, 2] + 0.05) * 100.0
        aperture_mm = data[:, 4] * 1e3
        order = np.argsort(x_left)
        x_left = x_left[order]
        x_right = x_right[order]
        aperture_mm = aperture_mm[order]
        x_step = np.empty(2 * len(aperture_mm))
        y_step = np.empty(2 * len(aperture_mm))
        x_step[0::2] = x_left
        x_step[1::2] = x_right
        y_step[0::2] = aperture_mm
        y_step[1::2] = aperture_mm
        label = "t = 0 s" if t_value == 0 else f"t = {t_value:.0e} s"
        ax.plot(x_step, y_step, lw=1.8, color=color, label=label)
    ax.axhline(APERTURE0 * 1e3, color="#666666", lw=1.0, ls="--")
    ax.set_xlim(0.0, 10.0)
    ax.set_xlabel("Position along fracture [cm]")
    ax.set_ylabel("Aperture [mm]")
    ax.set_title(r"Aperture-profile evolution for $t_{cross}=10^4$ s, $n=100$")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_aperture_evolution_t1e04_n100.png", dpi=220)
    plt.close(fig)


def compute_time_series(case_label: str) -> tuple[np.ndarray, np.ndarray]:
    raw_dir = RESULTS_DIR / "cases" / case_label / "raw_output"
    snapshot_paths = sorted(raw_dir.glob("DFN_step*.txt"))
    t_values = []
    dissolved_values = []
    for path in snapshot_paths:
        time_value, data = load_snapshot(path)
        lengths = np.sqrt((data[:, 2] - data[:, 0]) ** 2 + (data[:, 3] - data[:, 1]) ** 2)
        delta_b = data[:, 4] - APERTURE0
        dissolved_values.append(float(np.sum(2.0 * delta_b * lengths * THICKNESS)))
        t_values.append(time_value)
    return np.array(t_values), np.array(dissolved_values)


def plot_combined_evolution_figure() -> None:
    fig = plt.figure(figsize=(12.0, 8.2))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.25, 1.0], hspace=0.22, wspace=0.25)

    profile_times = PROFILE_TIMES
    profile_colors = ["#6c757d", "#9b2226", "#ca6702", "#0a9396"]
    volume_colors = {1.0e2: "#9b2226", 1.0e3: "#ca6702", 1.0e4: "#005f73"}

    for row_index, crossing_time in enumerate(COMBINED_TIMES):
        ax = fig.add_subplot(gs[row_index, 0])
        raw_dir = RESULTS_DIR / "cases" / f"t{crossing_time:.0e}_n100".replace("+", "") / "raw_output"
        for t_value, color in zip(profile_times, profile_colors):
            _, data = load_snapshot_nearest(raw_dir, t_value)
            x_left = (data[:, 0] + 0.05) * 100.0
            x_right = (data[:, 2] + 0.05) * 100.0
            aperture_mm = data[:, 4] * 1e3
            order = np.argsort(x_left)
            x_left = x_left[order]
            x_right = x_right[order]
            aperture_mm = aperture_mm[order]
            x_step = np.empty(2 * len(aperture_mm))
            y_step = np.empty(2 * len(aperture_mm))
            x_step[0::2] = x_left
            x_step[1::2] = x_right
            y_step[0::2] = aperture_mm
            y_step[1::2] = aperture_mm
            label = "0" if t_value == 0 else f"{t_value:.0e}"
            ax.plot(x_step, y_step, lw=1.5, color=color, label=label)
        ax.axhline(APERTURE0 * 1e3, color="#666666", lw=0.9, ls="--")
        ax.set_xlim(0.0, 10.0)
        ax.grid(True, axis="y", alpha=0.22)
        ax.set_ylabel("Aperture [mm]")
        ax.set_title(f"t_cross = {crossing_time:.0e} s, n = 100", fontsize=10)
        if row_index == len(COMBINED_TIMES) - 1:
            ax.set_xlabel("Position along fracture [cm]")
        if row_index == 0:
            ax.legend(title="time [s]", frameon=False, ncol=4, fontsize=8, title_fontsize=8, loc="upper right")

    ax_right = fig.add_subplot(gs[:, 1])
    for crossing_time, color in volume_colors.items():
        case_label = f"t{crossing_time:.0e}_n100".replace("+", "")
        t_values, dissolved_values = compute_time_series(case_label)
        ax_right.plot(t_values, dissolved_values, lw=2.0, color=color, label=f"t_cross={crossing_time:.0e} s")
    ax_right.set_xlabel("Time [s]")
    ax_right.set_ylabel("Total dissolved volume [m$^3$]")
    ax_right.set_title("Total dissolved-volume evolution")
    ax_right.grid(True, alpha=0.25)
    ax_right.legend(frameon=False)

    fig.suptitle("Aperture and Dissolution Evolution for the 5 mm Benchmark", y=0.995)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_combined_evolution_profiles_and_volume.png", dpi=220)
    plt.close(fig)


def write_notes(rows: list[dict[str, float | str]]) -> None:
    by_time: dict[float, list[dict[str, float | str]]] = defaultdict(list)
    for row in rows:
        by_time[float(row["travel_time_s"])].append(row)

    lines = [
        "# Single-Fracture Pure Gypsum Benchmark: 5 mm Aperture Variant",
        "",
        "This benchmark is identical to the current 10 cm single-fracture gypsum test except for one change:",
        "",
        "- initial fracture aperture is `5 mm` instead of `1 mm`.",
        "",
        "Shared setup:",
        "",
        "- domain size `10 cm x 10 cm`;",
        "- one centered horizontal fracture spanning the full domain length;",
        "- out-of-plane thickness `10 cm`;",
        "- chemistry kept as pure gypsum using only the `fast2` term;",
        "- all crossing-time cases share one gypsum-only volume law with `A2 = 2.50122e-6`, `k2 = 2.72730e-5`, and `Vref = 1e-3 m^3`;",
        "- an effective diffusion-height factor `m = sqrt(D * t_cross) / b` is applied with `D = 1e-9 m^2/s` and `b` approximated by the initial aperture in the analytical reference plots;",
        "- crossing times `100 s`, `1000 s`, `1e4 s`, `1e5 s`;",
        "- segment counts `n = 1, 5, 10, 50, 100`.",
        "",
        "Compared with the 1 mm benchmark, the larger aperture changes both flow rate and Reynolds number.",
        "",
        "Files:",
        "",
        "- `generated_inputs/`: copies of the exact input files used for this variant.",
        "- `results_snapshot/`: copied summaries and figures for this benchmark.",
        "- `Output/benchmark_gypsum_single_fracture_10cm_ap5mm/`: raw solver outputs and analysis products.",
    ]
    (BENCH_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    lines_cn = [
        "# 5 mm 裂缝宽度版本 Benchmark",
        "",
        "这组 benchmark 与当前 10 cm 单裂缝石膏 benchmark 的唯一差别是：",
        "",
        "- 初始裂缝宽度改为 `5 mm`，原测试组为 `1 mm`。",
        "",
        "其他设定保持一致：",
        "",
        "- 域大小 `10 cm x 10 cm`；",
        "- 一条居中水平裂缝贯穿整个域；",
        "- out-of-plane thickness 为 `10 cm`；",
        "- chemistry 仍然只保留纯 gypsum 的 `fast2` 项；",
        "- 所有 crossing time 都共用同一条 gypsum-only 体积公式，`A2 = 2.50122e-6`、`k2 = 2.72730e-5`，并统一 `Vref = 1e-3 m^3`；",
        "- 另外加入有效扩散高度因子 `m = sqrt(D * t_cross) / b`，其中 `D = 1e-9 m^2/s`，分析参考曲线中用初始 aperture 近似 `b`；",
        "- 穿越时间为 `100 s`、`1000 s`、`1e4 s`、`1e5 s`；",
        "- segment 数为 `n = 1, 5, 10, 50, 100`。",
        "",
        "相较于 1 mm 版本，5 mm 裂缝会改变流量、Re 数以及总注入水量，因此溶解体积和 aperture 演化都会系统变化。",
        "",
        "文件说明：",
        "",
        "- `generated_inputs/`：本组 benchmark 实际使用的输入文件副本；",
        "- `results_snapshot/`：本组 benchmark 的 summary 和图件副本；",
        "- `Output/benchmark_gypsum_single_fracture_10cm_ap5mm/`：原始求解器输出和分析结果。",
    ]
    (RESULTS_DIR / "README_CN.md").write_text("\n".join(lines_cn) + "\n", encoding="utf-8")
    (RESULTS_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_results_snapshot() -> None:
    LOCAL_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    LOCAL_FIG_DIR.mkdir(parents=True, exist_ok=True)
    for name in ["summary.csv", "case_index.csv", "README.md", "README_CN.md"]:
        source = RESULTS_DIR / name
        if source.exists():
            shutil.copy2(source, LOCAL_RESULTS_DIR / name)
    for figure in FIG_DIR.glob("*.png"):
        shutil.copy2(figure, LOCAL_FIG_DIR / figure.name)


def main() -> None:
    if not CASE_INDEX.exists():
        raise FileNotFoundError(f"Missing case index: {CASE_INDEX}")
    with CASE_INDEX.open("r", newline="", encoding="utf-8") as f:
        case_rows = list(csv.DictReader(f))
    rows = [compute_case_metrics(row) for row in case_rows]
    rows.sort(key=lambda item: (float(item["travel_time_s"]), int(item["segments"])))
    save_summary(rows)
    plot_reaction_law()
    plot_re_da_map(rows)
    plot_vs_segments(rows)
    plot_custom_t1e5(rows)
    plot_aperture_profile_t1e5()
    plot_aperture_evolution_t1e04_n100()
    plot_combined_evolution_figure()
    write_notes(rows)
    copy_results_snapshot()
    print(f"Wrote analysis products to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
