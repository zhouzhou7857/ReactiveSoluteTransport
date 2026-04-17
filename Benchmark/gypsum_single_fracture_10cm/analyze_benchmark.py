#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "Output" / "benchmark_gypsum_single_fracture_10cm"
CASE_INDEX = RESULTS_DIR / "case_index.csv"
FIG_DIR = RESULTS_DIR / "figures"

LENGTH = 0.1
APERTURE0 = 1.0e-3
THICKNESS = 0.1
A2 = 2.50122e-6
K2 = 2.72730e-5
VREF = 1.0e-3
RHO = 1.0e3
MU = 1.0e-3
GLOBAL_INJECTION_TIME = 1.0e5

def load_aperture_delta(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def compute_case_metrics(row: dict[str, str]) -> dict[str, float | str]:
    case_dir = Path(row["case_output_dir"])
    raw_dir = case_dir / "raw_output"
    aperture_delta_path = raw_dir / "DFN_aperture_delta.txt"
    if not aperture_delta_path.exists():
        raise FileNotFoundError(f"No DFN_aperture_delta.txt found in {raw_dir}")
    final = load_aperture_delta(aperture_delta_path)
    lengths = np.sqrt((final[:, 2] - final[:, 0]) ** 2 + (final[:, 3] - final[:, 1]) ** 2)
    delta_b = final[:, 6]
    dissolved_volume = float(np.sum(2.0 * delta_b * lengths * THICKNESS))
    travel_time = float(row["travel_time_s"])
    velocity = float(row["velocity_m_per_s"])
    reference_volume = (velocity * APERTURE0 * THICKNESS / VREF) * A2 * (
        GLOBAL_INJECTION_TIME - (1.0 - math.exp(-K2 * GLOBAL_INJECTION_TIME)) / K2
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
    y = A2 * (1.0 - np.exp(-K2 * t))
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(t, y * 1e6, color="#b13a1b", lw=2.2, label=r"$\Delta V_{ref}(t)=A_2(1-e^{-k_2 t})$")
    for tt in [1e2, 1e3, 1e4, 1e5]:
        ax.axvline(tt, color="#4a5568", lw=0.9, ls="--")
        ax.text(tt, A2 * (1.0 - math.exp(-K2 * tt)) * 1e6 * 1.03, f"{tt:.0e} s", fontsize=9, ha="left", va="bottom")
    ax.set_xscale("log")
    ax.set_xlabel("Residence Time [s]")
    ax.set_ylabel(r"Reference $\Delta V$ [$10^{-6}$ m$^3$]")
    ax.set_title("Pure Gypsum Fast2-Only Mineral Volume Law")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_reaction_law_and_targets.png", dpi=220)
    plt.close(fig)


def plot_re_da_map(rows: list[dict[str, float | str]]) -> None:
    unique = {}
    for row in rows:
        unique[row["travel_time_s"]] = (row["reynolds"], row["damkohler"])
    times = sorted(unique)
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    for tt in times:
        re, da = unique[tt]
        ax.scatter(re, da, s=70, color="#1f4e79")
        ax.text(re * 1.08, da * 1.06, f"t={tt:.0e} s", fontsize=9)
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

    colors = {
        1e2: "#9b2226",
        1e3: "#ca6702",
        1e4: "#005f73",
        1e5: "#0a9396",
    }

    fig, ax = plt.subplots(figsize=(7.8, 5.0))
    for tt in sorted(grouped):
        data = sorted(grouped[tt], key=lambda item: int(item["segments"]))
        x = [int(item["segments"]) for item in data]
        y = [float(item["dissolved_volume_m3"]) * 1e6 for item in data]
        ref = float(data[0]["reference_dissolved_volume_m3"]) * 1e6
        ax.plot(x, y, marker="o", lw=1.8, color=colors[tt], label=f"t={tt:.0e} s")
        ax.hlines(ref, x[0], x[-1], colors=colors[tt], linestyles="dashed", lw=1.0)
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
        by_seg = {}
        for item in rows_at_time:
            seg = int(item["segments"])
            value = float(item["dissolved_volume_m3"])
            if ref_value != 0.0:
                by_seg[seg] = abs(100.0 * (value - ref_value) / ref_value)
            else:
                by_seg[seg] = 0.0 if value == 0.0 else math.nan
        for j, seg in enumerate(segments):
            z[i, j] = by_seg[seg]
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


def write_note(rows: list[dict[str, float | str]]) -> None:
    by_time: dict[float, list[dict[str, float | str]]] = defaultdict(list)
    for row in rows:
        by_time[float(row["travel_time_s"])].append(row)

    lines = [
        "# 10 cm Single-Fracture Pure Gypsum Benchmark",
        "",
        "## Meaning Of This Benchmark",
        "",
        "- It isolates the numerical effect of fracture segmentation in the simplest possible geometry: one straight fracture crossing the full 10 cm domain.",
        "- Because the flow path is unique and matrix diffusion is suppressed, differences in dissolved gypsum volume mainly reflect chemistry-flow residence time and segment discretization, not DFN topology.",
        "- The pure gypsum law keeps only the code-defined `fast2` term: `DeltaV(t) = A2 * (1 - exp(-k2 * t))` with `A1 = 0` and `L = 0`.",
        f"- For every case, both injection time and total simulation time are fixed at `{GLOBAL_INJECTION_TIME:.0e} s`; only the fracture crossing time changes.",
        "",
        "## Main Findings",
        "",
    ]

    for tt in sorted(by_time):
        data = sorted(by_time[tt], key=lambda item: int(item["segments"]))
        best = min(data, key=lambda item: abs(float(item["relative_error"])))
        coarse = data[0]
        fine = data[-1]
        lines.extend(
            [
                f"- `t_cross = {tt:.0e} s`, `Re = {float(best['reynolds']):.3g}`, `Da = {float(best['damkohler']):.3g}`.",
                f"  Reference dissolved volume = {float(best['reference_dissolved_volume_m3']):.6e} m^3.",
                f"  `n = {int(coarse['segments'])}` gives error {100.0 * float(coarse['relative_error']):+.3f}%, while `n = {int(fine['segments'])}` gives {100.0 * float(fine['relative_error']):+.3f}%.",
                f"  Best case in this sweep is `n = {int(best['segments'])}` with error {100.0 * float(best['relative_error']):+.3f}%.",
            ]
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Larger `Da` means more of the gypsum reaction occurs inside the fracture during one crossing, so any segment-averaging error becomes more visible.",
            "- Increasing segment count reduces the spatial lumping error because mineral-volume increments are applied over shorter pieces of the fracture.",
            "- This benchmark is therefore a clean way to choose a minimum segment number before moving back to complex DFNs.",
            "",
            "## Files",
            "",
            "- `summary.csv`: numeric summary for all 20 cases.",
            "- `figures/`: benchmark plots.",
            "- `cases/<case>/raw_output/`: copied solver outputs for each run.",
        ]
    )

    (RESULTS_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    lines_cn = [
        "# 10 cm 单裂缝纯石膏 Benchmark",
        "",
        "## 这个 Benchmark 的意义",
        "",
        "- 它把问题压缩到最简单的几何条件：10 cm 直裂缝、唯一流路、无复杂 DFN 拓扑影响。",
        "- 因为这里只保留单一裂缝并把基质扩散压到近乎零，所以结果差异主要来自反应停留时间和 segment 离散，而不是网络连通性。",
        "- 当前活跃 chemistry 被改成仅保留代码中 `fast2` 项的纯石膏形式：`DeltaV(t) = A2 * (1 - exp(-k2 * t))`，并令 `A1 = 0`、`L = 0`。",
        f"- 所有 case 的注入时间和总模拟时间都固定为 `{GLOBAL_INJECTION_TIME:.0e} s`，只有裂缝穿越时间不同。",
        "- 这里采用的参考解对应程序实际输出的最后一个 snapshot，即 `t = t_injection = t_cross` 时系统内累计的总溶解体积。",
        "",
        "## 主要结论",
        "",
    ]

    for tt in sorted(by_time):
        data = sorted(by_time[tt], key=lambda item: int(item["segments"]))
        best = min(data, key=lambda item: abs(float(item["relative_error"])))
        coarse = data[0]
        fine = data[-1]
        lines_cn.extend(
            [
                f"- `t_cross = {tt:.0e} s`，`Re = {float(best['reynolds']):.3g}`，`Da = {float(best['damkohler']):.3g}`。",
                f"  参考溶解体积 = {float(best['reference_dissolved_volume_m3']):.6e} m^3。",
                f"  `n = {int(coarse['segments'])}` 时误差为 {100.0 * float(coarse['relative_error']):+.3f}% ，`n = {int(fine['segments'])}` 时为 {100.0 * float(fine['relative_error']):+.3f}% 。",
                f"  这一组里最佳离散是 `n = {int(best['segments'])}`，误差 {100.0 * float(best['relative_error']):+.3f}% 。",
            ]
        )

    lines_cn.extend(
        [
            "",
            "## 解释",
            "",
            "- `Da` 越大，粒子在一次穿越中经历的反应越充分，segment 平均化造成的误差也越容易显现。",
            "- 随着 segment 数增加，矿物体积变化被分配到更短的空间单元上，空间离散误差整体下降。",
            "- 因此，这个 benchmark 可以用来先确定单裂缝下所需的最小 segment 数，再回到复杂 DFN 计算。",
            "",
            "## 文件",
            "",
            "- `summary.csv`：20 组工况的数值汇总。",
            "- `figures/`：解释性图件。",
            "- `cases/<case>/raw_output/`：每组算例复制出来的原始求解器输出。",
        ]
    )

    (RESULTS_DIR / "README_CN.md").write_text("\n".join(lines_cn) + "\n", encoding="utf-8")


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
    write_note(rows)
    print(f"Wrote analysis products to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
