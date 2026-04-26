#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.ticker import FuncFormatter


REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_TAG = os.environ.get(
    "RST_RUN_TAG",
    "song_selected_10cm_phreeqc_pipeline_dp5000Pa_tinj5000",
)
RUN_ROOT = REPO_ROOT / "Output" / RUN_TAG / "rtm_runs"
OUTPUT_ROOT = REPO_ROOT / "Output" / RUN_TAG / "figures"

ALL_CASES = [
    ("p9_a1p5", 9, 1.5),
    ("p12_a2p0", 12, 2.0),
    ("p16_a2p5", 16, 2.5),
]
case_filter = [name.strip() for name in os.environ.get("RST_CASES", "").split(",") if name.strip()]
if case_filter:
    case_lookup = {name: (name, p, a) for name, p, a in ALL_CASES}
    CASES = [case_lookup[name] for name in case_filter]
else:
    CASES = ALL_CASES

COORD_DECIMALS = 8


@dataclass(frozen=True)
class SegmentKey:
    p1: tuple[float, float]
    p2: tuple[float, float]


def format_a(a_value: float) -> str:
    return str(int(a_value)) if float(a_value).is_integer() else str(a_value)


def build_key(row: np.ndarray) -> SegmentKey:
    p1 = tuple(np.round(row[:2], COORD_DECIMALS))
    p2 = tuple(np.round(row[2:4], COORD_DECIMALS))
    ordered = tuple(sorted((p1, p2)))
    return SegmentKey(ordered[0], ordered[1])


def load_dfn(path: Path, expected_cols: int) -> np.ndarray:
    data = np.loadtxt(path)
    data = np.atleast_2d(data)
    if data.shape[1] != expected_cols:
        raise ValueError(f"Unexpected shape in {path}: {data.shape}, expected {expected_cols} columns")
    return data


def to_plot_segments(rows: np.ndarray) -> np.ndarray:
    plotted = rows[:, :4].copy()
    plotted[:, [0, 2]] += 0.05
    plotted[:, [1, 3]] += 0.05
    plotted *= 100.0
    return plotted.reshape(-1, 2, 2)


def aperture_formatter(value: float, _pos: int) -> str:
    return f"{value * 1.0e6:.0f}"


def delta_formatter(value: float, _pos: int) -> str:
    value_um = value * 1.0e6
    if abs(value_um) < 5.0e-7:
        return "0"
    if abs(value_um - round(value_um)) < 5.0e-7:
        return f"{int(round(value_um))}"
    if value_um >= 1.0:
        return f"{value_um:.1f}"
    if value_um >= 0.01:
        return f"{value_um:.2f}"
    return f"{value_um:.4f}"


def add_line_collection(ax, rows: np.ndarray, values: np.ndarray, cmap, norm, linewidth: float) -> LineCollection:
    lc = LineCollection(to_plot_segments(rows), cmap=cmap, norm=norm, linewidths=linewidth, capstyle="round")
    lc.set_array(values)
    ax.add_collection(lc)
    ax.set_xlim(0.0, 10.0)
    ax.set_ylim(0.0, 10.0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([0.0, 5.0, 10.0])
    ax.set_yticks([0.0, 5.0, 10.0])
    ax.tick_params(axis="both", labelsize=11, length=3, width=0.8, pad=4)
    for spine in ax.spines.values():
        spine.set_color("#b9b9b9")
        spine.set_linewidth(0.8)
    return lc


def collect_case_data(case_name: str) -> dict[str, np.ndarray]:
    output_dir = RUN_ROOT / case_name / "output"
    init = load_dfn(output_dir / "DFN_init.txt", expected_cols=6)
    final = load_dfn(output_dir / "DFN.txt", expected_cols=6)
    delta = load_dfn(output_dir / "DFN_aperture_delta.txt", expected_cols=7)

    init_map = {build_key(row): row for row in init}
    final_map = {build_key(row): row for row in final}
    delta_map = {build_key(row): row for row in delta}

    removed_rows = []
    diff_rows = []
    diff_values = []
    for key, init_row in init_map.items():
        if key in final_map:
            diff_rows.append(final_map[key])
            diff_values.append(delta_map[key][6])
        else:
            diff_rows.append(init_row)
            diff_values.append(-init_row[4])
            removed_rows.append(init_row)

    return {
        "init": init,
        "final": final,
        "diff_rows": np.asarray(diff_rows, dtype=float),
        "diff_values": np.asarray(diff_values, dtype=float),
        "removed": np.asarray(removed_rows, dtype=float) if removed_rows else np.empty((0, 6), dtype=float),
    }


def save_summary(summary_rows: list[dict[str, float | int | str]]) -> None:
    csv_path = OUTPUT_ROOT / "song_selected_10cm_phreeqc_rtm_geometry_summary.csv"
    fieldnames = [
        "case",
        "p",
        "a",
        "initial_segments",
        "final_segments",
        "removed_segments",
        "initial_aperture_mean_um",
        "final_aperture_mean_um",
        "delta_mean_um",
        "delta_min_um",
        "delta_max_um",
    ]
    with csv_path.open("w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "stix",
            "axes.linewidth": 0.8,
            "xtick.direction": "out",
            "ytick.direction": "out",
        }
    )

    case_data = {case_name: collect_case_data(case_name) for case_name, _, _ in CASES}

    summary_rows: list[dict[str, float | int | str]] = []
    for case_name, p_value, a_value in CASES:
        data = case_data[case_name]
        fig, axes = plt.subplots(
            1,
            3,
            figsize=(12.8, 3.8),
            gridspec_kw={"width_ratios": [1.0, 1.0, 1.0]},
        )
        fig.subplots_adjust(left=0.055, right=0.905, top=0.84, bottom=0.18, wspace=0.42)

        init_norm = Normalize(vmin=0.0, vmax=20.0e-6)
        final_norm = Normalize(vmin=0.0, vmax=20.0e-6)
        diff_values = data["diff_values"]
        diff_vmax = float(np.max(diff_values))
        if diff_vmax <= 0.0:
            diff_vmax = 1.0e-12
        diff_vmax_um = float(np.ceil(diff_vmax * 1.0e6))
        if diff_vmax_um <= 0.0:
            diff_vmax_um = 1.0
        diff_norm = Normalize(vmin=0.0, vmax=diff_vmax_um * 1.0e-6)

        lc_init = add_line_collection(
            axes[0],
            data["init"],
            data["init"][:, 4],
            cmap=plt.get_cmap("viridis"),
            norm=init_norm,
            linewidth=1.05,
        )
        axes[0].set_title("Initial geometry", fontsize=13)

        lc_final = add_line_collection(
            axes[1],
            data["final"],
            data["final"][:, 4],
            cmap=plt.get_cmap("viridis"),
            norm=final_norm,
            linewidth=1.05,
        )
        axes[1].set_title("Reacted geometry", fontsize=13)

        lc_diff = add_line_collection(
            axes[2],
            data["diff_rows"],
            data["diff_values"],
            cmap=plt.get_cmap("viridis"),
            norm=diff_norm,
            linewidth=1.15,
        )
        axes[2].set_title("Aperture change", fontsize=13)

        for idx, ax in enumerate(axes):
            ax.set_xlabel("x (cm)", fontsize=12, labelpad=3)
            if idx == 0:
                ax.set_ylabel("y (cm)", fontsize=12, labelpad=5)

        fig.suptitle(rf"$(p,\, a)=({p_value},\, {format_a(a_value)})$", fontsize=16, y=0.96)

        cbar_init = fig.colorbar(lc_init, ax=axes[0], fraction=0.04, pad=0.02)
        cbar_init.ax.yaxis.set_major_formatter(FuncFormatter(aperture_formatter))
        cbar_init.set_label("Initial aperture (μm)", fontsize=12, labelpad=6)
        cbar_init.ax.tick_params(labelsize=11, length=3, width=0.8)

        cbar_ap = fig.colorbar(lc_final, ax=axes[1], fraction=0.04, pad=0.04)
        cbar_ap.ax.yaxis.set_major_formatter(FuncFormatter(aperture_formatter))
        cbar_ap.set_label("Reacted aperture (μm)", fontsize=12, labelpad=6)
        cbar_ap.ax.tick_params(labelsize=11, length=3, width=0.8)

        cbar_diff = fig.colorbar(lc_diff, ax=axes[2], fraction=0.04, pad=0.04)
        cbar_diff.ax.yaxis.set_major_formatter(FuncFormatter(delta_formatter))
        cbar_diff.set_label("Δ aperture (μm)", fontsize=12, labelpad=6)
        cbar_diff.ax.tick_params(labelsize=11, length=3, width=0.8)

        png_path = OUTPUT_ROOT / f"{case_name}_geometry_triptych.png"
        pdf_path = OUTPUT_ROOT / f"{case_name}_geometry_triptych.pdf"
        fig.savefig(png_path, dpi=600, facecolor="white", bbox_inches="tight")
        fig.savefig(pdf_path, facecolor="white", bbox_inches="tight")
        plt.close(fig)

        summary_rows.append(
            {
                "case": case_name,
                "p": p_value,
                "a": a_value,
                "initial_segments": int(data["init"].shape[0]),
                "final_segments": int(data["final"].shape[0]),
                "removed_segments": int(data["removed"].shape[0]),
                "initial_aperture_mean_um": float(np.mean(data["init"][:, 4]) * 1.0e6),
                "final_aperture_mean_um": float(np.mean(data["final"][:, 4]) * 1.0e6),
                "delta_mean_um": float(np.mean(data["diff_values"]) * 1.0e6),
                "delta_min_um": float(np.min(data["diff_values"]) * 1.0e6),
                "delta_max_um": float(np.max(data["diff_values"]) * 1.0e6),
            }
        )

    save_summary(summary_rows)


if __name__ == "__main__":
    main()
