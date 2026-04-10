#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, TwoSlopeNorm


ROOT = Path(__file__).resolve().parent
PARTICLE_PATTERN = re.compile(r"particle_positions_t(\d+)\.csv$")
GROUPS: dict[str, list[str]] = {
    "Np_group": ["Np1", "Np4", "Np2", "V4"],
    "V_group": ["V1", "V2", "V3", "V4"],
    "tg_group": ["tg1", "tg2", "V4", "tg4"],
}
CASE_LEGEND_LABELS: dict[str, str] = {
    "V1": "P=1e-1, Np=1e6",
    "V2": "P=1e-2, Np=1e6",
    "V3": "P=1e-3, Np=1e6",
    "V4": "P=1e-4, Np=1e6",
    "Np1": "P=1e-4, Np=1e4",
    "Np4": "P=1e-4, Np=5e4",
    "Np2": "P=1e-4, Np=1e5",
    "Test 1": "P=1e-4, Np=1e6",
    "tg1": "P=1e-4, tg=1 s",
    "tg2": "P=1e-4, tg=10 s",
    "tg4": "P=1e-4, tg=1000 s",
}

plt.style.use("seaborn-v0_8-whitegrid")


def load_table(file_path: Path, expected_cols: int) -> np.ndarray:
    data = np.loadtxt(file_path)
    if data.ndim != 2 or data.shape[1] != expected_cols:
        raise ValueError(f"{file_path} should have {expected_cols} columns, got {data.shape}")
    return data


def build_segments(data: np.ndarray) -> np.ndarray:
    return np.stack((data[:, 0:2], data[:, 2:4]), axis=1)


def segment_key(row: np.ndarray) -> tuple[float, float, float, float]:
    return tuple(np.round(row[:4], 12))


def reversed_segment_key(row: np.ndarray) -> tuple[float, float, float, float]:
    return tuple(np.round(row[[2, 3, 0, 1]], 12))


def describe(name: str, values: np.ndarray) -> list[str]:
    return [
        f"{name}:",
        f"  count = {values.size}",
        f"  min = {values.min():.12g}",
        f"  max = {values.max():.12g}",
        f"  mean = {values.mean():.12g}",
        f"  median = {np.median(values):.12g}",
        f"  std = {values.std():.12g}",
    ]


def add_segment_plot(
    ax: plt.Axes,
    segments: np.ndarray,
    values: np.ndarray,
    cmap: str,
    norm: Normalize,
    title: str,
    colorbar_label: str,
) -> None:
    collection = LineCollection(segments, cmap=cmap, norm=norm, linewidths=1.6)
    collection.set_array(values)
    ax.add_collection(collection)
    ax.autoscale()
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True, linestyle="--", alpha=0.25)
    plt.colorbar(collection, ax=ax, fraction=0.046, pad=0.04, label=colorbar_label)


def compute_corrected_delta(init_data: np.ndarray, delta_data: np.ndarray) -> np.ndarray:
    init_aperture_by_segment: dict[tuple[float, float, float, float], float] = {}
    for row in init_data:
        init_aperture_by_segment[segment_key(row)] = float(row[4])
        init_aperture_by_segment[reversed_segment_key(row)] = float(row[4])

    corrected = delta_data.copy()
    corrected_change = np.empty(delta_data.shape[0], dtype=float)

    for idx, row in enumerate(delta_data):
        key = segment_key(row)
        b_init = init_aperture_by_segment.get(key)
        if b_init is None:
            raise ValueError(f"Could not match segment by coordinates in row {idx} for case data")
        corrected_change[idx] = float(row[4]) - b_init

    corrected[:, 6] = corrected_change
    return corrected


def read_particle_frames(case_dir: Path) -> list[tuple[Path, np.ndarray]]:
    frames: list[tuple[Path, np.ndarray]] = []
    for path in sorted(case_dir.glob("particle_positions_t*.csv")):
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
        if not rows:
            points = np.empty((0, 2), dtype=float)
        else:
            points = np.array([[float(row["x"]), float(row["y"])] for row in rows], dtype=float)
        frames.append((path, points))
    return frames


def frame_label(path: Path) -> str:
    match = PARTICLE_PATTERN.search(path.name)
    return f"t = {int(match.group(1))}" if match else path.stem


def generate_case_outputs(case_dir: Path) -> None:
    init_file = case_dir / "DFN_init.txt"
    delta_file = case_dir / "DFN_aperture_delta.txt"
    if not init_file.exists() or not delta_file.exists():
        return

    init_data = load_table(init_file, expected_cols=6)
    raw_delta_data = load_table(delta_file, expected_cols=7)
    corrected_delta_data = compute_corrected_delta(init_data, raw_delta_data)
    corrected_delta_file = case_dir / "DFN_aperture_delta_corrected.txt"
    np.savetxt(corrected_delta_file, corrected_delta_data, fmt="%.12g")

    init_segments = build_segments(init_data)
    delta_segments = build_segments(corrected_delta_data)
    init_aperture = init_data[:, 4]
    aperture_change = corrected_delta_data[:, 6]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    add_segment_plot(
        axes[0],
        init_segments,
        init_aperture,
        cmap="viridis",
        norm=Normalize(vmin=float(init_aperture.min()), vmax=float(init_aperture.max())),
        title="Initial aperture",
        colorbar_label="aperture",
    )
    delta_abs_max = float(np.max(np.abs(aperture_change)))
    add_segment_plot(
        axes[1],
        delta_segments,
        aperture_change,
        cmap="coolwarm",
        norm=TwoSlopeNorm(vmin=-delta_abs_max, vcenter=0.0, vmax=delta_abs_max) if delta_abs_max > 0 else Normalize(vmin=0.0, vmax=1.0),
        title="Aperture change",
        colorbar_label="aperture_change",
    )
    fig.suptitle(f"{case_dir.name}: DFN aperture and aperture change")
    fig.savefig(case_dir / "dfn_aperture_and_change.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    (case_dir / "dfn_aperture_and_change_summary.txt").write_text(
        "\n".join(
            [
                f"Init file: {init_file.name}",
                f"Delta file: {delta_file.name}",
                f"Corrected delta file: {corrected_delta_file.name}",
                "Aperture change is recomputed by coordinate matching against DFN_init.txt.",
                "",
                *describe("Initial aperture", init_aperture),
                "",
                *describe("Aperture change", aperture_change),
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    init_velocity = init_data[:, 5]
    updated_velocity = corrected_delta_data[:, 5]
    vel_min = float(min(init_velocity.min(), updated_velocity.min()))
    vel_max = float(max(init_velocity.max(), updated_velocity.max()))
    vel_norm = Normalize(vmin=vel_min, vmax=vel_max)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    add_segment_plot(axes[0], init_segments, init_velocity, "plasma", vel_norm, "Initial velocity", "velocity")
    add_segment_plot(axes[1], delta_segments, updated_velocity, "plasma", vel_norm, "Updated velocity", "velocity")
    fig.suptitle(f"{case_dir.name}: DFN velocity comparison")
    fig.savefig(case_dir / "dfn_velocity_comparison.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    (case_dir / "dfn_velocity_comparison_summary.txt").write_text(
        "\n".join(
            [
                f"Initial file: {init_file.name}",
                f"Updated file: {corrected_delta_file.name}",
                "",
                *describe("Initial velocity", init_velocity),
                "",
                *describe("Updated velocity", updated_velocity),
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    frames = read_particle_frames(case_dir)
    if frames:
        fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
        collection = LineCollection(
            init_segments,
            cmap="viridis",
            norm=Normalize(vmin=float(init_aperture.min()), vmax=float(init_aperture.max())),
            linewidths=1.5,
        )
        collection.set_array(init_aperture)
        ax.add_collection(collection)
        ax.autoscale()
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.grid(True, linestyle="--", alpha=0.25)
        plt.colorbar(collection, ax=ax, fraction=0.046, pad=0.04, label="initial aperture")
        scatter = ax.scatter([], [], s=14, c="black", edgecolors="white", linewidths=0.35, zorder=3)
        title = ax.set_title("")

        def update(frame_index: int):
            frame_path, points = frames[frame_index]
            scatter.set_offsets(points if points.size else np.empty((0, 2), dtype=float))
            title.set_text(f"{case_dir.name}: particle positions on initial aperture ({frame_label(frame_path)}, n={len(points)})")
            return scatter, title

        animation = FuncAnimation(fig, update, frames=len(frames), interval=600, blit=False, repeat=True)
        animation.save(case_dir / "particle_positions_on_initial_aperture.gif", writer=PillowWriter(fps=2))
        plt.close(fig)

        counts = [len(points) for _, points in frames]
        (case_dir / "particle_positions_summary.txt").write_text(
            "\n".join(
                [
                    f"Frames = {len(frames)}",
                    f"Min particles/frame = {min(counts)}",
                    f"Max particles/frame = {max(counts)}",
                    f"Mean particles/frame = {np.mean(counts):.12g}",
                ]
            )
            + "\n",
            encoding="utf-8",
        )


def read_two_column_file(file_path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(file_path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected format in {file_path}")
    return data[:, 0], data[:, 1]


def max_abs_difference(curve_a: tuple[np.ndarray, np.ndarray], curve_b: tuple[np.ndarray, np.ndarray]) -> float:
    x_a, y_a = curve_a
    x_b, y_b = curve_b
    interp_b = np.interp(x_a, x_b, y_b)
    return float(np.max(np.abs(y_a - interp_b)))


def generate_group_comparison(group_name: str, case_names: list[str]) -> None:
    available_cases = [ROOT / case_name for case_name in case_names if (ROOT / case_name / "cdf.txt").exists() and (ROOT / case_name / "pdf.txt").exists()]
    missing_cases = [case_name for case_name in case_names if not (ROOT / case_name / "cdf.txt").exists() or not (ROOT / case_name / "pdf.txt").exists()]
    if not available_cases:
        return

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharex=True, constrained_layout=True)
    cdf_ax, pdf_ax = axes
    cdf_curves: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    pdf_curves: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    for case_dir in available_cases:
        cdf_curve = read_two_column_file(case_dir / "cdf.txt")
        pdf_curve = read_two_column_file(case_dir / "pdf.txt")
        cdf_curves[case_dir.name] = cdf_curve
        pdf_curves[case_dir.name] = pdf_curve
        cdf_ax.plot(cdf_curve[0], cdf_curve[1], linewidth=2, label=case_dir.name)
        pdf_ax.plot(pdf_curve[0], pdf_curve[1], linewidth=2, label=case_dir.name)

    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("Residence time / x")
        ax.grid(True, which="both", linestyle="--", alpha=0.35)

    cdf_ax.set_ylabel("CDF")
    cdf_ax.set_title(f"CDF comparison: {group_name}")
    cdf_ax.legend(title="Case")
    pdf_ax.set_ylabel("PDF")
    pdf_ax.set_title(f"PDF comparison: {group_name}")
    pdf_ax.legend(title="Case")

    fig.savefig(ROOT / f"cdf_pdf_comparison_{group_name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    lines = ["metric\tcase_a\tcase_b\tmax_abs_diff"]
    case_labels = list(cdf_curves.keys())
    for i, case_a in enumerate(case_labels):
        for case_b in case_labels[i + 1 :]:
            lines.append(f"CDF\t{case_a}\t{case_b}\t{max_abs_difference(cdf_curves[case_a], cdf_curves[case_b]):.12g}")
            lines.append(f"PDF\t{case_a}\t{case_b}\t{max_abs_difference(pdf_curves[case_a], pdf_curves[case_b]):.12g}")
    if missing_cases:
        lines.append("")
        lines.append("missing_cases")
        lines.extend(missing_cases)
    (ROOT / f"cdf_pdf_diff_{group_name}.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_delta_aperture_histogram(group_name: str, case_names: list[str]) -> None:
    available_cases = []
    missing_cases = []
    for case_name in case_names:
        corrected_path = ROOT / case_name / "DFN_aperture_delta_corrected.txt"
        if corrected_path.exists():
            available_cases.append((case_name, corrected_path))
        else:
            missing_cases.append(case_name)

    if not available_cases:
        return

    labels_map = DELTA_GROUP_PLOT_LABELS.get(group_name, {})
    labels: list[str] = []
    means: list[float] = []

    for case_name, corrected_path in available_cases:
        data = np.loadtxt(corrected_path)
        labels.append(labels_map.get(case_name, case_name))
        means.append(float(np.mean(data[:, 6])))

    fig, ax = plt.subplots(figsize=(8.2, 5.2), constrained_layout=True)
    bars = ax.bar(labels, means, color=["#3a7d44", "#669bbc", "#c17c74", "#9d4edd"][: len(labels)], edgecolor="black", linewidth=0.8)
    ax.set_ylabel("Mean delta aperture")
    ax.set_title(f"Mean delta aperture comparison: {group_name}")
    ax.grid(True, axis="y", linestyle="--", alpha=0.35)

    ymax = max(means) if means else 0.0
    for bar, mean in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + ymax * 0.02 if ymax > 0 else 0.0, f"{mean:.3e}", ha="center", va="bottom", fontsize=9)

    if missing_cases:
        ax.text(0.99, 0.98, f"Missing: {', '.join(missing_cases)}", transform=ax.transAxes, ha="right", va="top", fontsize=9)

    fig.savefig(ROOT / f"mean_delta_aperture_hist_{group_name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    lines = ["plot_label\tcase_dir\tmean_delta_aperture"]
    for label, (case_name, corrected_path) in zip(labels, available_cases):
        data = np.loadtxt(corrected_path)
        lines.append(f"{label}\t{case_name}\t{float(np.mean(data[:, 6])):.12g}")
    if missing_cases:
        lines.append("")
        lines.append("missing_cases")
        lines.extend(missing_cases)
    (ROOT / f"mean_delta_aperture_hist_{group_name}.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")


def get_group_case_label(group_name: str, case_name: str) -> str:
    if group_name == "Np_group" and case_name == "V4":
        return "P=1e-4, Np=1e6"
    if group_name == "tg_group" and case_name == "V4":
        return "P=1e-4, tg=100 s"
    return CASE_LEGEND_LABELS.get(case_name, case_name)


def load_case_geometry_metrics(case_dir: Path) -> dict[str, np.ndarray]:
    init_data = load_table(case_dir / "DFN_init.txt", expected_cols=6)
    corrected_data = load_table(case_dir / "DFN_aperture_delta_corrected.txt", expected_cols=7)
    lengths = np.sqrt((corrected_data[:, 2] - corrected_data[:, 0]) ** 2 + (corrected_data[:, 3] - corrected_data[:, 1]) ** 2)
    return {
        "delta_b": corrected_data[:, 6],
        "initial_velocity_abs": np.abs(init_data[: corrected_data.shape[0], 5]),
        "length": lengths,
    }


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    total_weight = float(np.sum(weights))
    if total_weight <= 0.0:
        return float("nan")
    return float(np.sum(values * weights) / total_weight)


def weighted_fraction(mask: np.ndarray, weights: np.ndarray) -> float:
    total_weight = float(np.sum(weights))
    if total_weight <= 0.0:
        return float("nan")
    return float(np.sum(weights[mask]) / total_weight)


def top_contribution_share(values: np.ndarray, weights: np.ndarray, top_fraction: float = 0.1) -> float:
    if values.size == 0:
        return float("nan")
    contribution = values * weights
    total = float(np.sum(contribution))
    if total <= 0.0:
        return 0.0
    count = max(1, int(math.ceil(values.size * top_fraction)))
    order = np.argsort(values)[::-1]
    top_total = float(np.sum(contribution[order[:count]]))
    return top_total / total


def cumulative_contribution_curve(values: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if values.size == 0:
        return np.array([0.0]), np.array([0.0])
    contribution = values * weights
    total = float(np.sum(contribution))
    if total <= 0.0:
        x = np.linspace(0.0, 1.0, values.size + 1)
        return x, np.zeros(values.size + 1, dtype=float)
    order = np.argsort(values)[::-1]
    sorted_contribution = contribution[order]
    cumulative = np.cumsum(sorted_contribution) / total
    x = np.arange(1, values.size + 1, dtype=float) / values.size
    return np.concatenate(([0.0], x)), np.concatenate(([0.0], cumulative))


def generate_delta_distribution_comparison(group_name: str, case_names: list[str]) -> None:
    available_cases = [ROOT / case_name for case_name in case_names if (ROOT / case_name / "DFN_aperture_delta_corrected.txt").exists()]
    missing_cases = [case_name for case_name in case_names if not (ROOT / case_name / "DFN_aperture_delta_corrected.txt").exists()]
    if not available_cases:
        return

    datasets: dict[str, np.ndarray] = {}
    global_max = 0.0
    for case_dir in available_cases:
        delta_b = load_table(case_dir / "DFN_aperture_delta_corrected.txt", expected_cols=7)[:, 6]
        datasets[case_dir.name] = delta_b
        global_max = max(global_max, float(np.max(delta_b)))

    bins = np.linspace(0.0, global_max if global_max > 0 else 1.0, 80)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)
    cdf_ax, pdf_ax = axes

    for case_name, delta_b in datasets.items():
        label = get_group_case_label(group_name, case_name)
        sorted_delta = np.sort(delta_b)
        cdf = np.arange(1, sorted_delta.size + 1, dtype=float) / sorted_delta.size
        counts, edges = np.histogram(delta_b, bins=bins, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])

        cdf_ax.plot(sorted_delta, cdf, linewidth=2, label=label)
        pdf_ax.plot(centers, counts, linewidth=2, label=label)

    for ax in axes:
        ax.set_xlabel("Delta aperture")
        ax.grid(True, linestyle="--", alpha=0.35)

    cdf_ax.set_ylabel("CDF")
    cdf_ax.set_title(f"Delta aperture CDF: {group_name}")
    cdf_ax.legend(title="Case")
    pdf_ax.set_ylabel("PDF")
    pdf_ax.set_title(f"Delta aperture PDF: {group_name}")
    pdf_ax.legend(title="Case")

    if missing_cases:
        pdf_ax.text(0.99, 0.98, f"Missing: {', '.join(missing_cases)}", transform=pdf_ax.transAxes, ha="right", va="top", fontsize=9)

    fig.savefig(ROOT / f"delta_aperture_cdf_pdf_{group_name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def generate_delta_velocity_scatter(group_name: str, case_names: list[str]) -> None:
    available_cases = [ROOT / case_name for case_name in case_names if (ROOT / case_name / "DFN_aperture_delta_corrected.txt").exists() and (ROOT / case_name / "DFN_init.txt").exists()]
    if not available_cases:
        return

    ncols = 2
    nrows = int(math.ceil(len(available_cases) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.4 * ncols, 4.8 * nrows), constrained_layout=True, squeeze=False)

    x_max = 0.0
    y_max = 0.0
    metrics_by_case: dict[str, dict[str, np.ndarray]] = {}
    for case_dir in available_cases:
        metrics = load_case_geometry_metrics(case_dir)
        metrics_by_case[case_dir.name] = metrics
        x_max = max(x_max, float(np.max(metrics["initial_velocity_abs"])))
        y_max = max(y_max, float(np.max(metrics["delta_b"])))

    for ax in axes.flat:
        ax.set_visible(False)

    for idx, case_dir in enumerate(available_cases):
        ax = axes[idx // ncols][idx % ncols]
        ax.set_visible(True)
        metrics = metrics_by_case[case_dir.name]
        ax.scatter(metrics["initial_velocity_abs"], metrics["delta_b"], s=8, alpha=0.28, color="#1d3557", edgecolors="none")
        ax.set_title(get_group_case_label(group_name, case_dir.name))
        ax.set_xlabel("|Initial velocity|")
        ax.set_ylabel("Delta aperture")
        ax.set_xlim(0.0, x_max * 1.02 if x_max > 0 else 1.0)
        ax.set_ylim(0.0, y_max * 1.02 if y_max > 0 else 1.0)
        ax.grid(True, linestyle="--", alpha=0.3)

    fig.savefig(ROOT / f"delta_aperture_vs_initial_velocity_{group_name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def generate_delta_geometry_summary(group_name: str, case_names: list[str]) -> None:
    lines = [
        "plot_label\tcase_dir\tcount\tmean_delta_b\tlength_weighted_mean_delta_b\tmedian_delta_b\tp90_delta_b\tmax_delta_b\tchanged_fraction\tlength_weighted_changed_fraction\ttop10pct_contribution_share"
    ]
    missing_cases: list[str] = []

    for case_name in case_names:
        case_dir = ROOT / case_name
        corrected = case_dir / "DFN_aperture_delta_corrected.txt"
        init_file = case_dir / "DFN_init.txt"
        if not corrected.exists() or not init_file.exists():
            missing_cases.append(case_name)
            continue

        metrics = load_case_geometry_metrics(case_dir)
        delta_b = metrics["delta_b"]
        lengths = metrics["length"]
        changed = delta_b > 0.0

        lines.append(
            "\t".join(
                [
                    get_group_case_label(group_name, case_name),
                    case_name,
                    str(delta_b.size),
                    f"{float(np.mean(delta_b)):.12g}",
                    f"{weighted_mean(delta_b, lengths):.12g}",
                    f"{float(np.median(delta_b)):.12g}",
                    f"{float(np.quantile(delta_b, 0.9)):.12g}",
                    f"{float(np.max(delta_b)):.12g}",
                    f"{float(np.mean(changed)):.12g}",
                    f"{weighted_fraction(changed, lengths):.12g}",
                    f"{top_contribution_share(delta_b, lengths, top_fraction=0.1):.12g}",
                ]
            )
        )

    if missing_cases:
        lines.append("")
        lines.append("missing_cases")
        lines.extend(missing_cases)

    (ROOT / f"delta_aperture_geometry_summary_{group_name}.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_cumulative_contribution_plot(group_name: str, case_names: list[str]) -> None:
    available_cases = [ROOT / case_name for case_name in case_names if (ROOT / case_name / "DFN_aperture_delta_corrected.txt").exists()]
    missing_cases = [case_name for case_name in case_names if not (ROOT / case_name / "DFN_aperture_delta_corrected.txt").exists()]
    if not available_cases:
        return

    fig, ax = plt.subplots(figsize=(8.4, 5.4), constrained_layout=True)
    lines = ["plot_label\tcase_dir\tx_fraction_segments\tcumulative_contribution_share"]

    for case_dir in available_cases:
        metrics = load_case_geometry_metrics(case_dir)
        x, y = cumulative_contribution_curve(metrics["delta_b"], metrics["length"])
        label = get_group_case_label(group_name, case_dir.name)
        ax.plot(x * 100.0, y, linewidth=2, label=label)
        sample_points = np.linspace(0, len(x) - 1, min(101, len(x)), dtype=int)
        for idx in np.unique(sample_points):
            lines.append(f"{label}\t{case_dir.name}\t{x[idx]:.12g}\t{y[idx]:.12g}")

    ax.set_xlabel("Top-ranked fracture segments (%)")
    ax.set_ylabel("Cumulative share of total widening")
    ax.set_title(f"Length-weighted cumulative contribution: {group_name}")
    ax.set_xlim(0.0, 100.0)
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(title="Case")

    if missing_cases:
        ax.text(0.99, 0.02, f"Missing: {', '.join(missing_cases)}", transform=ax.transAxes, ha="right", va="bottom", fontsize=9)
        lines.append("")
        lines.append("missing_cases")
        lines.extend(missing_cases)

    fig.savefig(ROOT / f"delta_aperture_cumulative_contribution_{group_name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    (ROOT / f"delta_aperture_cumulative_contribution_{group_name}.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_delta_aperture_difference_map(case_a: str, case_b: str) -> None:
    file_a = ROOT / case_a / "DFN_aperture_delta_corrected.txt"
    file_b = ROOT / case_b / "DFN_aperture_delta_corrected.txt"
    if not file_a.exists() or not file_b.exists():
        return

    data_a = load_table(file_a, expected_cols=7)
    data_b = load_table(file_b, expected_cols=7)

    map_b: dict[tuple[float, float, float, float], np.ndarray] = {}
    for row in data_b:
        map_b[segment_key(row)] = row
        map_b[reversed_segment_key(row)] = row

    common_rows: list[np.ndarray] = []
    diff_values: list[float] = []
    missing_in_b = 0
    for row_a in data_a:
        row_b = map_b.get(segment_key(row_a))
        if row_b is None:
            missing_in_b += 1
            continue
        common_rows.append(row_a)
        diff_values.append(float(row_a[6] - row_b[6]))

    if not common_rows:
        return

    common_data = np.vstack(common_rows)
    diff_array = np.asarray(diff_values, dtype=float)
    segments = build_segments(common_data)
    diff_abs_max = float(np.max(np.abs(diff_array)))

    fig, ax = plt.subplots(figsize=(8.4, 6.4), constrained_layout=True)
    add_segment_plot(
        ax,
        segments,
        diff_array,
        cmap="coolwarm",
        norm=TwoSlopeNorm(vmin=-diff_abs_max, vcenter=0.0, vmax=diff_abs_max) if diff_abs_max > 0 else Normalize(vmin=0.0, vmax=1.0),
        title=f"Delta aperture difference: {case_a} - {case_b}",
        colorbar_label=f"delta_b({case_a}) - delta_b({case_b})",
    )
    fig.savefig(ROOT / f"delta_aperture_difference_{case_a}_minus_{case_b}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    summary_lines = [
        f"Case A: {case_a}",
        f"Case B: {case_b}",
        f"Common segments compared = {diff_array.size}",
        f"Segments in {case_a} without match in {case_b} = {missing_in_b}",
        "",
        *describe(f"delta_b({case_a}) - delta_b({case_b})", diff_array),
    ]
    (ROOT / f"delta_aperture_difference_{case_a}_minus_{case_b}_summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


def matched_delta_difference(case_a: str, case_b: str) -> tuple[np.ndarray, np.ndarray, int]:
    file_a = ROOT / case_a / "DFN_aperture_delta_corrected.txt"
    file_b = ROOT / case_b / "DFN_aperture_delta_corrected.txt"
    data_a = load_table(file_a, expected_cols=7)
    data_b = load_table(file_b, expected_cols=7)

    map_b: dict[tuple[float, float, float, float], np.ndarray] = {}
    for row in data_b:
        map_b[segment_key(row)] = row
        map_b[reversed_segment_key(row)] = row

    common_rows: list[np.ndarray] = []
    diff_values: list[float] = []
    missing_in_b = 0
    for row_a in data_a:
        row_b = map_b.get(segment_key(row_a))
        if row_b is None:
            missing_in_b += 1
            continue
        common_rows.append(row_a)
        diff_values.append(float(row_a[6] - row_b[6]))

    if not common_rows:
        return np.empty((0, 7), dtype=float), np.empty(0, dtype=float), missing_in_b
    return np.vstack(common_rows), np.asarray(diff_values, dtype=float), missing_in_b


def generate_tg_convergence_analysis(reference_case: str = "tg1", compare_cases: list[str] | None = None) -> None:
    if compare_cases is None:
        compare_cases = ["tg2", "V4", "tg4"]
    if not (ROOT / reference_case / "DFN_aperture_delta_corrected.txt").exists():
        return

    lines = ["plot_label\tcase_dir\tcommon_segments\tmean_abs_diff\tmax_abs_diff\trms_diff\tmissing_segments_vs_ref"]
    x_vals: list[float] = [1.0]
    mean_abs_vals: list[float] = [0.0]
    rms_vals: list[float] = [0.0]
    max_abs_vals: list[float] = [0.0]

    tg_value_map = {"tg1": 1.0, "tg2": 10.0, "V4": 100.0, "tg4": 1000.0}
    ref_label = get_group_case_label("tg_group", reference_case)
    ref_count = load_table(ROOT / reference_case / "DFN_aperture_delta_corrected.txt", expected_cols=7).shape[0]
    lines.append(f"{ref_label}\t{reference_case}\t{ref_count}\t0\t0\t0\t0")

    for case_name in compare_cases:
        case_file = ROOT / case_name / "DFN_aperture_delta_corrected.txt"
        if not case_file.exists():
            continue
        _, diff, missing = matched_delta_difference(case_name, reference_case)
        if diff.size == 0:
            continue
        abs_diff = np.abs(diff)
        mean_abs = float(np.mean(abs_diff))
        max_abs = float(np.max(abs_diff))
        rms = float(np.sqrt(np.mean(diff * diff)))
        label = get_group_case_label("tg_group", case_name)
        lines.append(f"{label}\t{case_name}\t{diff.size}\t{mean_abs:.12g}\t{max_abs:.12g}\t{rms:.12g}\t{missing}")
        x_vals.append(tg_value_map[case_name])
        mean_abs_vals.append(mean_abs)
        rms_vals.append(rms)
        max_abs_vals.append(max_abs)

    if not x_vals:
        return

    fig, ax = plt.subplots(figsize=(8.4, 5.2), constrained_layout=True)
    ax.plot(x_vals, mean_abs_vals, marker="o", linewidth=2, label="mean |Δb - Δb_ref| (ref: tg=1 s)")
    ax.plot(x_vals, rms_vals, marker="s", linewidth=2, label="RMS error (ref: tg=1 s)")
    ax.plot(x_vals, max_abs_vals, marker="^", linewidth=2, label="max |Δb - Δb_ref| (ref: tg=1 s)")
    ax.set_xscale("log")
    ax.set_xlabel("Geometry update interval tg (s)")
    ax.set_ylabel("Difference relative to tg = 1 s")
    ax.set_title("Geometry update interval convergence")
    ax.grid(True, which="both", linestyle="--", alpha=0.35)
    ax.legend()
    fig.savefig(ROOT / "tg_group_delta_aperture_convergence.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    (ROOT / "tg_group_delta_aperture_convergence.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_tg_cdf_error_plot(reference_case: str = "tg1", compare_cases: list[str] | None = None) -> None:
    if compare_cases is None:
        compare_cases = ["tg2", "V4", "tg4"]
    ref_file = ROOT / reference_case / "cdf.txt"
    if not ref_file.exists():
        return

    ref_x, ref_y = read_two_column_file(ref_file)
    tg_value_map = {"tg1": 1.0, "tg2": 10.0, "V4": 100.0, "tg4": 1000.0}

    x_vals = [tg_value_map[reference_case]]
    max_abs_vals = [0.0]
    mean_abs_vals = [0.0]
    rms_vals = [0.0]
    lines = ["plot_label\tcase_dir\tmax_abs_cdf_diff\tmean_abs_cdf_diff\trms_cdf_diff"]
    lines.append(f"{get_group_case_label('tg_group', reference_case)}\t{reference_case}\t0\t0\t0")

    for case_name in compare_cases:
        case_file = ROOT / case_name / "cdf.txt"
        if not case_file.exists():
            continue
        x, y = read_two_column_file(case_file)
        interp = np.interp(ref_x, x, y)
        diff = interp - ref_y
        max_abs = float(np.max(np.abs(diff)))
        mean_abs = float(np.mean(np.abs(diff)))
        rms = float(np.sqrt(np.mean(diff * diff)))
        x_vals.append(tg_value_map[case_name])
        max_abs_vals.append(max_abs)
        mean_abs_vals.append(mean_abs)
        rms_vals.append(rms)
        lines.append(f"{get_group_case_label('tg_group', case_name)}\t{case_name}\t{max_abs:.12g}\t{mean_abs:.12g}\t{rms:.12g}")

    fig, ax = plt.subplots(figsize=(8.4, 5.2), constrained_layout=True)
    ax.plot(x_vals, max_abs_vals, marker="o", linewidth=2, label="max |CDF - CDF_ref| (ref: tg=1 s)")
    ax.plot(x_vals, mean_abs_vals, marker="s", linewidth=2, label="mean |CDF - CDF_ref| (ref: tg=1 s)")
    ax.plot(x_vals, rms_vals, marker="^", linewidth=2, label="RMS CDF error (ref: tg=1 s)")
    ax.set_xscale("log")
    ax.set_xlabel("Geometry update interval tg (s)")
    ax.set_ylabel("CDF error relative to tg = 1 s")
    ax.set_title("CDF error vs geometry update interval")
    ax.grid(True, which="both", linestyle="--", alpha=0.35)
    ax.legend()
    fig.savefig(ROOT / "tg_group_cdf_error_vs_tg.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    (ROOT / "tg_group_cdf_error_vs_tg.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_tg_difference_maps(reference_case: str = "tg1", compare_cases: list[str] | None = None) -> None:
    if compare_cases is None:
        compare_cases = ["tg2", "V4", "tg4"]
    if not (ROOT / reference_case / "DFN_aperture_delta_corrected.txt").exists():
        return

    for case_name in compare_cases:
        case_file = ROOT / case_name / "DFN_aperture_delta_corrected.txt"
        if not case_file.exists():
            continue
        common_data, diff_array, missing = matched_delta_difference(case_name, reference_case)
        if diff_array.size == 0:
            continue
        segments = build_segments(common_data)
        diff_abs_max = float(np.max(np.abs(diff_array)))

        fig, ax = plt.subplots(figsize=(8.4, 6.4), constrained_layout=True)
        add_segment_plot(
            ax,
            segments,
            diff_array,
            cmap="coolwarm",
            norm=TwoSlopeNorm(vmin=-diff_abs_max, vcenter=0.0, vmax=diff_abs_max) if diff_abs_max > 0 else Normalize(vmin=0.0, vmax=1.0),
            title=f"Delta aperture difference: {case_name} - {reference_case}",
            colorbar_label=f"delta_b({case_name}) - delta_b({reference_case})",
        )
        fig.savefig(ROOT / f"delta_aperture_difference_{case_name}_minus_{reference_case}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

        summary_lines = [
            f"Case A: {case_name}",
            f"Case B: {reference_case}",
            f"Common segments compared = {diff_array.size}",
            f"Segments in {case_name} without match in {reference_case} = {missing}",
            "",
            *describe(f"delta_b({case_name}) - delta_b({reference_case})", diff_array),
        ]
        (ROOT / f"delta_aperture_difference_{case_name}_minus_{reference_case}_summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


def main() -> None:
    for case_dir in sorted(path for path in ROOT.iterdir() if path.is_dir()):
        generate_case_outputs(case_dir)

    for group_name, case_names in GROUPS.items():
        generate_group_comparison(group_name, case_names)
        generate_delta_distribution_comparison(group_name, case_names)
        generate_delta_velocity_scatter(group_name, case_names)
        generate_delta_geometry_summary(group_name, case_names)
        generate_cumulative_contribution_plot(group_name, case_names)

    generate_delta_aperture_difference_map("V1", "V4")
    generate_tg_convergence_analysis(reference_case="tg1", compare_cases=["tg2", "V4", "tg4"])
    generate_tg_cdf_error_plot(reference_case="tg1", compare_cases=["tg2", "V4", "tg4"])
    generate_tg_difference_maps(reference_case="tg1", compare_cases=["tg2", "V4", "tg4"])


if __name__ == "__main__":
    main()
