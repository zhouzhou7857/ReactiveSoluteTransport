#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, TwoSlopeNorm


ROOT = Path(__file__).resolve().parent
INIT_FILE = ROOT / "DFN_init.txt"
DELTA_FILE = ROOT / "DFN_aperture_delta.txt"
CORRECTED_DELTA_FILE = ROOT / "DFN_aperture_delta_corrected.txt"
OUTPUT_FIG = ROOT / "dfn_aperture_and_change.png"
OUTPUT_SUMMARY = ROOT / "dfn_aperture_and_change_summary.txt"


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


def compute_corrected_delta(init_data: np.ndarray, delta_data: np.ndarray) -> np.ndarray:
    init_aperture_by_segment: dict[tuple[float, float, float, float], float] = {}
    for row in init_data:
        init_aperture_by_segment[segment_key(row)] = float(row[4])
        init_aperture_by_segment[reversed_segment_key(row)] = float(row[4])

    corrected = delta_data.copy()
    corrected_change = np.empty(delta_data.shape[0], dtype=float)
    missing_segments: list[tuple[float, float, float, float]] = []

    for idx, row in enumerate(delta_data):
        key = segment_key(row)
        b_init = init_aperture_by_segment.get(key)
        if b_init is None:
            missing_segments.append(key)
            b_init = float("nan")
        corrected_change[idx] = float(row[4]) - b_init

    if missing_segments:
        raise ValueError(f"Could not match {len(missing_segments)} segments by coordinates")

    corrected[:, 6] = corrected_change
    np.savetxt(CORRECTED_DELTA_FILE, corrected, fmt="%.12g")
    return corrected


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


def main() -> None:
    init_data = load_table(INIT_FILE, expected_cols=6)
    raw_delta_data = load_table(DELTA_FILE, expected_cols=7)
    delta_data = compute_corrected_delta(init_data, raw_delta_data)

    init_segments = build_segments(init_data)
    delta_segments = build_segments(delta_data)
    init_aperture = init_data[:, 4]
    aperture_change = delta_data[:, 6]

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
        norm=TwoSlopeNorm(vmin=-delta_abs_max, vcenter=0.0, vmax=delta_abs_max),
        title="Aperture change",
        colorbar_label="aperture_change",
    )

    fig.suptitle("DFN aperture and aperture change")
    fig.savefig(OUTPUT_FIG, dpi=300, bbox_inches="tight")
    plt.close(fig)

    summary_lines = [
        f"Init file: {INIT_FILE.name}",
        f"Delta file: {DELTA_FILE.name}",
        f"Corrected delta file: {CORRECTED_DELTA_FILE.name}",
        "Aperture change is recomputed by coordinate matching against DFN_init.txt.",
        "",
        *describe("Initial aperture", init_aperture),
        "",
        *describe("Aperture change", aperture_change),
    ]
    OUTPUT_SUMMARY.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
