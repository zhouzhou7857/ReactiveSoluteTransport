#!/usr/bin/env python3

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize


ROOT = Path(__file__).resolve().parent
INIT_FILE = ROOT / "DFN_init.txt"
DELTA_FILE = ROOT / "DFN_aperture_delta.txt"
VELOCITY_FIG = ROOT / "dfn_velocity_comparison.png"
VELOCITY_SUMMARY = ROOT / "dfn_velocity_comparison_summary.txt"
PARTICLE_GIF = ROOT / "particle_positions_on_initial_aperture.gif"
PARTICLE_SUMMARY = ROOT / "particle_positions_summary.txt"
PARTICLE_PATTERN = re.compile(r"particle_positions_t(\d+)\.csv$")


def load_table(file_path: Path, expected_cols: int) -> np.ndarray:
    data = np.loadtxt(file_path)
    if data.ndim != 2 or data.shape[1] != expected_cols:
        raise ValueError(f"{file_path} should have {expected_cols} columns, got {data.shape}")
    return data


def build_segments(data: np.ndarray) -> np.ndarray:
    return np.stack((data[:, 0:2], data[:, 2:4]), axis=1)


def add_segment_plot(
    ax: plt.Axes,
    segments: np.ndarray,
    values: np.ndarray,
    norm: Normalize,
    title: str,
    colorbar_label: str,
) -> None:
    collection = LineCollection(segments, cmap="plasma", norm=norm, linewidths=1.6)
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


def plot_velocity() -> None:
    init_data = load_table(INIT_FILE, expected_cols=6)
    delta_data = load_table(DELTA_FILE, expected_cols=7)

    init_segments = build_segments(init_data)
    delta_segments = build_segments(delta_data)
    init_velocity = init_data[:, 5]
    updated_velocity = delta_data[:, 5]

    vel_min = float(min(init_velocity.min(), updated_velocity.min()))
    vel_max = float(max(init_velocity.max(), updated_velocity.max()))
    norm = Normalize(vmin=vel_min, vmax=vel_max)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    add_segment_plot(axes[0], init_segments, init_velocity, norm, "Initial velocity", "velocity")
    add_segment_plot(axes[1], delta_segments, updated_velocity, norm, "Updated velocity", "velocity")
    fig.suptitle("DFN velocity comparison")
    fig.savefig(VELOCITY_FIG, dpi=300, bbox_inches="tight")
    plt.close(fig)

    summary_lines = [
        f"Initial file: {INIT_FILE.name}",
        f"Updated file: {DELTA_FILE.name}",
        "",
        *describe("Initial velocity", init_velocity),
        "",
        *describe("Updated velocity", updated_velocity),
    ]
    VELOCITY_SUMMARY.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


def read_particle_frames() -> list[tuple[Path, np.ndarray]]:
    frames: list[tuple[Path, np.ndarray]] = []
    for path in sorted(ROOT.glob("particle_positions_t*.csv")):
        data = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
        if getattr(data, "dtype", None) is None or data.dtype.names is None:
            points = np.empty((0, 2), dtype=float)
        elif data.shape == ():
            if np.isnan(data["x"]) or np.isnan(data["y"]):
                points = np.empty((0, 2), dtype=float)
            else:
                points = np.array([[float(data["x"]), float(data["y"])]], dtype=float)
        elif data.size == 0:
            points = np.empty((0, 2), dtype=float)
        else:
            points = np.column_stack((data["x"], data["y"])).astype(float, copy=False)
        frames.append((path, points))
    return frames


def frame_label(path: Path) -> str:
    match = PARTICLE_PATTERN.search(path.name)
    if not match:
        return path.stem
    return f"t = {int(match.group(1))}"


def make_particle_gif() -> None:
    init_data = load_table(INIT_FILE, expected_cols=6)
    segments = build_segments(init_data)
    apertures = init_data[:, 4]
    frames = read_particle_frames()
    if not frames:
        raise ValueError(f"No particle position files found in {ROOT}")

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    collection = LineCollection(
        segments,
        cmap="viridis",
        norm=Normalize(vmin=float(apertures.min()), vmax=float(apertures.max())),
        linewidths=1.5,
    )
    collection.set_array(apertures)
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
        if points.size == 0:
            scatter.set_offsets(np.empty((0, 2), dtype=float))
        else:
            scatter.set_offsets(points)
        title.set_text(f"Particle positions on initial aperture background ({frame_label(frame_path)}, n={len(points)})")
        return scatter, title

    animation = FuncAnimation(fig, update, frames=len(frames), interval=600, blit=False, repeat=True)
    animation.save(PARTICLE_GIF, writer=PillowWriter(fps=2))
    plt.close(fig)

    counts = [len(points) for _, points in frames]
    summary_lines = [
        f"Frames = {len(frames)}",
        f"Min particles/frame = {min(counts)}",
        f"Max particles/frame = {max(counts)}",
        f"Mean particles/frame = {np.mean(counts):.12g}",
    ]
    PARTICLE_SUMMARY.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


def main() -> None:
    plot_velocity()
    make_particle_gif()


if __name__ == "__main__":
    main()
