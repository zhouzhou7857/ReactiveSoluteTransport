import csv
import glob
import math
import os

import matplotlib.pyplot as plt


BASE_DIR = os.path.dirname(__file__)
APERTURE_FILE = os.path.join(BASE_DIR, "DFN_aperture_delta.txt")
PARTICLE_GLOB = os.path.join(BASE_DIR, "particle_positions_t*.csv")


def read_aperture_delta(path):
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("x1 y1 x2 y2"):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            x1, y1, x2, y2, aperture, velocity, delta = map(float, parts[:7])
            rows.append(
                {
                    "x1": x1,
                    "x2": x2,
                    "center_x": 0.5 * (x1 + x2),
                    "length": math.hypot(x2 - x1, y2 - y1),
                    "delta": delta,
                    "aperture": aperture,
                }
            )
    rows.sort(key=lambda r: r["center_x"])
    return rows


def read_particle_snapshots(pattern):
    snapshots = []
    for path in sorted(glob.glob(pattern)):
        with open(path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = []
            for row in reader:
                rows.append(
                    {
                        "time": float(row["time"]),
                        "particle_id": int(row["particle_id"]),
                        "x": float(row["x"]),
                    }
                )
        if rows:
            snapshots.append((rows[0]["time"], rows))
    return snapshots


def main():
    aperture_rows = read_aperture_delta(APERTURE_FILE)
    snapshots = read_particle_snapshots(PARTICLE_GLOB)
    if not aperture_rows:
        raise RuntimeError("No readable aperture-delta rows found")
    if not snapshots:
        raise RuntimeError("No particle snapshots found")

    segment_centers = [row["center_x"] for row in aperture_rows]
    segment_lengths = [row["length"] for row in aperture_rows]
    aperture_delta = [row["delta"] for row in aperture_rows]

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10, 7.5),
        sharex=True,
        gridspec_kw={"height_ratios": [1.1, 1.0]},
        constrained_layout=True,
    )

    ax = axes[0]
    ax.bar(
        segment_centers,
        aperture_delta,
        width=[0.9 * length for length in segment_lengths],
        color="#c85c43",
        edgecolor="black",
        linewidth=0.4,
    )
    ax.plot(segment_centers, aperture_delta, color="#7a2e1f", linewidth=1.0, alpha=0.9)
    ax.set_ylabel("Aperture change [m]")
    ax.set_title("Aperture change and particle positions along the single fracture")
    ax.grid(True, axis="y", alpha=0.3)

    ax = axes[1]
    colors = ["#33658a", "#86bbd8", "#758e4f", "#f6ae2d", "#f26419"]
    y_levels = list(range(len(snapshots)))
    for idx, (time_value, rows) in enumerate(snapshots):
        xs = [row["x"] for row in rows]
        ys = [y_levels[idx]] * len(rows)
        ax.scatter(xs, ys, s=28, color=colors[idx % len(colors)], alpha=0.85, label=f"t = {time_value:g} s")
    ax.set_yticks(y_levels)
    ax.set_yticklabels([f"t = {time_value:g} s" for time_value, _ in snapshots])
    ax.set_xlabel("x [m]")
    ax.set_ylabel("Particle snapshots")
    ax.grid(True, axis="x", alpha=0.3)
    ax.legend(loc="upper right")

    x_min = min(row["x1"] for row in aperture_rows)
    x_max = max(row["x2"] for row in aperture_rows)
    ax.set_xlim(x_min - 0.0002, x_max + 0.0002)

    output_path = os.path.join(BASE_DIR, "single_fracture_particle_aperture_combined.png")
    fig.savefig(output_path, dpi=220)

    summary_path = os.path.join(BASE_DIR, "single_fracture_particle_aperture_combined_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("Combined particle/aperture plotting summary\n")
        f.write(f"Number of segments: {len(aperture_rows)}\n")
        f.write(f"Snapshot count: {len(snapshots)}\n")
        f.write("Snapshot times: " + ", ".join(str(time_value) for time_value, _ in snapshots) + "\n")
        f.write(f"Max aperture change: {max(aperture_delta):.6e}\n")
        f.write(f"Min aperture change: {min(aperture_delta):.6e}\n")
        f.write(f"x-range: {x_min} to {x_max}\n")


if __name__ == "__main__":
    main()
