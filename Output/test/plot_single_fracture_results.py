import glob
import math
import os
import re

import matplotlib.pyplot as plt


BASE_DIR = os.path.dirname(__file__)
INIT_FILE = os.path.join(BASE_DIR, "DFN_init.txt")
DELTA_FILE = os.path.join(BASE_DIR, "DFN_aperture_delta.txt")
STEP_GLOB = os.path.join(BASE_DIR, "DFN_step*.txt")


def read_dfn_table(path):
    rows = []
    time_value = None
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("# time="):
                try:
                    time_value = float(line.split("=", 1)[1])
                except ValueError:
                    time_value = None
                continue
            if line.startswith("x1 y1 x2 y2"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            x1, y1, x2, y2, aperture, velocity = map(float, parts[:6])
            delta = float(parts[6]) if len(parts) >= 7 else None
            rows.append(
                {
                    "x1": x1,
                    "y1": y1,
                    "x2": x2,
                    "y2": y2,
                    "aperture": aperture,
                    "velocity": velocity,
                    "delta": delta,
                    "center_x": 0.5 * (x1 + x2),
                    "length": math.hypot(x2 - x1, y2 - y1),
                }
            )
    rows.sort(key=lambda r: r["center_x"])
    return time_value, rows


def same_geometry(rows_a, rows_b, tol=1e-12):
    if len(rows_a) != len(rows_b):
        return False
    for a, b in zip(rows_a, rows_b):
        for key in ("x1", "y1", "x2", "y2"):
            if abs(a[key] - b[key]) > tol:
                return False
    return True


def ensure_nonempty(rows, label):
    if not rows:
        raise RuntimeError(f"{label} has no readable DFN rows")


def main():
    _, init_rows = read_dfn_table(INIT_FILE)
    ensure_nonempty(init_rows, INIT_FILE)

    _, delta_rows = read_dfn_table(DELTA_FILE)
    ensure_nonempty(delta_rows, DELTA_FILE)

    valid_steps = []
    skipped_steps = []
    for step_path in sorted(glob.glob(STEP_GLOB)):
        time_value, rows = read_dfn_table(step_path)
        if same_geometry(init_rows, rows):
            valid_steps.append((step_path, time_value, rows))
        else:
            skipped_steps.append(os.path.basename(step_path))

    x = [row["center_x"] for row in init_rows]
    segment_index = list(range(1, len(init_rows) + 1))
    init_aperture = [row["aperture"] for row in init_rows]
    final_aperture = [row["aperture"] for row in delta_rows]
    delta_aperture = [row["delta"] if row["delta"] is not None else row["aperture"] - init_rows[i]["aperture"] for i, row in enumerate(delta_rows)]
    final_velocity = [row["velocity"] for row in delta_rows]

    plt.figure(figsize=(8, 4.8))
    plt.plot(segment_index, init_aperture, marker="o", label="Initial aperture")
    plt.plot(segment_index, final_aperture, marker="s", label="Final aperture")
    plt.xlabel("Segment index")
    plt.ylabel("Aperture [m]")
    plt.title("Single-fracture aperture profile")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(BASE_DIR, "single_fracture_aperture_profile.png"), dpi=200)
    plt.close()

    plt.figure(figsize=(8, 4.8))
    plt.bar(segment_index, delta_aperture, color="#c85c43")
    plt.xlabel("Segment index")
    plt.ylabel("Aperture change [m]")
    plt.title("Single-fracture aperture change by segment")
    plt.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(BASE_DIR, "single_fracture_aperture_delta.png"), dpi=200)
    plt.close()

    plt.figure(figsize=(8, 4.8))
    plt.plot(segment_index, final_velocity, marker="d", color="#33658a")
    plt.xlabel("Segment index")
    plt.ylabel("Velocity [m/s]")
    plt.title("Single-fracture final segment velocity")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(BASE_DIR, "single_fracture_velocity_profile.png"), dpi=200)
    plt.close()

    if valid_steps:
        plt.figure(figsize=(8, 4.8))
        for _, time_value, rows in valid_steps:
            label = f"t={time_value:g}" if time_value is not None else "t=unknown"
            plt.plot(segment_index, [row["aperture"] for row in rows], marker="o", linewidth=1.2, label=label)
        plt.xlabel("Segment index")
        plt.ylabel("Aperture [m]")
        plt.title("Aperture evolution on the valid single-fracture snapshots")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(BASE_DIR, "single_fracture_aperture_evolution.png"), dpi=200)
        plt.close()

    summary_path = os.path.join(BASE_DIR, "single_fracture_plot_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("Single-fracture plotting summary\n")
        f.write(f"Segments in current geometry: {len(init_rows)}\n")
        f.write(f"Valid step files used: {len(valid_steps)}\n")
        for step_path, time_value, _ in valid_steps:
            f.write(f"  used: {os.path.basename(step_path)} time={time_value}\n")
        f.write(f"Skipped step files: {len(skipped_steps)}\n")
        for name in skipped_steps:
            f.write(f"  skipped: {name}\n")
        f.write(f"Max aperture change: {max(delta_aperture):.6e}\n")
        f.write(f"Min aperture change: {min(delta_aperture):.6e}\n")
        f.write(f"Mean aperture change: {sum(delta_aperture)/len(delta_aperture):.6e}\n")


if __name__ == "__main__":
    main()
