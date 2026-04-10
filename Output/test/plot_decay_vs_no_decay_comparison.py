import csv
import math
import os
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np


BASE_DIR = os.path.dirname(__file__)
NO_DECAY_DIR = os.path.join(BASE_DIR, "no decay")
WITH_DECAY_DIR = os.path.join(BASE_DIR, "with decay")
CHEM_NO_DECAY = os.path.abspath(os.path.join(BASE_DIR, "..", "..", "Input", "Chemistry_files", "Chemistry_V3_match_no_decay.txt"))
CHEM_WITH_DECAY = os.path.abspath(os.path.join(BASE_DIR, "..", "..", "Input", "Chemistry_files", "Chemistry_V3_match_decay.txt"))


def read_cdf(path):
    times = []
    values = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                times.append(float(parts[0]))
                values.append(float(parts[1]))
    return times, values


def read_particle_snapshots(case_dir):
    snapshots = []
    for name in sorted(os.listdir(case_dir)):
        if not name.startswith("particle_positions_t") or not name.endswith(".csv"):
            continue
        path = os.path.join(case_dir, name)
        xs = []
        ys = []
        time_value = None
        with open(path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if time_value is None:
                    time_value = float(row["time"])
                xs.append(float(row["x"]))
                ys.append(float(row["y"]))
        snapshots.append(
            {
                "name": name,
                "time": time_value,
                "xs": xs,
                "ys": ys,
                "count": len(xs),
                "mean_x": sum(xs) / len(xs) if xs else None,
            }
        )
    return snapshots


def normalize_segment_key(x1, y1, x2, y2, ndigits=10):
    p1 = (round(x1, ndigits), round(y1, ndigits))
    p2 = (round(x2, ndigits), round(y2, ndigits))
    return tuple(sorted((p1, p2)))


def read_aperture_delta(path):
    rows = OrderedDict()
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 7:
                continue
            x1, y1, x2, y2, aperture, velocity, delta = map(float, parts[:7])
            key = normalize_segment_key(x1, y1, x2, y2)
            rows[key] = {
                "x1": x1,
                "y1": y1,
                "x2": x2,
                "y2": y2,
                "center_x": 0.5 * (x1 + x2),
                "center_y": 0.5 * (y1 + y2),
                "aperture": aperture,
                "velocity": velocity,
                "delta": delta,
                "length": math.hypot(x2 - x1, y2 - y1),
            }
    return rows


def first_crossing(times, values, threshold):
    for t, v in zip(times, values):
        if v >= threshold:
            return t
    return None


def read_chemistry_file(path):
    with open(path, "r", encoding="utf-8") as f:
        values = [float(line.strip()) for line in f if line.strip()]
    return {
        "C0": values[0],
        "k": values[1],
        "stoich": values[2],
        "molar_volume": values[3],
        "thickness": values[4],
    }


def draw_comparison_figure(
    no_decay_cdf_t,
    no_decay_cdf_v,
    with_decay_cdf_t,
    with_decay_cdf_v,
    no_decay_snapshots,
    with_decay_snapshots,
    no_decay_common,
    x_no,
    x_with,
    diff,
    chem_no_decay,
    chem_with_decay,
    nrows,
    ncols,
    figsize,
    out_path,
):
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, constrained_layout=True)
    flat_axes = np.asarray(axes, dtype=object).reshape(-1)

    ax = flat_axes[0]
    ax.plot(no_decay_cdf_t, no_decay_cdf_v, linewidth=2.0, label="No decay")
    ax.plot(with_decay_cdf_t, with_decay_cdf_v, linewidth=2.0, label="With decay")
    ax.set_xscale("log")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("CDF")
    ax.set_title("Breakthrough comparison")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = flat_axes[1]
    no_decay_times = [snap["time"] for snap in no_decay_snapshots if snap["mean_x"] is not None]
    no_decay_mean_x = [snap["mean_x"] for snap in no_decay_snapshots if snap["mean_x"] is not None]
    with_decay_times = [snap["time"] for snap in with_decay_snapshots if snap["mean_x"] is not None]
    with_decay_mean_x = [snap["mean_x"] for snap in with_decay_snapshots if snap["mean_x"] is not None]
    ax.plot(no_decay_times, no_decay_mean_x, marker="o", linewidth=1.8, label="No decay")
    ax.plot(with_decay_times, with_decay_mean_x, marker="s", linewidth=1.8, label="With decay")
    ax.set_xlabel("Snapshot time [s]")
    ax.set_ylabel("Mean particle x-position [m]")
    ax.set_title("Particle-cloud mean position")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = flat_axes[2]
    min_val = min(min(x_no), min(x_with))
    max_val = max(max(x_no), max(x_with))
    ax.scatter(x_no, x_with, s=10, alpha=0.35, color="#33658a")
    ax.plot([min_val, max_val], [min_val, max_val], linestyle="--", color="black", linewidth=1.0)
    ax.set_xlabel("Aperture change without decay [m]")
    ax.set_ylabel("Aperture change with decay [m]")
    ax.set_title("Per-segment geometry response")
    ax.grid(True, alpha=0.3)

    ax = flat_axes[3]
    centers_x = [row["center_x"] for row in no_decay_common]
    centers_y = [row["center_y"] for row in no_decay_common]
    sc = ax.scatter(centers_x, centers_y, c=diff, cmap="coolwarm", s=8)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("With decay minus no decay aperture change")
    ax.grid(True, alpha=0.2)
    fig.colorbar(sc, ax=ax, label="Delta difference [m]")

    ax = flat_axes[4]
    bins = 80
    ax.hist(x_no, bins=bins, alpha=0.55, label="No decay", color="#33658a", density=True)
    ax.hist(x_with, bins=bins, alpha=0.55, label="With decay", color="#f26419", density=True)
    ax.set_xlabel("Aperture change [m]")
    ax.set_ylabel("Density")
    ax.set_title("Aperture-change histogram comparison")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = flat_axes[5]
    max_time = max(no_decay_cdf_t[-1], with_decay_cdf_t[-1])
    time_grid = np.logspace(
        math.log10(min(no_decay_cdf_t[0], with_decay_cdf_t[0])),
        math.log10(max_time),
        400,
    )
    c_no = chem_no_decay["C0"] * np.exp(-chem_no_decay["k"] * time_grid)
    c_with = chem_with_decay["C0"] * np.exp(-chem_with_decay["k"] * time_grid)
    ax.plot(time_grid, c_no, linewidth=2.0, label="No decay chemistry")
    ax.plot(time_grid, c_with, linewidth=2.0, label="With decay chemistry")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Residence time [s]")
    ax.set_ylabel("Reactive concentration [mol/m^3]")
    ax.set_title("Model-imposed particle concentration decay")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def main():
    no_decay_cdf_t, no_decay_cdf_v = read_cdf(os.path.join(NO_DECAY_DIR, "cdf.txt"))
    with_decay_cdf_t, with_decay_cdf_v = read_cdf(os.path.join(WITH_DECAY_DIR, "cdf.txt"))

    no_decay_snapshots = read_particle_snapshots(NO_DECAY_DIR)
    with_decay_snapshots = read_particle_snapshots(WITH_DECAY_DIR)

    no_decay_delta = read_aperture_delta(os.path.join(NO_DECAY_DIR, "DFN_aperture_delta.txt"))
    with_decay_delta = read_aperture_delta(os.path.join(WITH_DECAY_DIR, "DFN_aperture_delta.txt"))
    chem_no_decay = read_chemistry_file(CHEM_NO_DECAY)
    chem_with_decay = read_chemistry_file(CHEM_WITH_DECAY)

    common_keys = sorted(set(no_decay_delta.keys()) & set(with_decay_delta.keys()))
    no_decay_common = [no_decay_delta[key] for key in common_keys]
    with_decay_common = [with_decay_delta[key] for key in common_keys]

    x_no = [row["delta"] for row in no_decay_common]
    x_with = [row["delta"] for row in with_decay_common]
    diff = [b - a for a, b in zip(x_no, x_with)]

    draw_comparison_figure(
        no_decay_cdf_t,
        no_decay_cdf_v,
        with_decay_cdf_t,
        with_decay_cdf_v,
        no_decay_snapshots,
        with_decay_snapshots,
        no_decay_common,
        x_no,
        x_with,
        diff,
        chem_no_decay,
        chem_with_decay,
        nrows=3,
        ncols=2,
        figsize=(14, 14),
        out_path=os.path.join(BASE_DIR, "decay_vs_no_decay_comparison.png"),
    )

    draw_comparison_figure(
        no_decay_cdf_t,
        no_decay_cdf_v,
        with_decay_cdf_t,
        with_decay_cdf_v,
        no_decay_snapshots,
        with_decay_snapshots,
        no_decay_common,
        x_no,
        x_with,
        diff,
        chem_no_decay,
        chem_with_decay,
        nrows=2,
        ncols=3,
        figsize=(21, 10),
        out_path=os.path.join(BASE_DIR, "decay_vs_no_decay_comparison_wide.png"),
    )

    summary_path = os.path.join(BASE_DIR, "decay_vs_no_decay_comparison_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("Decay vs no-decay comparison summary\n")
        f.write(f"Matched segments: {len(common_keys)}\n")
        f.write(f"No-decay aperture sum: {sum(x_no):.12g}\n")
        f.write(f"With-decay aperture sum: {sum(x_with):.12g}\n")
        if abs(sum(x_no)) > 0.0:
            f.write(f"Net-sum ratio (with/no): {sum(x_with)/sum(x_no):.12g}\n")
        f.write(f"Mean absolute difference per matched segment: {sum(abs(v) for v in diff)/len(diff):.12g}\n")
        f.write(f"Min difference (with-no): {min(diff):.12g}\n")
        f.write(f"Max difference (with-no): {max(diff):.12g}\n")
        f.write(f"No-decay first crossing CDF>=0.5: {first_crossing(no_decay_cdf_t, no_decay_cdf_v, 0.5)}\n")
        f.write(f"With-decay first crossing CDF>=0.5: {first_crossing(with_decay_cdf_t, with_decay_cdf_v, 0.5)}\n")
        f.write(f"No-decay first crossing CDF>=0.9: {first_crossing(no_decay_cdf_t, no_decay_cdf_v, 0.9)}\n")
        f.write(f"With-decay first crossing CDF>=0.9: {first_crossing(with_decay_cdf_t, with_decay_cdf_v, 0.9)}\n")
        f.write(f"No-decay chemistry C0: {chem_no_decay['C0']}\n")
        f.write(f"No-decay chemistry k: {chem_no_decay['k']}\n")
        f.write(f"With-decay chemistry C0: {chem_with_decay['C0']}\n")
        f.write(f"With-decay chemistry k: {chem_with_decay['k']}\n")
        f.write("No-decay mean x by snapshot:\n")
        for snap in no_decay_snapshots:
            if snap["mean_x"] is not None:
                f.write(f"  t={snap['time']}: mean_x={snap['mean_x']} count={snap['count']}\n")
        f.write("With-decay mean x by snapshot:\n")
        for snap in with_decay_snapshots:
            if snap["mean_x"] is not None:
                f.write(f"  t={snap['time']}: mean_x={snap['mean_x']} count={snap['count']}\n")
        f.write("Concentration-decay interpretation:\n")
        f.write("  The concentration curves are drawn from the active chemistry law C(t)=C0*exp(-k*t),\n")
        f.write("  using the two chemistry input files rather than particle-by-particle exported concentrations.\n")


if __name__ == "__main__":
    main()
