import csv
import math
import os
import re
import statistics

import matplotlib.pyplot as plt
import numpy as np


BASE_DIR = os.path.dirname(__file__)
CASE_DIRS = [
    "Case 1-C0.01",
    "Case 2-C0.03",
    "Case 3-C0.06",
]
CASE_COLORS = {
    "Case 1-C0.01": "#33658a",
    "Case 2-C0.03": "#f6ae2d",
    "Case 3-C0.06": "#f26419",
}


def case_c0(case_name):
    match = re.search(r"C0\.([0-9]+)", case_name)
    if match:
        return float("0." + match.group(1))
    return None


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


def read_particle_snapshot_means(case_dir):
    times = []
    mean_x = []
    counts = []
    for name in sorted(os.listdir(case_dir)):
        if not name.startswith("particle_positions_t") or not name.endswith(".csv"):
            continue
        path = os.path.join(case_dir, name)
        xs = []
        time_value = None
        with open(path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                time_value = float(row["time"])
                xs.append(float(row["x"]))
        if xs:
            times.append(time_value)
            mean_x.append(sum(xs) / len(xs))
            counts.append(len(xs))
    return times, mean_x, counts


def normalize_segment_key(x1, y1, x2, y2, ndigits=10):
    p1 = (round(x1, ndigits), round(y1, ndigits))
    p2 = (round(x2, ndigits), round(y2, ndigits))
    return tuple(sorted((p1, p2)))


def read_aperture_delta(path):
    rows = {}
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
                "delta": delta,
                "velocity": velocity,
            }
    return rows


def crossing_time(times, values, threshold):
    for t, v in zip(times, values):
        if v >= threshold:
            return t
    return None


def load_case(case_name):
    case_dir = os.path.join(BASE_DIR, case_name)
    cdf_t, cdf_v = read_cdf(os.path.join(case_dir, "cdf.txt"))
    ap_rows = read_aperture_delta(os.path.join(case_dir, "DFN_aperture_delta.txt"))
    times, mean_x, counts = read_particle_snapshot_means(case_dir)
    deltas = [row["delta"] for row in ap_rows.values()]
    pos = [d for d in deltas if d > 0]
    neg = [d for d in deltas if d < 0]
    return {
        "name": case_name,
        "c0": case_c0(case_name),
        "cdf_t": cdf_t,
        "cdf_v": cdf_v,
        "ap_rows": ap_rows,
        "deltas": deltas,
        "pos": pos,
        "neg": neg,
        "snapshot_t": times,
        "snapshot_mean_x": mean_x,
        "snapshot_counts": counts,
    }


def main():
    cases = [load_case(name) for name in CASE_DIRS]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

    ax = axes[0][0]
    for case in cases:
        ax.plot(case["cdf_t"], case["cdf_v"], linewidth=2.0,
                color=CASE_COLORS[case["name"]],
                label=f"{case['name']}")
    ax.set_xscale("log")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("CDF")
    ax.set_title("Breakthrough comparison across initial concentrations")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    ax = axes[0][1]
    for case in cases:
        ax.plot(case["snapshot_t"], case["snapshot_mean_x"], marker="o", linewidth=1.8,
                color=CASE_COLORS[case["name"]],
                label=f"{case['name']}")
    ax.set_xlabel("Snapshot time [s]")
    ax.set_ylabel("Mean particle x-position [m]")
    ax.set_title("Mean particle-cloud position")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    ax = axes[1][0]
    bins = np.linspace(
        min(min(case["deltas"]) for case in cases),
        max(max(case["deltas"]) for case in cases),
        120,
    )
    for case in cases:
        ax.hist(case["deltas"], bins=bins, density=True, histtype="step", linewidth=2.0,
                color=CASE_COLORS[case["name"]], label=case["name"])
    ax.set_xlabel("Aperture change [m]")
    ax.set_ylabel("Density")
    ax.set_title("Aperture-change distribution")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    ax = axes[1][1]
    x = np.arange(len(cases))
    width = 0.24
    net = [sum(case["deltas"]) for case in cases]
    positive = [sum(case["pos"]) for case in cases]
    absolute = [sum(abs(v) for v in case["deltas"]) for case in cases]
    ax.bar(x - width, net, width=width, label="Net sum", color="#33658a")
    ax.bar(x, positive, width=width, label="Positive sum", color="#86bbd8")
    ax.bar(x + width, absolute, width=width, label="Absolute sum", color="#f26419")
    ax.set_xticks(x)
    ax.set_xticklabels([f"C0={case['c0']}" for case in cases])
    ax.set_ylabel("Integrated aperture change [m]")
    ax.set_title("Integrated geometry-response metrics")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(fontsize=9)

    out_path = os.path.join(BASE_DIR, "case123_reactivity_comparison.png")
    fig.savefig(out_path, dpi=220)

    fig2, axes2 = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    vmax = max(max(abs(row["delta"]) for row in case["ap_rows"].values()) for case in cases)
    for ax, case in zip(axes2, cases):
        xs = [row["center_x"] for row in case["ap_rows"].values()]
        ys = [row["center_y"] for row in case["ap_rows"].values()]
        cs = [row["delta"] for row in case["ap_rows"].values()]
        sc = ax.scatter(xs, ys, c=cs, cmap="coolwarm", s=8, vmin=-vmax, vmax=vmax)
        ax.set_title(f"{case['name']}")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.grid(True, alpha=0.2)
    fig2.colorbar(sc, ax=axes2, label="Aperture change [m]")
    out_path2 = os.path.join(BASE_DIR, "case123_reactivity_spatial_maps.png")
    fig2.savefig(out_path2, dpi=220)

    summary_path = os.path.join(BASE_DIR, "case123_reactivity_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("Case 1/2/3 reactivity comparison summary\n")
        for case in cases:
            f.write(f"\n{case['name']}\n")
            f.write(f"  Initial concentration proxy C0: {case['c0']}\n")
            f.write(f"  Segment count: {len(case['deltas'])}\n")
            f.write(f"  Net aperture-change sum: {sum(case['deltas']):.12g}\n")
            f.write(f"  Positive aperture-change sum: {sum(case['pos']):.12g}\n")
            f.write(f"  Negative aperture-change sum: {sum(case['neg']):.12g}\n")
            f.write(f"  Absolute aperture-change sum: {sum(abs(v) for v in case['deltas']):.12g}\n")
            f.write(f"  Mean aperture change: {statistics.mean(case['deltas']):.12g}\n")
            f.write(f"  Median aperture change: {statistics.median(case['deltas']):.12g}\n")
            f.write(f"  Max aperture change: {max(case['deltas']):.12g}\n")
            f.write(f"  Min aperture change: {min(case['deltas']):.12g}\n")
            f.write(f"  First nonzero arrival: {crossing_time(case['cdf_t'], case['cdf_v'], 1e-12)}\n")
            f.write(f"  t50: {crossing_time(case['cdf_t'], case['cdf_v'], 0.5)}\n")
            f.write(f"  t90: {crossing_time(case['cdf_t'], case['cdf_v'], 0.9)}\n")
            if case["snapshot_t"]:
                f.write("  Mean x by snapshots:\n")
                for t, mx, n in zip(case["snapshot_t"], case["snapshot_mean_x"], case["snapshot_counts"]):
                    f.write(f"    t={t:g}: mean_x={mx:.12g}, count={n}\n")

        f.write("\nCross-case interpretation\n")
        net = [sum(case["deltas"]) for case in cases]
        abs_sum = [sum(abs(v) for v in case["deltas"]) for case in cases]
        f.write(f"  Net-sum ratios: case2/case1={net[1]/net[0]:.12g}, case3/case1={net[2]/net[0]:.12g}\n")
        f.write(f"  Abs-sum ratios: case2/case1={abs_sum[1]/abs_sum[0]:.12g}, case3/case1={abs_sum[2]/abs_sum[0]:.12g}\n")
        f.write("  Note: the breakthrough curves remain effectively identical in the current outputs,\n")
        f.write("  so the clearest effect of increasing initial concentration is currently on geometry-response magnitude,\n")
        f.write("  not on the bulk CDF timing.\n")


if __name__ == "__main__":
    main()
