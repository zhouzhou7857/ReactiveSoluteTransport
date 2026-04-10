import csv
import glob
import os
import re

import matplotlib.pyplot as plt


BASE_DIR = os.path.dirname(__file__)
SIM_FILE = os.path.abspath(os.path.join(BASE_DIR, "..", "..", "Input", "Simulation_files", "Simulation_single_test.txt"))
DFN_INIT_FILE = os.path.join(BASE_DIR, "DFN_init.txt")
PARTICLE_GLOB = os.path.join(BASE_DIR, "particle_positions_t*.csv")


def read_simulation_parameters(path):
    with open(path, "r", encoding="utf-8") as f:
        values = [line.strip() for line in f if line.strip()]
    nb_part = int(float(values[0]))
    t_injection = float(values[7]) if len(values) >= 8 else 0.0
    return nb_part, t_injection


def read_fracture_extent(path):
    xs = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("x1 y1 x2 y2"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                xs.append(float(parts[0]))
                xs.append(float(parts[2]))
    if not xs:
        raise RuntimeError("No readable DFN rows in DFN_init.txt")
    return min(xs), max(xs)


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
                        "y": float(row["y"]),
                        "mesh_index": int(row["mesh_index"]),
                    }
                )
        if rows:
            snapshots.append((path, rows[0]["time"], rows))
    return snapshots


def representative_particle_ids(nb_part):
    if nb_part <= 5:
        return list(range(nb_part))
    return sorted({0, nb_part // 4, nb_part // 2, (3 * nb_part) // 4, nb_part - 1})


def main():
    nb_part, t_injection = read_simulation_parameters(SIM_FILE)
    inlet_x, outlet_x = read_fracture_extent(DFN_INIT_FILE)
    snapshots = read_particle_snapshots(PARTICLE_GLOB)
    if not snapshots:
        raise RuntimeError("No particle position snapshots found")

    dt_inject = t_injection / float(nb_part - 1) if nb_part > 1 and t_injection > 0.0 else 0.0
    selected_ids = representative_particle_ids(nb_part)

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), constrained_layout=True)

    ax = axes[0]
    colors = ["#33658a", "#86bbd8", "#758e4f", "#f6ae2d", "#f26419"]
    for idx, (_, time_value, rows) in enumerate(snapshots):
        xs = [row["x"] for row in rows]
        ys = [0.0 for _ in rows]
        ax.scatter(xs, ys, s=28, alpha=0.85, color=colors[idx % len(colors)], label=f"t = {time_value:g} s")
    ax.plot([inlet_x, outlet_x], [0.0, 0.0], color="black", linewidth=1.0, alpha=0.5)
    ax.set_title("Particle cloud snapshots along the single fracture")
    ax.set_xlabel("x [m]")
    ax.set_yticks([])
    ax.grid(True, axis="x", alpha=0.3)
    ax.legend()

    ax = axes[1]
    snapshot_maps = []
    for _, _, rows in snapshots:
        snapshot_maps.append({row["particle_id"]: row for row in rows})

    for particle_id in selected_ids:
        times = []
        xs = []
        inject_time = particle_id * dt_inject
        times.append(inject_time)
        xs.append(inlet_x)
        for (_, time_value, _), snapshot_map in zip(snapshots, snapshot_maps):
            row = snapshot_map.get(particle_id)
            if row is None:
                continue
            if time_value >= inject_time:
                times.append(time_value)
                xs.append(row["x"])
        ax.plot(times, xs, marker="o", linewidth=1.5, label=f"particle {particle_id}")

    ax.axhline(outlet_x, color="gray", linestyle="--", linewidth=1.0, label="outlet position")
    ax.set_title("Representative particle motion in time")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("x position [m]")
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=2, fontsize=9)

    output_path = os.path.join(BASE_DIR, "single_fracture_particle_motion.png")
    fig.savefig(output_path, dpi=220)

    summary_path = os.path.join(BASE_DIR, "single_fracture_particle_motion_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("Particle motion plotting summary\n")
        f.write(f"Number of particles in simulation: {nb_part}\n")
        f.write(f"Injection duration: {t_injection}\n")
        f.write(f"Inter-particle injection interval: {dt_inject}\n")
        f.write("Selected representative particle ids: " + ", ".join(str(pid) for pid in selected_ids) + "\n")
        f.write("Snapshot times used: " + ", ".join(str(time_value) for _, time_value, _ in snapshots) + "\n")
        f.write(f"Inlet x: {inlet_x}\n")
        f.write(f"Outlet x: {outlet_x}\n")


if __name__ == "__main__":
    main()
