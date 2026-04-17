#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "Output" / "benchmark_particle_count_sensitivity"
CASE_INDEX = RESULTS_DIR / "case_index.csv"
FIG_DIR = RESULTS_DIR / "figures"
THICKNESS = 1.0

def load_aperture_delta(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def compute_case_metrics(row: dict[str, str]) -> dict[str, float | str]:
    case_dir = Path(row["case_output_dir"])
    raw_dir = case_dir / "raw_output"
    final = load_aperture_delta(raw_dir / "DFN_aperture_delta.txt")
    lengths = np.sqrt((final[:, 2] - final[:, 0]) ** 2 + (final[:, 3] - final[:, 1]) ** 2)
    delta_b = final[:, 6]
    dissolved_volume = float(np.sum(2.0 * delta_b * lengths * THICKNESS))
    return {
        "case": row["case"],
        "travel_time_s": float(row["travel_time_s"]),
        "nb_part": int(row["nb_part"]),
        "dissolved_volume_m3": dissolved_volume,
    }


def save_summary(rows: list[dict[str, float | str]]) -> None:
    with (RESULTS_DIR / "summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def make_plots(rows: list[dict[str, float | str]]) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    grouped: dict[float, list[dict[str, float | str]]] = defaultdict(list)
    for row in rows:
        grouped[float(row["travel_time_s"])].append(row)

    colors = {1e2: "#9b2226", 1e3: "#ca6702", 1e4: "#005f73", 1e5: "#0a9396"}
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for tt in sorted(grouped):
        data = sorted(grouped[tt], key=lambda item: int(item["nb_part"]))
        ax.plot(
            [int(item["nb_part"]) for item in data],
            [float(item["dissolved_volume_m3"]) for item in data],
            marker="o",
            lw=2.0,
            color=colors[tt],
            label=f"t={tt:.0e} s",
        )
    ax.set_xscale("log")
    ax.set_xlabel("Injected Particle Number")
    ax.set_ylabel("Dissolved Volume [m^3]")
    ax.set_title("Particle Count Sensitivity at n = 100")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_particle_count_sensitivity_volume.png", dpi=220)
    plt.close(fig)

    particle_counts = sorted({int(row["nb_part"]) for row in rows})
    travel_times = sorted(grouped)
    z = np.zeros((len(travel_times), len(particle_counts)))
    for i, tt in enumerate(travel_times):
        data = {int(item["nb_part"]): float(item["dissolved_volume_m3"]) for item in grouped[tt]}
        ref = data[100000]
        for j, npart in enumerate(particle_counts):
            z[i, j] = 0.0 if ref == 0.0 else 100.0 * (data[npart] - ref) / ref
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    im = ax.imshow(z, cmap="coolwarm", aspect="auto")
    ax.set_xticks(range(len(particle_counts)), [str(v) for v in particle_counts])
    ax.set_yticks(range(len(travel_times)), [f"{tt:.0e}" for tt in travel_times])
    ax.set_xlabel("Injected Particle Number")
    ax.set_ylabel("Crossing Time [s]")
    ax.set_title("Relative Difference vs 100000 Particles [%]")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("(V_N - V_100000) / V_100000 [%]")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_particle_count_sensitivity_heatmap.png", dpi=220)
    plt.close(fig)


def write_note(rows: list[dict[str, float | str]]) -> None:
    grouped: dict[float, list[dict[str, float | str]]] = defaultdict(list)
    for row in rows:
        grouped[float(row["travel_time_s"])].append(row)
    lines = [
        "# Particle Count Sensitivity",
        "",
        "- Geometry and chemistry are fixed to the current single-fracture benchmark with `n=100` and `gypsum_fast2` mode.",
        "- Injection time and total simulation time are both `1e5 s`.",
        "- Only injected particle number changes: `100`, `1000`, `10000`, `100000`.",
        "",
        "## Findings",
        "",
    ]
    for tt in sorted(grouped):
        data = sorted(grouped[tt], key=lambda item: int(item["nb_part"]))
        ref = float(data[-1]["dissolved_volume_m3"])
        low = float(data[0]["dissolved_volume_m3"])
        rel = 0.0 if ref == 0.0 else 100.0 * (low - ref) / ref
        lines.append(
            f"- `t={tt:.0e} s`: `N=100` vs `N=100000` changes dissolved volume by {rel:+.3f}%."
        )
    (RESULTS_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    with CASE_INDEX.open("r", newline="", encoding="utf-8") as f:
        case_rows = list(csv.DictReader(f))
    rows = [compute_case_metrics(row) for row in case_rows]
    rows.sort(key=lambda item: (float(item["travel_time_s"]), int(item["nb_part"])))
    save_summary(rows)
    make_plots(rows)
    write_note(rows)
    print(f"Wrote analysis products to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
