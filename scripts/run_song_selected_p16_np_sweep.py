#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import shutil
import subprocess
from pathlib import Path


REPO_ROOT = Path("/home/zhouw/Documents/Codes/ReactiveSoluteTransport")
INPUT_DIR = REPO_ROOT / "Input" / "Simulation_files"
OUTPUT_DIR = REPO_ROOT / "Output"
PIPELINE_SCRIPT = REPO_ROOT / "scripts" / "run_song_selected_phreeqc_pipeline.py"
PLOT_SCRIPT = REPO_ROOT / "scripts" / "plot_song_selected_10cm_phreeqc_rtm_geometry_changes.py"

ROOT_TAG = "song_selected_10cm_phreeqc_pipeline_dp5000Pa_tinj5000_p16_a2p5_np_sweep_no_vp_scaling"
ROOT_DIR = OUTPUT_DIR / ROOT_TAG

BASE_SIM_NAME = "song_selected_10cm_rtm_vp5e-6_np1e5_tinj5000.txt"
CASE_NAME = "p16_a2p5"
PRESSURE_PA = 5000
TINJ_S = 5000

NP_FACTORS = [
    ("np1e4", 10_000, 0.1),
    ("np1e6", 1_000_000, 10.0),
    ("np1e7", 10_000_000, 100.0),
]


def build_sim_file(np_value: int, label: str) -> tuple[str, str]:
    src = INPUT_DIR / BASE_SIM_NAME
    lines = src.read_text(encoding="ascii").splitlines()
    if not lines:
        raise RuntimeError(f"Empty simulation file: {src}")
    lines[0] = str(np_value)
    sim_name = f"song_selected_10cm_rtm_vp5e-6_{label}_tinj5000.txt"
    precheck_name = f"song_selected_10cm_rtm_vp5e-6_{label}_tinj5000_precheck_autogen.txt"
    (INPUT_DIR / sim_name).write_text("\n".join(lines) + "\n", encoding="ascii")
    return sim_name, precheck_name


def read_one_row(csv_path: Path) -> dict[str, str]:
    with csv_path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 1:
        raise RuntimeError(f"Expected 1 row in {csv_path}, got {len(rows)}")
    return rows[0]


def copy_if_exists(src: Path, dst: Path) -> None:
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def main() -> None:
    ROOT_DIR.mkdir(parents=True, exist_ok=True)
    sweep_rows: list[dict[str, object]] = []

    for label, np_value, factor in NP_FACTORS:
        sim_name, precheck_name = build_sim_file(np_value, label)
        run_tag = f"{ROOT_TAG}/{label}"
        env = dict(os.environ)
        env.update(
            {
                "RST_RUN_TAG": run_tag,
                "RST_USE_VP_WIDTH_CORRECTION": "0",
                "RST_CASES": CASE_NAME,
                "RST_SIM_NORMAL_NAME": sim_name,
                "RST_SIM_PRECHECK_NAME": precheck_name,
            }
        )

        subprocess.run(
            ["python", str(PIPELINE_SCRIPT)],
            cwd=REPO_ROOT,
            env=env,
            check=True,
            text=True,
        )
        subprocess.run(
            ["python", str(PLOT_SCRIPT)],
            cwd=REPO_ROOT,
            env=env,
            check=True,
            text=True,
        )

        run_root = OUTPUT_DIR / run_tag
        case_row = read_one_row(run_root / "case_summary.csv")
        geom_row = read_one_row(run_root / "figures" / "song_selected_10cm_phreeqc_rtm_geometry_summary.csv")

        # Copy key images to the sweep root for quick access.
        copy_if_exists(
            run_root / "figures" / f"{CASE_NAME}_geometry_triptych.png",
            ROOT_DIR / "figures" / f"{label}_{CASE_NAME}_geometry_triptych.png",
        )
        copy_if_exists(
            run_root / "figures" / f"{CASE_NAME}_geometry_triptych.pdf",
            ROOT_DIR / "figures" / f"{label}_{CASE_NAME}_geometry_triptych.pdf",
        )

        sweep_rows.append(
            {
                "label": label,
                "np_factor": factor,
                "nb_part": int(case_row["nb_part"]),
                "simulation_file": case_row["simulation_file"],
                "run_root": str(run_root),
                "particle_time_share_s": float(case_row["particle_time_share_s"]),
                "particle_volume_m3": float(case_row["particle_volume_m3"]),
                "surface_area_m2": float(case_row["surface_area_m2"]),
                "A1_m3": float(case_row["A1_m3"]),
                "k1_s^-1": float(case_row["k1_s^-1"]),
                "A2_m3": float(case_row["A2_m3"]),
                "k2_s^-1": float(case_row["k2_s^-1"]),
                "c_m3_s^-1": float(case_row["c_m3_s^-1"]),
                "fit_rmse_m3": float(case_row["fit_rmse_m3"]),
                "initial_aperture_mean_um": float(geom_row["initial_aperture_mean_um"]),
                "final_aperture_mean_um": float(geom_row["final_aperture_mean_um"]),
                "delta_mean_um": float(geom_row["delta_mean_um"]),
                "delta_min_um": float(geom_row["delta_min_um"]),
                "delta_max_um": float(geom_row["delta_max_um"]),
            }
        )

    summary_csv = ROOT_DIR / "particle_number_sweep_summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "label",
                "np_factor",
                "nb_part",
                "simulation_file",
                "run_root",
                "particle_time_share_s",
                "particle_volume_m3",
                "surface_area_m2",
                "A1_m3",
                "k1_s^-1",
                "A2_m3",
                "k2_s^-1",
                "c_m3_s^-1",
                "fit_rmse_m3",
                "initial_aperture_mean_um",
                "final_aperture_mean_um",
                "delta_mean_um",
                "delta_min_um",
                "delta_max_um",
            ],
        )
        writer.writeheader()
        writer.writerows(sweep_rows)

    readme = ROOT_DIR / "README.txt"
    lines = [
        f"Root run directory: {ROOT_DIR}",
        f"Case: {CASE_NAME}",
        f"Pressure: {PRESSURE_PA} Pa",
        f"Injection time: {TINJ_S} s",
        "",
        "Subruns:",
    ]
    for row in sweep_rows:
        lines.append(
            f"  {row['label']}: nb_part={row['nb_part']}, "
            f"delta_mean_um={row['delta_mean_um']}, delta_max_um={row['delta_max_um']}, "
            f"run_root={row['run_root']}"
        )
    readme.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(ROOT_DIR)
    print(summary_csv)


if __name__ == "__main__":
    main()
