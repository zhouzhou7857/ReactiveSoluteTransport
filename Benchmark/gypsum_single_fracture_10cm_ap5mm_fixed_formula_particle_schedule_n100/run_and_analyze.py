#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import os
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = REPO_ROOT / "Input"
OUTPUT_DIR = REPO_ROOT / "Output"
BENCH_DIR = REPO_ROOT / "Benchmark" / "gypsum_single_fracture_10cm_ap5mm_fixed_formula_particle_schedule_n100"
FILE_NAMES_DIR = BENCH_DIR / "file_names"
LOCAL_INPUT_DIR = BENCH_DIR / "generated_inputs"
RESULTS_DIR = OUTPUT_DIR / "benchmark_gypsum_single_fracture_10cm_ap5mm_fixed_formula_particle_schedule_n100"
CASE_DIR = RESULTS_DIR / "cases"
FIG_DIR = RESULTS_DIR / "figures"
EXECUTABLE = REPO_ROOT / "Code" / "Release" / "ReactiveTransportPart"
EXEC_CWD = REPO_ROOT / "Code" / "Release"
CHEM_MODE = "custom_delta_v"

LENGTH = 0.1
APERTURE = 5.0e-3
THICKNESS = 0.1
POROSITY = 0.05
DM = 1.0e-30
SIMU_OPTION = 0
T_MIN = 1.0e-3
NT = 200
SEED = 1
GLOBAL_INJECTION_TIME = 1.0e5
GLOBAL_TOTAL_TIME = 1.0e5
SEGMENTS = 100

RHO = 1.0e3
G = 9.8
MU = 1.0e-3

TRAVEL_TIMES = [1.0e2, 1.0e3, 1.0e4, 1.0e5]
PARTICLE_COUNTS = {
    1.0e2: 100000,
    1.0e3: 10000,
    1.0e4: 1000,
    1.0e5: 100,
}

# Use the Vp = 5e-6 m^3 gypsum-only law for every case.
FIXED_A2 = 1.25061420e-08
FIXED_K2 = 1.59497791e-04
FIXED_L = 0.0
FIXED_VREF = 5.0e-06


def ensure_dirs() -> None:
    for path in [
        INPUT_DIR / "Domain_files",
        INPUT_DIR / "Simulation_files",
        INPUT_DIR / "DFN_files",
        FILE_NAMES_DIR,
        LOCAL_INPUT_DIR / "Domain_files",
        LOCAL_INPUT_DIR / "Simulation_files",
        LOCAL_INPUT_DIR / "DFN_files",
        CASE_DIR,
        FIG_DIR,
    ]:
        path.mkdir(parents=True, exist_ok=True)


def mirror_input_file(source: Path, relative_subdir: str) -> None:
    shutil.copy2(source, LOCAL_INPUT_DIR / relative_subdir / source.name)


def head_drop_for_travel_time(travel_time: float) -> float:
    velocity = LENGTH / travel_time
    return velocity * LENGTH * 12.0 * MU / (RHO * G * APERTURE * APERTURE)


def write_domain_file(label: str, travel_time: float) -> str:
    name = f"benchmark_gsf10_ap5mm_fixedvp_domain_{label}.txt"
    path = INPUT_DIR / "Domain_files" / name
    head_drop = head_drop_for_travel_time(travel_time)
    path.write_text(
        f"{LENGTH} {LENGTH}\n"
        f"{DM:.6e} {POROSITY}\n"
        f"{head_drop:.12e} 0.0\n",
        encoding="ascii",
    )
    mirror_input_file(path, "Domain_files")
    return name


def write_simulation_file(label: str, travel_time: float, nb_part: int) -> str:
    name = f"benchmark_gsf10_ap5mm_fixedvp_sim_{label}.txt"
    path = INPUT_DIR / "Simulation_files" / name
    output_interval = min(1.0e4, max(travel_time / 10.0, 1.0))
    reaction_dt = travel_time / 100.0
    path.write_text(
        f"{nb_part}\n"
        f"{0.05}\n"
        f"{SIMU_OPTION}\n"
        f"{T_MIN:.6e}\n"
        f"{GLOBAL_TOTAL_TIME:.6e}\n"
        f"{NT}\n"
        f"{SEED}\n"
        f"{GLOBAL_INJECTION_TIME:.6e}\n"
        f"{output_interval:.6e}\n"
        f"{reaction_dt:.6e}\n"
        f"{FIXED_VREF:.6e}\n"
        f"{THICKNESS:.6e}\n",
        encoding="ascii",
    )
    mirror_input_file(path, "Simulation_files")
    return name


def write_dfn_file() -> str:
    name = "benchmark_gsf10_ap5mm_fixedvp_dfn_n100.txt"
    path = INPUT_DIR / "DFN_files" / name
    dx = LENGTH / SEGMENTS
    x0 = -0.5 * LENGTH
    lines = ["file", str(SEGMENTS)]
    for idx in range(SEGMENTS):
        x1 = x0 + idx * dx
        x2 = x1 + dx
        lines.append(f"{x1:.12e} 0.0 {x2:.12e} 0.0 {APERTURE:.12e}")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    mirror_input_file(path, "DFN_files")
    return name


def write_file_names(case_label: str, domain_name: str, sim_name: str, dfn_name: str) -> Path:
    path = FILE_NAMES_DIR / f"{case_label}.txt"
    path.write_text("\n".join([domain_name, sim_name, dfn_name]) + "\n", encoding="ascii")
    return path


def clean_main_output() -> None:
    for pattern in [
        "DFN_step*.txt",
        "particle_positions_t*.csv",
        "DFN.txt",
        "DFN_aperture_delta.txt",
        "DFN_init.txt",
        "DFN_raw.txt",
        "cdf.txt",
        "pdf.txt",
    ]:
        for path in OUTPUT_DIR.glob(pattern):
            path.unlink()


def copy_main_output(target_dir: Path) -> None:
    if target_dir.exists():
        shutil.rmtree(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    for pattern in [
        "DFN_step*.txt",
        "particle_positions_t*.csv",
        "DFN.txt",
        "DFN_aperture_delta.txt",
        "DFN_init.txt",
        "DFN_raw.txt",
        "cdf.txt",
        "pdf.txt",
    ]:
        for path in OUTPUT_DIR.glob(pattern):
            shutil.copy2(path, target_dir / path.name)


def run_case(case_label: str, file_names_path: Path, case_output_dir: Path) -> None:
    clean_main_output()
    if case_output_dir.exists():
        shutil.rmtree(case_output_dir)
    case_output_dir.mkdir(parents=True, exist_ok=True)
    run_env = {
        **os.environ,
        "RST_CHEM_MODE": CHEM_MODE,
        "RST_DELTA_V_A1": "0.0",
        "RST_DELTA_V_K1": "0.0",
        "RST_DELTA_V_A2": f"{FIXED_A2:.12e}",
        "RST_DELTA_V_K2": f"{FIXED_K2:.12e}",
        "RST_DELTA_V_L": f"{FIXED_L:.12e}",
        "RST_DELTA_V_VREF": f"{FIXED_VREF:.12e}",
        "RST_FRACTURE_THICKNESS": f"{THICKNESS:.12e}",
        "RST_USE_EFFECTIVE_DIFFUSION_HEIGHT_FACTOR": "0",
        "RST_USE_VP_WIDTH_CORRECTION": "0",
    }
    completed = subprocess.run(
        [str(EXECUTABLE), str(file_names_path)],
        cwd=EXEC_CWD,
        text=True,
        capture_output=True,
        check=False,
        env=run_env,
    )
    (case_output_dir / "run.log").write_text(
        completed.stdout + "\n\n[stderr]\n" + completed.stderr,
        encoding="utf-8",
    )
    if completed.returncode != 0:
        raise RuntimeError(f"Case {case_label} failed with exit code {completed.returncode}")
    copy_main_output(case_output_dir / "raw_output")


def compute_metrics(case_dir: Path, travel_time: float, nb_part: int) -> dict[str, float]:
    path = case_dir / "raw_output" / "DFN_aperture_delta.txt"
    dissolved_volume = 0.0
    max_delta_b = 0.0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.split()
            x1, x2 = float(parts[0]), float(parts[2])
            delta_b = float(parts[-1])
            seg_len = abs(x2 - x1)
            dissolved_volume += 2.0 * THICKNESS * seg_len * delta_b
            max_delta_b = max(max_delta_b, delta_b)
    velocity = LENGTH / travel_time
    dt_particle = GLOBAL_INJECTION_TIME / (nb_part - 1) if nb_part > 1 else GLOBAL_INJECTION_TIME
    actual_vp = velocity * APERTURE * THICKNESS * dt_particle
    return {
        "velocity_m_per_s": velocity,
        "actual_vp_m3": actual_vp,
        "dissolved_volume_m3": dissolved_volume,
        "max_delta_aperture_m": max_delta_b,
        "max_delta_aperture_mm": max_delta_b * 1.0e3,
    }


def load_baseline_summary() -> dict[float, dict[str, float]]:
    baseline = {}
    path = OUTPUT_DIR / "benchmark_gypsum_single_fracture_10cm_ap5mm" / "summary.csv"
    with path.open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["segments"]) != 100:
                continue
            baseline[float(row["travel_time_s"])] = {
                "dissolved_volume_m3": float(row["dissolved_volume_m3"]),
            }
    for tt in TRAVEL_TIMES:
        case_name = f"t{tt:.0e}".replace("+", "") + "_n100"
        p = OUTPUT_DIR / "benchmark_gypsum_single_fracture_10cm_ap5mm" / "cases" / case_name / "raw_output" / "DFN_aperture_delta.txt"
        mx = 0.0
        with p.open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    mx = max(mx, float(line.split()[-1]))
        baseline[tt]["max_delta_aperture_mm"] = mx * 1.0e3
    return baseline


def plot_results(rows: list[dict[str, float]]) -> None:
    travel_times = [row["travel_time_s"] for row in rows]
    dissolved = [row["dissolved_volume_m3"] for row in rows]
    aperture = [row["max_delta_aperture_mm"] for row in rows]
    baseline_dissolved = [row["baseline_dissolved_volume_m3"] for row in rows]
    baseline_aperture = [row["baseline_max_delta_aperture_mm"] for row in rows]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5))

    axes[0].plot(travel_times, dissolved, marker="o", lw=1.8, color="#b85c00", label="Fixed Vp=5e-6 law")
    axes[0].plot(travel_times, baseline_dissolved, marker="s", lw=1.5, color="#005f73", label="Vp-matched law")
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("t_cross [s]")
    axes[0].set_ylabel("Dissolved Volume [m$^3$]")
    axes[0].set_title("Dissolved Volume")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend(frameon=False)

    axes[1].plot(travel_times, aperture, marker="o", lw=1.8, color="#b85c00", label="Fixed Vp=5e-6 law")
    axes[1].plot(travel_times, baseline_aperture, marker="s", lw=1.5, color="#005f73", label="Vp-matched law")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("t_cross [s]")
    axes[1].set_ylabel("Max Aperture Increase [mm]")
    axes[1].set_title("Max Aperture Increase")
    axes[1].grid(True, which="both", alpha=0.25)
    axes[1].legend(frameon=False)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_fixed_formula_vs_vp_matched.png", dpi=220)
    plt.close(fig)


def main() -> None:
    ensure_dirs()
    dfn_name = write_dfn_file()
    baseline = load_baseline_summary()
    rows = []

    for travel_time in TRAVEL_TIMES:
        label = f"t{travel_time:.0e}".replace("+", "")
        nb_part = PARTICLE_COUNTS[travel_time]
        domain_name = write_domain_file(label, travel_time)
        sim_name = write_simulation_file(label, travel_time, nb_part)
        case_label = f"{label}_n100_np{nb_part}"
        file_names_path = write_file_names(case_label, domain_name, sim_name, dfn_name)
        case_output_dir = CASE_DIR / case_label
        print(f"Running {case_label}")
        run_case(case_label, file_names_path, case_output_dir)
        metrics = compute_metrics(case_output_dir, travel_time, nb_part)
        rows.append(
            {
                "case": case_label,
                "travel_time_s": travel_time,
                "particle_count": nb_part,
                "velocity_m_per_s": metrics["velocity_m_per_s"],
                "actual_vp_m3": metrics["actual_vp_m3"],
                "formula_vref_m3": FIXED_VREF,
                "formula_a2_m3": FIXED_A2,
                "formula_k2_s_inv": FIXED_K2,
                "dissolved_volume_m3": metrics["dissolved_volume_m3"],
                "max_delta_aperture_mm": metrics["max_delta_aperture_mm"],
                "baseline_dissolved_volume_m3": baseline[travel_time]["dissolved_volume_m3"],
                "baseline_max_delta_aperture_mm": baseline[travel_time]["max_delta_aperture_mm"],
                "dissolved_volume_ratio_to_baseline": metrics["dissolved_volume_m3"] / baseline[travel_time]["dissolved_volume_m3"],
                "max_delta_ratio_to_baseline": metrics["max_delta_aperture_mm"] / baseline[travel_time]["max_delta_aperture_mm"],
                "domain_file": domain_name,
                "simulation_file": sim_name,
                "dfn_file": dfn_name,
                "file_names_path": str(file_names_path),
                "case_output_dir": str(case_output_dir),
            }
        )

    results_path = RESULTS_DIR / "summary.csv"
    with results_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    plot_results(rows)
    print(f"Wrote {results_path}")


if __name__ == "__main__":
    main()
