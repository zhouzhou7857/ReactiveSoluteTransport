#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import os
import shutil
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = REPO_ROOT / "Input"
OUTPUT_DIR = REPO_ROOT / "Output"
BENCH_DIR = REPO_ROOT / "Benchmark" / "gypsum_single_fracture_10cm_ap5mm"
FILE_NAMES_DIR = BENCH_DIR / "file_names"
LOCAL_INPUT_DIR = BENCH_DIR / "generated_inputs"
LOCAL_RESULTS_DIR = BENCH_DIR / "results_snapshot"
RESULTS_DIR = OUTPUT_DIR / "benchmark_gypsum_single_fracture_10cm_ap5mm"
CASE_DIR = RESULTS_DIR / "cases"
EXECUTABLE = REPO_ROOT / "Code" / "Release" / "ReactiveTransportPart"
EXEC_CWD = REPO_ROOT / "Code" / "Release"
CHEM_MODE = "custom_delta_v"

DOMAIN_NAME = "benchmark_gsf10_ap5mm_domain_{label}.txt"
SIM_NAME = "benchmark_gsf10_ap5mm_sim_{label}.txt"
DFN_NAME = "benchmark_gsf10_ap5mm_dfn_n{segments}.txt"

LENGTH = 0.1
APERTURE = 5.0e-3
THICKNESS = 0.1
POROSITY = 0.05
DM = 1.0e-30

NB_PART = 10000
PROBA_TRANSFER = 0.05
SIMU_OPTION = 0
T_MIN = 1.0e-3
NT = 200
SEED = 1
GLOBAL_INJECTION_TIME = 1.0e5
GLOBAL_TOTAL_TIME = 1.0e5

TRAVEL_TIMES = [1.0e2, 1.0e3, 1.0e4, 1.0e5]
SEGMENT_COUNTS = [1, 5, 10, 50, 100]

RHO = 1.0e3
G = 9.8
MU = 1.0e-3

CHEMISTRY_BY_TRAVEL_TIME = {
    1.0e2: {"A2": 1.25061420e-08, "K2": 1.59497791e-04, "L": 0.0, "VREF": 5.0e-06},
    1.0e3: {"A2": 1.25061420e-09, "K2": 3.43627573e-04, "L": 0.0, "VREF": 5.0e-07},
    1.0e4: {"A2": 1.25061420e-10, "K2": 7.40323164e-04, "L": 0.0, "VREF": 5.0e-08},
    1.0e5: {"A2": 1.25061420e-11, "K2": 1.59497791e-03, "L": 0.0, "VREF": 5.0e-09},
}


def ensure_dirs() -> None:
    for path in [
        INPUT_DIR / "Domain_files",
        INPUT_DIR / "Simulation_files",
        INPUT_DIR / "DFN_files",
        FILE_NAMES_DIR,
        CASE_DIR,
        LOCAL_INPUT_DIR / "Domain_files",
        LOCAL_INPUT_DIR / "Simulation_files",
        LOCAL_INPUT_DIR / "DFN_files",
        LOCAL_RESULTS_DIR,
    ]:
        path.mkdir(parents=True, exist_ok=True)


def mirror_input_file(source: Path, relative_subdir: str) -> None:
    target = LOCAL_INPUT_DIR / relative_subdir / source.name
    shutil.copy2(source, target)


def head_drop_for_travel_time(travel_time: float) -> float:
    velocity = LENGTH / travel_time
    return velocity * LENGTH * 12.0 * MU / (RHO * G * APERTURE * APERTURE)


def write_domain_file(label: str, travel_time: float) -> str:
    name = DOMAIN_NAME.format(label=label)
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


def chemistry_coefficients(travel_time: float) -> dict[str, float]:
    return dict(CHEMISTRY_BY_TRAVEL_TIME[travel_time])


def write_simulation_file(label: str, travel_time: float) -> str:
    name = SIM_NAME.format(label=label)
    path = INPUT_DIR / "Simulation_files" / name
    output_interval = min(1.0e4, max(travel_time / 10.0, 1.0))
    reaction_dt = travel_time / 100.0
    vref = chemistry_coefficients(travel_time)["VREF"]
    path.write_text(
        f"{NB_PART}\n"
        f"{PROBA_TRANSFER}\n"
        f"{SIMU_OPTION}\n"
        f"{T_MIN:.6e}\n"
        f"{GLOBAL_TOTAL_TIME:.6e}\n"
        f"{NT}\n"
        f"{SEED}\n"
        f"{GLOBAL_INJECTION_TIME:.6e}\n"
        f"{output_interval:.6e}\n"
        f"{reaction_dt:.6e}\n"
        f"{vref:.6e}\n"
        f"{THICKNESS:.6e}\n",
        encoding="ascii",
    )
    mirror_input_file(path, "Simulation_files")
    return name


def write_dfn_file(segments: int) -> str:
    name = DFN_NAME.format(segments=segments)
    path = INPUT_DIR / "DFN_files" / name
    dx = LENGTH / segments
    x0 = -0.5 * LENGTH
    lines = ["file", str(segments)]
    for idx in range(segments):
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


def run_case(case_label: str, file_names_path: Path, case_output_dir: Path, chemistry: dict[str, float], travel_time: float) -> None:
    clean_main_output()
    if case_output_dir.exists():
        shutil.rmtree(case_output_dir)
    case_output_dir.mkdir(parents=True, exist_ok=True)
    run_env = {
        **os.environ,
        "RST_CHEM_MODE": CHEM_MODE,
        "RST_DELTA_V_A1": "0.0",
        "RST_DELTA_V_K1": "0.0",
        "RST_DELTA_V_A2": f"{chemistry['A2']:.12e}",
        "RST_DELTA_V_K2": f"{chemistry['K2']:.12e}",
        "RST_DELTA_V_L": f"{chemistry['L']:.12e}",
        "RST_DELTA_V_VREF": f"{chemistry['VREF']:.12e}",
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


def analytic_reference(travel_time: float) -> tuple[float, float, float, float]:
    chemistry = chemistry_coefficients(travel_time)
    velocity = LENGTH / travel_time
    reynolds = RHO * velocity * APERTURE / MU
    damkohler = chemistry["K2"] * travel_time
    flow_rate = velocity * APERTURE * THICKNESS
    dissolved_volume = (flow_rate / chemistry["VREF"]) * chemistry["A2"] * (
        GLOBAL_INJECTION_TIME
        - (1.0 - math.exp(-chemistry["K2"] * GLOBAL_INJECTION_TIME)) / chemistry["K2"]
    )
    return velocity, reynolds, damkohler, dissolved_volume


def main() -> None:
    ensure_dirs()
    summary_rows = []

    for travel_time in TRAVEL_TIMES:
        label = f"t{travel_time:.0e}".replace("+", "")
        domain_name = write_domain_file(label, travel_time)
        sim_name = write_simulation_file(label, travel_time)
        chemistry = chemistry_coefficients(travel_time)
        velocity, reynolds, damkohler, dissolved_ref = analytic_reference(travel_time)
        for segments in SEGMENT_COUNTS:
            dfn_name = write_dfn_file(segments)
            case_label = f"{label}_n{segments:03d}"
            file_names_path = write_file_names(case_label, domain_name, sim_name, dfn_name)
            case_output_dir = CASE_DIR / case_label
            print(f"Running {case_label}")
            run_case(case_label, file_names_path, case_output_dir, chemistry, travel_time)
            summary_rows.append(
                {
                    "case": case_label,
                    "travel_time_s": travel_time,
                    "segments": segments,
                    "velocity_m_per_s": velocity,
                    "reynolds": reynolds,
                    "damkohler": damkohler,
                    "reference_dissolved_volume_m3": dissolved_ref,
                    "delta_v_a2_m3": chemistry["A2"],
                    "delta_v_k2_s_inv": chemistry["K2"],
                    "delta_v_linear_m3_per_s": chemistry["L"],
                    "vref_m3": chemistry["VREF"],
                    "effective_diffusion_height_factor": 1.0,
                    "domain_file": domain_name,
                    "simulation_file": sim_name,
                    "dfn_file": dfn_name,
                    "file_names_path": str(file_names_path),
                    "case_output_dir": str(case_output_dir),
                }
            )

    summary_path = RESULTS_DIR / "case_index.csv"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)
    shutil.copy2(summary_path, LOCAL_RESULTS_DIR / "case_index.csv")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
