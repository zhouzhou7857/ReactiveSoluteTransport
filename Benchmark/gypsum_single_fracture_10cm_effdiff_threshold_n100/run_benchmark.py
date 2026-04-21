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
BENCH_DIR = REPO_ROOT / "Benchmark" / "gypsum_single_fracture_10cm_effdiff_threshold_n100"
FILE_NAMES_DIR = BENCH_DIR / "file_names"
LOCAL_INPUT_DIR = BENCH_DIR / "generated_inputs"
RESULTS_DIR = OUTPUT_DIR / "benchmark_gypsum_single_fracture_10cm_effdiff_threshold_n100"
CASE_DIR = RESULTS_DIR / "cases"
EXECUTABLE = REPO_ROOT / "Code" / "Release" / "ReactiveTransportPart"
EXEC_CWD = REPO_ROOT / "Code" / "Release"
CHEM_MODE = "custom_delta_v"

LENGTH = 0.1
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
APERTURES = [1.0e-3, 5.0e-4]
SEGMENTS = 100

RHO = 1.0e3
G = 9.8
MU = 1.0e-3
D_EFFECTIVE = 1.0e-9

A2 = 2.50122e-6
K2 = 2.72730e-5
VREF = 1.0e-3


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
    ]:
        path.mkdir(parents=True, exist_ok=True)


def mirror_input_file(source: Path, relative_subdir: str) -> None:
    shutil.copy2(source, LOCAL_INPUT_DIR / relative_subdir / source.name)


def aperture_label(aperture: float) -> str:
    mm = aperture * 1.0e3
    if abs(mm - round(mm)) < 1.0e-12:
        return f"ap{int(round(mm))}mm"
    return f"ap{str(mm).replace('.', 'p')}mm"


def head_drop_for_travel_time(travel_time: float, aperture: float) -> float:
    velocity = LENGTH / travel_time
    return velocity * LENGTH * 12.0 * MU / (RHO * G * aperture * aperture)


def write_domain_file(label: str, travel_time: float, aperture: float) -> str:
    name = f"benchmark_gsf10_{label}_domain.txt"
    path = INPUT_DIR / "Domain_files" / name
    head_drop = head_drop_for_travel_time(travel_time, aperture)
    path.write_text(
        f"{LENGTH} {LENGTH}\n"
        f"{DM:.6e} {POROSITY}\n"
        f"{head_drop:.12e} 0.0\n",
        encoding="ascii",
    )
    mirror_input_file(path, "Domain_files")
    return name


def write_simulation_file(label: str, travel_time: float) -> str:
    name = f"benchmark_gsf10_{label}_sim.txt"
    path = INPUT_DIR / "Simulation_files" / name
    output_interval = min(1.0e4, max(travel_time / 10.0, 1.0))
    reaction_dt = travel_time / 100.0
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
        f"{VREF:.6e}\n"
        f"{THICKNESS:.6e}\n",
        encoding="ascii",
    )
    mirror_input_file(path, "Simulation_files")
    return name


def write_dfn_file(label: str, aperture: float) -> str:
    name = f"benchmark_gsf10_{label}_dfn_n{SEGMENTS}.txt"
    path = INPUT_DIR / "DFN_files" / name
    dx = LENGTH / SEGMENTS
    x0 = -0.5 * LENGTH
    lines = ["file", str(SEGMENTS)]
    for idx in range(SEGMENTS):
        x1 = x0 + idx * dx
        x2 = x1 + dx
        lines.append(f"{x1:.12e} 0.0 {x2:.12e} 0.0 {aperture:.12e}")
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


def effective_diffusion_height_factor(travel_time: float, aperture: float) -> float:
    if aperture <= 0.0:
        return 0.0
    return math.sqrt(D_EFFECTIVE * travel_time) / aperture


def run_case(case_label: str, file_names_path: Path, case_output_dir: Path, travel_time: float) -> None:
    clean_main_output()
    if case_output_dir.exists():
        shutil.rmtree(case_output_dir)
    case_output_dir.mkdir(parents=True, exist_ok=True)
    run_env = {
        **os.environ,
        "RST_CHEM_MODE": CHEM_MODE,
        "RST_DELTA_V_A1": "0.0",
        "RST_DELTA_V_K1": "0.0",
        "RST_DELTA_V_A2": f"{A2:.12e}",
        "RST_DELTA_V_K2": f"{K2:.12e}",
        "RST_DELTA_V_L": "0.0",
        "RST_DELTA_V_VREF": f"{VREF:.12e}",
        "RST_FRACTURE_THICKNESS": f"{THICKNESS:.12e}",
        "RST_USE_EFFECTIVE_DIFFUSION_HEIGHT_FACTOR": "1",
        "RST_EFFECTIVE_DIFFUSION_COEFFICIENT": f"{D_EFFECTIVE:.12e}",
        "RST_EFFECTIVE_DIFFUSION_TIME": f"{travel_time:.12e}",
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


def dissolved_volume_from_aperture_delta(path: Path, thickness: float) -> tuple[float, float]:
    total = 0.0
    max_delta = 0.0
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.split()
            x1 = float(parts[0])
            x2 = float(parts[2])
            delta_b = float(parts[-1])
            seg_len = abs(x2 - x1)
            total += 2.0 * thickness * seg_len * delta_b
            max_delta = max(max_delta, delta_b)
    return total, max_delta


def main() -> None:
    ensure_dirs()
    rows = []
    for aperture in APERTURES:
        alabel = aperture_label(aperture)
        for travel_time in TRAVEL_TIMES:
            tlabel = f"t{travel_time:.0e}".replace("+", "")
            label = f"{alabel}_{tlabel}"
            domain_name = write_domain_file(label, travel_time, aperture)
            sim_name = write_simulation_file(label, travel_time)
            dfn_name = write_dfn_file(label, aperture)
            case_label = f"{label}_n100"
            file_names_path = write_file_names(case_label, domain_name, sim_name, dfn_name)
            case_output_dir = CASE_DIR / case_label
            print(f"Running {case_label}")
            run_case(case_label, file_names_path, case_output_dir, travel_time)
            dissolved_volume, max_delta_b = dissolved_volume_from_aperture_delta(
                case_output_dir / "raw_output" / "DFN_aperture_delta.txt",
                THICKNESS,
            )
            rows.append(
                {
                    "case": case_label,
                    "aperture_m": aperture,
                    "aperture_mm": aperture * 1.0e3,
                    "travel_time_s": travel_time,
                    "effective_diffusion_height_m": math.sqrt(D_EFFECTIVE * travel_time),
                    "effective_diffusion_factor": effective_diffusion_height_factor(travel_time, aperture),
                    "dissolved_volume_m3": dissolved_volume,
                    "max_delta_aperture_m": max_delta_b,
                    "max_delta_aperture_mm": max_delta_b * 1.0e3,
                }
            )

    results_path = RESULTS_DIR / "summary.csv"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with results_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {results_path}")


if __name__ == "__main__":
    main()
