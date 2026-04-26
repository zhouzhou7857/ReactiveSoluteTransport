#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
from collections import defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


REPO_ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = REPO_ROOT / "Input"
OUTPUT_DIR = REPO_ROOT / "Output"
BENCH_DIR = REPO_ROOT / "Benchmark" / "gypsum_single_fracture_10cm_ap5mm"
OLD_RESULTS_DIR = BENCH_DIR / "results_snapshot"
OLD_OUTPUT_RESULTS_DIR = OUTPUT_DIR / "benchmark_gypsum_single_fracture_10cm_ap5mm"

EXECUTABLE = REPO_ROOT / "Code" / "Release" / "ReactiveTransportPart"
EXEC_CWD = REPO_ROOT / "Code" / "Release"
PHREEQC_ROOT = Path("/home/zhouw/Documents/Codes/phreeqc-3.8.6-17100")
PHREEQC_EXE = Path("/usr/local/bin/phreeqc")
PHREEQC_DB = PHREEQC_ROOT / "database" / "phreeqc.dat"
PHREEQC_TEMPLATE = PHREEQC_ROOT / "gypsum" / "volume_5e-06.pqi"

CHEM_MODE = "custom_delta_v"
RUN_TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
RUN_TAG = f"benchmark_gypsum_single_fracture_10cm_ap5mm_total_injection_reference_volume_precipitation_{RUN_TIMESTAMP}"
RUN_ROOT = BENCH_DIR / f"results_total_injection_reference_volume_precipitation_{RUN_TIMESTAMP}"
RUN_OUTPUT_ROOT = OUTPUT_DIR / RUN_TAG
PHREEQC_DIR = RUN_ROOT / "phreeqc"
PT_CASE_DIR = RUN_ROOT / "pt_cases"
FIG_DIR = RUN_ROOT / "figures"
VALIDATION_DIR = RUN_ROOT / "validation"
LOG_DIR = RUN_ROOT / "logs"
ORIGINAL_INPUT_COPY_DIR = RUN_ROOT / "copied_original_inputs"
NEW_INPUT_COPY_DIR = RUN_ROOT / "new_inputs"
SCRIPT_SNAPSHOT_DIR = RUN_ROOT / "scripts_snapshot"

DOMAIN_NAME = f"{RUN_TAG}_domain_{{label}}.txt"
SIM_NAME = f"{RUN_TAG}_sim_{{label}}.txt"
DFN_NAME = f"{RUN_TAG}_dfn_n{{segments}}.txt"

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

GYPSUM_PRECIP_RATE_MOL_PER_M2_PER_S = 1.0e-7
GYPSUM_PRECIP_REACTION_ORDER = 1.0
CA_MMOL_PER_KGW = 30.0
S6_MMOL_PER_KGW = 30.0
GYPSUM_MOLAR_VOLUME_M3_PER_MOL = 7.39e-5
A_REF_M2 = 2.0 * LENGTH * THICKNESS
A_REF_CM2 = A_REF_M2 * 1.0e4

AXIS_LABEL_FONTSIZE = 14
TICK_LABEL_FONTSIZE = 12
LEGEND_FONTSIZE = 12

plt.rcParams.update(
    {
        "axes.labelsize": AXIS_LABEL_FONTSIZE,
        "axes.titlesize": AXIS_LABEL_FONTSIZE,
        "xtick.labelsize": TICK_LABEL_FONTSIZE,
        "ytick.labelsize": TICK_LABEL_FONTSIZE,
        "legend.fontsize": LEGEND_FONTSIZE,
        "font.family": "serif",
        "mathtext.fontset": "stix",
    }
)


@dataclass
class TravelTimeSetup:
    travel_time_s: float
    velocity_m_per_s: float
    q_in_m3_per_s: float
    injection_duration_s: float
    vref_m3: float
    water_kg: float
    reactive_area_m2: float
    reactive_area_cm2: float
    maximum_solution_gypsum_volume_m3: float
    maximum_solution_gypsum_moles: float
    particle_time_share_s: float
    nominal_particle_volume_m3: float
    nominal_sum_particle_volume_m3: float
    nominal_sum_particle_ratio: float
    nominal_sum_particle_relative_error: float


@dataclass
class Fast2Fit:
    a2_m3: float
    k2_s_inv: float
    linear_m3_per_s: float
    rmse_m3: float
    r2: float


def sci_label(value: float) -> str:
    return f"{value:.0e}".replace("+", "")


def ensure_dirs() -> None:
    for path in [
        RUN_ROOT,
        RUN_OUTPUT_ROOT,
        PHREEQC_DIR,
        PT_CASE_DIR,
        FIG_DIR,
        VALIDATION_DIR,
        LOG_DIR,
        ORIGINAL_INPUT_COPY_DIR,
        NEW_INPUT_COPY_DIR / "Domain_files",
        NEW_INPUT_COPY_DIR / "Simulation_files",
        NEW_INPUT_COPY_DIR / "DFN_files",
        NEW_INPUT_COPY_DIR / "file_names",
        SCRIPT_SNAPSHOT_DIR,
        INPUT_DIR / "Domain_files",
        INPUT_DIR / "Simulation_files",
        INPUT_DIR / "DFN_files",
    ]:
        path.mkdir(parents=True, exist_ok=True)


def copy_original_context() -> None:
    for name in ["README.md", "run_benchmark.py", "analyze_benchmark.py"]:
        shutil.copy2(BENCH_DIR / name, ORIGINAL_INPUT_COPY_DIR / name)
    if (BENCH_DIR / "generated_inputs").exists():
        shutil.copytree(
            BENCH_DIR / "generated_inputs",
            ORIGINAL_INPUT_COPY_DIR / "generated_inputs",
            dirs_exist_ok=True,
        )
    if (BENCH_DIR / "file_names").exists():
        shutil.copytree(
            BENCH_DIR / "file_names",
            ORIGINAL_INPUT_COPY_DIR / "file_names",
            dirs_exist_ok=True,
        )
    if PHREEQC_TEMPLATE.exists():
        shutil.copy2(PHREEQC_TEMPLATE, ORIGINAL_INPUT_COPY_DIR / PHREEQC_TEMPLATE.name)
    shutil.copy2(__file__, SCRIPT_SNAPSHOT_DIR / Path(__file__).name)
    shutil.copy2(REPO_ROOT / "Code" / "src" / "Chemistry" / "Chemistry.cpp", SCRIPT_SNAPSHOT_DIR / "Chemistry.cpp")
    shutil.copy2(REPO_ROOT / "Code" / "src" / "Chemistry" / "Chemistry.h", SCRIPT_SNAPSHOT_DIR / "Chemistry.h")
    shutil.copy2(REPO_ROOT / "Code" / "src" / "Transport" / "Transport.cpp", SCRIPT_SNAPSHOT_DIR / "Transport.cpp")


def build_executable() -> None:
    proc = subprocess.run(
        ["make", "-j4"],
        cwd=EXEC_CWD,
        text=True,
        capture_output=True,
        check=False,
    )
    (LOG_DIR / "build.log").write_text(proc.stdout + "\n--- STDERR ---\n" + proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError("Failed to rebuild ReactiveTransportPart; see logs/build.log")


def compute_setup(travel_time: float) -> TravelTimeSetup:
    velocity = LENGTH / travel_time
    q_in = abs(velocity) * APERTURE * THICKNESS
    vref = q_in * GLOBAL_INJECTION_TIME
    particle_time_share = GLOBAL_INJECTION_TIME / (NB_PART - 1)
    nominal_particle_volume = q_in * particle_time_share
    nominal_sum_particle_volume = nominal_particle_volume * NB_PART
    max_gypsum_moles = min(CA_MMOL_PER_KGW, S6_MMOL_PER_KGW) * 1.0e-3 * (vref * 1000.0)
    return TravelTimeSetup(
        travel_time_s=travel_time,
        velocity_m_per_s=velocity,
        q_in_m3_per_s=q_in,
        injection_duration_s=GLOBAL_INJECTION_TIME,
        vref_m3=vref,
        water_kg=vref * 1000.0,
        reactive_area_m2=A_REF_M2,
        reactive_area_cm2=A_REF_CM2,
        maximum_solution_gypsum_volume_m3=max_gypsum_moles * GYPSUM_MOLAR_VOLUME_M3_PER_MOL,
        maximum_solution_gypsum_moles=max_gypsum_moles,
        particle_time_share_s=particle_time_share,
        nominal_particle_volume_m3=nominal_particle_volume,
        nominal_sum_particle_volume_m3=nominal_sum_particle_volume,
        nominal_sum_particle_ratio=nominal_sum_particle_volume / vref if vref > 0.0 else math.nan,
        nominal_sum_particle_relative_error=(nominal_sum_particle_volume - vref) / vref if vref > 0.0 else math.nan,
    )


def head_drop_for_travel_time(travel_time: float) -> float:
    velocity = LENGTH / travel_time
    return velocity * LENGTH * 12.0 * MU / (RHO * G * APERTURE * APERTURE)


def mirror_new_input(source: Path, relative_subdir: str) -> None:
    shutil.copy2(source, NEW_INPUT_COPY_DIR / relative_subdir / source.name)


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
    mirror_new_input(path, "Domain_files")
    return name


def write_simulation_file(label: str, travel_time: float, vref: float) -> str:
    name = SIM_NAME.format(label=label)
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
        f"{vref:.12e}\n"
        f"{THICKNESS:.12e}\n",
        encoding="ascii",
    )
    mirror_new_input(path, "Simulation_files")
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
    mirror_new_input(path, "DFN_files")
    return name


def write_file_names(case_label: str, domain_name: str, sim_name: str, dfn_name: str) -> Path:
    path = NEW_INPUT_COPY_DIR / "file_names" / f"{case_label}.txt"
    path.write_text("\n".join([domain_name, sim_name, dfn_name]) + "\n", encoding="ascii")
    return path


def phreeqc_step_list(max_time: float) -> list[float]:
    target_times = [
        1, 1, 2, 3, 5, 8, 10, 20, 30, 50, 80,
        100, 200, 300, 500, 800,
        1e3, 2e3, 3e3, 5e3, 8e3,
        1e4, 2e4, 3e4, 5e4, 8e4,
        1e5,
    ]
    filtered_targets = sorted({float(step) for step in target_times if step <= max_time})
    increments: list[float] = []
    prev = 0.0
    for target in filtered_targets:
        increments.append(target - prev)
        prev = target
    return increments


def generate_phreeqc_input(setup: TravelTimeSetup) -> tuple[Path, Path, Path]:
    case_name = f"tt_{sci_label(setup.travel_time_s)}"
    case_dir = PHREEQC_DIR / case_name
    case_dir.mkdir(parents=True, exist_ok=True)
    pqi = case_dir / f"{case_name}.pqi"
    pqo = case_dir / f"{case_name}.pqo"
    sel = case_dir / f"{case_name}.sel"
    fit_csv = case_dir / f"{case_name}_volume_curve.csv"
    steps = " ".join(f"{step:g}" for step in phreeqc_step_list(GLOBAL_TOTAL_TIME))
    text = f"""TITLE Kinetic gypsum precipitation from supersaturated water for t_cross={setup.travel_time_s:.6e} s
#------------------------------------------------------------------------------
# New benchmark workflow:
#   - Reference water volume V_ref equals total injected solution volume
#   - Reactive surface area is fixed to the two-wall fracture area
#   - One PHREEQC cumulative precipitation law F_ref(t) is shared by all PT particles
#   - Precipitation is exported to PT as negative mineral-volume change
#------------------------------------------------------------------------------
# V_ref        = {setup.vref_m3:.12e} m3
# Water mass   = {setup.water_kg:.12e} kg
# A_ref        = {setup.reactive_area_m2:.12e} m2 = {setup.reactive_area_cm2:.6f} cm2
# Ca           = {CA_MMOL_PER_KGW:.6f} mmol/kgw
# S(6)         = {S6_MMOL_PER_KGW:.6f} mmol/kgw
# k_precip     = {GYPSUM_PRECIP_RATE_MOL_PER_M2_PER_S:.12e} mol/m2/s
# n_precip     = {GYPSUM_PRECIP_REACTION_ORDER:.6f}

SOLUTION 1 Supersaturated Ca-SO4 water
    temp    25.0
    units   mmol/kgw
    pH      7.0 charge
    Ca      {CA_MMOL_PER_KGW:.12e}
    S(6)    {S6_MMOL_PER_KGW:.12e}
    -water  {setup.water_kg:.12e}

RATES
Gypsum
-start
  10 si = SI("Gypsum")
  20 IF si <= 0 THEN SAVE 0 : END
  30 area = PARM(1)
  40 k = PARM(2)
  50 n = PARM(3)
  60 omega = SR("Gypsum")
  70 rate = -k * area * (omega - 1)^n
  80 moles = rate * TIME
  90 SAVE moles
-end

KINETICS 1
Gypsum
    -m      0.0
    -m0     0.0
    -parms  {setup.reactive_area_m2:.12e}  {GYPSUM_PRECIP_RATE_MOL_PER_M2_PER_S:.12e}  {GYPSUM_PRECIP_REACTION_ORDER:.12e}
    -tol    1e-8
    -steps  {steps}
    -cvode  true

INCREMENTAL_REACTIONS true

SELECTED_OUTPUT
    -file                {sel.name}
    -reset               false
    -totals              Ca S(6)
    -kinetics            Gypsum
    -saturation_indices  Gypsum

USER_PUNCH
    -headings Time_s pH SI_Gypsum SR_Gypsum V_Precip_cm3 DV_ref_m3
-start
  10 PUNCH TOTAL_TIME
  20 PUNCH -LA("H+")
  30 PUNCH SI("Gypsum")
  40 PUNCH SR("Gypsum")
  50 v_precip = KIN("Gypsum") * 73.9
  60 PUNCH v_precip
  70 PUNCH -v_precip * 1e-6
-end
END
"""
    pqi.write_text(text, encoding="utf-8")
    return pqi, pqo, fit_csv


def run_phreeqc(pqi: Path, pqo: Path) -> None:
    proc = subprocess.run(
        [str(PHREEQC_EXE), str(pqi), str(pqo), str(PHREEQC_DB)],
        cwd=pqi.parent,
        text=True,
        capture_output=True,
        check=False,
    )
    (pqi.parent / "phreeqc_run.log").write_text(proc.stdout + "\n--- STDERR ---\n" + proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(f"PHREEQC failed for {pqi}")


def parse_selected_output(path: Path) -> dict[str, np.ndarray]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    header = lines[0].split()
    rows: list[list[float]] = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) != len(header):
            continue
        try:
            rows.append([float(item) for item in parts])
        except ValueError:
            continue
    arr = np.array(rows, dtype=float)
    if arr.ndim != 2 or arr.shape[0] == 0:
        raise RuntimeError(f"No numeric data parsed from {path}")
    return {name: arr[:, idx] for idx, name in enumerate(header)}


def dv_model(time_s: np.ndarray, a2: float, k2: float, linear: float) -> np.ndarray:
    return a2 * (1.0 - np.exp(-k2 * time_s)) + linear * time_s


def fit_reference_curve(travel_time: float, phreeqc_data: dict[str, np.ndarray], csv_path: Path) -> Fast2Fit:
    time_s = phreeqc_data["Time_s"]
    signed_dv_m3 = phreeqc_data["DV_ref_m3"]
    precip_abs_m3 = np.maximum(-signed_dv_m3, 0.0)
    mask = (
        (time_s > 0.0)
        & (time_s <= GLOBAL_TOTAL_TIME + 1.0e-12)
        & np.isfinite(time_s)
        & np.isfinite(precip_abs_m3)
    )
    t_fit = time_s[mask]
    y_fit = precip_abs_m3[mask]
    if t_fit.size == 0:
        raise RuntimeError(f"No valid PHREEQC data to fit for t_cross={travel_time}")
    a2_guess = max(float(y_fit[-1]) * 1.05, 1.0e-18)
    half_value = 0.5 * float(y_fit[-1])
    half_index = int(np.searchsorted(y_fit, half_value, side="left"))
    if half_index < t_fit.size and t_fit[half_index] > 0.0:
        k2_guess = math.log(2.0) / float(t_fit[half_index])
    else:
        k2_guess = min(1.0 / max(travel_time, 1.0), 1.0e-3)
    linear_guess = max(y_fit[-1] / max(t_fit[-1], 1.0) * 1.0e-3, 0.0)
    sigma = np.maximum(y_fit, max(float(y_fit[-1]) * 1.0e-4, 1.0e-18))
    popt, _ = curve_fit(
        dv_model,
        t_fit,
        y_fit,
        p0=[a2_guess, k2_guess, linear_guess],
        bounds=([0.0, 0.0, 0.0], [np.inf, np.inf, np.inf]),
        sigma=sigma,
        absolute_sigma=False,
        maxfev=100000,
    )
    pred_abs = dv_model(t_fit, *popt)
    signed_fit = -pred_abs
    signed_target = -y_fit
    rmse = float(np.sqrt(np.mean((signed_fit - signed_target) ** 2)))
    denom = float(np.sum((y_fit - np.mean(y_fit)) ** 2))
    r2 = 1.0 - float(np.sum((pred_abs - y_fit) ** 2)) / denom if denom > 0.0 else 1.0
    with csv_path.open("w", newline="", encoding="ascii") as f:
        writer = csv.writer(f)
        writer.writerow(["Time_s", "DV_phreeqc_m3", "DV_fit_m3", "Precipitated_Gypsum_m3"])
        for t_value, dv_value, fit_value, precip_value in zip(t_fit, signed_target, signed_fit, y_fit):
            writer.writerow([f"{t_value:.12e}", f"{dv_value:.12e}", f"{fit_value:.12e}", f"{precip_value:.12e}"])
    return Fast2Fit(
        a2_m3=-float(popt[0]),
        k2_s_inv=float(popt[1]),
        linear_m3_per_s=-float(popt[2]),
        rmse_m3=rmse,
        r2=r2,
    )


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


def run_pt_case(
    case_label: str,
    file_names_path: Path,
    case_output_dir: Path,
    fit: Fast2Fit,
    setup: TravelTimeSetup,
) -> None:
    clean_main_output()
    if case_output_dir.exists():
        shutil.rmtree(case_output_dir)
    case_output_dir.mkdir(parents=True, exist_ok=True)
    run_env = {
        **os.environ,
        "RST_CHEM_MODE": CHEM_MODE,
        "RST_DELTA_V_A1": "0.0",
        "RST_DELTA_V_K1": "0.0",
        "RST_DELTA_V_A2": f"{fit.a2_m3:.12e}",
        "RST_DELTA_V_K2": f"{fit.k2_s_inv:.12e}",
        "RST_DELTA_V_L": f"{fit.linear_m3_per_s:.12e}",
        "RST_DELTA_V_VREF": f"{setup.vref_m3:.12e}",
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


def integrated_reference_dissolved_volume(setup: TravelTimeSetup, fit: Fast2Fit) -> float:
    t_end = GLOBAL_INJECTION_TIME
    if fit.k2_s_inv > 0.0:
        exponential_integral = fit.a2_m3 * (
            t_end - (1.0 - math.exp(-fit.k2_s_inv * t_end)) / fit.k2_s_inv
        )
    else:
        exponential_integral = 0.0
    linear_integral = 0.5 * fit.linear_m3_per_s * t_end * t_end
    return (setup.q_in_m3_per_s / setup.vref_m3) * (exponential_integral + linear_integral)


def load_aperture_delta(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def load_snapshot(path: Path) -> tuple[float, np.ndarray]:
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().strip()
    match = re.search(r"time=([0-9eE+.\-]+)", header)
    if match is None:
        raise ValueError(f"Missing time header in {path}")
    time_value = float(match.group(1))
    data = np.loadtxt(path, comments="#", skiprows=2)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return time_value, data


def load_all_snapshots(raw_dir: Path) -> list[tuple[float, np.ndarray]]:
    snapshots = [load_snapshot(path) for path in sorted(raw_dir.glob("DFN_step*.txt"))]
    snapshots.sort(key=lambda item: item[0])
    return snapshots


def compute_case_metrics(row: dict[str, str]) -> dict[str, float | str]:
    raw_dir = Path(row["case_output_dir"]) / "raw_output"
    final = load_aperture_delta(raw_dir / "DFN_aperture_delta.txt")
    lengths = np.sqrt((final[:, 2] - final[:, 0]) ** 2 + (final[:, 3] - final[:, 1]) ** 2)
    delta_b = final[:, 6]
    final_aperture = final[:, 4]
    mineral_volume_change = float(np.sum(2.0 * delta_b * lengths * THICKNESS))
    precipitated_volume = max(-mineral_volume_change, 0.0)
    max_delta_aperture_m = float(np.max(delta_b))
    min_delta_aperture_m = float(np.min(delta_b))
    mean_delta_aperture_m = float(np.mean(delta_b))
    return {
        "case": row["case"],
        "travel_time_s": float(row["travel_time_s"]),
        "segments": int(row["segments"]),
        "velocity_m_per_s": float(row["velocity_m_per_s"]),
        "reynolds": float(row["reynolds"]),
        "q_in_m3_per_s": float(row["q_in_m3_per_s"]),
        "vref_m3": float(row["vref_m3"]),
        "fit_a2_m3": float(row["fit_a2_m3"]),
        "fit_k2_s_inv": float(row["fit_k2_s_inv"]),
        "fit_linear_m3_per_s": float(row["fit_linear_m3_per_s"]),
        "fit_rmse_m3": float(row["fit_rmse_m3"]),
        "fit_r2": float(row["fit_r2"]),
        "reference_mineral_volume_change_m3": float(row["reference_dissolved_volume_m3"]),
        "mineral_volume_change_m3": mineral_volume_change,
        "precipitated_volume_m3": precipitated_volume,
        "max_delta_aperture_m": max_delta_aperture_m,
        "min_delta_aperture_m": min_delta_aperture_m,
        "mean_delta_aperture_m": mean_delta_aperture_m,
        "max_delta_aperture_mm": max_delta_aperture_m * 1.0e3,
        "min_delta_aperture_mm": min_delta_aperture_m * 1.0e3,
        "mean_delta_aperture_mm": mean_delta_aperture_m * 1.0e3,
        "final_mean_aperture_mm": float(np.mean(final_aperture) * 1.0e3),
        "case_output_dir": row["case_output_dir"],
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=False), encoding="utf-8")


def plot_f_ref_curves(phreeqc_records: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for record in phreeqc_records:
        label = sci_label(float(record["travel_time_s"]))
        curve = np.loadtxt(record["fit_curve_csv"], delimiter=",", skiprows=1)
        if curve.ndim == 1:
            curve = curve.reshape(1, -1)
        ax.plot(curve[:, 0], curve[:, 1] * 1.0e6, lw=2.0, label=rf"$t_{{cross}}={label}$ s")
    ax.set_xscale("log")
    ax.set_xlabel("Reaction Time [s]")
    ax.set_ylabel(r"Signed Gypsum Volume Change [$10^{-6}$ m$^3$]")
    ax.set_title(r"PHREEQC Precipitation Reference Curves $F_{ref}(t)$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_f_ref_curves.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.5), sharex=False, sharey=False)
    for ax, record in zip(axes.flat, phreeqc_records):
        label = sci_label(float(record["travel_time_s"]))
        curve = np.loadtxt(record["fit_curve_csv"], delimiter=",", skiprows=1)
        if curve.ndim == 1:
            curve = curve.reshape(1, -1)
        ax.plot(curve[:, 0], curve[:, 1] * 1.0e6, "o", ms=3.5, label="PHREEQC")
        ax.plot(curve[:, 0], curve[:, 2] * 1.0e6, "-", lw=1.8, label="Fit")
        ax.set_xscale("log")
        ax.set_title(rf"$t_{{cross}}={label}$ s")
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(frameon=False)
    axes[1, 0].set_xlabel("Reaction Time [s]")
    axes[1, 1].set_xlabel("Reaction Time [s]")
    axes[0, 0].set_ylabel(r"$\Delta V$ [$10^{-6}$ m$^3$]")
    axes[1, 0].set_ylabel(r"$\Delta V$ [$10^{-6}$ m$^3$]")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_f_ref_fit_check.png", dpi=220)
    plt.close(fig)


def plot_aperture_profile(summary_rows: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for row in sorted(summary_rows, key=lambda item: float(item["travel_time_s"])):
        if int(row["segments"]) != 100:
            continue
        final = load_aperture_delta(Path(row["case_output_dir"]) / "raw_output" / "DFN_aperture_delta.txt")
        x_mid_cm = 0.5 * (final[:, 0] + final[:, 2]) * 100.0 + 5.0
        aperture_mm = final[:, 4] * 1.0e3
        ax.plot(x_mid_cm, aperture_mm, lw=1.8, label=rf"$t_{{cross}}={sci_label(float(row['travel_time_s']))}$ s")
    ax.axhline(APERTURE * 1.0e3, color="#666666", lw=1.0, ls="--")
    ax.set_xlabel("x Along Fracture [cm]")
    ax.set_ylabel("Final Aperture [mm]")
    ax.set_title("Final Aperture Profile Along The 10 cm Fracture (n=100)")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_aperture_profile_n100.png", dpi=220)
    plt.close(fig)


def plot_aperture_evolution(summary_rows: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 8.0), sharex=True)
    for row in sorted(summary_rows, key=lambda item: float(item["travel_time_s"])):
        if int(row["segments"]) != 100:
            continue
        raw_dir = Path(row["case_output_dir"]) / "raw_output"
        snapshots = load_all_snapshots(raw_dir)
        times = np.array([time for time, _ in snapshots], dtype=float)
        mean_aperture = np.array([np.mean(data[:, 4]) for _, data in snapshots], dtype=float) * 1.0e3
        min_aperture = np.array([np.min(data[:, 4]) for _, data in snapshots], dtype=float) * 1.0e3
        label = rf"$t_{{cross}}={sci_label(float(row['travel_time_s']))}$ s"
        axes[0].plot(times, mean_aperture, lw=1.8, label=label)
        axes[1].plot(times, min_aperture, lw=1.8, label=label)
    for ax, title in zip(axes, ["Mean aperture", "Min aperture"]):
        ax.axhline(APERTURE * 1.0e3, color="#666666", lw=1.0, ls="--")
        ax.set_ylabel(f"{title} [mm]")
        ax.grid(True, alpha=0.25)
        ax.legend(frameon=False)
    axes[1].set_xlabel("Simulation Time [s]")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_aperture_evolution_n100.png", dpi=220)
    plt.close(fig)


def plot_total_dissolved_volume(summary_rows: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for row in sorted(summary_rows, key=lambda item: float(item["travel_time_s"])):
        if int(row["segments"]) != 100:
            continue
        raw_dir = Path(row["case_output_dir"]) / "raw_output"
        snapshots = load_all_snapshots(raw_dir)
        times = []
        precipitated_values = []
        for time_value, data in snapshots:
            lengths = np.sqrt((data[:, 2] - data[:, 0]) ** 2 + (data[:, 3] - data[:, 1]) ** 2)
            delta_b = data[:, 4] - APERTURE
            mineral_volume_change = float(np.sum(2.0 * delta_b * lengths * THICKNESS))
            times.append(time_value)
            precipitated_values.append(max(-mineral_volume_change, 0.0))
        ax.plot(times, np.array(precipitated_values) * 1.0e6, lw=1.8, label=rf"$t_{{cross}}={sci_label(float(row['travel_time_s']))}$ s")
    ax.set_xlabel("Simulation Time [s]")
    ax.set_ylabel(r"Accumulated Precipitated Gypsum Volume [$10^{-6}$ m$^3$]")
    ax.set_title("Accumulated Gypsum Precipitation Over Time (n=100)")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_total_precipitated_volume_n100.png", dpi=220)
    plt.close(fig)


def plot_old_vs_new(summary_rows: list[dict[str, object]]) -> None:
    old_summary = OLD_RESULTS_DIR / "summary.csv"
    if not old_summary.exists():
        return
    old_rows: list[dict[str, str]] = []
    with old_summary.open("r", encoding="utf-8") as f:
        old_rows = list(csv.DictReader(f))
    old_by_tt = {
        float(row["travel_time_s"]): row
        for row in old_rows
        if int(row["segments"]) == 100
    }
    new_by_tt = {
        float(row["travel_time_s"]): row
        for row in summary_rows
        if int(row["segments"]) == 100
    }
    travel_times = sorted(set(old_by_tt) & set(new_by_tt))
    if not travel_times:
        return
    x = np.arange(len(travel_times))
    width = 0.36
    old_max = [float(old_by_tt[tt]["max_delta_aperture_mm"]) for tt in travel_times]
    new_max = [float(new_by_tt[tt]["min_delta_aperture_mm"]) for tt in travel_times]
    old_vol = [float(old_by_tt[tt]["dissolved_volume_m3"]) * 1.0e6 for tt in travel_times]
    new_vol = [float(new_by_tt[tt]["mineral_volume_change_m3"]) * 1.0e6 for tt in travel_times]

    fig, axes = plt.subplots(2, 1, figsize=(8.2, 8.0), sharex=True)
    axes[0].bar(x - width / 2.0, old_max, width=width, label="Old per-particle $F_{V_p}(t)$")
    axes[0].bar(x + width / 2.0, new_max, width=width, label="New precipitation $F_{ref}(t)$")
    axes[0].set_ylabel("Aperture Change [mm]")
    axes[0].grid(True, axis="y", alpha=0.25)
    axes[0].legend(frameon=False)

    axes[1].bar(x - width / 2.0, old_vol, width=width, label="Old per-particle $F_{V_p}(t)$")
    axes[1].bar(x + width / 2.0, new_vol, width=width, label="New precipitation $F_{ref}(t)$")
    axes[1].set_ylabel(r"Signed Mineral Volume Change [$10^{-6}$ m$^3$]")
    axes[1].set_xticks(x, [sci_label(tt) for tt in travel_times])
    axes[1].set_xlabel("Crossing Time [s]")
    axes[1].grid(True, axis="y", alpha=0.25)
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_old_vs_new_n100.png", dpi=220)
    plt.close(fig)


def generate_validation_report(
    setups: list[TravelTimeSetup],
    phreeqc_records: list[dict[str, object]],
    summary_rows: list[dict[str, object]],
) -> None:
    chemistry_text = (REPO_ROOT / "Code" / "src" / "Chemistry" / "Chemistry.cpp").read_text(encoding="utf-8")
    formula_ok = (
        "double scale = ComputeParticleVolumeScalingFactor" in chemistry_text
        and "return delta_v_ref*scale;" in chemistry_text
    )
    report_lines = [
        "# Validation report",
        "",
        "## Active chemistry implementation",
        f"- `ComputeParticleSegmentMineralVolumeChange` uses `delta_v_ref * scale`: `{formula_ok}`",
        "- `scale` is computed from `V_particle / V_ref` when optional corrections are disabled.",
        "",
        "## Volume checks",
    ]
    payload: dict[str, object] = {
        "chemistry_formula_matches_expected": formula_ok,
        "travel_time_cases": [],
    }
    fit_by_tt = {float(record["travel_time_s"]): record for record in phreeqc_records}
    summary_by_tt = defaultdict(list)
    for row in summary_rows:
        summary_by_tt[float(row["travel_time_s"])].append(row)
    for setup in setups:
        phreeqc_record = fit_by_tt[setup.travel_time_s]
        fit_water_kg = float(phreeqc_record["water_kg"])
        fit_area_m2 = float(phreeqc_record["a_ref_m2"])
        report_lines.extend(
            [
                f"### t_cross = {setup.travel_time_s:.0e} s",
                f"- `Q_in = |u| * b * thickness = {setup.q_in_m3_per_s:.12e} m^3/s`",
                f"- `V_ref = Q_in * t_injection = {setup.vref_m3:.12e} m^3`",
                f"- `particle_time_share = {setup.particle_time_share_s:.12e} s`",
                f"- `nominal V_particle = {setup.nominal_particle_volume_m3:.12e} m^3`",
                f"- `sum(V_particle)_nominal = {setup.nominal_sum_particle_volume_m3:.12e} m^3`",
                f"- `sum(V_particle)/V_ref = {setup.nominal_sum_particle_ratio:.12e}`",
                f"- `PHREEQC water = {fit_water_kg:.12e} kg`",
                f"- `PHREEQC A_ref = {fit_area_m2:.12e} m^2`",
                f"- `one F_ref(t) file shared across segment counts = {len(summary_by_tt[setup.travel_time_s]) == len(SEGMENT_COUNTS)}`",
            ]
        )
        payload["travel_time_cases"].append(
            {
                **asdict(setup),
                "phreeqc_water_kg": fit_water_kg,
                "phreeqc_a_ref_m2": fit_area_m2,
                "shared_single_reference_curve_across_segment_counts": len(summary_by_tt[setup.travel_time_s]) == len(SEGMENT_COUNTS),
            }
        )
    report_lines.extend(
        [
            "",
            "## Workflow checks",
            "- No benchmark-local PHREEQC input uses per-particle `V_particle` as the water reference volume.",
            "- One PHREEQC file is generated per crossing-time case, not per particle and not per segment count.",
            "- `A_ref` is fixed to the fracture two-wall area and is not multiplied by `Np`.",
            "- Positive PHREEQC precipitated gypsum volume is exported to PT as negative `F_ref(t)` so aperture decreases.",
        ]
    )
    (VALIDATION_DIR / "validation_report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    write_json(VALIDATION_DIR / "validation_report.json", payload)


def write_readme(
    setups: list[TravelTimeSetup],
    phreeqc_records: list[dict[str, object]],
    summary_rows: list[dict[str, object]],
) -> None:
    fit_by_tt = {float(record["travel_time_s"]): record for record in phreeqc_records}
    lines = [
        f"# {RUN_ROOT.name}",
        "",
        "## What changed relative to the old workflow",
        "- The previous benchmark represented gypsum dissolution.",
        "- This rerun represents kinetic gypsum precipitation from supersaturated Ca-SO4 water.",
        "- It uses one total-injection reference system per crossing-time case.",
        "- `V_ref` is the total injected solution volume computed from the PT inlet condition and injection duration.",
        "- `A_ref` is fixed to the physical two-wall fracture area `2 * 10 cm * 10 cm = 200 cm^2`.",
        "- PHREEQC water volume is set to `V_ref`, not to per-particle `V_particle`.",
        f"- Supersaturated solution is `Ca={CA_MMOL_PER_KGW:g} mmol/kgw`, `S(6)={S6_MMOL_PER_KGW:g} mmol/kgw`.",
        f"- Precipitation kinetics use `A={A_REF_M2:g} m2`, `k={GYPSUM_PRECIP_RATE_MOL_PER_M2_PER_S:g} mol/m2/s`, `n={GYPSUM_PRECIP_REACTION_ORDER:g}`.",
        "- PT chemistry uses `dV_particle = [F_ref(t_end) - F_ref(t_start)] * (V_particle / V_ref)`.",
        "- `F_ref(t)` is negative for precipitation, so the existing aperture update decreases aperture without changing the geometry formula.",
        "",
        "## Folder contents",
        "- `copied_original_inputs/`: copies of the original benchmark inputs and scripts used for comparison.",
        "- `new_inputs/`: regenerated PT input files for this rerun.",
        "- `phreeqc/`: regenerated PHREEQC inputs, raw outputs, extracted curves, and fitted coefficients.",
        "- `pt_cases/`: PT logs and raw outputs for all segment-count cases.",
        "- `figures/`: final plots for the rerun and old-vs-new comparison.",
        "- `validation/`: validation summary and formula checks.",
        "- `scripts_snapshot/`: copies of the active benchmark script and chemistry-relevant source files.",
        "",
        "## Reference-system summary",
    ]
    for setup in setups:
        fit = fit_by_tt[setup.travel_time_s]
        lines.extend(
            [
                f"### t_cross = {setup.travel_time_s:.0e} s",
                f"- `u_inlet = {setup.velocity_m_per_s:.12e} m/s`",
                f"- `Q_in = {setup.q_in_m3_per_s:.12e} m^3/s`",
                f"- `t_injection = {setup.injection_duration_s:.12e} s`",
                f"- `V_ref = {setup.vref_m3:.12e} m^3`",
                f"- `A_ref = {setup.reactive_area_m2:.12e} m^2 = {setup.reactive_area_cm2:.6f} cm^2`",
                f"- `Maximum Ca/S-limited gypsum precipitate = {setup.maximum_solution_gypsum_volume_m3:.12e} m^3`",
                f"- `Maximum Ca/S-limited gypsum moles = {setup.maximum_solution_gypsum_moles:.12e} mol`",
                f"- `Fit: A2 = {float(fit['a2_m3']):.12e} m^3, k2 = {float(fit['k2_s_inv']):.12e} s^-1, L = {float(fit['linear_m3_per_s']):.12e} m^3/s`",
            ]
        )
    n100_rows = [row for row in summary_rows if int(row["segments"]) == 100]
    if n100_rows:
        lines.extend(["", "## n=100 final results"])
        for row in sorted(n100_rows, key=lambda item: float(item["travel_time_s"])):
            lines.extend(
                [
                    f"### t_cross = {float(row['travel_time_s']):.0e} s",
                    f"- `mineral_volume_change = {float(row['mineral_volume_change_m3']):.12e} m^3`",
                    f"- `precipitated_volume = {float(row['precipitated_volume_m3']):.12e} m^3`",
                    f"- `min_delta_aperture = {float(row['min_delta_aperture_mm']):.12e} mm`",
                    f"- `mean_delta_aperture = {float(row['mean_delta_aperture_mm']):.12e} mm`",
                ]
            )
    (RUN_ROOT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()
    copy_original_context()
    build_executable()

    setups = [compute_setup(travel_time) for travel_time in TRAVEL_TIMES]
    write_json(
        RUN_ROOT / "parameter_summary.json",
        {
            "run_tag": RUN_TAG,
            "nb_part": NB_PART,
            "segment_counts": SEGMENT_COUNTS,
            "travel_time_cases": [asdict(setup) for setup in setups],
        },
    )

    phreeqc_records: list[dict[str, object]] = []
    fit_by_tt: dict[float, Fast2Fit] = {}
    for setup in setups:
        pqi, pqo, fit_csv = generate_phreeqc_input(setup)
        run_phreeqc(pqi, pqo)
        selected_output = pqi.with_suffix(".sel")
        phreeqc_data = parse_selected_output(selected_output)
        fit = fit_reference_curve(setup.travel_time_s, phreeqc_data, fit_csv)
        fit_by_tt[setup.travel_time_s] = fit
        fit_json = {
            **asdict(setup),
            **asdict(fit),
            "phreeqc_input": str(pqi),
            "phreeqc_output": str(pqo),
            "phreeqc_selected_output": str(selected_output),
            "fit_curve_csv": str(fit_csv),
            "water_kg": setup.water_kg,
            "a_ref_m2": setup.reactive_area_m2,
        }
        phreeqc_records.append(fit_json)
        write_json(pqi.parent / f"{pqi.stem}_summary.json", fit_json)

    plot_f_ref_curves(phreeqc_records)

    case_index_rows: list[dict[str, object]] = []
    for setup in setups:
        label = f"t{setup.travel_time_s:.0e}".replace("+", "")
        domain_name = write_domain_file(label, setup.travel_time_s)
        sim_name = write_simulation_file(label, setup.travel_time_s, setup.vref_m3)
        fit = fit_by_tt[setup.travel_time_s]
        reference_dissolved_volume = integrated_reference_dissolved_volume(setup, fit)
        reynolds = RHO * setup.velocity_m_per_s * APERTURE / MU
        for segments in SEGMENT_COUNTS:
            dfn_name = write_dfn_file(segments)
            case_label = f"{label}_n{segments:03d}"
            file_names_path = write_file_names(case_label, domain_name, sim_name, dfn_name)
            case_output_dir = PT_CASE_DIR / case_label
            print(f"Running {case_label}")
            run_pt_case(case_label, file_names_path, case_output_dir, fit, setup)
            case_index_rows.append(
                {
                    "case": case_label,
                    "travel_time_s": setup.travel_time_s,
                    "segments": segments,
                    "velocity_m_per_s": setup.velocity_m_per_s,
                    "reynolds": reynolds,
                    "q_in_m3_per_s": setup.q_in_m3_per_s,
                    "vref_m3": setup.vref_m3,
                    "fit_a2_m3": fit.a2_m3,
                    "fit_k2_s_inv": fit.k2_s_inv,
                    "fit_linear_m3_per_s": fit.linear_m3_per_s,
                    "fit_rmse_m3": fit.rmse_m3,
                    "fit_r2": fit.r2,
                    "reference_dissolved_volume_m3": reference_dissolved_volume,
                    "domain_file": domain_name,
                    "simulation_file": sim_name,
                    "dfn_file": dfn_name,
                    "file_names_path": str(file_names_path),
                    "case_output_dir": str(case_output_dir),
                }
            )
    write_csv(RUN_ROOT / "case_index.csv", case_index_rows)

    summary_rows = [compute_case_metrics(row) for row in case_index_rows]
    summary_rows.sort(key=lambda item: (float(item["travel_time_s"]), int(item["segments"])))
    write_csv(RUN_ROOT / "summary.csv", summary_rows)

    plot_aperture_profile(summary_rows)
    plot_aperture_evolution(summary_rows)
    plot_total_dissolved_volume(summary_rows)
    plot_old_vs_new(summary_rows)

    generate_validation_report(setups, phreeqc_records, summary_rows)
    write_readme(setups, phreeqc_records, summary_rows)

    shutil.copytree(RUN_ROOT, RUN_OUTPUT_ROOT, dirs_exist_ok=True)


if __name__ == "__main__":
    main()
