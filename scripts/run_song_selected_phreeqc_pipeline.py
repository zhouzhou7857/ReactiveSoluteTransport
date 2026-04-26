#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit


REPO_ROOT = Path("/home/zhouw/Documents/Codes/ReactiveSoluteTransport")
INPUT_DIR = REPO_ROOT / "Input"
OUTPUT_DIR = REPO_ROOT / "Output"
RELEASE_DIR = REPO_ROOT / "Code" / "Release"
EXECUTABLE = RELEASE_DIR / "ReactiveTransportPart"

PHREEQC_ROOT = Path("/home/zhouw/Documents/Codes/phreeqc-3.8.6-17100")
PHREEQC_EXE = Path("/usr/local/bin/phreeqc")
PHREEQC_DB = PHREEQC_ROOT / "database" / "phreeqc.dat"
PHREEQC_TEMPLATE = PHREEQC_ROOT / "gypsum" / "volume_5e-06.pqi"

DOMAIN_NAME = "song_selected_10cm_dp5000Pa_domain.txt"
SIM_NORMAL_NAME = os.environ.get(
    "RST_SIM_NORMAL_NAME",
    "song_selected_10cm_rtm_vp5e-6_np1e5_tinj5000.txt",
)
SIM_PRECHECK_NAME = os.environ.get(
    "RST_SIM_PRECHECK_NAME",
    "song_selected_10cm_rtm_vp5e-6_np1e5_tinj5000_precheck_autogen.txt",
)

RUN_TAG = os.environ.get(
    "RST_RUN_TAG",
    "song_selected_10cm_phreeqc_pipeline_dp5000Pa_tinj5000_no_vp_width_correction",
)
RUN_ROOT = OUTPUT_DIR / RUN_TAG
USE_VP_WIDTH_CORRECTION = os.environ.get("RST_USE_VP_WIDTH_CORRECTION", "0")

ALL_CASES = {
    "p9_a1p5": "song_selected_10cm/p9_a1p5/p9_a1p5_raw_filemode.txt",
    "p12_a2p0": "song_selected_10cm/p12_a2p0/p12_a2p0_raw_filemode.txt",
    "p16_a2p5": "song_selected_10cm/p16_a2p5/p16_a2p5_raw_filemode.txt",
}
case_filter = [name.strip() for name in os.environ.get("RST_CASES", "").split(",") if name.strip()]
if case_filter:
    CASES = {name: ALL_CASES[name] for name in case_filter}
else:
    CASES = ALL_CASES

REF_SOLID_VOLUME_M3 = 5.0e-6
REF_GYPSUM_MOL = 0.04805
REF_DOLOMITE_MOL = 0.02095
REF_SYLVITE_MOL = 0.00265
GYPSUM_VM_CM3_PER_MOL = 73.9
DOLOMITE_VM_CM3_PER_MOL = 64.5
SYLVITE_VM_CM3_PER_MOL = 37.5
CALCITE_VM_CM3_PER_MOL = 36.9

FIT_A1_REF = 2.02191810e-5
FIT_K1_REF = 2.14128542e-3
FIT_A2_REF = 2.50122841e-6
FIT_K2_REF = 2.72737385e-5
FIT_C_REF = 2.48279978e-16
FIT_VP_REF = 1.0e-3


@dataclass
class Diagnostics:
    nb_part: int
    t_injection: float
    particle_time_share: float
    fracture_volume: float
    particle_volume: float
    estimated_particle_residence_time: float


@dataclass
class FormulaFit:
    A1: float
    k1: float
    A2: float
    k2: float
    c: float
    rmse: float


def ensure_dirs() -> None:
    for path in [
        RUN_ROOT,
        RUN_ROOT / "file_names",
        RUN_ROOT / "generated_inputs",
        RUN_ROOT / "generated_inputs" / "Simulation_files",
        RUN_ROOT / "generated_inputs" / "Domain_files",
        RUN_ROOT / "generated_inputs" / "DFN_files",
        RUN_ROOT / "phreeqc",
        RUN_ROOT / "rtm_runs",
    ]:
        path.mkdir(parents=True, exist_ok=True)


def clean_main_output() -> None:
    patterns = [
        "DFN_step*.txt",
        "particle_positions_t*.csv",
        "DFN.txt",
        "DFN_aperture_delta.txt",
        "DFN_init.txt",
        "DFN_raw.txt",
        "cdf.txt",
        "pdf.txt",
        "initial_transport_diagnostics.txt",
    ]
    for pattern in patterns:
        for path in OUTPUT_DIR.glob(pattern):
            path.unlink()


def create_precheck_simulation_file() -> Path:
    normal = INPUT_DIR / "Simulation_files" / SIM_NORMAL_NAME
    precheck = INPUT_DIR / "Simulation_files" / SIM_PRECHECK_NAME
    lines = normal.read_text(encoding="ascii").rstrip().splitlines()
    if lines and lines[-1].strip() in {"1", "2"}:
        lines[-1] = "1"
    else:
        lines.append("1")
    precheck.write_text("\n".join(lines) + "\n", encoding="ascii")
    shutil.copy2(precheck, RUN_ROOT / "generated_inputs" / "Simulation_files" / precheck.name)
    return precheck


def write_file_names(case_name: str, simulation_name: str) -> Path:
    path = RUN_ROOT / "file_names" / f"{case_name}_{simulation_name}.txt"
    path.write_text(
        "\n".join([DOMAIN_NAME, simulation_name, CASES[case_name]]) + "\n",
        encoding="ascii",
    )
    return path


def run_rtm(file_names_path: Path, env_updates: dict[str, str], log_path: Path) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    env.update(env_updates)
    clean_main_output()
    proc = subprocess.run(
        [str(EXECUTABLE), str(file_names_path)],
        cwd=RELEASE_DIR,
        env=env,
        text=True,
        capture_output=True,
    )
    log_path.write_text(proc.stdout + "\n--- STDERR ---\n" + proc.stderr, encoding="utf-8")
    return proc


def parse_initial_diagnostics(path: Path) -> Diagnostics:
    values: dict[str, float] = {}
    for line in path.read_text(encoding="ascii").splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0] != "initial_transport_diagnostics":
            values[parts[0]] = float(parts[1])
    return Diagnostics(
        nb_part=int(values["nb_part"]),
        t_injection=values["t_injection"],
        particle_time_share=values["particle_time_share"],
        fracture_volume=values["fracture_volume"],
        particle_volume=values["particle_volume"],
        estimated_particle_residence_time=values["estimated_particle_residence_time"],
    )


def snapshot_precheck_inputs(case_name: str, sim_name: str) -> None:
    case_dir = RUN_ROOT / "precheck" / case_name / "input"
    case_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(INPUT_DIR / "Domain_files" / DOMAIN_NAME, case_dir / DOMAIN_NAME)
    shutil.copy2(INPUT_DIR / "Simulation_files" / sim_name, case_dir / sim_name)
    dfn_src = INPUT_DIR / "DFN_files" / CASES[case_name]
    shutil.copy2(dfn_src, case_dir / dfn_src.name)


def reference_guess(vp: float) -> FormulaFit:
    s = vp / FIT_VP_REF
    return FormulaFit(
        A1=FIT_A1_REF * s,
        k1=FIT_K1_REF * s ** (-1.0 / 3.0),
        A2=FIT_A2_REF * s,
        k2=FIT_K2_REF * s ** (-1.0 / 3.0),
        c=FIT_C_REF * s ** (2.0 / 3.0),
        rmse=0.0,
    )


def generate_phreeqc_input(case_name: str, vp: float) -> tuple[Path, Path, Path]:
    case_dir = RUN_ROOT / "phreeqc" / case_name
    case_dir.mkdir(parents=True, exist_ok=True)
    scale = vp / REF_SOLID_VOLUME_M3
    water_kg = vp * 1000.0
    area_m2 = 2.0 * vp ** (2.0 / 3.0)
    gypsum_mol = REF_GYPSUM_MOL * scale
    dolomite_mol = REF_DOLOMITE_MOL * scale
    sylvite_mol = REF_SYLVITE_MOL * scale
    kinetic_sel_name = f"{case_name}.sel"
    equilibrium_sel_name = f"{case_name}_equilibrium.sel"
    lines = PHREEQC_TEMPLATE.read_text(encoding="utf-8").splitlines()
    current_section = ""
    current_species = ""
    selected_output_index = 0
    out_lines: list[str] = []

    for line in lines:
        stripped = line.strip()
        if line.startswith("TITLE "):
            if "Part B" in line:
                out_lines.append(f"TITLE Part B -- Equilibrium reference, scaled solid volume {vp:.12e} m3, R:W=1:1")
            else:
                out_lines.append(f"TITLE Mineral volume evolution: scaled solid volume {vp:.12e} m3, R:W=1:1")
            continue
        if stripped.startswith("SOLUTION 1"):
            current_section = "solution1"
            current_species = ""
        elif stripped.startswith("KINETICS 1"):
            current_section = "kinetics1"
            current_species = ""
        elif stripped.startswith("EQUILIBRIUM_PHASES 1"):
            current_section = "equilibrium1"
            current_species = ""
        elif stripped.startswith("SELECTED_OUTPUT"):
            current_section = "selected_output"
            current_species = ""
            selected_output_index += 1
        elif stripped.startswith("SOLUTION 2"):
            current_section = "solution2"
            current_species = ""
        elif stripped.startswith("EQUILIBRIUM_PHASES 2"):
            current_section = "equilibrium2"
            current_species = ""
        elif stripped.startswith("USER_PUNCH") or stripped.startswith("USER_GRAPH") or stripped == "END":
            current_species = ""

        if current_section == "solution1" and stripped.startswith("-water"):
            out_lines.append(re.sub(r"(-water\s+)[0-9.eE+-]+", rf"\g<1>{water_kg:.12e}", line))
            continue

        if current_section == "solution2" and stripped.startswith("-water"):
            out_lines.append(re.sub(r"(-water\s+)[0-9.eE+-]+", rf"\g<1>{water_kg:.12e}", line))
            continue

        if current_section == "kinetics1":
            if stripped in {"Gypsum", "Dolomite", "Sylvite"}:
                current_species = stripped
            if stripped.startswith("-parms"):
                rate_map = {"Gypsum": "5e-5", "Dolomite": "5e-10", "Sylvite": "5e-2"}
                rate = rate_map[current_species]
                prefix = re.match(r"^(\s*-parms\s+)", line).group(1)
                comment = ""
                if "#" in line:
                    comment = "  " + line.split("#", 1)[1]
                    comment = f"# {comment.strip('# ').strip()}"
                out_lines.append(f"{prefix}{rate:<6} {area_m2:.12e}  {comment}".rstrip())
                continue
            if stripped.startswith("-m0"):
                mol = {"Gypsum": gypsum_mol, "Dolomite": dolomite_mol, "Sylvite": sylvite_mol}[current_species]
                out_lines.append(re.sub(r"(-m0\s+)[0-9.eE+-]+", rf"\g<1>{mol:.12e}", line))
                continue
            if stripped.startswith("-m"):
                mol = {"Gypsum": gypsum_mol, "Dolomite": dolomite_mol, "Sylvite": sylvite_mol}[current_species]
                out_lines.append(re.sub(r"(-m\s+)[0-9.eE+-]+", rf"\g<1>{mol:.12e}", line))
                continue

        if current_section == "equilibrium2":
            if stripped.startswith("Gypsum"):
                out_lines.append(re.sub(r"(Gypsum\s+0\.0\s+)[0-9.eE+-]+", rf"\g<1>{gypsum_mol:.12e}", line))
                continue
            if stripped.startswith("Dolomite"):
                out_lines.append(re.sub(r"(Dolomite\s+0\.0\s+)[0-9.eE+-]+", rf"\g<1>{dolomite_mol:.12e}", line))
                continue
            if stripped.startswith("Sylvite"):
                out_lines.append(re.sub(r"(Sylvite\s+0\.0\s+)[0-9.eE+-]+", rf"\g<1>{sylvite_mol:.12e}", line))
                continue

        if current_section == "selected_output" and stripped.startswith("-file"):
            sel_name = kinetic_sel_name if selected_output_index == 1 else equilibrium_sel_name
            out_lines.append(re.sub(r"(-file\s+)\S+", rf"\g<1>{sel_name}", line))
            continue

        out_lines.append(line)

    text = "\n".join(out_lines) + "\n"
    pqi = case_dir / f"{case_name}.pqi"
    pqo = case_dir / f"{case_name}.pqo"
    sel = case_dir / kinetic_sel_name
    pqi.write_text(text, encoding="utf-8")
    return pqi, pqo, sel


def run_phreeqc(pqi: Path, pqo: Path) -> None:
    proc = subprocess.run(
        [str(PHREEQC_EXE), str(pqi), str(pqo), str(PHREEQC_DB)],
        cwd=pqi.parent,
        text=True,
        capture_output=True,
    )
    (pqi.parent / "phreeqc_run.log").write_text(proc.stdout + "\n--- STDERR ---\n" + proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(f"PHREEQC failed for {pqi}")


def parse_selected_output(sel: Path) -> dict[str, np.ndarray]:
    lines = [line.strip() for line in sel.read_text(encoding="utf-8").splitlines() if line.strip()]
    header = lines[0].split()
    rows = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) != len(header):
            continue
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            continue
    arr = np.array(rows, dtype=float)
    data = {name: arr[:, idx] for idx, name in enumerate(header)}
    return data


def dv_model(t: np.ndarray, A1: float, k1: float, A2: float, k2: float, c: float) -> np.ndarray:
    return A1 * (1.0 - np.exp(-k1 * t)) + A2 * (1.0 - np.exp(-k2 * t)) + c * t


def fit_formula(case_name: str, vp: float, sel: Path) -> FormulaFit:
    data = parse_selected_output(sel)
    t = data["Time_s"]
    total_cm3 = data["V_total_cm3"]
    initial_total_cm3 = (
        REF_GYPSUM_MOL * (vp / REF_SOLID_VOLUME_M3) * GYPSUM_VM_CM3_PER_MOL
        + REF_DOLOMITE_MOL * (vp / REF_SOLID_VOLUME_M3) * DOLOMITE_VM_CM3_PER_MOL
        + REF_SYLVITE_MOL * (vp / REF_SOLID_VOLUME_M3) * SYLVITE_VM_CM3_PER_MOL
    )
    mask = (t > 0.0) & np.isfinite(total_cm3) & (total_cm3 > 0.0)
    t_fit = t[mask]
    dv = np.maximum((initial_total_cm3 - total_cm3[mask]) * 1.0e-6, 0.0)
    if t_fit.size == 0:
        raise RuntimeError(f"No valid PHREEQC data for {case_name}")
    guess = reference_guess(vp)
    p0 = [guess.A1, guess.k1, guess.A2, guess.k2, max(guess.c, 1.0e-24)]
    lower = [0.0, 0.0, 0.0, 0.0, 0.0]
    upper = [np.inf, np.inf, np.inf, np.inf, np.inf]
    popt, _ = curve_fit(dv_model, t_fit, dv, p0=p0, bounds=(lower, upper), maxfev=50000)
    dv_pred = dv_model(t_fit, *popt)
    rmse = float(np.sqrt(np.mean((dv_pred - dv) ** 2)))
    fit = FormulaFit(A1=float(popt[0]), k1=float(popt[1]), A2=float(popt[2]), k2=float(popt[3]), c=float(popt[4]), rmse=rmse)
    fit_csv = RUN_ROOT / "phreeqc" / case_name / f"{case_name}_fit_curve.csv"
    with fit_csv.open("w", newline="", encoding="ascii") as f:
        writer = csv.writer(f)
        writer.writerow(["Time_s", "DV_phreeqc_m3", "DV_fit_m3"])
        writer.writerows(zip(t_fit, dv, dv_pred))
    return fit


def copy_runtime_outputs(target_dir: Path) -> None:
    target_dir.mkdir(parents=True, exist_ok=True)
    patterns = [
        "DFN_step*.txt",
        "particle_positions_t*.csv",
        "DFN.txt",
        "DFN_aperture_delta.txt",
        "DFN_init.txt",
        "DFN_raw.txt",
        "cdf.txt",
        "pdf.txt",
        "initial_transport_diagnostics.txt",
    ]
    for pattern in patterns:
        for path in OUTPUT_DIR.glob(pattern):
            shutil.copy2(path, target_dir / path.name)


def write_formula_csv(rows: list[dict[str, object]]) -> Path:
    path = RUN_ROOT / "derived_volume_formula_coefficients.csv"
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "case",
                "particle_volume_m3",
                "surface_area_m2",
                "A1_m3",
                "k1_s^-1",
                "A2_m3",
                "k2_s^-1",
                "c_m3_s^-1",
                "rmse_m3",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_summary_csv(rows: list[dict[str, object]]) -> Path:
    path = RUN_ROOT / "case_summary.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "case",
                "domain_file",
                "simulation_file",
                "dfn_file",
                "pressure_pa",
                "t_injection_s",
                "nb_part",
                "particle_time_share_s",
                "fracture_volume_m3",
                "particle_volume_m3",
                "estimated_particle_residence_time_s",
                "surface_area_m2",
                "A1_m3",
                "k1_s^-1",
                "A2_m3",
                "k2_s^-1",
                "c_m3_s^-1",
                "fit_rmse_m3",
                "rtm_returncode",
                "rtm_log",
                "phreeqc_input",
                "phreeqc_output",
                "phreeqc_selected_output",
                "precheck_diagnostics",
                "rtm_output_dir",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


def convert_csv_to_xlsx(csv_path: Path) -> Path | None:
    proc = subprocess.run(
        [
            "soffice",
            "--headless",
            "--convert-to",
            "xlsx",
            "--outdir",
            str(csv_path.parent),
            str(csv_path),
        ],
        text=True,
        capture_output=True,
    )
    log_path = csv_path.parent / "xlsx_conversion.log"
    log_path.write_text(proc.stdout + "\n--- STDERR ---\n" + proc.stderr, encoding="utf-8")
    xlsx_path = csv_path.with_suffix(".xlsx")
    return xlsx_path if xlsx_path.exists() else None


def main() -> None:
    ensure_dirs()
    shutil.copy2(INPUT_DIR / "Domain_files" / DOMAIN_NAME, RUN_ROOT / "generated_inputs" / "Domain_files" / DOMAIN_NAME)
    shutil.copy2(INPUT_DIR / "Simulation_files" / SIM_NORMAL_NAME, RUN_ROOT / "generated_inputs" / "Simulation_files" / SIM_NORMAL_NAME)
    precheck_sim = create_precheck_simulation_file()
    formula_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for case_name in CASES:
        snapshot_precheck_inputs(case_name, precheck_sim.name)
        precheck_file_names = write_file_names(case_name, precheck_sim.name)
        precheck_log = RUN_ROOT / "precheck" / case_name / "precheck.log"
        precheck_proc = run_rtm(
            precheck_file_names,
            {"RST_USE_VP_WIDTH_CORRECTION": USE_VP_WIDTH_CORRECTION},
            precheck_log,
        )
        if precheck_proc.returncode != 0:
            raise RuntimeError(f"Precheck failed for {case_name}")
        diag_src = OUTPUT_DIR / "initial_transport_diagnostics.txt"
        precheck_diag_path = RUN_ROOT / "precheck" / case_name / "initial_transport_diagnostics.txt"
        precheck_diag_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(diag_src, precheck_diag_path)
        diagnostics = parse_initial_diagnostics(precheck_diag_path)

        pqi, pqo, sel = generate_phreeqc_input(case_name, diagnostics.particle_volume)
        run_phreeqc(pqi, pqo)
        fit = fit_formula(case_name, diagnostics.particle_volume, sel)
        surface_area = 2.0 * diagnostics.particle_volume ** (2.0 / 3.0)
        formula_rows.append(
            {
                "case": case_name,
                "particle_volume_m3": diagnostics.particle_volume,
                "surface_area_m2": surface_area,
                "A1_m3": fit.A1,
                "k1_s^-1": fit.k1,
                "A2_m3": fit.A2,
                "k2_s^-1": fit.k2,
                "c_m3_s^-1": fit.c,
                "rmse_m3": fit.rmse,
            }
        )

        normal_file_names = write_file_names(case_name, SIM_NORMAL_NAME)
        rtm_case_dir = RUN_ROOT / "rtm_runs" / case_name
        rtm_case_dir.mkdir(parents=True, exist_ok=True)
        rtm_log = rtm_case_dir / "run.log"
        env = {
            "RST_USE_VP_WIDTH_CORRECTION": USE_VP_WIDTH_CORRECTION,
            "RST_CHEM_MODE": "custom_delta_v",
            "RST_DELTA_V_A1": f"{fit.A1:.16e}",
            "RST_DELTA_V_K1": f"{fit.k1:.16e}",
            "RST_DELTA_V_A2": f"{fit.A2:.16e}",
            "RST_DELTA_V_K2": f"{fit.k2:.16e}",
            "RST_DELTA_V_L": f"{fit.c:.16e}",
        }
        rtm_proc = run_rtm(normal_file_names, env, rtm_log)
        copy_runtime_outputs(rtm_case_dir / "output")
        input_dir = rtm_case_dir / "input"
        input_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(INPUT_DIR / "Domain_files" / DOMAIN_NAME, input_dir / DOMAIN_NAME)
        shutil.copy2(INPUT_DIR / "Simulation_files" / SIM_NORMAL_NAME, input_dir / SIM_NORMAL_NAME)
        dfn_src = INPUT_DIR / "DFN_files" / CASES[case_name]
        shutil.copy2(dfn_src, input_dir / dfn_src.name)
        shutil.copy2(normal_file_names, input_dir / "File_names.txt")

        summary_rows.append(
            {
                "case": case_name,
                "domain_file": str(INPUT_DIR / "Domain_files" / DOMAIN_NAME),
                "simulation_file": str(INPUT_DIR / "Simulation_files" / SIM_NORMAL_NAME),
                "dfn_file": str(INPUT_DIR / "DFN_files" / CASES[case_name]),
                "pressure_pa": 5000.0,
                "t_injection_s": diagnostics.t_injection,
                "nb_part": diagnostics.nb_part,
                "particle_time_share_s": diagnostics.particle_time_share,
                "fracture_volume_m3": diagnostics.fracture_volume,
                "particle_volume_m3": diagnostics.particle_volume,
                "estimated_particle_residence_time_s": diagnostics.estimated_particle_residence_time,
                "surface_area_m2": surface_area,
                "A1_m3": fit.A1,
                "k1_s^-1": fit.k1,
                "A2_m3": fit.A2,
                "k2_s^-1": fit.k2,
                "c_m3_s^-1": fit.c,
                "fit_rmse_m3": fit.rmse,
                "rtm_returncode": rtm_proc.returncode,
                "rtm_log": str(rtm_log),
                "phreeqc_input": str(pqi),
                "phreeqc_output": str(pqo),
                "phreeqc_selected_output": str(sel),
                "precheck_diagnostics": str(precheck_diag_path),
                "rtm_output_dir": str(rtm_case_dir / "output"),
            }
        )

    formula_csv = write_formula_csv(formula_rows)
    summary_csv = write_summary_csv(summary_rows)
    summary_xlsx = convert_csv_to_xlsx(summary_csv)

    readme = RUN_ROOT / "README.txt"
    lines = [
        f"Run root: {RUN_ROOT}",
        f"Formula CSV: {formula_csv}",
        f"Summary CSV: {summary_csv}",
        f"Summary XLSX: {summary_xlsx if summary_xlsx else 'conversion failed'}",
        "",
        "Case directories:",
    ]
    for row in summary_rows:
        lines.append(f"{row['case']}:")
        lines.append(f"  precheck diagnostics: {row['precheck_diagnostics']}")
        lines.append(f"  PHREEQC input: {row['phreeqc_input']}")
        lines.append(f"  PHREEQC output: {row['phreeqc_output']}")
        lines.append(f"  PHREEQC selected output: {row['phreeqc_selected_output']}")
        lines.append(f"  RTM log: {row['rtm_log']}")
        lines.append(f"  RTM output dir: {row['rtm_output_dir']}")
    readme.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(RUN_ROOT)
    print(formula_csv)
    print(summary_csv)
    if summary_xlsx:
        print(summary_xlsx)


if __name__ == "__main__":
    main()
