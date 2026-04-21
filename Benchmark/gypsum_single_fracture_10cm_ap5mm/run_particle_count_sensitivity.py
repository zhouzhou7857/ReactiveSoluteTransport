#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import shutil
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = REPO_ROOT / "Input"
OUTPUT_DIR = REPO_ROOT / "Output"
BENCH_DIR = REPO_ROOT / "Benchmark" / "gypsum_single_fracture_10cm"
FILE_NAMES_DIR = BENCH_DIR / "particle_count_file_names"
RESULTS_DIR = OUTPUT_DIR / "benchmark_particle_count_sensitivity"
CASE_DIR = RESULTS_DIR / "cases"
EXECUTABLE = REPO_ROOT / "Code" / "Release" / "ReactiveTransportPart"
EXEC_CWD = REPO_ROOT / "Code" / "Release"

CHEM_MODE = "gypsum_fast2"
DFN_NAME = "benchmark_gsf10_dfn_n100.txt"
TRAVEL_TIMES = [1.0e2, 1.0e3, 1.0e4, 1.0e5]
PARTICLE_COUNTS = [100, 1000, 10000, 100000]
VREF = 1.0e-3
THICKNESS = 1.0


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


def write_simulation_file(travel_time: float, nb_part: int) -> str:
    label = f"pcsens_t{travel_time:.0e}_np{nb_part}".replace("+", "")
    name = f"{label}.txt"
    path = INPUT_DIR / "Simulation_files" / name
    output_interval = min(1.0e4, max(travel_time / 10.0, 1.0))
    reaction_dt = travel_time / 100.0
    path.write_text(
        "\n".join(
            [
                str(nb_part),
                "0.05",
                "0",
                "1.000000e-03",
                "1.000000e+05",
                "200",
                "1",
                "1.000000e+05",
                f"{output_interval:.6e}",
                f"{reaction_dt:.6e}",
                f"{VREF:.6e}",
                f"{THICKNESS:.6e}",
            ]
        )
        + "\n",
        encoding="ascii",
    )
    return name


def write_file_names(case_label: str, domain_name: str, sim_name: str) -> Path:
    path = FILE_NAMES_DIR / f"{case_label}.txt"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join([domain_name, sim_name, DFN_NAME]) + "\n",
        encoding="ascii",
    )
    return path


def run_case(case_label: str, file_names_path: Path, case_output_dir: Path) -> None:
    clean_main_output()
    if case_output_dir.exists():
        shutil.rmtree(case_output_dir)
    case_output_dir.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        [str(EXECUTABLE), str(file_names_path)],
        cwd=EXEC_CWD,
        text=True,
        capture_output=True,
        check=False,
        env={**os.environ, "RST_CHEM_MODE": CHEM_MODE},
    )
    (case_output_dir / "run.log").write_text(completed.stdout + "\n\n[stderr]\n" + completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"{case_label} failed with exit code {completed.returncode}")
    copy_main_output(case_output_dir / "raw_output")


def main() -> None:
    rows = []
    CASE_DIR.mkdir(parents=True, exist_ok=True)
    for travel_time in TRAVEL_TIMES:
        domain_name = f"benchmark_gsf10_domain_t{travel_time:.0e}.txt".replace("+", "")
        for nb_part in PARTICLE_COUNTS:
            sim_name = write_simulation_file(travel_time, nb_part)
            case_label = f"t{travel_time:.0e}_np{nb_part}".replace("+", "")
            file_names_path = write_file_names(case_label, domain_name, sim_name)
            case_output_dir = CASE_DIR / case_label
            print(f"Running {case_label}")
            run_case(case_label, file_names_path, case_output_dir)
            rows.append(
                {
                    "case": case_label,
                    "travel_time_s": travel_time,
                    "nb_part": nb_part,
                    "domain_file": domain_name,
                    "simulation_file": sim_name,
                    "dfn_file": DFN_NAME,
                    "case_output_dir": str(case_output_dir),
                }
            )
    summary_path = RESULTS_DIR / "case_index.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
