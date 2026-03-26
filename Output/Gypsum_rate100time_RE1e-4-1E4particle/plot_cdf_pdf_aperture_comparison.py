#!/usr/bin/env python3

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "plots"
CASE_PATTERN = re.compile(r"^n=(\d+)$")
TIME_PATTERN = re.compile(r"#\s*time\s*=\s*([0-9.eE+-]+)")

plt.style.use("seaborn-v0_8-whitegrid")


def discover_cases(root: Path) -> list[Path]:
    cases: list[tuple[int, Path]] = []
    for path in root.iterdir():
        if not path.is_dir():
            continue
        match = CASE_PATTERN.match(path.name)
        if match:
            cases.append((int(match.group(1)), path))
    return [path for _, path in sorted(cases, key=lambda item: item[0])]


def read_two_column_file(file_path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(file_path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected format in {file_path}")
    return data[:, 0], data[:, 1]


def read_aperture_history(case_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    times: list[float] = []
    aperture_means: list[float] = []
    aperture_mins: list[float] = []
    aperture_maxs: list[float] = []

    for step_file in sorted(case_dir.glob("DFN_step*.txt")):
        with step_file.open("r", encoding="utf-8") as handle:
            lines = [line.strip() for line in handle if line.strip()]

        if len(lines) < 3:
            continue

        time_match = TIME_PATTERN.match(lines[0])
        if not time_match:
            raise ValueError(f"Cannot parse time from {step_file}")

        apertures = np.array([float(line.split()[4]) for line in lines[2:]], dtype=float)
        times.append(float(time_match.group(1)))
        aperture_means.append(float(apertures.mean()))
        aperture_mins.append(float(apertures.min()))
        aperture_maxs.append(float(apertures.max()))

    if not times:
        raise ValueError(f"No aperture history found in {case_dir}")

    order = np.argsort(times)
    times_array = np.asarray(times)[order]
    mean_array = np.asarray(aperture_means)[order]
    min_array = np.asarray(aperture_mins)[order]
    max_array = np.asarray(aperture_maxs)[order]
    return times_array, mean_array, min_array, max_array


def read_aperture_profile(case_dir: Path, target_time: float = 1.0e4) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    closest_time: float | None = None
    closest_left: np.ndarray | None = None
    closest_right: np.ndarray | None = None
    closest_aperture: np.ndarray | None = None

    for step_file in sorted(case_dir.glob("DFN_step*.txt")):
        with step_file.open("r", encoding="utf-8") as handle:
            lines = [line.strip() for line in handle if line.strip()]

        if len(lines) < 3:
            continue

        time_match = TIME_PATTERN.match(lines[0])
        if not time_match:
            raise ValueError(f"Cannot parse time from {step_file}")

        time_value = float(time_match.group(1))
        segment_rows = [line.split() for line in lines[2:]]
        x_left = np.array([float(fields[0]) for fields in segment_rows], dtype=float)
        x_right = np.array([float(fields[2]) for fields in segment_rows], dtype=float)
        apertures = np.array([float(fields[4]) for fields in segment_rows], dtype=float)

        if closest_time is None or abs(time_value - target_time) < abs(closest_time - target_time):
            closest_time = time_value
            closest_left = x_left
            closest_right = x_right
            closest_aperture = apertures

    if closest_time is None or closest_left is None or closest_right is None or closest_aperture is None:
        raise ValueError(f"No aperture profile found in {case_dir}")

    return closest_time, closest_left, closest_right, closest_aperture


def plot_cdf_pdf(cases: list[Path]) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharex=True)
    cdf_ax, pdf_ax = axes

    for case_dir in cases:
        x_cdf, y_cdf = read_two_column_file(case_dir / "cdf.txt")
        x_pdf, y_pdf = read_two_column_file(case_dir / "pdf.txt")
        label = case_dir.name
        cdf_ax.plot(x_cdf, y_cdf, linewidth=2, label=label)
        pdf_ax.plot(x_pdf, y_pdf, linewidth=2, label=label)

    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("Residence time / x")
        ax.grid(True, which="both", linestyle="--", alpha=0.35)

    cdf_ax.set_ylabel("CDF")
    cdf_ax.set_title("CDF comparison across cases")
    cdf_ax.legend(title="Case")

    pdf_ax.set_ylabel("PDF")
    pdf_ax.set_title("PDF comparison across cases")
    pdf_ax.legend(title="Case")

    fig.suptitle("Gypsum_rate100time_RE1e-4-1E4particle: CDF / PDF comparison", fontsize=14)
    fig.tight_layout()

    output_path = OUTPUT_DIR / "cdf_pdf_comparison.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_aperture(cases: list[Path]) -> Path:
    fig, ax = plt.subplots(figsize=(8.8, 5.4))

    for case_dir in cases:
        times, aperture_mean, aperture_min, aperture_max = read_aperture_history(case_dir)
        ax.plot(times, aperture_mean, marker="o", linewidth=2, label=case_dir.name)
        ax.fill_between(times, aperture_min, aperture_max, alpha=0.18)

    ax.set_xlabel("Time")
    ax.set_ylabel("Aperture")
    ax.set_title("Average aperture evolution across cases")
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(title="Case")
    fig.tight_layout()

    output_path = OUTPUT_DIR / "aperture_vs_time_comparison.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_aperture_profile(cases: list[Path], target_time: float = 1.0e4) -> Path:
    fig, ax = plt.subplots(figsize=(8.8, 5.4))

    for case_dir in cases:
        matched_time, x_left, x_right, apertures = read_aperture_profile(case_dir, target_time=target_time)
        x_edges = np.concatenate(([x_left[0]], x_right))
        y_steps = np.concatenate((apertures, [apertures[-1]]))
        x_centers = 0.5 * (x_left + x_right)

        ax.step(x_edges, y_steps, where="post", linewidth=2, label=f"{case_dir.name} (t={matched_time:.0f})")
        ax.plot(x_centers, apertures, linestyle="none", marker="o", markersize=4)

    ax.set_xlabel("x coordinate")
    ax.set_ylabel("Aperture")
    ax.set_title("Fracture aperture profile at t = 1E4")
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(title="Case")
    fig.tight_layout()

    output_path = OUTPUT_DIR / "aperture_profile_t1e4.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def write_summary(cases: list[Path]) -> Path:
    reference_case = cases[0]
    ref_x_cdf, ref_y_cdf = read_two_column_file(reference_case / "cdf.txt")
    ref_x_pdf, ref_y_pdf = read_two_column_file(reference_case / "pdf.txt")

    lines = [
        "Case\tMaxAbsDiff_CDF_vs_n=1\tMaxAbsDiff_PDF_vs_n=1\tInitialMeanAperture\tFinalMeanAperture\tApertureGrowthFactor"
    ]

    for case_dir in cases:
        x_cdf, y_cdf = read_two_column_file(case_dir / "cdf.txt")
        x_pdf, y_pdf = read_two_column_file(case_dir / "pdf.txt")
        times, aperture_mean, _, _ = read_aperture_history(case_dir)

        interp_cdf = np.interp(ref_x_cdf, x_cdf, y_cdf)
        interp_pdf = np.interp(ref_x_pdf, x_pdf, y_pdf)
        cdf_diff = float(np.max(np.abs(interp_cdf - ref_y_cdf)))
        pdf_diff = float(np.max(np.abs(interp_pdf - ref_y_pdf)))
        growth_factor = float(aperture_mean[-1] / aperture_mean[0]) if aperture_mean[0] != 0 else float("inf")

        lines.append(
            "\t".join(
                [
                    case_dir.name,
                    f"{cdf_diff:.6g}",
                    f"{pdf_diff:.6g}",
                    f"{aperture_mean[0]:.6g}",
                    f"{aperture_mean[-1]:.6g}",
                    f"{growth_factor:.6g}",
                ]
            )
        )

    output_path = OUTPUT_DIR / "comparison_summary.tsv"
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_path


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    cases = discover_cases(ROOT)
    if not cases:
        raise SystemExit("No case directories named like 'n=*' were found.")

    cdf_pdf_path = plot_cdf_pdf(cases)
    aperture_path = plot_aperture(cases)
    aperture_profile_path = plot_aperture_profile(cases)
    summary_path = write_summary(cases)

    print(f"Saved {cdf_pdf_path}")
    print(f"Saved {aperture_path}")
    print(f"Saved {aperture_profile_path}")
    print(f"Saved {summary_path}")


if __name__ == "__main__":
    main()
