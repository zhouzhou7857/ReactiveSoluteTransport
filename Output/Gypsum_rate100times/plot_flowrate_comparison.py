#!/usr/bin/env python3

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
CASE_PATTERN = re.compile(r"^Re1E-\d+$")
TIME_PATTERN = re.compile(r"#\s*time\s*=\s*([0-9.eE+-]+)")


def discover_cases(root: Path) -> list[Path]:
    cases = [path for path in root.iterdir() if path.is_dir() and CASE_PATTERN.match(path.name)]
    return sorted(cases, key=lambda path: float(path.name.replace("Re1E", "")), reverse=True)


def read_cdf(case_dir: Path) -> tuple[list[float], list[float]]:
    x_values: list[float] = []
    cdf_values: list[float] = []
    with (case_dir / "cdf.txt").open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            x_value, cdf_value = stripped.split()[:2]
            x_values.append(float(x_value))
            cdf_values.append(float(cdf_value))
    return x_values, cdf_values


def read_aperture_history(case_dir: Path) -> dict[int, list[tuple[float, float]]]:
    history: dict[int, list[tuple[float, float]]] = {}
    for step_file in sorted(case_dir.glob("DFN_step*.txt")):
        with step_file.open("r", encoding="utf-8") as handle:
            lines = [line.strip() for line in handle if line.strip()]

        if not lines:
            continue
        first_line = lines[0]
        match = TIME_PATTERN.match(first_line)
        if not match:
            raise ValueError(f"Cannot parse time from {step_file}")

        time_value = float(match.group(1))
        if len(lines) <= 2:
            continue

        for segment_idx, line in enumerate(lines[2:]):
            fields = line.split()
            aperture = float(fields[4])
            history.setdefault(segment_idx, []).append((time_value, aperture))

    if not history:
        raise ValueError(f"No aperture history found in {case_dir}")

    return history


def plot_cdf(cases: list[Path]) -> Path:
    fig, ax = plt.subplots(figsize=(8, 5))

    for case_dir in cases:
        x_values, cdf_values = read_cdf(case_dir)
        ax.plot(x_values, cdf_values, linewidth=2, label=case_dir.name)

    ax.set_xscale("log")
    ax.set_xlabel("Value")
    ax.set_ylabel("CDF")
    ax.set_title("CDF Comparison Across Flow Rates")
    ax.grid(True, which="both", linestyle="--", alpha=0.35)
    ax.legend(title="Case")
    fig.tight_layout()

    output = ROOT / "cdf_comparison.png"
    fig.savefig(output, dpi=300)
    plt.close(fig)
    return output


def plot_aperture(cases: list[Path]) -> Path:
    fig, ax = plt.subplots(figsize=(8, 5))

    for case_dir in cases:
        history = read_aperture_history(case_dir)

        if len(history) == 1:
            segment_id = next(iter(history))
            series = sorted(history[segment_id], key=lambda item: item[0])
            times = [item[0] for item in series]
            apertures = [item[1] for item in series]
            ax.plot(times, apertures, marker="o", linewidth=2, label=case_dir.name)
            continue

        for segment_id, series in sorted(history.items()):
            series = sorted(series, key=lambda item: item[0])
            times = [item[0] for item in series]
            apertures = [item[1] for item in series]
            ax.plot(
                times,
                apertures,
                marker="o",
                linewidth=1.6,
                label=f"{case_dir.name} seg{segment_id}",
            )

    ax.set_xlabel("Time")
    ax.set_ylabel("Aperture")
    ax.set_title("Fracture Segment Aperture Growth Over Time")
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(title="Case")
    fig.tight_layout()

    output = ROOT / "aperture_vs_time_comparison.png"
    fig.savefig(output, dpi=300)
    plt.close(fig)
    return output


def main() -> None:
    cases = discover_cases(ROOT)
    if not cases:
        raise SystemExit("No case directories matching 'Re1E-*' were found.")

    cdf_output = plot_cdf(cases)
    aperture_output = plot_aperture(cases)
    print(f"Saved {cdf_output}")
    print(f"Saved {aperture_output}")


if __name__ == "__main__":
    main()
