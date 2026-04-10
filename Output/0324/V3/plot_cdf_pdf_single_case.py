#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
CDF_FILE = ROOT / "cdf.txt"
PDF_FILE = ROOT / "pdf.txt"
OUTPUT_FIG = ROOT / "cdf_pdf_single_case_V3.png"

plt.style.use("seaborn-v0_8-whitegrid")


def read_two_column_file(file_path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(file_path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected format in {file_path}")
    return data[:, 0], data[:, 1]


def main() -> None:
    x_cdf, y_cdf = read_two_column_file(CDF_FILE)
    x_pdf, y_pdf = read_two_column_file(PDF_FILE)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0), constrained_layout=True)
    cdf_ax, pdf_ax = axes

    cdf_ax.plot(x_cdf, y_cdf, linewidth=2.2, color="#1d3557")
    cdf_ax.set_xscale("log")
    cdf_ax.set_xlabel("Residence time / x")
    cdf_ax.set_ylabel("CDF")
    cdf_ax.set_title("V3: CDF")
    cdf_ax.grid(True, which="both", linestyle="--", alpha=0.35)

    pdf_ax.plot(x_pdf, y_pdf, linewidth=2.2, color="#c1121f")
    pdf_ax.set_xscale("log")
    pdf_ax.set_xlabel("Residence time / x")
    pdf_ax.set_ylabel("PDF")
    pdf_ax.set_title("V3: PDF")
    pdf_ax.grid(True, which="both", linestyle="--", alpha=0.35)

    fig.suptitle("V3: CDF and PDF")
    fig.savefig(OUTPUT_FIG, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
