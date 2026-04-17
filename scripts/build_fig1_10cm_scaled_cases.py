#!/usr/bin/env python3
from pathlib import Path
import shutil

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_RUNS_DIR = REPO_ROOT / "Output" / "fig1_dfn_runs"
CASE_ROOT = REPO_ROOT / "Output" / "fig1_10cm_scaled_cases"
SCALE = 0.01  # 10 m -> 0.1 m = 10 cm

P_VALUES = [9, 12, 16]
A_VALUES = [1.5, 2.0, 2.5]


def case_name(p: int, a: float) -> str:
    return f"p{p}_a{str(a).replace('.', 'p')}"


def format_a(a: float) -> str:
    return str(int(a)) if float(a).is_integer() else str(a)


def write_file_mode_dfn(dst: Path, data: np.ndarray) -> None:
    with dst.open("w", encoding="ascii") as f:
        f.write("file\n")
        f.write(f"{len(data)}\n")
        for row in data:
            f.write(
                f"{row[0]:.12g} {row[1]:.12g} {row[2]:.12g} {row[3]:.12g} {row[4]:.12g}\n"
            )


def load_scaled_cases():
    scaled = {}
    all_apertures = []
    for p in P_VALUES:
        for a in A_VALUES:
            src_raw = SRC_RUNS_DIR / case_name(p, a) / "DFN_raw.txt"
            data = np.loadtxt(src_raw)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            data_scaled = data.copy()
            data_scaled[:, :4] *= SCALE
            data_scaled[:, 4] *= SCALE
            scaled[(p, a)] = data_scaled
            all_apertures.append(data_scaled[:, 4])
    return scaled, np.concatenate(all_apertures)


def build_case_folders():
    scaled_cases, all_apertures = load_scaled_cases()
    CASE_ROOT.mkdir(parents=True, exist_ok=True)

    summary_lines = []
    for p in P_VALUES:
        for a in A_VALUES:
            name = case_name(p, a)
            case_dir = CASE_ROOT / name
            input_dir = case_dir / "input"
            output_dir = case_dir / "output"
            input_dir.mkdir(parents=True, exist_ok=True)
            output_dir.mkdir(parents=True, exist_ok=True)

            data_scaled = scaled_cases[(p, a)]
            write_file_mode_dfn(input_dir / "DFN_scaled_10cm.txt", data_scaled)

            with (input_dir / "Domain_10cm_square.txt").open("w", encoding="ascii") as f:
                f.write("0.1 0.1\n")
                f.write("1e-10 0.05\n")
                f.write("1.0 0.0\n")

            with (input_dir / "README.txt").open("w", encoding="ascii") as f:
                f.write(f"Scaled 10 cm case for (p,a)=({p},{format_a(a)})\n")
                f.write("Source: Output/fig1_dfn_runs original 10 m DFN_raw.txt\n")
                f.write("Scaling applied uniformly to x, y, and aperture with factor 0.01\n")
                f.write("DFN file format: file + segment count + x1 y1 x2 y2 aperture\n")

            shutil.copy2(input_dir / "DFN_scaled_10cm.txt", output_dir / "DFN_raw.txt")
            shutil.copy2(input_dir / "DFN_scaled_10cm.txt", output_dir / "DFN.txt")
            shutil.copy2(input_dir / "DFN_scaled_10cm.txt", output_dir / "DFN_init.txt")

            lengths = np.sqrt(
                (data_scaled[:, 2] - data_scaled[:, 0]) ** 2
                + (data_scaled[:, 3] - data_scaled[:, 1]) ** 2
            )
            summary_lines.append(
                f"{name}: fractures={len(data_scaled)} mean_length_m={lengths.mean():.6g} "
                f"mean_aperture_m={data_scaled[:,4].mean():.6g}"
            )

    with (CASE_ROOT / "summary.txt").open("w", encoding="ascii") as f:
        f.write("Scaled 10 cm versions of the 9 Song Fig. 1 DFN cases\n")
        f.write("Scaling factor applied to geometry and aperture: 0.01\n\n")
        for line in summary_lines:
            f.write(line + "\n")

    return scaled_cases, all_apertures


def plot_scaled_cases(scaled_cases, all_apertures):
    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "stix",
            "axes.linewidth": 0.8,
            "xtick.direction": "out",
            "ytick.direction": "out",
        }
    )

    norm = Normalize(vmin=float(all_apertures.min()), vmax=float(all_apertures.max()))
    cmap = LinearSegmentedColormap.from_list("black_red", ["#111111", "#ff1a1a"])

    fig, axes = plt.subplots(len(P_VALUES), len(A_VALUES), figsize=(8.5, 9.6))
    fig.subplots_adjust(left=0.07, right=0.985, top=0.985, bottom=0.08, wspace=0.18, hspace=0.34)

    for i, p in enumerate(P_VALUES):
        for j, a in enumerate(A_VALUES):
            ax = axes[i, j]
            data = scaled_cases[(p, a)]
            data_plot = data.copy()
            # Display the scaled DFN in a 10 cm x 10 cm frame with origin at
            # the lower-left corner instead of the domain center.
            data_plot[:, [0, 2]] += 0.05
            data_plot[:, [1, 3]] += 0.05
            data_plot[:, :4] *= 100.0  # m -> cm for plotting only

            segments = data_plot[:, :4].reshape(-1, 2, 2)
            apertures = data[:, 4]
            lc = LineCollection(
                segments,
                cmap=cmap,
                norm=norm,
                linewidths=1.1,
                capstyle="round",
            )
            lc.set_array(apertures)
            ax.add_collection(lc)
            ax.set_xlim(0.0, 10.0)
            ax.set_ylim(0.0, 10.0)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xticks([0.0, 5.0, 10.0])
            ax.set_yticks([0.0, 5.0, 10.0])
            ax.tick_params(axis="both", labelsize=8, length=3, width=0.8, pad=6)
            for spine in ax.spines.values():
                spine.set_color("#b9b9b9")
                spine.set_linewidth(0.8)
            if i == len(P_VALUES) - 1:
                ax.set_xlabel("x (cm)", fontsize=8, labelpad=4)
            if j == 0:
                ax.set_ylabel("y (cm)", fontsize=8, labelpad=8)
            ax.text(
                0.5,
                1.04,
                rf"$(p,\, a)=({p},\, {format_a(a)})$",
                transform=ax.transAxes,
                ha="center",
                va="bottom",
                fontsize=11,
            )

    out_png = CASE_ROOT / "fig1_like_10cm.png"
    out_pdf = CASE_ROOT / "fig1_like_10cm.pdf"
    fig.savefig(out_png, dpi=300, facecolor="white", bbox_inches="tight")
    fig.savefig(out_pdf, facecolor="white", bbox_inches="tight")
    print(f"Saved {out_png}")
    print(f"Saved {out_pdf}")


def main():
    scaled_cases, all_apertures = build_case_folders()
    plot_scaled_cases(scaled_cases, all_apertures)
    print(f"Saved scaled cases under {CASE_ROOT}")


if __name__ == "__main__":
    main()
