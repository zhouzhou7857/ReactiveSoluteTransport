#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize


REPO_ROOT = Path(__file__).resolve().parents[1]
RUNS_DIR = REPO_ROOT / "Output" / "fig1_10cm_direct_generated"
OUT_PNG = RUNS_DIR / "fig1_10cm_direct_label_top.png"
OUT_PDF = RUNS_DIR / "fig1_10cm_direct_label_top.pdf"

P_VALUES = [9, 12, 16]
A_VALUES = [1.5, 2.0, 2.5]


def case_dir(p: int, a: float) -> Path:
    return RUNS_DIR / f"p{p}_a{str(a).replace('.', 'p')}"


def format_a(a: float) -> str:
    return str(int(a)) if float(a).is_integer() else str(a)


def load_case(path: Path):
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def main():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "stix",
            "axes.linewidth": 0.8,
            "xtick.direction": "out",
            "ytick.direction": "out",
        }
    )

    loaded = {}
    all_apertures = []
    for p in P_VALUES:
        for a in A_VALUES:
            data = load_case(case_dir(p, a) / "output" / "DFN_raw.txt")
            loaded[(p, a)] = data
            all_apertures.append(data[:, 4])

    all_apertures = np.concatenate(all_apertures)
    norm = Normalize(vmin=float(all_apertures.min()), vmax=float(all_apertures.max()))
    cmap = LinearSegmentedColormap.from_list("black_red", ["#111111", "#ff1a1a"])

    fig, axes = plt.subplots(len(P_VALUES), len(A_VALUES), figsize=(8.5, 9.6))
    fig.subplots_adjust(left=0.07, right=0.985, top=0.985, bottom=0.08, wspace=0.18, hspace=0.34)

    for i, p in enumerate(P_VALUES):
        for j, a in enumerate(A_VALUES):
            ax = axes[i, j]
            data = loaded[(p, a)].copy()
            data[:, [0, 2]] += 0.05
            data[:, [1, 3]] += 0.05
            data[:, :4] *= 100.0
            segments = data[:, :4].reshape(-1, 2, 2)
            apertures = loaded[(p, a)][:, 4]

            lc = LineCollection(segments, cmap=cmap, norm=norm, linewidths=1.1, capstyle="round")
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

    fig.savefig(OUT_PNG, dpi=300, facecolor="white", bbox_inches="tight")
    fig.savefig(OUT_PDF, facecolor="white", bbox_inches="tight")
    print(f"Saved {OUT_PNG}")
    print(f"Saved {OUT_PDF}")


if __name__ == "__main__":
    main()
