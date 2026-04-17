#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize


REPO_ROOT = Path(__file__).resolve().parents[1]
RUNS_DIR = REPO_ROOT / "Output" / "fig1_dfn_runs"
OUT_PNG = RUNS_DIR / "fig1_like.png"
OUT_PDF = RUNS_DIR / "fig1_like.pdf"

P_VALUES = [9, 12, 16]
A_VALUES = [1.5, 2.0, 2.5]


def case_dir(p: int, a: float) -> Path:
    return RUNS_DIR / f"p{p}_a{str(a).replace('.', 'p')}"


def load_segments(path: Path):
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    segments = data[:, :4].reshape(-1, 2, 2)
    apertures = data[:, 4]
    return segments, apertures


def format_a(a: float) -> str:
    return str(int(a)) if float(a).is_integer() else str(a)


def main():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "stix",
            "axes.linewidth": 0.8,
        }
    )

    loaded = {}
    all_apertures = []
    for p in P_VALUES:
        for a in A_VALUES:
            raw_path = case_dir(p, a) / "DFN_raw.txt"
            segments, apertures = load_segments(raw_path)
            loaded[(p, a)] = (segments, apertures)
            all_apertures.append(apertures)

    all_apertures = np.concatenate(all_apertures)
    norm = Normalize(vmin=float(all_apertures.min()), vmax=float(all_apertures.max()))
    cmap = LinearSegmentedColormap.from_list("black_red", ["#111111", "#ff1a1a"])

    fig, axes = plt.subplots(len(P_VALUES), len(A_VALUES), figsize=(8.1, 9.3))
    fig.subplots_adjust(
        left=0.06,
        right=0.985,
        top=0.985,
        bottom=0.05,
        wspace=0.10,
        hspace=0.24,
    )

    for i, p in enumerate(P_VALUES):
        for j, a in enumerate(A_VALUES):
            ax = axes[i, j]
            segments, apertures = loaded[(p, a)]
            lc = LineCollection(
                segments,
                cmap=cmap,
                norm=norm,
                linewidths=1.1,
                capstyle="round",
            )
            lc.set_array(apertures)
            ax.add_collection(lc)
            ax.set_xlim(-5, 5)
            ax.set_ylim(-5, 5)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_color("#b9b9b9")
                spine.set_linewidth(0.8)
            ax.text(
                0.5,
                -0.045,
                rf"$(p,\, a)=({p},\, {format_a(a)})$",
                transform=ax.transAxes,
                ha="center",
                va="top",
                fontsize=11.5,
            )

    fig.savefig(OUT_PNG, dpi=300, facecolor="white", bbox_inches="tight")
    fig.savefig(OUT_PDF, facecolor="white", bbox_inches="tight")
    print(f"Saved {OUT_PNG}")
    print(f"Saved {OUT_PDF}")


if __name__ == "__main__":
    main()
