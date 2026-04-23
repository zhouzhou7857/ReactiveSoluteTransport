#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.ticker import FuncFormatter


REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = REPO_ROOT / "Input" / "DFN_files" / "song_selected_10cm"
OUTPUT_ROOT = REPO_ROOT / "Output" / "song_selected_10cm_triptych"
OUTPUT_PNG = OUTPUT_ROOT / "song_selected_10cm_triptych_viridis.png"
OUTPUT_PDF = OUTPUT_ROOT / "song_selected_10cm_triptych_viridis.pdf"

CASES = [
    ("p9_a1p5", 9, 1.5),
    ("p12_a2p0", 12, 2.0),
    ("p16_a2p5", 16, 2.5),
]


def format_a(a: float) -> str:
    return str(int(a)) if float(a).is_integer() else str(a)


def load_filemode_dfn(path: Path) -> np.ndarray:
    with path.open("r", encoding="ascii") as f:
        mode = f.readline().strip()
        if mode != "file":
            raise ValueError(f"Unsupported DFN mode in {path}: {mode}")
        count = int(f.readline().strip())

    data = np.loadtxt(path, skiprows=2)
    if count == 0:
        return np.empty((0, 5), dtype=float)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape != (count, 5):
        raise ValueError(f"Unexpected DFN shape in {path}: {data.shape}, expected ({count}, 5)")
    return data


def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

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
    for case_name, _, _ in CASES:
        path = INPUT_ROOT / case_name / f"{case_name}_raw_filemode.txt"
        loaded[case_name] = load_filemode_dfn(path)

    norm = Normalize(vmin=0.0, vmax=3.0e-5)
    cmap = plt.get_cmap("viridis")

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.85))
    fig.subplots_adjust(left=0.08, right=0.90, top=0.83, bottom=0.22, wspace=0.18)

    last_lc = None
    for idx, (ax, (case_name, p, a)) in enumerate(zip(axes, CASES)):
        data = loaded[case_name].copy()
        data[:, [0, 2]] += 0.05
        data[:, [1, 3]] += 0.05
        data[:, :4] *= 100.0
        segments = data[:, :4].reshape(-1, 2, 2)

        lc = LineCollection(segments, cmap=cmap, norm=norm, linewidths=1.1, capstyle="round")
        lc.set_array(loaded[case_name][:, 4])
        ax.add_collection(lc)
        last_lc = lc

        ax.set_xlim(0.0, 10.0)
        ax.set_ylim(0.0, 10.0)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([0.0, 5.0, 10.0])
        ax.set_yticks([0.0, 5.0, 10.0])
        ax.tick_params(axis="both", labelsize=10, length=3, width=0.8, pad=6)
        ax.set_xlabel("x (cm)", fontsize=10, labelpad=3)
        if idx == 0:
            ax.set_ylabel("y (cm)", fontsize=10, labelpad=6)

        for spine in ax.spines.values():
            spine.set_color("#b9b9b9")
            spine.set_linewidth(0.8)

        ax.text(
            0.5,
            1.04,
            rf"$(p,\, a)=({p},\, {format_a(a)})$",
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=13,
        )

    cbar = fig.colorbar(last_lc, ax=axes, fraction=0.035, pad=0.02)
    cbar.set_ticks([0.0e-5, 1.0e-5, 2.0e-5, 3.0e-5])
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x*1.0e6:.0f}"))
    cbar.set_label("Aperture (μm)", fontsize=10, labelpad=6)
    cbar.ax.tick_params(labelsize=10, length=3, width=0.8)

    fig.savefig(OUTPUT_PNG, dpi=600, facecolor="white", bbox_inches="tight")
    fig.savefig(OUTPUT_PDF, facecolor="white", bbox_inches="tight")
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
