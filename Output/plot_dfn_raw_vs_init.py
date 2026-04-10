import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_PATH = os.path.join(BASE_DIR, "DFN_raw.txt")
INIT_PATH = os.path.join(BASE_DIR, "DFN_init.txt")
OUT_PATH = os.path.join(BASE_DIR, "dfn_raw_vs_init.png")


def read_dfn(path):
    segments = []
    apertures = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            x1, y1, x2, y2, aperture = map(float, parts[:5])
            segments.append([(x1, y1), (x2, y2)])
            apertures.append(aperture)
    return np.asarray(segments, dtype=float), np.asarray(apertures, dtype=float)


def build_linewidths(apertures, aperture_min, aperture_max, width_min=1.0, width_max=6.0):
    if aperture_max <= aperture_min:
        return np.full_like(apertures, 2.5, dtype=float)
    scaled = (apertures - aperture_min) / (aperture_max - aperture_min)
    return width_min + scaled * (width_max - width_min)


def style_axis(ax, title, xlim, ylim):
    ax.set_title(title)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.grid(True, alpha=0.20)


def main():
    raw_segments, raw_apertures = read_dfn(RAW_PATH)
    init_segments, init_apertures = read_dfn(INIT_PATH)

    all_segments = np.concatenate([raw_segments, init_segments], axis=0)
    all_apertures = np.concatenate([raw_apertures, init_apertures], axis=0)

    xs = all_segments[:, :, 0].ravel()
    ys = all_segments[:, :, 1].ravel()
    pad_x = 0.03 * max(xs.max() - xs.min(), 1e-6)
    pad_y = 0.03 * max(ys.max() - ys.min(), 1e-6)
    xlim = (xs.min() - pad_x, xs.max() + pad_x)
    ylim = (ys.min() - pad_y, ys.max() + pad_y)

    aperture_min = float(all_apertures.min())
    aperture_max = float(all_apertures.max())

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True, sharex=True, sharey=True)
    cmap = plt.get_cmap("viridis")

    raw_collection = LineCollection(
        raw_segments,
        linewidths=build_linewidths(raw_apertures, aperture_min, aperture_max),
        array=raw_apertures,
        cmap=cmap,
        alpha=0.95,
    )
    axes[0].add_collection(raw_collection)
    style_axis(axes[0], f"DFN raw ({len(raw_segments)} segments)", xlim, ylim)

    init_collection = LineCollection(
        init_segments,
        linewidths=build_linewidths(init_apertures, aperture_min, aperture_max),
        array=init_apertures,
        cmap=cmap,
        alpha=0.95,
    )
    axes[1].add_collection(init_collection)
    style_axis(axes[1], f"DFN init ({len(init_segments)} segments)", xlim, ylim)

    mappable = plt.cm.ScalarMappable(cmap=cmap)
    mappable.set_array(all_apertures)
    mappable.set_clim(aperture_min, aperture_max)
    fig.colorbar(mappable, ax=axes, label="Aperture [m]", shrink=0.95)

    fig.suptitle("DFN Structure Comparison: raw vs init", fontsize=14)
    fig.savefig(OUT_PATH, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
