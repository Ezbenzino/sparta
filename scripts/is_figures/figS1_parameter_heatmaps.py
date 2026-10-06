"""Fig. S1 -- xi0/beta/lambda parameter sensitivity heatmaps."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

import figstyle as S  # noqa: E402
from isdata import j  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/figures/is/supplement"


def median_grids(data: dict, grid_name: str):
    slides = list(data["slides"].values())
    grid = slides[0][grid_name]["grid"]
    crosslink = np.median(np.stack([np.asarray(s[grid_name]["crosslink_pct"], float)
                                    for s in slides]), axis=0)
    potential = np.median(np.stack([np.asarray(s[grid_name]["potential"], float)
                                    for s in slides]), axis=0)
    log_potential = np.log10(np.clip(potential, 1e-6, None))
    return grid, crosslink, log_potential


def mark_default(ax, xs, ys, x_default, y_default):
    col = np.argmin(np.abs(np.asarray(xs, float) - x_default))
    row = np.argmin(np.abs(np.asarray(ys, float) - y_default))
    ax.add_patch(Rectangle((col - 0.5, row - 0.5), 1, 1, fill=False,
                           edgecolor=S.INK, lw=1.2))


def heatmap(ax, values, xs, ys, xlabel, ylabel, cmap, vmin, vmax, cbar_label,
            x_default=None, y_default=None):
    im = ax.imshow(values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                   origin="upper")
    ax.set_xticks(range(len(xs)))
    ax.set_xticklabels([f"{x:g}" for x in xs], fontsize=6.2)
    ax.set_yticks(range(len(ys)))
    ax.set_yticklabels([f"{y:g}" for y in ys], fontsize=6.2)
    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel(ylabel, fontsize=7)
    if x_default is not None:
        mark_default(ax, xs, ys, x_default, y_default)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.ax.tick_params(labelsize=6)
    cb.set_label(cbar_label, fontsize=6.5)
    return im


def main():
    S.apply()
    data = j("bmab_sensitivity.json")
    lam_grid, lam_cross, lam_pot = median_grids(data, "lam_x_beta")
    xi_grid, xi_cross, xi_pot = median_grids(data, "xi0_x_beta")

    fig, axs = plt.subplots(2, 2, figsize=(S.FULL_W, 118 * S.MM))
    cmap_cross = S.mpl.colormaps["RdBu_r"].resampled(256)
    cmap_pot = S.mpl.colormaps["viridis"].resampled(256)

    heatmap(axs[0, 0], lam_cross, lam_grid["xs"], lam_grid["ys"],
            r"$\lambda$ (ECM decay)", r"$\beta$ (mesh contraction)",
            cmap_cross, 0, 100, "Crosslink contribution (%)", 3.0, 3.0)
    heatmap(axs[0, 1], lam_pot, lam_grid["xs"], lam_grid["ys"],
            r"$\lambda$ (ECM decay)", r"$\beta$ (mesh contraction)",
            cmap_pot, None, None, r"$\log_{10}$(dissociation potential)")
    heatmap(axs[1, 0], xi_cross, xi_grid["xs"], xi_grid["ys"],
            r"$\xi_0$ (baseline mesh, nm)", r"$\beta$ (mesh contraction)",
            cmap_cross, 0, 100, "Crosslink contribution (%)", 20.0, 3.0)
    heatmap(axs[1, 1], xi_pot, xi_grid["xs"], xi_grid["ys"],
            r"$\xi_0$ (baseline mesh, nm)", r"$\beta$ (mesh contraction)",
            cmap_pot, None, None, r"$\log_{10}$(dissociation potential)")

    for ax, label in zip(axs.ravel(), "abcd"):
        S.panel(ax, label, x=-0.20, y=1.04)
    fig.tight_layout(h_pad=1.8, w_pad=1.2)
    S.save(fig, OUT, "FigS1_parameter_sensitivity", eps=False)


if __name__ == "__main__":
    main()
