"""Fig. 1 -- SPARTA overview: workflow, the two operators on a toy lattice, and the median real section."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from scipy.stats import rankdata  # noqa: E402
import scipy.sparse as sp  # noqa: E402

import figstyle as S  # noqa: E402
from isdata import j  # noqa: E402
from sparta.barrier import compute_b_cell, compute_b_cell_field, compute_b_mab  # noqa: E402
from sparta.io_ import Paths, load_config, load_graph  # noqa: E402
from sparta.node_tables import coords_from_nodes, load_nodes, scores_from_nodes  # noqa: E402

OUT = Path(__file__).resolve().parents[2] / "results" / "figures" / "is"
CMAP_CELL = LinearSegmentedColormap.from_list("cell", ["#eef5fd", "#b7d3f6", "#5598e7", "#1c5cab", "#0d366b"])
CMAP_MAB = LinearSegmentedColormap.from_list("mab", ["#fdf1ea", "#f6b896", "#eb6834", "#b8461b", "#6e2a10"])
CMAP_GREY = LinearSegmentedColormap.from_list("grey", ["#f4f3f0", "#c3c2b7", "#898781", "#52514e", "#1f1f1e"])
EXAMPLE = "CSCC04"     # the section whose association equals the primary-cohort median


def toy_lattice(nx_=14, ny_=9, h=1.0):
    pts = []
    for jj in range(ny_):
        for i in range(nx_):
            pts.append((i * h + (0.5 * h if jj % 2 else 0.0), jj * h * np.sqrt(3) / 2))
    xy = np.array(pts)
    pairs = np.array(sorted(cKDTree(xy).query_pairs(1.05 * h)))
    n = len(xy)
    A = sp.coo_matrix((np.ones(2 * len(pairs)), (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])),
                      shape=(n, n)).tocsr()
    return xy, A, pairs


def box(ax, x, y, w, h, text, fc="#ffffff", fs=6.8, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.0,rounding_size=0.02", fc=fc, ec=S.INK2,
                                lw=0.6, transform=ax.transAxes))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=S.INK,
            transform=ax.transAxes, fontweight="bold" if bold else "normal", linespacing=1.2)


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), transform=ax.transAxes, arrowstyle="-|>",
                                 mutation_scale=6, lw=0.6, color=S.INK2, shrinkA=0, shrinkB=0))


def fence(ax, xy, edges, half=0.42, lw=0.7, color=None):
    color = color or S.INK
    segs = []
    for (u, v) in edges:
        mid = xy[[u, v]].mean(axis=0)
        d = xy[v] - xy[u]
        nrm = np.array([-d[1], d[0]]) / (np.linalg.norm(d) + 1e-12) * half * np.linalg.norm(d)
        segs.append([[mid[0] - nrm[0], mid[1] - nrm[1]], [mid[0] + nrm[0], mid[1] + nrm[1]]])
    from matplotlib.collections import LineCollection
    ax.add_collection(LineCollection(segs, colors=color, linewidths=lw, zorder=4, capstyle="round"))


def spot_area(ax, xy, pitch, frac=0.95):
    fig = ax.figure
    fig.canvas.draw()
    bbox = ax.get_window_extent()
    x0, x1 = ax.get_xlim()
    pts_per_um = bbox.width / fig.dpi * 72 / (x1 - x0)
    d = pitch * pts_per_um * frac
    return d ** 2


def main():
    S.apply()
    cfg = load_config(None)
    P = Paths(cfg)
    fig = plt.figure(figsize=(S.FULL_W, 142 * S.MM))

    # ---------------- a: workflow ----------------
    ax0 = fig.add_axes([0.0, 0.865, 1.0, 0.125])
    ax0.axis("off")
    labels = ["Spot × gene counts\n+ spot coordinates",
              "Signature scores,\nrank-normalised\nwithin the section",
              "Radius graph;\nvessels, immune-entry\nsources, tumour-core sinks",
              "Minimum cut  (B$_{\\rm cell}$)\nScreened diffusion–\nabsorption  (B$_{\\rm mAb}$)",
              "Field association,\nsurrogate and structural\nnulls, patient model"]
    w, gap = 0.172, 0.0275
    for k, t in enumerate(labels):
        x = 0.012 + k * (w + gap)
        box(ax0, x, 0.08, w, 0.80, t, fc="#f6f6f4" if k == 3 else "#ffffff", bold=False)
        if k < 4:
            arrow(ax0, x + w + 0.003, 0.48, x + w + gap - 0.003, 0.48)
    S.panel(ax0, "a", x=0.0, y=0.93)

    # ---------------- b: minimum cut on a toy lattice ----------------
    xy, A, pairs = toy_lattice()
    n = len(xy)
    rng = np.random.default_rng(3)
    cx, cy = xy[:, 0], xy[:, 1]
    ecm = np.full(n, 0.10) + 0.06 * rng.random(n)
    band = (cx > 5.0) & (cx < 7.1)
    ecm[band] = 0.92
    ecm[band & (cy > 5.1)] = 0.20                       # a gap in the band
    caf = ecm.copy()
    source = np.flatnonzero(cx < 0.6)
    sink = np.flatnonzero(cx > 11.4)
    out = compute_b_cell(A, ecm, caf, source, sink, a=3.0, b_ecm=8.0, c_caf=4.0, return_capacity=True)
    axb = fig.add_axes([0.035, 0.475, 0.43, 0.36])
    axc = fig.add_axes([0.535, 0.475, 0.43, 0.36])
    for ax in (axb, axc):
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_xlim(-0.7, 14.2)
        ax.set_ylim(-2.6, 8.9)
    for (u, v), c in zip(out["edge_pairs"], out["edge_capacity"]):
        axb.plot(xy[[u, v], 0], xy[[u, v], 1], color=CMAP_CELL(0.18 + 0.82 * c), lw=0.3 + 1.7 * c,
                 solid_capstyle="round", zorder=1)
    axb.scatter(cx, cy, s=10, c=["#52514e" if e > 0.5 else "#ffffff" for e in ecm], edgecolors=S.INK2,
                linewidths=0.4, zorder=2)
    fence(axb, xy, out["cut_edges"], half=0.45, lw=1.1)
    axb.scatter(xy[source, 0], xy[source, 1], s=18, c=S.CELL, edgecolors=S.INK, linewidths=0.4, zorder=5)
    axb.scatter(xy[sink, 0], xy[sink, 1], s=18, marker="s", c=S.INK, edgecolors="white", linewidths=0.4, zorder=5)
    axb.text(0.0, 7.55, "sources", fontsize=6.8, color=S.INK2, ha="left", va="bottom")
    axb.text(6.05, 7.55, "matrix-rich band (grey)\nwith a gap", fontsize=6.8, color=S.INK2, ha="center", va="bottom",
             linespacing=1.0)
    axb.text(13.6, 7.55, "sinks", fontsize=6.8, color=S.INK2, ha="right", va="bottom")
    axb.text(6.6, -0.95, r"$c_{uv}=\sigma(a-b\,\bar E_{uv}-c\,\bar F_{uv})$,   "
                         r"$B_{\rm cell}=1/(\mathrm{max\ flow}+\epsilon)$",
             fontsize=7.6, ha="center", va="top", color=S.INK)
    axb.text(6.6, -1.95, "edge width/shade = capacity;  black bars = minimum cut", fontsize=6.6,
             ha="center", va="top", color=S.INK2)
    S.panel(axb, "b", x=-0.05, y=0.93)

    # ---------------- c: screened diffusion-absorption on the same lattice ----------------
    crosslink = np.where(band & (cy <= 5.1), 0.95, 0.15 + 0.1 * rng.random(n))
    ag = np.full(n, 0.03)
    rim = (cx > 9.0) & (cx < 10.1)
    ag[rim] = 0.9
    vessel = source
    bm = compute_b_mab(A, ecm, crosslink, ag, vessel)["b_mab"]
    for (u, v) in pairs:
        axc.plot(xy[[u, v], 0], xy[[u, v], 1], color=S.HAIR, lw=0.6, zorder=1)
    vmax = float(np.nanpercentile(bm, 99))
    sc = axc.scatter(cx, cy, s=24, c=bm, cmap=CMAP_MAB, vmin=0, vmax=vmax, edgecolors="none", zorder=2)
    absorb = np.flatnonzero(rim)
    axc.scatter(xy[absorb, 0], xy[absorb, 1], s=34, facecolors="none", edgecolors=S.INK, linewidths=0.6, zorder=3)
    axc.scatter(xy[vessel, 0], xy[vessel, 1], s=18, c=S.MAB, edgecolors=S.INK, linewidths=0.4, zorder=4)
    axc.text(0.0, 7.55, "vessels (φ = 1)", fontsize=6.8, color=S.INK2, ha="left", va="bottom")
    axc.text(6.05, 7.55, "size-excluding band\nwith a gap", fontsize=6.8, color=S.INK2, ha="center", va="bottom",
             linespacing=1.0)
    axc.text(9.55, 7.55 - 0.0, "absorbing\nrim", fontsize=6.8, color=S.INK2, ha="left", va="bottom", linespacing=1.0)
    axc.text(6.6, -0.95, r"$(L_g+\mathrm{diag}\,\kappa)\,\varphi=0$,   "
                         r"$g_{uv}=g_0\,e^{-\lambda\bar E_{uv}}\,(1-r/\xi_{uv})_+^{2}$,   $B_{\rm mAb}=-\log\varphi$",
             fontsize=7.6, ha="center", va="top", color=S.INK)
    axc.text(6.6, -1.95, "fill = B$_{\\rm mAb}$ (light: reached, dark: depleted);  rings = absorbing spots",
             fontsize=6.6, ha="center", va="top", color=S.INK2)
    S.panel(axc, "c", x=-0.05, y=0.93)

    # ---------------- d-g: the median real section ----------------
    sid = EXAMPLE
    nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
    A2, D2, src2, snk2, ves2, _ = load_graph(P.graph(sid))
    Sx = scores_from_nodes(nodes)
    xy2 = coords_from_nodes(nodes)
    xy2 = xy2 - xy2.min(axis=0)
    xy2[:, 1] = xy2[:, 1].max() - xy2[:, 1]
    bcf = compute_b_cell_field(A2, Sx["ecm"], Sx["caf"], src2, **cfg["barrier"]["b_cell"])["b_cell_field"]
    bmf = compute_b_mab(A2, Sx["ecm"], Sx["crosslink"], Sx["ag_target"], ves2, **cfg["barrier"]["b_mab"])["b_mab"]
    cut2 = json.load(open(P.mincut(sid), encoding="utf-8"))["cut_edges"]
    rho = j("spatial_null_check.json")["per_slide"][sid]["real_rho_partial"]
    comp = np.full(len(xy2), 0)
    comp[ves2] = 1
    comp[src2] = 2
    comp[snk2] = 3
    pad = 90
    titles = ["Compartments", "ECM score (rank)", "B$_{\\rm cell}$ field (rank) + minimum cut",
              "B$_{\\rm mAb}$ field (rank)"]
    W = 0.235
    for k in range(4):
        ax = fig.add_axes([0.01 + k * 0.2475, 0.02, W, 0.36])
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_xlim(xy2[:, 0].min() - pad, xy2[:, 0].max() + pad)
        ax.set_ylim(xy2[:, 1].min() - pad - 330, xy2[:, 1].max() + pad)
        s = spot_area(ax, xy2, 100.0, frac=0.92)
        if k == 0:
            cols = np.array(["#e9e8e3", S.MAB, S.CELL, S.INK])[comp]
            ax.scatter(xy2[:, 0], xy2[:, 1], c=cols, s=s, edgecolors="none", rasterized=True)
            hs = [Line2D([0], [0], marker="o", ls="none", ms=3.6, mfc=c, mec="none", label=l)
                  for c, l in [(S.MAB, "vessel"), (S.CELL, "source"), (S.INK, "sink"), ("#e9e8e3", "other")]]
            ax.legend(handles=hs, loc="lower left", bbox_to_anchor=(-0.02, -0.04), ncol=4, fontsize=6.2,
                      handletextpad=0.05, columnspacing=0.5, frameon=True, facecolor="white",
                      edgecolor="none", framealpha=1.0, borderpad=0.15)
        else:
            val = [None, Sx["ecm"], rankdata(bcf) / len(bcf), rankdata(bmf) / len(bmf)][k]
            cm = [None, CMAP_GREY, CMAP_CELL, CMAP_MAB][k]
            sc2 = ax.scatter(xy2[:, 0], xy2[:, 1], c=val, s=s, cmap=cm, vmin=0, vmax=1, edgecolors="none",
                             rasterized=True)
            if k == 2:
                fence(ax, xy2, cut2, half=0.55, lw=0.9)
            cax = ax.inset_axes([0.03, 0.045, 0.36, 0.026])
            cb = fig.colorbar(sc2, cax=cax, orientation="horizontal", ticks=[0, 1])
            cb.ax.set_xticklabels(["low", "high"], fontsize=6.2)
            cb.outline.set_linewidth(0.4)
            cb.ax.tick_params(length=1.5, width=0.4, pad=1)
        ax.text(0.0, 1.0, titles[k], transform=ax.transAxes, fontsize=7, ha="left", va="bottom", color=S.INK)
        if k > 0:
            xm = xy2[:, 0].max()
            ym = xy2[:, 1].min() - pad - 120
            ax.plot([xm - 1000, xm], [ym, ym], color=S.INK, lw=1.0)
            ax.text(xm - 500, ym - 30, "1 mm", fontsize=6.2, ha="center", va="top")
        S.panel(ax, "defg"[k], x=-0.04, y=1.07)
    fig.text(0.01, 0.405, f"Primary-cohort section with the median field association ({sid}: cSCC, Visium; "
             f"partial ρ = {rho:.2f})", fontsize=7, color=S.INK, ha="left", va="bottom", style="italic")
    S.save(fig, OUT, "Fig1_overview")


if __name__ == "__main__":
    main()
