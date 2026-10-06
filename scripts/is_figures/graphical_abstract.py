"""SPARTA graphical abstract -- publication version, built from real section data.

Canvas 13.0 x 5.2 cm (w:h = 2.5:1), the size at which Springer displays a graphical
abstract, so every label is set at its final size and nothing is downscaled.

Reading order: one tissue graph (matrix score) -> the two operators computed on that same
graph (minimum-cut cellular barrier; screened diffusion-absorption field) -> the three
levels of evidence in the bottom strip. Colour carries one meaning throughout: grey is the
common input, blue is the cellular barrier, orange is the molecular barrier.

Output: results/figures/graphical_abstract.{pdf,png,tif}
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "is_figures"))

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from scipy.stats import rankdata  # noqa: E402

import figstyle as S  # noqa: E402
from sparta.barrier import compute_b_cell_field, compute_b_mab  # noqa: E402
from sparta.io_ import Paths, load_config, load_graph  # noqa: E402
from sparta.node_tables import coords_from_nodes, load_nodes, scores_from_nodes  # noqa: E402

OUT = ROOT / "results" / "figures"
W_MM, H_MM, DPI = 130.0, 52.0, 500
SECTION = "CSCC05"

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
HAIR, BG, PANEL = "#dfe6ec", "#f5f8fb", "#ffffff"
CELL, MAB, GREEN, MATRIX = S.CELL, S.MAB, "#168b70", "#52514e"
CMAP_CELL = LinearSegmentedColormap.from_list("cell", ["#eef5fd", "#b7d3f6", "#5598e7", "#1c5cab", "#0d366b"])
CMAP_MAB = LinearSegmentedColormap.from_list("mab", ["#fdf1ea", "#f6b896", "#eb6834", "#b8461b", "#6e2a10"])
CMAP_GREY = LinearSegmentedColormap.from_list("grey", ["#f4f3f0", "#c3c2b7", "#898781", "#52514e", "#1f1f1e"])


def txt(ax, x, y, s, size=5.2, color=INK2, weight="normal", ha="left", va="baseline", z=8):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha, va=va, zorder=z)


def rounded(ax, x, y, w, h, fc=PANEL, ec=HAIR, r=1.6, lw=0.6, z=1):
    from matplotlib.patches import FancyBboxPatch
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                fc=fc, ec=ec, lw=lw, zorder=z))


def load_section(sid):
    cfg = load_config(None)
    P = Paths(cfg)
    nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
    A, D, src, snk, ves, _ = load_graph(P.graph(sid))
    Sx = scores_from_nodes(nodes)
    xy = coords_from_nodes(nodes)
    xy = xy - xy.min(axis=0)
    xy[:, 1] = xy[:, 1].max() - xy[:, 1]
    bcf = compute_b_cell_field(A, Sx["ecm"], Sx["caf"], src, **cfg["barrier"]["b_cell"])["b_cell_field"]
    bmf = compute_b_mab(A, Sx["ecm"], Sx["crosslink"], Sx["ag_target"], ves, **cfg["barrier"]["b_mab"])["b_mab"]
    cut = json.load(open(P.mincut(sid), encoding="utf-8"))["cut_edges"]
    pitch = 200.0 if (sid.startswith("CSCC") and sid not in ("CSCC01", "CSCC02", "CSCC03", "CSCC04")) else 100.0
    radius = 150.0 if pitch == 100.0 else 300.0     # manuscript Table: Visium 150 um, first-generation ST 300 um
    return dict(xy=xy, ecm=Sx["ecm"], bcf=bcf, bmf=bmf, cut=cut, A=A, pitch=pitch, radius=radius, n=len(xy))


def draw_map(ax, sec, kind, x0, y0, size):
    """Draw one section map into the mm-square [x0, x0+size] x [y0, y0+size]."""
    xy, pitch = sec["xy"], sec["pitch"]
    span = max(np.ptp(xy[:, 0]), np.ptp(xy[:, 1])) + 2 * pitch
    cx, cy = 0.5 * (xy[:, 0].min() + xy[:, 0].max()), 0.5 * (xy[:, 1].min() + xy[:, 1].max())
    sc = size / span
    X = x0 + size / 2 + (xy[:, 0] - cx) * sc
    Y = y0 + size / 2 + (xy[:, 1] - cy) * sc
    dot_mm = pitch * sc
    area = (dot_mm * 72 / 25.4) ** 2
    if kind == "matrix":
        vals, cm = sec["ecm"], CMAP_GREY
    elif kind == "cell":
        vals, cm = rankdata(sec["bcf"]) / len(sec["bcf"]), CMAP_CELL
    else:
        vals, cm = rankdata(sec["bmf"]) / len(sec["bmf"]), CMAP_MAB
    if kind == "matrix":                       # the mesh itself: this panel is the graph
        A = sec["A"].tocoo()
        segs0 = [[(X[i], Y[i]), (X[j], Y[j])] for i, j in zip(A.row, A.col) if i < j]
        ax.add_collection(LineCollection(segs0, colors="#9aa4ad", linewidths=0.14, alpha=0.55, zorder=2))
    ax.scatter(X, Y, c=vals, s=area, cmap=cm, vmin=0.0, vmax=1.0, edgecolors="none",
               linewidths=0, rasterized=True, zorder=3)
    if kind == "cell":
        segs = []
        for (u, v) in sec["cut"]:
            p, q = np.array([X[u], Y[u]]), np.array([X[v], Y[v]])
            mid, d = (p + q) / 2, q - p
            nrm = np.array([-d[1], d[0]]) / (np.linalg.norm(d) + 1e-12) * (dot_mm / 2)
            segs.append([mid - nrm, mid + nrm])
        ax.add_collection(LineCollection(segs, colors="#0b0b0b", linewidths=0.55, zorder=5))


def main():
    S.apply()
    sec = load_section(SECTION)
    fig = plt.figure(figsize=(W_MM / 25.4, H_MM / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_MM); ax.set_ylim(0, H_MM); ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(plt.Rectangle((0, 0), W_MM, H_MM, fc=BG, ec="none", zorder=0))

    # ------------------------------------------------------------------ header
    txt(ax, 4, 49.6, "SPARTA", 6.2, CELL, "bold")
    txt(ax, 126, 49.6, "structural summaries of tissue architecture \u2014 not clinical predictions",
        4.9, MUTED, ha="right")
    txt(ax, 4, 45.2, "One tissue graph, two transport barriers", 9.4, INK, "bold")
    ax.plot([4, 126], [43.4, 43.4], color=HAIR, lw=0.6, zorder=2)

    # ------------------------------------------------------- three real maps
    size = 28.0
    panels = [
        (4.0, "matrix", MATRIX, "one tissue graph", None,
         "%s · %d spots · %d µm graph" % (SECTION, sec["n"], sec["radius"])),
        (45.5, "cell", CELL, "cellular barrier", r"$B_{\rm cell}=1/F_{\rm max}$", "minimum cut"),
        (87.0, "molecular", MAB, "molecular barrier", r"$B_{\rm mAb}=-\log\varphi$", "screened field"),
    ]
    for x0, kind, col, title, formula, note in panels:
        rounded(ax, x0 - 1.0, 8.3, 34.0, 34.4, fc="white", ec=HAIR)
        draw_map(ax, sec, kind, x0 + 2.0, 11.2, size)
        txt(ax, x0 - 1.0, 41.0, title, 5.6, col, "bold")
        if formula is None:
            txt(ax, x0 - 1.0, 9.9, note, 4.9, MUTED)
        else:
            txt(ax, x0 - 1.0, 9.9, formula, 6.4, col)
            txt(ax, x0 + 17.0, 9.9, note, 4.9, MUTED)

    # ------------------------------------------------------------ evidence row
    ax.plot([4, 126], [7.4, 7.4], color=HAIR, lw=0.6, zorder=2)
    rows = [
        (4.0, CELL, "540", "simulated tissues", "\u03c1 = 0.87 with agent-based access"),
        (45.5, GREEN, "114", "cores / 35 patients", "\u03c1 = \u22120.53 against held-out CD8+ cells"),
        (87.0, MAB, "77", "transcriptomic sections", "66/69 positive \u00b7 Visium, ST, Slide-seqV2"),
    ]
    for x0, col, num, lab, extra in rows:
        ax.plot([x0, x0], [1.3, 6.4], color=col, lw=1.6, solid_capstyle="round", zorder=4)
        txt(ax, x0 + 2.6, 4.9, num, 11.0, col, "bold", va="center")
        txt(ax, x0 + 12.2, 4.9, lab, 5.2, INK, va="center")
        txt(ax, x0 + 2.6, 2.0, extra, 4.8, INK2, va="center")

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "graphical_abstract.pdf")
    fig.savefig(OUT / "graphical_abstract.png", dpi=DPI)
    fig.savefig(OUT / "graphical_abstract.tif", dpi=DPI, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("wrote", OUT / "graphical_abstract.pdf")


if __name__ == "__main__":
    main()
