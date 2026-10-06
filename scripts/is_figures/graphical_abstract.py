"""SPARTA graphical abstract -- conceptual, mechanism-first, minimal text.

Graphical-abstract conventions applied here:
  * it shows the LOGIC of the method, not results -- no counts, no correlations, no
    formulas, no operator symbols, no negative scope statements;
  * words are kept to a handful of short labels; the diagrams carry the meaning;
  * one visual grammar throughout: grey = matrix, blue = cell-scale transport,
    orange = molecular-scale transport.

Canvas 13.0 x 5.2 cm (w:h = 2.5:1), the size at which the abstract is displayed, so every
label is set at its final size.

Output: results/figures/graphical_abstract.{pdf,png,tif}
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "is_figures"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import Circle, FancyBboxPatch  # noqa: E402

import figstyle as S  # noqa: E402

OUT = ROOT / "results" / "figures"
W_MM, H_MM, DPI = 130.0, 52.0, 500

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8b959e"
HAIR, BG, PANEL = "#dfe6ec", "#f5f8fb", "#ffffff"
CELL, MAB = S.CELL, S.MAB
STROMA, MATRIX, LIGHT = "#b9c2ca", "#5d6874", "#d7dee4"
MESH = "#ccd5dd"


def txt(ax, x, y, s, size=5.4, color=INK2, weight="normal", ha="left", va="baseline", z=9):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha, va=va, zorder=z)


def card(ax, x, y, w, h, r=1.6, ec=HAIR, fc=PANEL, lw=0.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                fc=fc, ec=ec, lw=lw, zorder=1))


def dot(ax, x, y, r, fc, ec="white", lw=0.35, z=4):
    ax.add_patch(Circle((x, y), r, fc=fc, ec=ec, lw=lw, zorder=z))


def blend(c1, c2, t):
    a = np.array([int(c1[i:i + 2], 16) for i in (1, 3, 5)])
    b = np.array([int(c2[i:i + 2], 16) for i in (1, 3, 5)])
    return "#%02x%02x%02x" % tuple((a + (b - a) * t).astype(int))


def tissue(ax, cx, cy, half):
    """A schematic section on a spot lattice: stroma, dense matrix capsule, tumour core, vessels."""
    dx, dy = 1.0, np.sqrt(3) / 2
    pts = [(i * dx + (0.5 * dx if j % 2 else 0.0), j * dy)
           for j in range(-9, 10) for i in range(-10, 11)]
    P = np.array(pts, float)
    th = np.arctan2(P[:, 1], P[:, 0]); r = np.linalg.norm(P, axis=1)
    keep = r <= 7.2 + 0.65 * np.sin(3 * th) + 0.45 * np.cos(5 * th + 1.0)
    P, r = P[keep], r[keep]
    sc = half / 7.8
    X, Y = cx + P[:, 0] * sc, cy + P[:, 1] * sc
    segs = [[(X[i], Y[i]), (X[j], Y[j])] for i in range(len(X)) for j in range(i + 1, len(X))
            if (P[i, 0] - P[j, 0]) ** 2 + (P[i, 1] - P[j, 1]) ** 2 <= (1.06 * dx) ** 2]
    ax.add_collection(LineCollection(segs, colors=MESH, linewidths=0.32, zorder=2))
    order = np.argsort(-r)
    vessels = set(order[:9])
    for i in range(len(X)):
        if i in vessels:
            col, rad = CELL, 0.72
        elif r[i] < 2.3:
            col, rad = MAB, 0.62
        elif 3.2 < r[i] < 4.3:
            col, rad = MATRIX, 0.56
        else:
            col, rad = STROMA, 0.50
        dot(ax, X[i], Y[i], rad, col)


def channel(ax, x0, y0, blocked):
    """A 5 x 3 lattice: source left, dense matrix in the middle, sink right."""
    dx, dy = 4.4, 3.4
    pts = [(x0 + i * dx, y0 + j * dy) for j in range(3) for i in range(5)]
    segs = [[p, q] for i, p in enumerate(pts) for q in pts[i + 1:]
            if abs(p[0] - q[0]) <= dx + 1e-9 and abs(p[1] - q[1]) <= dy + 1e-9]
    ax.add_collection(LineCollection(segs, colors=MESH, linewidths=0.32, zorder=2))
    for (x, y) in pts:
        col = 0 if x < x0 + dx else (4 if x > x0 + 3 * dx else 2)
        if blocked:
            c = CELL if col == 0 else (MAB if col == 4 else MATRIX)
        else:
            c = blend("#bcd8f6", MAB, 0.12 + 0.88 * (x - x0) / (4 * dx))
        dot(ax, x, y, 0.92, c)
    mid, ymid = x0 + 1.5 * dx, y0 + dy
    if blocked:
        ax.plot([mid, mid], [y0 - 1.4, y0 + 2 * dy + 1.4], color=INK, lw=0.85,
                ls=(0, (1.5, 1.3)), zorder=6)
        ax.annotate("", xy=(mid - 0.8, ymid), xytext=(x0 + 0.9, ymid),
                    arrowprops=dict(arrowstyle="-|>", color=CELL, lw=1.0,
                                    shrinkA=0, shrinkB=0, mutation_scale=7), zorder=6)
    else:
        ax.annotate("", xy=(x0 + 4 * dx - 0.3, ymid), xytext=(x0 + 0.9, ymid),
                    arrowprops=dict(arrowstyle="-|>", color=MAB, lw=1.0,
                                    shrinkA=0, shrinkB=0, mutation_scale=7), zorder=6)
        dot(ax, x0 + 4 * dx, ymid, 1.9, "none", ec=MAB, lw=0.85, z=7)


def mini_graph(ax, x0, ys, colours, dx=5.2, r=0.95):
    """A 4 x 3 lattice of the same node values, drawn at one of two spatial arrangements."""
    pts = [(x0 + i * dx, y) for y in ys for i in range(4)]
    segs = [[p, q] for i, p in enumerate(pts) for q in pts[i + 1:]
            if abs(p[0] - q[0]) <= dx + 1e-9 and abs(p[1] - q[1]) <= 4.1 + 1e-9]
    ax.add_collection(LineCollection(segs, colors=MESH, linewidths=0.32, zorder=2))
    for (x, y), c in zip(pts, colours):
        dot(ax, x, y, r, c)


def main():
    S.apply()
    fig = plt.figure(figsize=(W_MM / 25.4, H_MM / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_MM); ax.set_ylim(0, H_MM); ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(plt.Rectangle((0, 0), W_MM, H_MM, fc=BG, ec="none", zorder=0))

    txt(ax, 4, 46.4, "One tissue graph, two transport barriers", 9.2, INK, "bold")
    ax.plot([4, 126], [43.4, 43.4], color=HAIR, lw=0.6, zorder=2)

    # 1 ---------------------------------------------------- the tissue graph
    card(ax, 4, 8.6, 36.0, 32.2)
    tissue(ax, 22.0, 24.0, 12.6)
    txt(ax, 5.2, 38.6, "tissue graph", 5.4, INK, "bold")

    # 2 ---------------------------------------------------- the two barriers
    card(ax, 43.0, 8.6, 46.0, 32.2)
    txt(ax, 44.2, 38.6, "two barriers", 5.4, INK, "bold")
    card(ax, 44.6, 23.8, 42.8, 12.4, r=1.2, ec="#cfe1f5")
    channel(ax, 47.6, 26.4, blocked=True)
    txt(ax, 86.2, 29.8, "cells", 5.0, CELL, "bold", ha="right", va="center")
    card(ax, 44.6, 10.6, 42.8, 11.6, r=1.2, ec="#f8d8c9")
    channel(ax, 47.6, 13.2, blocked=False)
    txt(ax, 86.2, 16.4, "molecules", 5.0, MAB, "bold", ha="right", va="center")

    # 3 ---------------------------------------------------- the null models
    card(ax, 92.0, 8.6, 34.0, 32.2)
    txt(ax, 93.2, 38.6, "null models", 5.4, INK, "bold")
    # the same node values in their observed positions, then in permuted positions
    greys = ["#5d6874", "#8b959e", "#b9c2ca", "#d7dee4"]
    observed = greys * 3
    permuted = [greys[i] for i in (2, 0, 3, 1, 1, 3, 0, 2, 3, 1, 2, 0)]
    mini_graph(ax, 105.4, (35.0, 31.0, 27.0), observed)
    txt(ax, 93.6, 31.0, "observed", 5.0, MUTED, va="center")
    mini_graph(ax, 105.4, (19.0, 15.0, 11.0), permuted)
    txt(ax, 93.6, 15.0, "resampled", 5.0, MUTED, va="center")
    ax.annotate("", xy=(119.0, 22.6), xytext=(106.2, 22.6),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.9,
                                connectionstyle="arc3,rad=-0.34", shrinkA=0, shrinkB=0,
                                mutation_scale=8), zorder=6)

    # left-to-right flow between the three panels
    for xa in (40.6, 89.6):
        ax.annotate("", xy=(xa + 1.9, 24.0), xytext=(xa + 0.2, 24.0),
                    arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.9,
                                    shrinkA=0, shrinkB=0, mutation_scale=8), zorder=5)

    txt(ax, 65, 3.6, "Tissue structure separated from model construction", 5.6, INK, ha="center")

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "graphical_abstract.pdf")
    fig.savefig(OUT / "graphical_abstract.png", dpi=DPI)
    fig.savefig(OUT / "graphical_abstract.tif", dpi=DPI, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("wrote", OUT / "graphical_abstract.pdf")


if __name__ == "__main__":
    main()
