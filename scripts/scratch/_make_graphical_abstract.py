# -*- coding: utf-8 -*-
"""Graphical abstract — vertical portrait, matplotlib, zero-text-error.
Layout (top->bottom): title / spatial graph / B_cell box / B_mAb box /
one-substrate note / findings list. Output: results/figures/graphical_abstract_portrait.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

BLUE = "#2c7fb8"
ORANGE = "#d95f0e"
GREY = "#b8b8b8"
INK = "#222222"

fig = plt.figure(figsize=(4.6, 11.2), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

# ---- Title (italic) ----
ax.text(0.5, 0.955, "Two transport operators, one substrate:",
        ha="center", va="top", fontsize=13, style="italic", color=INK)
ax.text(0.5, 0.928, "graph-transport modelling of cell and antibody\ndelivery barriers",
        ha="center", va="top", fontsize=11, style="italic", color=INK, linespacing=1.3)

# ---- Spatial graph ----
ax.text(0.36, 0.855, "spatial graph", ha="center", va="center",
        fontsize=12, fontweight="bold", color=INK)

rng = np.random.default_rng(7)
# immune-entry (blue) cluster on the left
blue = np.column_stack([0.10 + 0.10*rng.random(5), 0.66 + 0.07*rng.random(5)])
# tumour core (orange) on the right
orange = np.column_stack([0.60 + 0.10*rng.random(5), 0.67 + 0.06*rng.random(5)])
# matrix (grey) scattered in between
grey = np.column_stack([0.22 + 0.30*rng.random(16), 0.64 + 0.10*rng.random(16)])

# edges: connect grey-grey neighbours, blue-grey, grey-orange
pts = np.vstack([blue, grey, orange])
for i in range(len(pts)):
    for j in range(i+1, len(pts)):
        d = np.hypot(pts[i,0]-pts[j,0], pts[i,1]-pts[j,1])
        if d < 0.16:
            ax.plot([pts[i,0], pts[j,0]], [pts[i,1], pts[j,1]],
                    color=GREY, lw=0.7, zorder=1, alpha=0.6)

ax.scatter(blue[:,0], blue[:,1], s=70, color=BLUE, zorder=3, edgecolors="white", linewidths=0.6)
ax.scatter(grey[:,0], grey[:,1], s=55, color="#cfcfcf", zorder=2, edgecolors="white", linewidths=0.5)
ax.scatter(orange[:,0], orange[:,1], s=70, color=ORANGE, zorder=3, edgecolors="white", linewidths=0.6)

# ---- B_cell box (blue) ----
ax.add_patch(FancyBboxPatch((0.26, 0.47), 0.40, 0.11,
             boxstyle="round,pad=0.012,rounding_size=0.018",
             linewidth=2.2, edgecolor=BLUE, facecolor="white"))
ax.text(0.46, 0.552, r"$B_{cell}$", ha="center", va="center",
        fontsize=15, color=BLUE, style="italic")
ax.text(0.46, 0.505, "min-cut barrier", ha="center", va="center",
        fontsize=10.5, color=INK)
# arrow from graph down to B_cell (left side)
ax.add_patch(FancyArrowPatch((0.28, 0.66), (0.36, 0.585),
             arrowstyle="-|>", mutation_scale=14, lw=1.8, color=BLUE))

# ---- B_mAb box (orange) ----
ax.add_patch(FancyBboxPatch((0.26, 0.31), 0.40, 0.11,
             boxstyle="round,pad=0.012,rounding_size=0.018",
             linewidth=2.2, edgecolor=ORANGE, facecolor="white"))
ax.text(0.46, 0.392, r"$B_{mAb}$", ha="center", va="center",
        fontsize=15, color=ORANGE, style="italic")
ax.text(0.46, 0.345, "screened Poisson, size exclusion", ha="center", va="center",
        fontsize=10.5, color=INK)
ax.add_patch(FancyArrowPatch((0.74, 0.66), (0.58, 0.425),
             arrowstyle="-|>", mutation_scale=14, lw=1.8, color=ORANGE))

# ---- one substrate ----
ax.text(0.46, 0.255, "one substrate: extracellular matrix",
        ha="center", va="center", fontsize=10.5, color="#555555", style="italic")

# ---- findings ----
ax.text(0.20, 0.185, "findings", ha="left", va="center",
        fontsize=12, fontweight="bold", color=INK)
rows = [
    (ORANGE, "percolation-limited"),
    (BLUE,   "barriers co-located"),
    ("#555555", "deterministic, < 0.2 s"),
]
y = 0.145
for c, lab in rows:
    ax.scatter([0.16], [y], s=45, color=c, zorder=3)
    ax.text(0.20, y, lab, ha="left", va="center", fontsize=10.5, color=INK)
    y -= 0.045

out = r"D:\sparta\results\figures\graphical_abstract_portrait.png"
fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.12)
from PIL import Image
print("saved", out, Image.open(out).size)
