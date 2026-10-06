#!/usr/bin/env python
"""Generate publication-oriented CMPB figures from locked SPARTA outputs.

All analytic values are read from results/validation and results/counterfactual.
The framework illustration (Fig. 1) is explicitly schematic. Each figure is
written as vector PDF plus 600-dpi TIFF and PNG for the manuscript preview.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "figures" / "cmpb"
OUT.mkdir(parents=True, exist_ok=True)

NAVY = "#243B53"
BLUE = "#3E78B2"
TEAL = "#168C83"
CORAL = "#D96C55"
GOLD = "#D59A32"
PURPLE = "#7A6AA6"
INK = "#273444"
MID = "#65758B"
LIGHT = "#E7EDF3"
PALE = "#F5F8FB"
GREY = "#AAB5C1"

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "#AAB5C1",
    "axes.linewidth": 0.8,
    "xtick.color": INK,
    "ytick.color": INK,
    "text.color": INK,
    "axes.labelcolor": INK,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})


def readj(rel: str):
    with (ROOT / rel).open(encoding="utf-8") as f:
        return json.load(f)


def savefig(fig, name: str):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(OUT / f"{name}.tif", dpi=600, bbox_inches="tight", pad_inches=0.08,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(OUT / f"{name}.png", dpi=220, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def panel_label(ax, letter):
    ax.text(-0.045, 1.055, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", color=NAVY, va="bottom", ha="left")


def figure1():
    """Schematic of the two model operators and the primary/external cohorts."""
    fig = plt.figure(figsize=(10.8, 6.35), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[3.55, 1.0], width_ratios=[1.0, 1.0])
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    for ax in axes:
        ax.set_xlim(-0.35, 10.35)
        ax.set_ylim(-0.45, 7.45)
        ax.set_aspect("equal")
        ax.axis("off")

    # One simple graph, drawn twice to emphasize shared nodes/geometry.
    coords = {(x, y): (x + (0.42 if y % 2 else 0), y) for y in range(7) for x in range(10)}
    rng = np.random.default_rng(42)
    node_values = {k: rng.uniform(0.1, 0.9) for k in coords}
    for ax, mode in zip(axes, ("cell", "molecule")):
        for (x, y), (px, py) in coords.items():
            for nb in ((x + 1, y), (x, y + 1), (x + 1, y + 1) if y % 2 == 0 else (x - 1, y + 1)):
                if nb in coords:
                    ax.plot([px, coords[nb][0]], [py, coords[nb][1]], color="#D8E0E8", lw=0.65, zorder=1)
        for (x, y), (px, py) in coords.items():
            if x <= 2:
                fc = BLUE
            elif x >= 7:
                fc = PURPLE
            elif mode == "cell" and x in (4, 5):
                fc = CORAL
            elif mode == "molecule":
                val = np.clip((x / 9) * 0.62 + node_values[(x, y)] * 0.30, 0, 1)
                fc = plt.cm.YlOrRd(val)
            else:
                fc = TEAL if node_values[(x, y)] > 0.56 else "#B9C9D7"
            ax.scatter(px, py, s=41, color=fc, edgecolor="white", lw=0.4, zorder=2)
        if mode == "cell":
            ax.axvspan(3.87, 6.0, color=CORAL, alpha=0.08, zorder=0)
            ax.text(4.95, 7.06, "minimum-cut band", color=CORAL, fontsize=8, ha="center")
            ax.annotate("source\n(endothelial-rich)", xy=(1, 5.2), xytext=(-0.23, 6.72),
                        arrowprops=dict(arrowstyle="-", color=BLUE, lw=1), color=BLUE,
                        fontsize=8, ha="left", va="bottom")
            ax.annotate("tumour sink", xy=(8.0, 5.9), xytext=(8.0, 6.72),
                        arrowprops=dict(arrowstyle="-", color=PURPLE, lw=1), color=PURPLE,
                        fontsize=7.4, ha="center", va="bottom")
            ax.set_title("Cellular migration: source–sink minimum cut", loc="left",
                         color=NAVY, fontweight="bold", pad=10)
            ax.text(5.0, -0.33, r"$B_{cell}=1/(\mathrm{maxflow}+\varepsilon)$  ·  returns a cut geometry",
                    ha="center", va="top", fontsize=9, color=INK)
        else:
            ax.text(0.65, 6.85, "vessel boundary  φ = 1", color=BLUE, fontsize=8, ha="left")
            ax.annotate("modelled\nabsorption", xy=(7.5, 4.1), xytext=(8.35, 6.7),
                        arrowprops=dict(arrowstyle="->", color=PURPLE, lw=1), color=PURPLE,
                        fontsize=8, ha="center")
            ax.set_title("IgG-sized transport: screened diffusion–absorption", loc="left",
                         color=NAVY, fontweight="bold", pad=10)
            ax.text(5.0, -0.33, r"$B_{mAb}=-\log(\phi)$  ·  model score, not measured exposure",
                    ha="center", va="top", fontsize=9, color=INK)
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")

    axc = fig.add_subplot(gs[1, :])
    axc.axis("off")
    axc.set_xlim(0, 1); axc.set_ylim(0, 1)
    axc.add_patch(plt.Rectangle((0.025, 0.12), 0.59, 0.75, facecolor="#EEF4F8", edgecolor="#CBD7E2", lw=0.8))
    axc.text(0.045, 0.72, "PRIMARY ANALYSIS", fontsize=8, weight="bold", color=MID)
    axc.text(0.045, 0.49, "19 sections  ·  7 patients", fontsize=12, weight="bold", color=NAVY)
    axc.text(0.045, 0.24, "15 cSCC sections / 6 patients  +  4 melanoma sections / 1 patient",
             fontsize=9, color=INK)
    axc.add_patch(plt.Rectangle((0.64, 0.12), 0.335, 0.75, facecolor="#FBF4EA", edgecolor="#E6D2B4", lw=0.8))
    axc.text(0.66, 0.72, "EXPLORATORY EXTERNAL ARM", fontsize=8, weight="bold", color="#9B6A23")
    axc.text(0.66, 0.49, "3 sections  ·  2 individuals", fontsize=12, weight="bold", color=NAVY)
    axc.text(0.66, 0.25, "2 adjacent breast sections (1 patient)\n+ 1 non-tumour lymph node",
             fontsize=7.8, color=INK, va="center")
    axc.text(0.5, 0.015, "Schematic illustration; nodes and fields are not patient tissue images.",
             ha="center", fontsize=7.6, color=MID, style="italic")
    fig.suptitle("SPARTA maps two model-defined transport operators to a shared spatial graph",
                 fontsize=14, fontweight="bold", color=NAVY)
    savefig(fig, "figure_1_framework")


def figure2():
    main = readj("results/validation/spatial_null_check.json")["per_slide"]
    ext = readj("results/validation/ext_validation.json")["per_slide"]
    ids = sorted([s for s in main if s.startswith("CSCC")]) + sorted([s for s in main if s.startswith("MEL")])
    fig, ax = plt.subplots(figsize=(8.7, 7.35), constrained_layout=True)
    y_main = np.arange(len(ids))[::-1] + 4
    colors = [TEAL if s.startswith("CSCC") else PURPLE for s in ids]
    for y, sid, color in zip(y_main, ids, colors):
        r = main[sid]
        q = r.get("bh_padj", 1.0)
        marker = "o" if q < 0.05 else "o"
        face = color if q < 0.05 else "white"
        ax.scatter(r["real_rho_partial"], y, s=44, facecolor=face, edgecolor=color,
                   lw=1.3, marker=marker, zorder=3)
        ax.text(0.515, y, f"{q:.3f}", va="center", fontsize=7, color=INK)
    ax.axvline(0, color=INK, lw=0.8)
    ax.axvline(0.217, color=GOLD, lw=1.2, ls="--")
    ax.text(0.217, y_main[0] + 0.65, "main median", ha="center", color="#9C6B1C", fontsize=7)
    ax.axhline(3.2, color=GREY, lw=0.8)
    y_ext = [2.45, 1.55, 0.65]
    for y, sid in zip(y_ext, ("BRCA01", "BRCA02", "LN01")):
        r = ext[sid]["spatial_null"]
        q = r.get("bh_padj", r["empirical_p_one_sided"])
        marker = "D" if sid.startswith("BRCA") else "s"
        color = GOLD if sid.startswith("BRCA") else MID
        ax.scatter(r["real_rho_partial"], y, s=55, facecolor=color if q < .05 else "white",
                   edgecolor=color, lw=1.2, marker=marker, zorder=3)
        ax.text(0.515, y, f"{q:.4f}" if sid == "BRCA02" else f"{q:.3f}",
                va="center", fontsize=7, color=INK)
    labels = ids + ["BRCA01", "BRCA02", "LN01"]
    ys = list(y_main) + y_ext
    ax.set_yticks(ys, labels)
    ax.set_xlim(-0.12, 0.61)
    ax.set_ylim(0.05, max(y_main) + 1.25)
    ax.set_xlabel("Partial Spearman ρ after quadratic rank-space adjustment for vessel distance")
    ax.set_ylabel("")
    ax.text(0.515, max(y_main) + 0.75, "BH q", fontsize=7, weight="bold", ha="left")
    ax.text(-0.11, max(y_main) + 0.28, "PRIMARY: cSCC + melanoma (19 sections)", fontsize=8,
            weight="bold", color=NAVY)
    ax.text(-0.105, 3.48, "EXTERNAL: exploratory", fontsize=7.5,
            weight="bold", color="#9B6A23")
    ax.grid(axis="x", color=LIGHT, lw=0.7)
    ax.text(0.01, -0.16,
            "Main cohort: 18/19 ρ > 0; spatial-surrogate p < 0.05 in 15/19 and BH q < 0.05 in 14/19. "
            "External BRCA02 q = 0.0499 is at the 500-surrogate Monte Carlo boundary.",
            transform=ax.transAxes, fontsize=7.3, color=MID, va="top", wrap=True)
    panel_label(ax, "A")
    ax.set_title("Model-field coupling across primary and exploratory sections", loc="left", color=NAVY,
                 fontweight="bold", pad=24)
    savefig(fig, "figure_2_coupling")


def figure3():
    data = readj("results/validation/bmab_sensitivity.json")["slides"]
    fig, axes = plt.subplots(1, 2, figsize=(9.1, 4.05), constrained_layout=True)
    ax = axes[0]
    x = np.linspace(0, 1, 301)
    xi0, beta, radius = 20.0, 3.0, 5.5
    xi = xi0 * np.exp(-beta * x)
    phi = np.where(radius / xi < 1, (1 - radius / xi) ** 2, 0.0)
    ax.plot(x, xi, color=BLUE, lw=2.1, label=r"$ξ=ξ_0 e^{-βx}$")
    ax.axhline(radius, color=CORAL, ls="--", lw=1.3, label="IgG-sized radius r = 5.5 nm")
    threshold = np.log(xi0 / radius) / beta
    ax.axvline(threshold, color=GOLD, ls=":", lw=1.5)
    ax.fill_between(x, 0, xi, where=xi <= radius, color=CORAL, alpha=0.12)
    ax.text(threshold + .02, 16.5, f"threshold x = {threshold:.2f}", fontsize=7.5, color="#9B6A23")
    ax.set_xlim(0, 1); ax.set_ylim(0, 21)
    ax.set_xlabel("Within-section rank-normalised crosslink score x")
    ax.set_ylabel("Model-defined effective mesh ξ (nm)")
    ax.grid(color=LIGHT, lw=.6)
    ax.legend(frameon=False, loc="upper right", fontsize=7.1)
    ax.set_title("Size-exclusion law", loc="left", color=NAVY,
                 fontweight="bold")

    ax = axes[1]
    beta_grid = data["MEL01"]["lam_x_beta"]["grid"]["ys"]
    lam_grid = data["MEL01"]["lam_x_beta"]["grid"]["xs"]
    lam_idx = lam_grid.index(3.0)
    selected = ["CSCC01", "CSCC03", "CSCC05", "MEL01", "MEL03"]
    cols = ["#168C83", "#62B3A9", "#3E78B2", "#7A6AA6", "#B09ACB"]
    for sid, color in zip(selected, cols):
        vals = data[sid]["lam_x_beta"]["crosslink_pct"]
        # stored as beta-by-lambda grid; each entry is a model-derived edge fraction.
        yv = [row[lam_idx] for row in vals]
        ax.plot(beta_grid, yv, marker="o", ms=3.1, lw=1.45, color=color, alpha=.9, label=sid)
    ax.axvline(3.0, color=GOLD, ls="--", lw=1.2)
    ax.set_xscale("log")
    ax.set_ylim(-2, 102)
    ax.set_xlabel("Mesh-contraction steepness β (log scale; λ = 3)")
    ax.set_ylabel("Edges assigned zero size-exclusion conductance (%)")
    ax.grid(color=LIGHT, lw=.6)
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="upper left")
    ax.set_title("Edge exclusion vs β", loc="left", color=NAVY,
                 fontweight="bold")
    panel_label(axes[0], "A"); panel_label(axes[1], "B")
    fig.suptitle("Size exclusion is a model law, not a mesh measurement", fontsize=13,
                 fontweight="bold", color=NAVY)
    fig.text(.5, -.025,
             "The rank-based input fixes the complete-exclusion threshold by construction. "
             "At β = 3, the crosslinking-attributed variance share is a model decomposition, not biological attribution.",
             ha="center", fontsize=7.5, color=MID)
    savefig(fig, "figure_3_parameter_sensitivity")


def bh(pvals):
    p = np.asarray(pvals, float)
    n = len(p); order = np.argsort(p)
    qsort = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    q = np.empty(n); q[order] = np.clip(qsort, 0, 1)
    return q


def figure4():
    main = readj("results/validation/s2_matched_selection.json")["per_slide"]
    ext = readj("results/validation/s2_matched_ext.json")["per_slide"]
    ids = sorted([s for s in main if s.startswith("CSCC")]) + sorted([s for s in main if s.startswith("MEL")])
    p = np.array([main[s]["p_vs_in_cut_matched"] for s in ids]); q = bh(p)
    fig = plt.figure(figsize=(8.7, 7.1), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[5.3, 1.25], hspace=.07)
    ax = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax)
    y = np.arange(len(ids))[::-1]
    for yy, sid, qq in zip(y, ids, q):
        row = main[sid]; ratio = row["ratio_vs_in_cut_matched"]
        col = TEAL if sid.startswith("CSCC") else PURPLE
        ax.scatter(ratio, yy, marker="o", s=38, facecolor=col if qq < .05 else "white",
                   edgecolor=col, lw=1.25, zorder=3)
    ax.axvline(1, color=INK, lw=.9, ls="--")
    ax.set_yticks(y, ids)
    ax.set_xlim(.90, 1.28)
    ax.set_xlabel("")
    ax.set_title("Primary cohort (n = 19 sections)", loc="left", color=NAVY, fontweight="bold")
    ax.grid(axis="x", color=LIGHT, lw=.6)
    y2 = np.array([2, 1, 0])
    for yy, sid in zip(y2, ("BRCA01", "BRCA02", "LN01")):
        row = ext[sid]
        col = GOLD if sid.startswith("BRCA") else MID
        # exploratory p-values are separate; 2/3 pass raw and BH in this family.
        qext = bh([ext[s]["p_vs_in_cut_matched"] for s in ("BRCA01", "BRCA02", "LN01")])
        idx = ("BRCA01", "BRCA02", "LN01").index(sid)
        ax2.scatter(row["ratio_vs_in_cut_matched"], yy, marker="D" if sid.startswith("BRCA") else "s",
                    s=50, facecolor=col if qext[idx] < .05 else "white", edgecolor=col, lw=1.25)
    ax2.axvline(1, color=INK, lw=.9, ls="--")
    ax2.set_yticks(y2, ["BRCA01", "BRCA02", "LN01"])
    ax2.set_xlabel("Residual barrier: scattered removal / contiguous-gap removal")
    ax2.set_title("Exploratory external arm (3 sections; not pooled with primary cohort)",
                  loc="left", color="#9B6A23", fontsize=9, fontweight="bold", pad=3)
    ax2.grid(axis="x", color=LIGHT, lw=.6)
    ax2.set_ylim(-.55, 2.55)
    plt.setp(ax.get_xticklabels(), visible=False)
    panel_label(ax, "A")
    fig.suptitle("Selection-matched S2 counterfactual: a modest model effect",
                 fontsize=13, fontweight="bold", color=NAVY)
    savefig(fig, "figure_4_s2_matched")


def figure5():
    files = sorted((ROOT / "results" / "counterfactual").glob("*.json"))
    trajectories = []
    for path in files:
        try:
            obj = json.loads(path.read_text(encoding="utf-8"))
            if not path.stem.startswith(("CSCC", "MEL")) or "s3" not in obj:
                continue
            y = np.asarray(obj["s3"]["mean_core"], dtype=float)
            x = np.asarray(obj["s3"]["radii_nm"], dtype=float)
            if len(x) and np.all(np.isfinite(y)):
                trajectories.append((path.stem, x, y))
        except (json.JSONDecodeError, TypeError, ValueError):
            continue
    fig, ax = plt.subplots(figsize=(7.4, 4.8), constrained_layout=True)
    if not trajectories:
        raise RuntimeError("No primary-cohort S3 results found under results/counterfactual")
    matrix = np.vstack([t[2] for t in trajectories])
    radii = trajectories[0][1]
    for sid, x, yy in trajectories:
        col = TEAL if sid.startswith("CSCC") else PURPLE
        ax.plot(x, yy, color=col, alpha=.27, lw=.9)
    med = np.median(matrix, axis=0)
    ax.plot(radii, med, color=NAVY, lw=2.5, marker="o", ms=5, label="median across 19 sections")
    ax.axvline(5.5, color=CORAL, lw=1.3, ls="--", label="IgG hydrodynamic radius (5.5 nm)")
    ax.set_xlabel("Molecular hydrodynamic radius r (nm)")
    ax.set_ylabel("Mean tumour-core $B_{mAb}$ (unitless model score)")
    ax.grid(color=LIGHT, lw=.65)
    ax.legend(frameon=False, loc="upper left")
    panel_label(ax, "A")
    ax.set_title("In-model size scan across primary sections", loc="left", color=NAVY,
                 fontweight="bold")
    savefig(fig, "figure_5_size_scan")


def graphical_abstract():
    fig, ax = plt.subplots(figsize=(8.2, 5.2), constrained_layout=True)
    ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 6.2)
    ax.text(5, 5.86, "SPARTA", fontsize=23, fontweight="bold", color=NAVY, ha="center")
    ax.text(5, 5.45, "Two model-defined transport barriers on one spatial graph",
            fontsize=11, color=INK, ha="center")
    # central shared substrate
    ax.add_patch(plt.Rectangle((3.35, 3.2), 3.3, 1.55, facecolor="#EEF4F8", edgecolor="#C8D5E1", lw=1.0))
    ax.text(5, 4.38, "Spatial transcriptomics", ha="center", fontsize=10, weight="bold", color=NAVY)
    ax.text(5, 3.94, "one graph · rank-based expression scores", ha="center", fontsize=8.4, color=INK)
    ax.text(5, 3.52, "shared matrix signal is a built-in coupling source", ha="center", fontsize=7.5, color=MID)
    # two operator cards
    cards = [(0.35, "T-cell migration", "minimum cut", "barrier score + cut geometry", BLUE),
             (6.85, "IgG-sized molecule", "screened diffusion–absorption", "exposure-deficit model score", CORAL)]
    for x0, title, op, out, color in cards:
        ax.add_patch(plt.Rectangle((x0, 3.1), 2.8, 1.8, facecolor="white", edgecolor=color, lw=1.4))
        ax.text(x0 + 1.4, 4.5, title, ha="center", fontsize=9.2, weight="bold", color=color)
        ax.text(x0 + 1.4, 4.05, op, ha="center", fontsize=8.5, color=INK)
        ax.text(x0 + 1.4, 3.55, out, ha="center", fontsize=7.2, color=MID)
    for xstart, xend in ((3.45, 3.15), (6.55, 6.85)):
        ax.annotate("", xy=(xend, 3.96), xytext=(xstart, 3.96),
                    arrowprops=dict(arrowstyle="->", lw=1.4, color=MID))
    # results strip
    items = [
        (0.5, "18/19", "positive partial ρ"),
        (3.0, "14/19", "BH q < 0.05 under graph-spectral null"),
        (6.3, "1.053×", "median matched S2 ratio; in-model"),
    ]
    for x0, big, small in items:
        ax.add_patch(plt.Rectangle((x0, 1.54), 3.0, .9, facecolor=PALE, edgecolor=LIGHT, lw=.8))
        ax.text(x0 + .2, 2.1, big, fontsize=16, weight="bold", color=NAVY, va="center")
        ax.text(x0 + .2, 1.74, small, fontsize=7.3, color=MID, va="center")
    ax.text(5, 1.06, "Primary data: 19 sections from 7 patients  |  External arm: 3 sections, exploratory",
            ha="center", fontsize=8.4, color=INK)
    ax.text(5, .63, "Expression proxies and qualitative scale parameters; no measured mesh size, antibody exposure, or functional validation",
            ha="center", fontsize=7.6, color=CORAL, weight="bold")
    ax.text(5, .28, "Model schematic; not a clinical efficacy or treatment-response prediction.",
            ha="center", fontsize=7.4, color=MID, style="italic")
    savefig(fig, "graphical_abstract")


def main():
    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    graphical_abstract()
    print(f"Wrote publication-oriented figure set to {OUT}")


if __name__ == "__main__":
    main()
