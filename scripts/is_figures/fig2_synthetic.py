"""Fig. 2 -- simulation benchmark with planted barriers and an agent-based ground truth."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

import figstyle as S  # noqa: E402
from isdata import VAL, j  # noqa: E402
import run_38_synthetic_benchmark as B  # noqa: E402
from sparta.barrier import compute_b_cell  # noqa: E402

OUT = Path(__file__).resolve().parents[2] / "results" / "figures" / "is"
GEOM_LABEL = {"closed": "Closed", "gap05": "Gap 5%", "gap10": "Gap 10%", "gap20": "Gap 20%", "gap40": "Gap 40%",
              "closed_far": "Distant\nclosed", "band": "Band", "patches": "Patches", "scattered": "Scatter"}
SUMMARY_LABEL = {"sparta_bcell": "SPARTA B$_{\\rm cell}$ (min cut)",
                 "sparta_field_core": "SPARTA field at core",
                 "nhood_enrichment_z": "Neighbourhood enrichment",
                 "domain_boundary_coverage": "Domain boundary coverage",
                 "ecm_peritumoural_mean": "Peritumoural ECM",
                 "ripley_L_200um": "Ripley's L",
                 "ecm_global_mean": "Global ECM",
                 "vessel_boundary_distance": "Vessel–nest distance"}
GREYS = LinearSegmentedColormap.from_list("g", ["#ffffff", "#e1e0d9", "#c3c2b7", "#898781", "#52514e", "#1f1f1e"])


def fence(ax, xy, edges, half=0.5, lw=0.9):
    segs = []
    for (u, v) in edges:
        mid = xy[[u, v]].mean(axis=0)
        d = xy[v] - xy[u]
        nrm = np.array([-d[1], d[0]]) * half
        segs.append([mid - nrm, mid + nrm])
    ax.add_collection(LineCollection(segs, colors=S.INK, linewidths=lw, zorder=4, capstyle="round"))


def main():
    S.apply()
    sb = j("synthetic_benchmark.json")
    with open(VAL / "synthetic_benchmark_tissues.csv", encoding="utf-8") as f:
        tissues = list(csv.DictReader(f))
    geoms = sb["geoms"]
    fig = plt.figure(figsize=(S.FULL_W, 128 * S.MM))

    # ---------------- a: example tissues ----------------
    xy, A = B.hex_lattice()
    rng = np.random.default_rng(44)
    show = ["closed", "gap10", "closed_far", "band"]
    W = 0.235
    for k, g in enumerate(show):
        ax = fig.add_axes([0.005 + k * 0.248, 0.555, W, 0.40])
        ax.set_aspect("equal")
        ax.axis("off")
        t = B.make_tissue(xy, A, g, rng)
        e, c = B.observe(t, rng, 0.15)
        out = compute_b_cell(A, e, c, t.vessel, t.core, **B.CELL)
        acc, _ = B.access_ground_truth(t, rng, W=2000, T=1500)
        s = 9.5
        col = np.full(len(xy), "#f1f0ec", dtype=object)
        col[t.tumour] = "#c3c2b7"
        col[t.fib] = "#52514e"
        ax.scatter(xy[:, 0], xy[:, 1], s=s, c=list(col), edgecolors="none", zorder=1)
        ax.scatter(xy[t.core, 0], xy[t.core, 1], s=s, marker="s", c=S.INK, edgecolors="none", zorder=2)
        ax.scatter(xy[t.vessel, 0], xy[t.vessel, 1], s=16, c=S.CELL, edgecolors="white", linewidths=0.4, zorder=3)
        fence(ax, xy, out["cut_edges"], half=0.42, lw=0.8)
        ax.set_xlim(-1560, 1560)
        ax.set_ylim(-1620, 1560)
        ax.text(0.5, 1.0, GEOM_LABEL[g].replace("\n", " ") + " capsule" if g in ("closed", "closed_far") else
                ("Capsule with 10% gap" if g == "gap10" else "Peritumoural band"),
                transform=ax.transAxes, ha="center", va="bottom", fontsize=7.2)
        ax.text(0.5, -0.01, f"access {acc:.2f}   B$_{{\\rm cell}}$ {out['b_cell']:.2f}", transform=ax.transAxes,
                ha="center", va="top", fontsize=6.6, color=S.INK2)
        if k == 0:
            S.panel(ax, "a", x=0.0, y=1.0)
    # legend for a
    axl = fig.add_axes([0.0, 0.505, 1.0, 0.03])
    axl.axis("off")
    items = [("o", S.CELL, "vessel (source)"), ("s", S.INK, "tumour core (sink)"), ("o", "#52514e", "matrix-rich spot"),
             ("o", "#c3c2b7", "tumour nest"), ("_", S.INK, "minimum cut")]
    x0 = 0.12
    for m, cc, lab in items:
        if m == "_":
            axl.plot([x0 - 0.006, x0 + 0.006], [0.5, 0.5], color=cc, lw=1.0, transform=axl.transAxes)
        else:
            axl.plot([x0], [0.5], marker=m, color=cc, ms=4.2, transform=axl.transAxes, ls="none")
        axl.text(x0 + 0.012, 0.5, lab, va="center", ha="left", fontsize=6.8, transform=axl.transAxes)
        x0 += 0.165

    # ---------------- b: ground-truth access by geometry ----------------
    axb = fig.add_axes([0.095, 0.085, 0.275, 0.34])
    for i, g in enumerate(geoms):
        v = np.array([float(r["access"]) for r in tissues if r["geom"] == g])
        jit = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.5
        axb.plot(i + jit, v, "o", ms=2.0, mfc="white", mec=S.MUTED, mew=0.5, ls="none", zorder=1)
        axb.plot([i - 0.3, i + 0.3], [v.mean(), v.mean()], color=S.INK, lw=1.4, zorder=2)
    axb.set_xticks(range(len(geoms)))
    axb.set_xticklabels([GEOM_LABEL[g] for g in geoms], fontsize=6.4, rotation=90)
    axb.set_ylabel("Ground-truth access\n(fraction of agents reaching the core)")
    axb.set_ylim(0, 1)
    S.panel(axb, "b", x=-0.30, y=1.03)

    # ---------------- c: Spearman with lost access ----------------
    order = ["sparta_bcell", "sparta_field_core", "nhood_enrichment_z", "domain_boundary_coverage",
             "ecm_peritumoural_mean", "ripley_L_200um", "ecm_global_mean", "vessel_boundary_distance"]
    axc = fig.add_axes([0.525, 0.085, 0.16, 0.34])
    for i, k in enumerate(order):
        pooled = sb["pooled"][k]["spearman_vs_lost_access"]
        colr = S.CELL if k.startswith("sparta") else S.INK2
        for nz, v in sb["per_noise"].items():
            axc.plot([v[k]["spearman_vs_lost_access"]], [i], "|", color=S.MUTED, ms=5, mew=0.8)
        axc.plot([pooled], [i], "o", color=colr, ms=4.4)
    axc.axvline(0, color=S.INK2, lw=0.6)
    axc.set_yticks(range(len(order)))
    axc.set_yticklabels([SUMMARY_LABEL[k] for k in order], fontsize=6.6)
    axc.set_ylim(len(order) - 0.5, -0.5)
    axc.set_xlim(-0.2, 1.0)
    axc.set_xlabel("Spearman ρ with lost access")
    axc.tick_params(axis="y", length=0)
    axc.spines["left"].set_visible(False)
    S.panel(axc, "c", x=-1.25, y=1.03)

    # ---------------- d: AUC matrix ----------------
    cols = [("auc_closed_vs_gap05", "Closed vs\n5% gap"), ("auc_closed_vs_band", "Closed vs\nband"),
            ("auc_closed_vs_scattered", "Closed vs\nscatter"), ("auc_closedfar_vs_scattered", "Distant vs\nscatter")]
    axd = fig.add_axes([0.715, 0.085, 0.27, 0.34])
    M = np.array([[sb["pooled"][k].get(c[0]) if sb["pooled"][k].get(c[0]) is not None else np.nan for c in cols]
                  for k in order])
    axd.imshow(M, cmap=GREYS, vmin=0.5, vmax=1.0, aspect="auto")
    for i in range(M.shape[0]):
        for jj in range(M.shape[1]):
            v = M[i, jj]
            if np.isfinite(v):
                axd.text(jj, i, f"{v:.2f}", ha="center", va="center", fontsize=6.6,
                         color="white" if v > 0.82 else S.INK)
    axd.set_xticks(range(len(cols)))
    axd.set_xticklabels([c[1] for c in cols], fontsize=6.4)
    axd.xaxis.tick_top()
    axd.set_yticks([])
    for sp in axd.spines.values():
        sp.set_visible(False)
    axd.tick_params(length=0)
    axd.text(1.5, len(order) - 0.35, "AUC (0.5 = chance)", ha="center", va="top", fontsize=6.6, color=S.INK2)
    S.panel(axd, "d", x=-0.05, y=1.13)
    S.save(fig, OUT, "Fig2_synthetic")


if __name__ == "__main__":
    main()
