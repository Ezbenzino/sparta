"""Fig. 2 -- simulation benchmark with planted barriers and process-level ground truths."""
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
RADII = [(0.5, "0.5", "o", "#f7c9a8"), (2.0, "2", "s", S.MAB), (5.5, "5.5", "^", "#9c3413")]
RHO_CMAP = LinearSegmentedColormap.from_list("r", [S.INK2, S.WASH_GREY, S.MAB])
MAB_ROWS = [("b_mab_core", "$B_{\\rm mAb}$ at the core"),
            ("b_mab_reach_mean", "$B_{\\rm mAb}$, mean reachable"),
            ("b_cell_field_core", "$B_{\\rm cell}$ field at the core"),
            ("ecm_peritumoural_mean", "Peritumoural ECM"),
            ("ecm_global_mean", "Global ECM"),
            ("vessel_boundary_distance", "Vessel–nest distance")]


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
    mg = j("mab_ground_truth.json")
    with open(VAL / "synthetic_benchmark_tissues.csv", encoding="utf-8") as f:
        tissues = list(csv.DictReader(f))
    with open(VAL / "mab_ground_truth_tissues.csv", encoding="utf-8") as f:
        mtis = list(csv.DictReader(f))
    geoms = sb["geoms"]
    fig = plt.figure(figsize=(S.FULL_W, 190 * S.MM))

    # ---------------- a: example tissues ----------------
    xy, A = B.hex_lattice()
    rng = np.random.default_rng(44)
    show = ["closed", "gap10", "closed_far", "band"]
    W = 0.235
    for k, g in enumerate(show):
        ax = fig.add_axes([0.005 + k * 0.248, 0.700, W, 0.27])
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
    axl = fig.add_axes([0.0, 0.660, 1.0, 0.03])
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
    axb = fig.add_axes([0.095, 0.355, 0.275, 0.235])
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
    axc = fig.add_axes([0.525, 0.355, 0.16, 0.235])
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
    axd = fig.add_axes([0.715, 0.355, 0.27, 0.235])
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

    # ---------------- e: ground-truth particle delivery by geometry and radius ----------------
    axe = fig.add_axes([0.095, 0.045, 0.275, 0.235])
    T_fin = str(mg["checkpoints"][-1])
    for kr, (r, rl, mk, cc) in enumerate(RADII):
        off = (kr - 1) * 0.27
        for i, g in enumerate(geoms):
            v = np.array([float(row["delivery"]) for row in mtis
                          if row["geom"] == g and float(row["radius"]) == r and row["T"] == T_fin])
            jit = (np.random.default_rng(100 + kr * 9 + i).random(len(v)) - 0.5) * 0.12
            axe.plot(i + off + jit, v, mk, ms=2.2, mfc=cc, mec=cc, ls="none", zorder=1)
            m = mg["delivery_by_radius"][rl]["by_geom"][g]
            axe.plot([i + off - 0.11, i + off + 0.11], [m, m], color=S.INK, lw=1.1, zorder=2)
    axe.set_xticks(range(len(geoms)))
    axe.set_xticklabels([GEOM_LABEL[g] for g in geoms], fontsize=6.4, rotation=90)
    axe.set_ylabel("Ground-truth core delivery\n(particles at 2,500 steps)")
    axe.set_ylim(-0.002, None)
    hnd = [plt.Line2D([], [], marker=mk, ls="none", color=cc, ms=4, label=f"{r:g} nm")
           for (r, rl, mk, cc) in RADII]
    axe.legend(handles=hnd, loc="upper right", handletextpad=0.2, borderpad=0.2, fontsize=6.6)
    S.panel(axe, "e", x=-0.30, y=1.03)

    # ---------------- f: Spearman with lost delivery by radius ----------------
    axf = fig.add_axes([0.50, 0.045, 0.135, 0.235])
    R = np.array([[mg["per_radius"][rl][k]["spearman_vs_lost_delivery"] for (_, rl, _, _) in RADII]
                  for (k, _) in MAB_ROWS])
    axf.imshow(R, cmap=RHO_CMAP, vmin=-0.7, vmax=0.7, aspect="auto")
    for i in range(R.shape[0]):
        for jj in range(R.shape[1]):
            v = R[i, jj]
            axf.text(jj, i, f"{v:.2f}", ha="center", va="center", fontsize=6.6,
                     color="white" if abs(v) > 0.5 else S.INK)
    axf.set_xticks(range(len(RADII)))
    axf.set_xticklabels([f"{r:g}" for (r, _, _, _) in RADII], fontsize=6.4)
    axf.xaxis.tick_top()
    axf.set_yticks(range(len(MAB_ROWS)))
    axf.set_yticklabels([lab for (_, lab) in MAB_ROWS], fontsize=6.6)
    axf.tick_params(axis="y", length=0)
    axf.tick_params(axis="x", length=0)
    for sp in axf.spines.values():
        sp.set_visible(False)
    axf.text(1.0, len(MAB_ROWS) - 0.35, "Spearman ρ with lost delivery", ha="center", va="top",
             fontsize=6.6, color=S.INK2)
    axf.text(1.0, len(MAB_ROWS) + 0.75, "probe radius (nm)", ha="center", va="bottom",
             fontsize=6.6, color=S.INK2)
    S.panel(axf, "f", x=-0.42, y=1.13)

    # ---------------- g: census-time dependence ----------------
    axg = fig.add_axes([0.72, 0.045, 0.245, 0.235])
    cps = mg["checkpoints"]
    for (r, rl, mk, cc) in RADII:
        v = [mg["by_checkpoint"][str(T)][rl]["b_mab_core"] for T in cps]
        axg.plot(range(len(cps)), v, marker=mk, ms=3.6, lw=1.1, color=cc, label=f"{r:g} nm")
    axg.axhline(0, color=S.INK2, lw=0.6)
    axg.set_xticks(range(len(cps)))
    axg.set_xticklabels([f"{T:,}" for T in cps])
    axg.set_xlabel("Census (steps)")
    axg.set_ylabel("ρ($B_{\\rm mAb}$ core, lost delivery)")
    axg.set_ylim(-0.35, 1.0)
    axg.legend(loc="lower right", handletextpad=0.3, borderpad=0.2, fontsize=6.6)
    S.panel(axg, "g", x=-0.20, y=1.03)
    S.save(fig, OUT, "Fig2_synthetic")


if __name__ == "__main__":
    main()
