"""Fig. 8 -- what the two transport operators add over simple spatial summaries."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import figstyle as S  # noqa: E402
from isdata import EXTERNAL_ORDER, PRIMARY_ORDER, REPLICATION_ORDER, j  # noqa: E402

OUT = ROOT / "results" / "figures" / "is"
EXAMPLE = "CSCC04"
BASE_LABEL = {"stromal_density": "dens.", "dist_tumour_boundary": "dist.", "niche_z": "niche"}


def example_arrays(sid):
    from run_52_simple_baselines_spearman import compartments, signed_distance
    from sparta.barrier import compute_b_cell_field
    from sparta.io_ import Paths, load_config, load_graph
    from sparta.node_tables import coords_from_nodes, load_nodes, scores_from_nodes
    cfg = load_config(None)
    cfg["paths"]["root"] = str(ROOT)
    P = Paths(cfg)
    A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
    A = sp.csr_matrix(A)
    nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
    Sc = scores_from_nodes(nodes)
    xy = coords_from_nodes(nodes)
    fc = compute_b_cell_field(A, Sc["ecm"], Sc["caf"], source, **cfg["barrier"]["b_cell"])
    bz = np.load(P.barrier(sid), allow_pickle=True)
    reach = np.asarray(fc["reachable"], bool) & np.asarray(bz["reachable"], bool)
    f_cell = np.where(reach, fc["b_cell_field"], np.nan)
    f_mab = np.where(reach, np.asarray(bz["b_mab"], float), np.nan)
    cut = np.zeros(A.shape[0], bool)
    cut[np.asarray(bz["cut_nodes"], int)] = True
    lab, _ = compartments(nodes)
    Ah = A + sp.identity(A.shape[0], format="csr")
    dens = np.asarray(Ah @ (0.5 * (Sc["ecm"] + Sc["caf"]))).ravel() / np.asarray(Ah.sum(axis=1)).ravel()
    dist = signed_distance(xy, lab == 0)
    return f_cell, f_mab, dens, dist, cut


def rank01(v):
    from scipy.stats import rankdata
    return (rankdata(v) - 1) / max(len(v) - 1, 1)


def main():
    S.apply()
    sb = j("simple_baselines_spearman.json")
    per = sb["per_section"]
    syn = j("synthetic_benchmark.json")["pooled"]
    cx = j("codex_validation.json")["summary"]["primary_and_comparators"]

    fig = plt.figure(figsize=(S.FULL_W, 128 * S.MM))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 1.25], wspace=0.62, hspace=0.62,
                          left=0.075, right=0.985, top=0.95, bottom=0.10)

    # ---------------- a, b: one section, spot level ----------------
    f_cell, f_mab, dens, dist, cut = example_arrays(EXAMPLE)
    ok = np.isfinite(f_cell) & np.isfinite(f_mab)
    axa = fig.add_subplot(gs[0, 0])
    y = rank01(f_cell[ok])
    x = dens[ok]
    c = cut[ok]
    axa.scatter(x[~c], y[~c], s=2.2, color=S.LIGHT, linewidths=0, rasterized=True)
    axa.scatter(x[c], y[c], s=4.0, color=S.CELL, linewidths=0, rasterized=True)
    r = spearmanr(dens[ok], f_cell[ok])[0]
    axa.text(1.0, 1.02, f"{EXAMPLE}, ρ = {r:.2f}", transform=axa.transAxes, ha="right", va="bottom", fontsize=6.6)
    axa.set_xlabel("Local stromal density (ECM + CAF)")
    axa.set_ylabel("B$_{\\rm cell}$ field (rank)")
    axa.set_ylim(-0.02, 1.02)
    axa.legend(handles=[Line2D([0], [0], marker="o", color="none", mfc=S.CELL, mec="none", ms=3.2,
                               label="minimum-cut spots"),
                        Line2D([0], [0], marker="o", color="none", mfc=S.LIGHT, mec="none", ms=3.2,
                               label="other spots")],
               loc="upper left", fontsize=6.2, handletextpad=0.1)
    S.panel(axa, "a", x=-0.30, y=1.03)

    axb = fig.add_subplot(gs[0, 1])
    yb = rank01(f_mab[ok])
    axb.scatter(dist[ok] / 1000.0, yb, s=2.2, color=S.MAB, linewidths=0, rasterized=True)
    rb = spearmanr(dist[ok], f_mab[ok])[0]
    axb.axvline(0, color=S.HAIR, lw=0.6, zorder=0)
    axb.text(1.0, 1.02, f"{EXAMPLE}, ρ = {rb:.2f}", transform=axb.transAxes, ha="right", va="bottom", fontsize=6.6)
    axb.set_xlabel("Signed distance to tumour (mm)")
    axb.set_ylabel("B$_{\\rm mAb}$ field (rank)")
    axb.set_ylim(-0.02, 1.02)
    S.panel(axb, "b", x=-0.30, y=1.03)

    # ---------------- c: within-section rho, 30 sections ----------------
    axc = fig.add_subplot(gs[0, 2])
    sids = PRIMARY_ORDER + EXTERNAL_ORDER + REPLICATION_ORDER
    rng = np.random.default_rng(3)
    xt, xl = [], []
    xpos = 0
    for fld, col in (("b_cell_field", S.CELL), ("b_mab_field", S.MAB)):
        for b in ("stromal_density", "dist_tumour_boundary", "niche_z"):
            v = np.array([per[s]["spot"][f"{fld}__{b}"] for s in sids], float)
            for s, val in zip(sids, v):
                axc.plot(xpos + (rng.random() - 0.5) * 0.42, val, S.SHAPE[S.platform_key(s)], ms=2.4,
                         mfc="white", mec=S.MUTED, mew=0.5, ls="none", zorder=2)
            axc.plot([xpos - 0.32, xpos + 0.32], [np.median(v)] * 2, color=col, lw=1.8, zorder=3)
            xt.append(xpos)
            xl.append(BASE_LABEL[b])
            xpos += 1
        xpos += 0.6
    axc.axhline(0, color=S.HAIR, lw=0.6, zorder=0)
    axc.set_xticks(xt)
    axc.set_xticklabels(xl, fontsize=6.4)
    axc.set_ylabel("Within-section Spearman ρ")
    axc.set_ylim(-0.62, 0.78)
    axc.text(1, 0.74, "B$_{\\rm cell}$ field", ha="center", va="top", fontsize=7, color=S.CELL, fontweight="bold")
    axc.text(4.6, 0.74, "B$_{\\rm mAb}$ field", ha="center", va="top", fontsize=7, color=S.MAB, fontweight="bold")
    S.panel(axc, "c", x=-0.22, y=1.03)

    # ---------------- d: section level ----------------
    axd = fig.add_subplot(gs[1, 0])
    tum = [s for s in sids if s != "LN01"]
    xs = np.array([per[s]["section"]["stromal_density"] for s in tum])
    ys = np.array([per[s]["section"]["log_b_rel"] for s in tum])
    for s, xv, yv in zip(tum, xs, ys):
        axd.plot(xv, yv, S.SHAPE[S.platform_key(s)], ms=3.4, mfc=S.CELL, mec="white", mew=0.4, ls="none")
    sl = sb["section_level"]["log_b_rel__stromal_density"]
    axd.text(1.0, 1.02, f"ρ = {sl['rho']:.2f} (p = {sl['p']:.2f}), {sl['n_sections']} sections",
             transform=axd.transAxes, ha="right", va="bottom", fontsize=6.6)
    axd.set_xlabel("Peritumoural stromal density")
    axd.set_ylabel("log B$_{\\rm rel}$ (cut vs matrix-free)")
    S.panel(axd, "d", x=-0.30, y=1.03)

    # ---------------- e: simulation ----------------
    axe = fig.add_subplot(gs[1, 1])
    items = [("sparta_bcell", "SPARTA B$_{\\rm cell}$", S.CELL),
             ("ecm_peritumoural_mean", "stromal density", S.MUTED),
             ("vessel_boundary_distance", "distance to tumour", S.MUTED),
             ("nhood_enrichment_z", "niche enrichment", S.MUTED)]
    for i, (k, lab, col) in enumerate(items):
        auc = syn[k]["auc_closed_vs_gap05"]
        rho = syn[k]["spearman_vs_lost_access"]
        axe.plot(auc, i, "o", ms=4.6, mfc=col, mec=col, ls="none", zorder=3)
        axe.plot(rho, i, "D", ms=3.8, mfc="white", mec=col, mew=0.9, ls="none", zorder=3)
    axe.axvline(0.5, color=S.HAIR, lw=0.6, zorder=0)
    axe.axvline(0.0, color=S.HAIR, lw=0.6, zorder=0)
    axe.set_yticks(range(len(items)))
    axe.set_yticklabels([it[1] for it in items], fontsize=6.6)
    axe.invert_yaxis()
    axe.set_xlim(-0.35, 1.05)
    axe.set_xlabel("Simulation: AUC, closed vs 5% gap (filled)\nor ρ with lost access (open)")
    S.panel(axe, "e", x=-0.62, y=1.03)

    # ---------------- f: CODEX, measured CD8 ----------------
    axf = fig.add_subplot(gs[1, 2])
    ph = j("codex_validation.json")["summary"]["post_hoc"]
    rows = [("log_b_rel", "SPARTA log B$_{\\rm rel}$", None),
            ("peritumoural_matrix", "stromal density", "b_rel_given_peritumoural_matrix"),
            ("vessel_tumour_distance", "vessel–tumour dist.*", "b_rel_given_vessel_tumour_distance"),
            ("core_depth", "tumour-core depth*", "b_rel_given_core_depth"),
            ("contact_enrichment_z", "niche enrichment", "b_rel_given_contact_enrichment_z"),
            (None, "all four*", "b_rel_given_all_simple")]
    for i, (k, lab, pk) in enumerate(rows):
        if k is not None and k in cx:
            col = S.CELL if k == "log_b_rel" else S.MUTED
            axf.plot(cx[k]["ci_core"], [i - 0.12, i - 0.12], color=col, lw=1.1, zorder=2)
            axf.plot(cx[k]["rho_core"], i - 0.12, "o", ms=4.4, mfc=col, mec=col, ls="none", zorder=3)
        if pk is not None and pk in ph:
            axf.plot(ph[pk]["ci"], [i + 0.14, i + 0.14], color=S.CELL, lw=1.1, zorder=2)
            axf.plot(ph[pk]["rho"], i + 0.14, "o", ms=4.4, mfc="white", mec=S.CELL, mew=1.0, ls="none", zorder=3)
    axf.axvline(0, color=S.HAIR, lw=0.6, zorder=0)
    axf.set_yticks(range(len(rows)))
    axf.set_yticklabels([r[1] for r in rows], fontsize=6.4)
    axf.invert_yaxis()
    axf.set_xlim(-0.8, 0.85)
    axf.set_xlabel("CODEX: ρ with CD8$^+$ infiltration ratio\n(cores; bars, 95% CI)")
    axf.legend(handles=[Line2D([0], [0], marker="o", color="none", mfc=S.MUTED, mec=S.MUTED, ms=3.8,
                               label="summary alone"),
                        Line2D([0], [0], marker="o", color="none", mfc="white", mec=S.CELL, mew=1.0, ms=3.8,
                               label="SPARTA given summary")],
               loc="upper right", fontsize=5.6, handletextpad=0.1, borderaxespad=0.1)
    S.panel(axf, "f", x=-0.52, y=1.03)
    S.save(fig, OUT, "Fig8_baselines")


if __name__ == "__main__":
    main()
