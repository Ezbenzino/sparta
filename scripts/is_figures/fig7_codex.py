"""Fig. 7 -- the cellular barrier against measured CD8+ T-cell positions (CODEX, colorectal cancer)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy.stats import rankdata  # noqa: E402

import figstyle as S  # noqa: E402
from isdata import j  # noqa: E402

OUT = ROOT / "results" / "figures" / "is"
CSV = ROOT / "data" / "external" / "codex_crc" / "CRC_clusters_neighborhoods_markers.csv"


def pick_examples(prim):
    """Deterministic: cores closest (in rank space) to high-barrier/excluded and low-barrier/infiltrated."""
    rb = rankdata(prim["log_b_rel"]) / len(prim)
    ri = rankdata(prim["ir"]) / len(prim)
    big = prim["n_lcc"] >= prim["n_lcc"].median()
    d_hi = np.where(big, (rb - 0.85) ** 2 + (ri - 0.15) ** 2, np.inf)
    d_lo = np.where(big, (rb - 0.15) ** 2 + (ri - 0.85) ** 2, np.inf)
    return prim.iloc[int(np.argmin(d_hi))]["core"], prim.iloc[int(np.argmin(d_lo))]["core"]


def core_geometry(core):
    import importlib.util
    spec = importlib.util.spec_from_file_location("r47", ROOT / "scripts" / "run_47_codex_validation.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    from sparta.barrier import compute_b_cell
    from sparta.cellgraph import core_sinks, delaunay_graph, largest_component, rank01

    df = pd.read_csv(CSV, usecols=list(m.COL.values()))
    sub = df[df[m.COL["spot"]] == core]
    xy = sub[[m.COL["x"], m.COL["y"]]].to_numpy(float) * m.PX_UM
    lab = sub[m.COL["lab"]].astype(str).to_numpy()
    ok = lab != "dirt"
    scaf = np.isin(lab, m.SCAFFOLD) & ok
    Sx = np.flatnonzero(scaf)
    A, _ = delaunay_graph(xy[Sx], 50.0)
    keep = np.flatnonzero(largest_component(A))
    A = A[keep][:, keep].tocsr()
    Sx = Sx[keep]
    labL = lab[Sx]
    E = rank01(sub[m.COL["col4"]].to_numpy(float)[Sx])
    F = rank01(0.5 * (rank01(sub[m.COL["asma"]].to_numpy(float)[Sx]) + rank01(sub[m.COL["vim"]].to_numpy(float)[Sx])))
    tumour = labL == "tumor cells"
    source = np.flatnonzero(labL == "vasculature")
    sink = np.setdiff1d(core_sinks(A, tumour), source)
    res = compute_b_cell(A, E, F, source, sink, a=3.0, b_ecm=8.0, c_caf=4.0)
    cd8 = ok & (lab == "CD8+ T cells")
    return dict(xy=xy[Sx], tumour=tumour, source=source, sink=sink, E=E, cut=res["cut_edges"],
                cd8=xy[cd8])


def draw_core(ax, g, title):
    xy = g["xy"]
    m = (8 * g["E"]) / 8.0
    m_nt = m[~g["tumour"]]
    xy_nt = xy[~g["tumour"]]
    civ_high = m_nt >= np.median(m_nt)
    ax.scatter(xy_nt[~civ_high, 0], xy_nt[~civ_high, 1], s=1.2, color="#c8c8c8",
               linewidths=0, zorder=1, rasterized=True)
    ax.scatter(xy_nt[civ_high, 0], xy_nt[civ_high, 1], s=1.2, color="#5a5a5a",
               linewidths=0, zorder=1, rasterized=True)
    ax.scatter(xy[g["tumour"], 0], xy[g["tumour"], 1], s=1.2, color=S.LIGHT, linewidths=0, zorder=1,
               rasterized=True)
    ax.scatter(xy[g["sink"], 0], xy[g["sink"], 1], s=1.6, marker="s", color=S.INK2, linewidths=0, zorder=2,
               rasterized=True)
    ax.scatter(xy[g["source"], 0], xy[g["source"], 1], s=4, color=S.CELL, linewidths=0, zorder=3,
               rasterized=True)
    segs = [[xy[u], xy[v]] for u, v in g["cut"]]
    ax.add_collection(LineCollection(segs, colors=S.INK, linewidths=0.9, zorder=4))
    ax.scatter(g["cd8"][:, 0], g["cd8"][:, 1], s=5, marker="+", color=S.AQUA, linewidths=0.6, zorder=5)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s_ in ax.spines.values():
        s_.set_visible(False)
    allp = np.vstack([xy, g["cd8"]]) if len(g["cd8"]) else xy
    cxp = (allp[:, 0].min() + allp[:, 0].max()) / 2
    cyp = (allp[:, 1].min() + allp[:, 1].max()) / 2
    half = max(np.ptp(allp[:, 0]), np.ptp(allp[:, 1])) / 2 + 10
    ax.set_xlim(cxp - half, cxp + half)
    ax.set_ylim(cyp - half, cyp + half)
    ax.invert_yaxis()
    ax.set_title(title, fontsize=7, color=S.INK2, loc="left")
    sx, sy = cxp - half + 10, cyp + half - 12
    ax.plot([sx, sx + 100], [sy, sy], color=S.INK, lw=1.0)
    ax.text(sx + 50, sy - 6, "100 μm", ha="center", va="bottom", fontsize=6.4, color=S.INK2)


def main():
    S.apply()
    res = j("codex_validation.json")["summary"]
    cores = pd.read_csv(ROOT / "results" / "validation" / "codex_cores.csv")
    prim = cores[(cores.variant == "primary") & (cores.included)].reset_index(drop=True)
    hi, lo = pick_examples(prim)
    rhi = prim.set_index("core").loc[hi]
    rlo = prim.set_index("core").loc[lo]

    fig = plt.figure(figsize=(S.FULL_W, 128 * S.MM))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 1.05], height_ratios=[1.22, 1.0],
                          wspace=0.42, hspace=0.48, left=0.03, right=0.985, top=0.95, bottom=0.08)
    ax_hi = fig.add_subplot(gs[0, 0])
    ax_lo = fig.add_subplot(gs[0, 1])
    axb = fig.add_subplot(gs[0, 2])
    axc = fig.add_subplot(gs[1, 0:2])
    axd = fig.add_subplot(gs[1, 2])

    # ---------------- a: two example cores ----------------
    draw_core(ax_hi, core_geometry(hi),
              f"core {hi}: log B$_{{\\rm rel}}$ {rhi['log_b_rel']:.1f}, IR {rhi['ir']:+.1f}")
    draw_core(ax_lo, core_geometry(lo),
              f"core {lo}: log B$_{{\\rm rel}}$ {rlo['log_b_rel']:.1f}, IR {rlo['ir']:+.1f}")
    S.panel(ax_hi, "a", x=-0.12, y=1.06)
    handles = [Line2D([0], [0], marker="o", color="none", mfc=S.CELL, mec="none", ms=3.5, label="vessel (source)"),
               Line2D([0], [0], marker="s", color="none", mfc=S.INK2, mec="none", ms=3, label="tumour core (sink)"),
               Line2D([0], [0], marker="o", color="none", mfc=S.LIGHT, mec="none", ms=3, label="other tumour"),
               Line2D([0], [0], marker="o", color="none", mfc="#5a5a5a", mec="none", ms=3,
                      label="stroma (dark: collagen IV–high)"),
               Line2D([0], [0], color=S.INK, lw=1.0, label="minimum cut"),
               Line2D([0], [0], marker="+", color="none", mec=S.AQUA, ms=4, mew=0.8,
                      label="CD8$^+$ T cell (withheld)")]
    ax_hi.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.02), ncol=3, fontsize=6.3,
                 columnspacing=0.8, handletextpad=0.3)

    # ---------------- b: patient means ----------------
    g = prim.groupby("patient").agg(b=("log_b_rel", "mean"), ir=("ir", "mean"), grp=("group", "first"))
    for grp, mk in ((1, "o"), (2, "s")):
        sub = g[g.grp == grp]
        axb.plot(sub["b"], sub["ir"], marker=mk, ls="none", ms=3.6, mec=S.INK, mew=0.7,
                 mfc="white" if grp == 1 else S.INK, label="CLR" if grp == 1 else "DII")
    pr = res["primary_and_comparators"]["log_b_rel"]
    axb.axhline(0, color=S.HAIR, lw=0.6, zorder=0)
    axb.set_xlabel("Patient mean log B$_{\\rm rel}$ (scaffold only)")
    axb.set_ylabel("Patient mean CD8$^+$ core IR (log$_2$)")
    axb.text(0.98, 0.98, f"ρ = {pr['patient']['rho']:.2f}\none-sided p = {pr['patient']['p_one_sided_perm']:.4f}\n"
             f"{pr['patient']['n_patients']} patients", transform=axb.transAxes, ha="right", va="top",
             fontsize=6.6, color=S.INK2)
    axb.legend(loc="lower left", fontsize=6.4, handletextpad=0.2)
    S.panel(axb, "b", x=-0.28, y=1.02)

    # ---------------- c: comparators ----------------
    order = [("log_b_rel", "SPARTA B$_{\\rm rel}$ (min cut / matrix-free cut)", True),
             ("log_b_cell", "SPARTA B$_{\\rm cell}$ (absolute)", True),
             ("field_core", "SPARTA B$_{\\rm cell}$ field at the core", True),
             ("peritumoural_matrix", "Peritumoural matrix density", False),
             ("stromal_fraction", "Stromal fraction", False),
             ("ripley_l_matrix", "Ripley's L of matrix-rich cells", False),
             ("contact_enrichment_z", "Tumour–stroma contact enrichment", False),
             ("vessel_tumour_distance", "Vessel-to-tumour distance (post hoc)", False),
             ("core_depth", "Tumour-core depth (post hoc)", False)]
    comp = res["primary_and_comparators"]
    order = [o for o in order if o[0] in comp]
    for i, (k, lab, sparta) in enumerate(order):
        c = comp[k]
        col = S.CELL if sparta else S.INK2
        axc.plot(c["ci_core"], [i, i], color=S.LIGHT if not sparta else "#9ec5f4", lw=2.4,
                 solid_capstyle="butt", zorder=1)
        axc.plot([c["rho_core"]], [i], "o", color=col, ms=4.2, zorder=3)
        axc.plot([c["patient"]["rho"]], [i], marker="|", color=col, ms=7, mew=1.2, zorder=3)
        axc.text(1.22, i, f"{c['rho_core']:+.2f} / {c['patient']['rho']:+.2f}", va="center", ha="right",
                 fontsize=6.4, color=S.INK2)
    axc.text(1.22, -0.95, "core / patient ρ", va="center", ha="right", fontsize=6.4, color=S.INK2,
             fontweight="bold")
    axc.axvline(0, color=S.INK2, lw=0.6, zorder=0)
    axc.set_yticks(range(len(order)))
    axc.set_yticklabels([o[1] for o in order], fontsize=6.8)
    axc.set_ylim(len(order) - 0.4, -1.4)
    axc.set_xlim(-0.8, 1.23)
    axc.set_xticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
    axc.set_xlabel("Spearman ρ with CD8$^+$ core infiltration ratio (negative = barrier-like)")
    axc.tick_params(axis="y", length=0)
    axc.spines["left"].set_visible(False)
    S.panel(axc, "c", x=-0.42, y=1.0)

    # ---------------- d: CLR vs DII ----------------
    rng = np.random.default_rng(3)
    for i, (grp, name) in enumerate(((1, "CLR"), (2, "DII"))):
        v = g.loc[g.grp == grp, "b"].to_numpy()
        jit = (rng.random(len(v)) - 0.5) * 0.3
        axd.plot(i + jit, v, "o" if grp == 1 else "s", ms=3.2, mec=S.INK, mew=0.6,
                 mfc="white" if grp == 1 else S.INK, ls="none")
        axd.plot([i - 0.25, i + 0.25], [np.median(v)] * 2, color=S.CELL, lw=1.6)
    s3 = res["s3_clr_vs_dii"]
    axd.set_xticks([0, 1])
    axd.set_xticklabels([f"CLR\n(n={int((g.grp == 1).sum())})", f"DII\n(n={int((g.grp == 2).sum())})"])
    axd.set_xlim(-0.6, 1.6)
    axd.set_ylabel("Patient mean log B$_{\\rm rel}$")
    lo_, hi_ = float(g["b"].min()), float(g["b"].max())
    axd.set_ylim(lo_ - 0.15, hi_ + 0.75)
    axd.text(0.5, 0.98, f"Mann–Whitney p = {s3['mannwhitney_p_b_rel']:.3f}\n(exploratory)",
             transform=axd.transAxes, ha="center", va="top", fontsize=6.6, color=S.INK2)
    S.panel(axd, "d", x=-0.32, y=1.02)
    S.save(fig, OUT, "Fig7_codex")
    return hi, lo


if __name__ == "__main__":
    print(main())
