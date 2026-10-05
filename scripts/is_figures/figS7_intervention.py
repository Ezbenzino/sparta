"""Fig. S7 (Online Resource 1) -- modality-specific counterfactual intervention maps (in-model).

Moved from the main text to Online Resource 1 in v2.2; panel e (compartments of the high-impact
spots, run_53) was added.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
ROOT = Path(__file__).resolve().parents[2]
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import figstyle as S  # noqa: E402
from isdata import PRIMARY_ORDER, REPLICATION_ORDER, j  # noqa: E402

OUT = ROOT / "results" / "figures" / "is" / "supplement"
EXAMPLE = "CSCC04"     # the section with the cohort-median field association (as in Fig. 1d-g)


def ramp(hex_hi):
    return LinearSegmentedColormap.from_list("r", ["#f4f3ef", hex_hi])


def map_panel(ax, xy, d, cut, color, label, top_frac=0.05):
    v = np.clip(d, 0, None)
    hi = np.quantile(v[v > 0], 0.97) if np.any(v > 0) else 1.0
    order = np.argsort(v)
    ax.scatter(xy[order, 0], xy[order, 1], c=v[order], cmap=ramp(color), vmin=0, vmax=hi, s=4.2,
               marker="h", linewidths=0, rasterized=True)
    k = max(1, int(round(top_frac * len(v))))
    top = np.argsort(-d)[:k]
    ax.scatter(xy[top, 0], xy[top, 1], s=6, facecolors="none", edgecolors=S.INK, linewidths=0.45)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for s_ in ax.spines.values():
        s_.set_visible(False)
    ax.set_title(label, fontsize=7, color=S.INK2, loc="left")
    x0 = xy[:, 0].min()
    y0 = xy[:, 1].max() + 120
    ax.plot([x0, x0 + 1000], [y0, y0], color=S.INK, lw=1.0)
    ax.text(x0 + 500, y0 + 60, "1 mm", ha="center", va="top", fontsize=6.4, color=S.INK2)


def main():
    S.apply()
    it = j("intervention_targeting.json")
    per, summ = it["per_section"], it["summary"]["primary"]
    z = np.load(ROOT / "results" / "intervention" / f"{EXAMPLE}.intervention.npz", allow_pickle=True)
    xy, dc, dm, cut = z["coords_um"], z["delta_b_cell"], z["delta_b_mab"], z["cut_mask"].astype(bool)

    fig = plt.figure(figsize=(S.FULL_W, 160 * S.MM))
    gs = fig.add_gridspec(3, 3, width_ratios=[1.0, 1.0, 1.15], height_ratios=[1.0, 0.95, 0.95],
                          wspace=0.35, hspace=0.75, left=0.02, right=0.985, top=0.96, bottom=0.07)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[0, 2])
    axd = fig.add_subplot(gs[1, :])
    axe = fig.add_subplot(gs[2, :])

    map_panel(axa, xy, dc, cut, S.CELL, f"{EXAMPLE}: single-spot ΔB$_{{\\rm cell}}$")
    map_panel(axb, xy, dm, cut, S.MAB, f"{EXAMPLE}: single-spot ΔB$_{{\\rm mAb}}$")
    S.panel(axa, "a", x=-0.14, y=1.06)
    S.panel(axb, "b", x=-0.14, y=1.06)
    axa.legend(handles=[Line2D([0], [0], marker="o", color="none", mfc="none", mec=S.INK, ms=3.4, mew=0.6,
                               label="top 5% of spots")], loc="upper left", bbox_to_anchor=(0.0, 0.02),
               fontsize=6.4, handletextpad=0.2)

    # ---------------- c: joint ablation curves ----------------
    for key, col, lab in (("curve_cell", S.CELL, "B$_{\\rm cell}$ (of a full cut breach)"),
                          ("curve_mab", S.MAB, "B$_{\\rm mAb}$ (of the all-spot limit)")):
        grid = np.array([0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 1.0])
        curves = []
        for sid in PRIMARY_ORDER:
            pts = np.array(per[sid][key], float)
            xk = pts[:, 0] / per[sid]["n"]
            curves.append(np.interp(np.log(grid), np.log(xk), pts[:, 1], left=np.nan, right=pts[-1, 1]))
        C = np.array(curves)
        med = np.nanmedian(C, axis=0)
        q1, q3 = np.nanpercentile(C, 25, axis=0), np.nanpercentile(C, 75, axis=0)
        axc.fill_between(grid, q1, q3, color=S.WASH_CELL if col == S.CELL else S.WASH_MAB, lw=0)
        axc.plot(grid, med, color=col, lw=1.4, label=lab)
    axc.axhline(0.5, color=S.HAIR, lw=0.6, zorder=0)
    axc.set_xscale("log")
    axc.set_xlim(0.002, 1.0)
    axc.set_ylim(0, 1.1)
    axc.set_xlabel("Spots ablated jointly (fraction of section, ranked)")
    axc.set_ylabel("Fraction of reference reduction")
    axc.legend(loc="lower right", fontsize=6.3)
    S.panel(axc, "c", x=-0.30, y=1.04)

    # ---------------- d: what a budget of 20 spots buys ----------------
    groups = [("cell", [("cell_sparta", "SPARTA\nranking"), ("cell_random_cut", "random\ncut spots"),
                        ("cell_density", "densest\nmatrix"), ("cell_random", "random\nspots"),
                        ("cell_with_mab_ranking", "B$_{\\rm mAb}$\nranking")]),
              ("mab", [("mab_sparta", "SPARTA\nranking"), ("mab_density", "densest\nmatrix"),
                       ("mab_random", "random\nspots"), ("mab_with_cell_ranking", "B$_{\\rm cell}$\nranking")])]
    xpos, xt, xl = 0, [], []
    rng = np.random.default_rng(5)
    for op, items in groups:
        col = S.CELL if op == "cell" else S.MAB
        for k, lab in items:
            v = np.array([per[s]["k20"][k] for s in PRIMARY_ORDER], float)
            jit = (rng.random(len(v)) - 0.5) * 0.36
            axd.plot(xpos + jit, v, "o", ms=2.6, mfc="white", mec=S.MUTED, mew=0.6, ls="none", zorder=2)
            axd.plot([xpos - 0.3, xpos + 0.3], [np.median(v)] * 2, color=col, lw=1.8, zorder=3)
            xt.append(xpos)
            xl.append(lab)
            xpos += 1
        xpos += 0.8
    axd.set_xticks(xt)
    axd.set_xticklabels(xl, fontsize=6.6)
    axd.set_ylabel("Joint reduction with 20 spots\n(fraction of reference)")
    axd.set_ylim(-0.03, 1.02)
    axd.set_xlim(-0.7, xpos - 0.6)
    n_cell = len(groups[0][1])
    axd.text((n_cell - 1) / 2, 1.0, "cellular barrier B$_{\\rm cell}$", ha="center", va="bottom", fontsize=7,
             color=S.CELL, fontweight="bold")
    axd.text(n_cell + 0.8 + (len(groups[1][1]) - 1) / 2, 1.0, "IgG-sized transport B$_{\\rm mAb}$", ha="center",
             va="bottom", fontsize=7, color=S.MAB, fontweight="bold")
    S.panel(axd, "d", x=-0.06, y=1.04)

    # ---------------- e: compartments of the high-impact spots (run_53) ----------------
    de = j("intervention_domain_enrichment.json")["summary"]
    shades = {"tumour": "#52514e", "stroma": "#a9a79e", "immune": "#e4e2db"}
    bars = [("all spots", None, None),
            ("B$_{\\rm cell}$ top 1%", "b_cell_top1pct", S.CELL), ("B$_{\\rm cell}$ top 5%", "b_cell_top5pct", S.CELL),
            ("B$_{\\rm mAb}$ top 1%", "b_mab_top1pct", S.MAB), ("B$_{\\rm mAb}$ top 5%", "b_mab_top5pct", S.MAB)]
    for i, (lab, key, col) in enumerate(bars):
        left = 0.0
        for d in ("tumour", "stroma", "immune"):
            v = de["b_cell_top1pct"][d]["pooled_frac_all"] if key is None else de[key][d]["pooled_frac_top"]
            axe.barh(i, v, left=left, color=shades[d], edgecolor="white", lw=0.6, height=0.62)
            if key is not None and d == "immune":
                p = de[key][d]["patient_signflip_p_two_sided"]
                axe.text(1.01, i, f"immune: p = {p:.3f}", va="center", ha="left", fontsize=6.2, color=S.INK2)
            axe.text(left + v / 2, i, f"{100 * v:.0f}%", va="center", ha="center", fontsize=6.0,
                     color="white" if d == "tumour" else S.INK)
            left += v
    axe.set_yticks(range(len(bars)))
    axe.set_yticklabels([b[0] for b in bars], fontsize=6.6)
    for t, b in zip(axe.get_yticklabels(), bars):
        if b[2] is not None:
            t.set_color(b[2])
    axe.invert_yaxis()
    axe.set_xlim(0, 1.0)
    axe.set_xlabel("Share of spots (pooled over 29 tumour sections); p, patient-level sign-flip test")
    axe.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=shades[d]) for d in ("tumour", "stroma", "immune")],
               labels=["tumour", "stromal", "immune"], loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=3,
               fontsize=6.4, handlelength=0.9)
    for sp_ in ("top", "right"):
        axe.spines[sp_].set_visible(False)
    S.panel(axe, "e", x=-0.06, y=1.10)
    S.save(fig, OUT, "FigS7_intervention", eps=False)


if __name__ == "__main__":
    main()
