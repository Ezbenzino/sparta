"""Fig. 5 -- how much of the field coupling does the model produce by itself? (structural nulls)"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import figstyle as S  # noqa: E402
from isdata import EXTERNAL_ORDER, PATIENT_ORDER, PRIMARY_ORDER, j, ledger  # noqa: E402

OUT = Path(__file__).resolve().parents[2] / "results" / "figures" / "is"
G_GEOM = "#e1e0d9"
G_CONS = "#a8a79f"


def main():
    S.apply()
    geo = j("geometry_null.json")["per_slide"]
    abl = j("shared_input_spatial_null.json")["per_slide"]
    pl = j("patient_level_inference.json")
    led = ledger()

    fig = plt.figure(figsize=(S.FULL_W, 116 * S.MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1.0], height_ratios=[1.0, 1.0], wspace=0.45, hspace=0.62,
                          left=0.085, right=0.985, top=0.955, bottom=0.15)
    ax = fig.add_subplot(gs[:, 0])
    axb = fig.add_subplot(gs[0, 1])
    gsr = gs[1, 1].subgridspec(1, 2, wspace=0.55)
    axc = fig.add_subplot(gsr[0, 0])
    axd = fig.add_subplot(gsr[0, 1])

    # ---------------- a: observed vs geometry and construction nulls ----------------
    y = 0.0
    yt, yl, heads = [], [], []
    prev = None
    rows = []
    for sid in PRIMARY_ORDER + EXTERNAL_ORDER:
        pat = led[sid]["patient"] if led[sid]["cancer_type"] != "other" else "EXT"
        if pat != prev:
            if prev is not None:
                y += 0.35
            heads.append((y, S.PATIENT_LABEL.get(pat, "External sections (exploratory)")))
            y += 0.85
        rows.append((y, sid))
        yt.append(y)
        yl.append(sid)
        prev = pat
        y += 0.85
    for (yy, sid) in rows:
        r = geo[sid]
        ax.plot(r["geometry_null_ci95"], [yy + 0.17] * 2, color=G_GEOM, lw=2.2, solid_capstyle="butt", zorder=1)
        ax.plot(r["construction_null_ci95"], [yy - 0.17] * 2, color=G_CONS, lw=2.2, solid_capstyle="butt", zorder=1)
        pk = S.platform_key(sid)
        sig = r["p_vs_construction_null"] < 0.05
        ax.plot([r["rho_obs"]], [yy], marker=S.SHAPE[pk], ms=4.3, mec=S.INK, mew=0.8,
                mfc=S.INK if sig else "white", ls="none", zorder=3)
    for (yy, lab) in heads:
        ax.text(-0.335, yy, lab, va="center", ha="left", fontsize=6.6, color=S.INK2, style="italic")
    ax.axvline(0, color=S.INK2, lw=0.6, zorder=0)
    ax.set_yticks(yt)
    ax.set_yticklabels(yl, fontsize=6.6)
    ax.set_ylim(y - 0.2, -1.0)
    ax.set_xlim(-0.34, 0.52)
    ax.set_xlabel("Partial Spearman ρ")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    hs = [Line2D([0], [0], color=G_GEOM, lw=2.2, label="geometry null, 95%"),
          Line2D([0], [0], color=G_CONS, lw=2.2, label="construction null, 95%"),
          Line2D([0], [0], marker="o", color="none", mec=S.INK, mfc=S.INK, ms=4.3, label="observed, p < 0.05 vs construction"),
          Line2D([0], [0], marker="o", color="none", mec=S.INK, mfc="white", ms=4.3, label="observed, p ≥ 0.05")]
    ax.legend(handles=hs, loc="upper left", bbox_to_anchor=(-0.02, -0.095), ncol=2, fontsize=6.4,
              columnspacing=0.8, handletextpad=0.4)
    S.panel(ax, "a", x=-0.17, y=1.0)

    # ---------------- b: patient-level excess over the construction null ----------------
    ex = pl["construction_null_excess"]
    nm = ex["nested_model"]
    by = {}
    for sid in PRIMARY_ORDER:
        by.setdefault(led[sid]["patient"], []).append(geo[sid]["excess_over_construction_null"])
    for i, pat in enumerate(PATIENT_ORDER):
        v = np.array(by[pat])
        axb.plot([v.min(), v.max()], [i, i], color=S.LIGHT, lw=1.2, zorder=1)
        axb.plot(v, [i] * len(v), "|", color=S.MUTED, ms=5, mew=0.8, zorder=2)
        axb.plot([v.mean()], [i], "o", color=S.INK, ms=4.2, zorder=3)
    yp = len(PATIENT_ORDER) + 0.6
    axb.plot(nm["ci95"], [yp, yp], color=S.INK, lw=1.2)
    axb.plot([nm["mean"]], [yp], marker="D", color=S.INK, ms=5.0)
    axb.axvline(0, color=S.INK2, lw=0.6, zorder=0)
    axb.set_yticks(list(range(len(PATIENT_ORDER))) + [yp])
    axb.set_yticklabels([S.PATIENT_LABEL[p] for p in PATIENT_ORDER] + ["Pooled (nested model)"], fontsize=6.8)
    axb.set_ylim(yp + 0.8, -0.8)
    axb.set_xlabel("Observed ρ − construction-null mean")
    axb.tick_params(axis="y", length=0)
    axb.spines["left"].set_visible(False)
    axb.text(0.0, -0.31, f"Pooled excess {nm['mean']:.2f} (95% CI {nm['ci95'][0]:.2f}–{nm['ci95'][1]:.2f}); "
             f"{ex['sign_flip']['n_positive']}/{ex['sign_flip']['J']} patients > 0",
             transform=axb.transAxes, fontsize=6.6, color=S.INK2, va="top")
    S.panel(axb, "b", x=-0.50, y=1.02)

    # ---------------- c: ECM ablation, paired ----------------
    full = np.array([abl[s]["rho_partial_full_model"] for s in PRIMARY_ORDER])
    ab = np.array([abl[s]["real_rho_partial"] for s in PRIMARY_ORDER])
    for f_, a_ in zip(full, ab):
        axc.plot([0, 1], [f_, a_], color=S.LIGHT, lw=0.6, zorder=1)
    axc.plot(np.zeros_like(full), full, "o", ms=2.8, color=S.INK2, zorder=2)
    axc.plot(np.ones_like(ab), ab, "o", ms=2.8, color=S.INK2, zorder=2)
    axc.plot([0, 1], [np.median(full), np.median(ab)], color=S.INK, lw=1.5, zorder=3)
    axc.axhline(0, color=S.INK2, lw=0.6, zorder=0)
    axc.set_xticks([0, 1])
    axc.set_xticklabels(["Full\nmodel", "ECM term\nablated"], fontsize=6.6)
    axc.set_xlim(-0.35, 1.35)
    axc.set_ylabel("Partial Spearman ρ")
    S.panel(axc, "c", x=-0.55, y=1.04)

    # ---------------- d: input correlations ----------------
    keys = [("ecm_caf", "ECM\nCAF"), ("ecm_crosslink", "ECM\nXL"), ("caf_crosslink", "CAF\nXL"),
            ("ecm_ag", "ECM\nlig.")]
    for i, (k, lab) in enumerate(keys):
        v = np.array([geo[s]["input_correlations"][k] for s in PRIMARY_ORDER
                      if geo[s]["input_correlations"][k] is not None])
        jit = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.4
        axd.plot(i + jit, v, "o", ms=2.4, mfc="white", mec=S.MUTED, mew=0.6, ls="none")
        axd.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=S.INK, lw=1.4)
    axd.axhline(0, color=S.INK2, lw=0.6)
    axd.set_xticks(range(len(keys)))
    axd.set_xticklabels([k[1] for k in keys], fontsize=6.2)
    axd.set_ylim(-0.4, 1.0)
    axd.set_ylabel("Within-section Spearman ρ")
    S.panel(axd, "d", x=-0.55, y=1.04)
    S.save(fig, OUT, "Fig5_structure")


if __name__ == "__main__":
    main()
