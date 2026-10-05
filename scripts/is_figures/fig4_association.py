"""Fig. 4 -- field association within sections and across patients."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import figstyle as S  # noqa: E402
from isdata import (EXTERNAL_ORDER, PATIENT_ORDER, PRIMARY_ORDER, REPLICATION_ORDER,  # noqa: E402
                    REPLICATION_PATIENTS, j, ledger)

OUT = Path(__file__).resolve().parents[2] / "results" / "figures" / "is"


def main():
    S.apply()
    sn = j("spatial_null_check.json")["per_slide"]
    ext = j("ext_validation.json")["per_slide"]
    pl = j("patient_level_inference.json")
    adj = j("adjustment_robustness.json")["per_slide"]
    rep = j("replication_melanoma.json")
    rsn = rep["per_slide"]
    led = ledger()

    fig = plt.figure(figsize=(S.FULL_W, 150 * S.MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 1.0], height_ratios=[1.35, 1.0],
                          wspace=0.45, hspace=0.42, left=0.09, right=0.985, top=0.965, bottom=0.13)
    ax = fig.add_subplot(gs[:, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[1, 1])

    # ---------------- a: per-section forest with spatial-null interval ----------------
    rows, headers = [], []
    y = 0.0
    yt, ylab = [], []
    prev_pat = None
    for sid in PRIMARY_ORDER:
        pat = led[sid]["patient"]
        if pat != prev_pat:
            if prev_pat is not None:
                y += 0.35
            headers.append((y, S.PATIENT_LABEL[pat]))
            y += 0.85
        r = sn[sid]
        rows.append((y, sid, r["real_rho_partial"], r["null_ci95"], r["bh_padj"], S.platform_key(sid)))
        yt.append(y)
        ylab.append(sid)
        prev_pat = pat
        y += 0.85
    y += 0.35
    headers.append((y, "Replication cohort (melanoma lymph-node metastases)"))
    y += 0.85
    prev_pat = None
    for sid in REPLICATION_ORDER:
        pat = led[sid]["patient"]
        if pat != prev_pat and prev_pat is not None:
            y += 0.15
        r = rsn[sid]
        rows.append((y, sid, r["real_rho_partial"], r["null_ci95"], r["bh_padj"], S.platform_key(sid)))
        yt.append(y)
        ylab.append(sid.replace("MEL_THR", "LN P").replace("_rep", " r"))
        prev_pat = pat
        y += 0.85
    y += 0.35
    headers.append((y, "External sections (exploratory)"))
    y += 0.85
    for sid in EXTERNAL_ORDER:
        r = ext[sid]["spatial_null"]
        rows.append((y, sid, r["real_rho_partial"], r["null_ci95"], r["bh_padj"], "external"))
        yt.append(y)
        ylab.append(sid)
        y += 0.85

    for (yy, sid, rho, ci, q, pk) in rows:
        ax.plot(ci, [yy, yy], color=S.LIGHT, lw=2.4, solid_capstyle="butt", zorder=1)
        filled = q < 0.05
        ax.plot([rho], [yy], marker=S.SHAPE[pk], ms=4.4, mec=S.INK, mew=0.8,
                mfc=S.INK if filled else "white", ls="none", zorder=3)
        qtxt = f"{q:.4f}" if 0.045 <= q < 0.055 else f"{q:.3f}"   # keep borderline values unambiguous
        ax.text(0.60, yy, qtxt, va="center", ha="right", fontsize=6.8, color=S.INK2)
    for (yy, lab) in headers:
        ax.text(-0.355, yy, lab, va="center", ha="left", fontsize=6.8, color=S.INK2, style="italic")
    ax.text(0.60, -0.9, "BH q", va="center", ha="right", fontsize=6.8, color=S.INK2, fontweight="bold")
    ax.axvline(0, color=S.INK2, lw=0.6, zorder=0)
    ax.set_yticks(yt)
    ax.set_yticklabels(ylab, fontsize=6.8)
    ax.set_ylim(y - 0.2, -1.4)
    ax.set_xlim(-0.36, 0.61)
    ax.set_xticks([-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4])
    ax.set_xlabel("Partial Spearman ρ (B$_{\\rm cell}$ field vs B$_{\\rm mAb}$ field)")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    S.panel(ax, "a", x=-0.16, y=1.0)
    handles = [Line2D([0], [0], marker=S.SHAPE[k], color="none", mec=S.INK, mfc=S.INK, ms=4.4,
                      label=S.SHAPE_LABEL[k]) for k in ("visium_cscc", "legacy_st", "visium_mel",
                                                         "legacy_mel", "external")]
    handles += [Line2D([0], [0], marker="o", color="none", mec=S.INK, mfc="white", ms=4.4,
                       label="Open: BH q ≥ 0.05"),
                Line2D([0], [0], color=S.LIGHT, lw=2.4, label="95% surrogate null")]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.02, -0.075), ncol=3,
              fontsize=6.4, columnspacing=0.8, handletextpad=0.3)

    # ---------------- b: patient level ----------------
    prim = pl["primary_association"]
    nm = prim["nested_model"]
    sf = prim["sign_flip"]
    rnm = rep["patient_level"]["association"]["nested_model"]
    rsf = rep["patient_level"]["association"]["sign_flip"]
    pnm = rep["patient_level"]["pooled_primary_plus_replication"]["nested_model"]
    psf = rep["patient_level"]["pooled_primary_plus_replication"]["sign_flip"]
    by_pat = {}
    for sid in PRIMARY_ORDER:
        by_pat.setdefault(led[sid]["patient"], []).append(sn[sid]["real_rho_partial"])
    for sid in REPLICATION_ORDER:
        by_pat.setdefault(led[sid]["patient"], []).append(rsn[sid]["real_rho_partial"])
    ylist = []
    yy = 0.0
    for pat in PATIENT_ORDER:
        ylist.append((yy, pat))
        yy += 1
    yy += 0.4
    for pat in REPLICATION_PATIENTS:
        ylist.append((yy, pat))
        yy += 1
    for i, pat in ylist:
        v = np.array(by_pat[pat])
        axb.plot([v.min(), v.max()], [i, i], color=S.LIGHT, lw=1.2, zorder=1)
        axb.plot(v, [i] * len(v), "|", color=S.MUTED, ms=5, mew=0.8, zorder=2)
        axb.plot([v.mean()], [i], "o", color=S.INK, ms=4.2, zorder=3)
        axb.text(0.47, i, f"n={len(v)}", va="center", ha="right", fontsize=6.6, color=S.INK2)
    pools = [(yy + 0.5, nm, "Primary, pooled (7)"), (yy + 1.5, rnm, "Replication, pooled (4)"),
             (yy + 2.5, pnm, "All 11 patients")]
    for ypool, fit, _ in pools:
        lo, hi = fit["ci95"]
        axb.plot([lo, hi], [ypool, ypool], color=S.INK, lw=1.2)
        axb.plot([fit["mean"]], [ypool], marker="D", color=S.INK, ms=5.0)
    axb.axvline(0, color=S.INK2, lw=0.6, zorder=0)
    axb.set_yticks([i for i, _ in ylist] + [p[0] for p in pools])
    axb.set_yticklabels([S.PATIENT_LABEL[p] for _, p in ylist] + [p[2] for p in pools], fontsize=6.8)
    axb.set_ylim(pools[-1][0] + 0.8, -0.8)
    axb.set_xlim(-0.05, 0.48)
    axb.set_xlabel("Partial Spearman ρ")
    axb.tick_params(axis="y", length=0)
    axb.spines["left"].set_visible(False)
    lo, hi = nm["ci95"]
    axb.text(1.02, 1.0,
             f"Primary: ρ = {nm['mean']:.2f} ({lo:.2f}–{hi:.2f}); {sf['n_positive']}/{sf['J']} > 0, "
             f"sign-flip p = {sf['p_one_sided']:.4f}\n"
             f"Replication: ρ = {rnm['mean']:.2f} ({rnm['ci95'][0]:.2f}–{rnm['ci95'][1]:.2f}); "
             f"{rsf['n_positive']}/{rsf['J']} > 0, p = {rsf['p_one_sided']:.4f}\n"
             f"All 11: ρ = {pnm['mean']:.2f} ({pnm['ci95'][0]:.2f}–{pnm['ci95'][1]:.2f}); "
             f"{psf['n_positive']}/{psf['J']} > 0, p = {psf['p_one_sided']:.5f}",
             transform=axb.transAxes, fontsize=6.0, color=S.INK2, va="top", ha="left", linespacing=1.35,
             rotation=0, clip_on=False, visible=False)
    S.panel(axb, "b", x=-0.47, y=1.02)

    # ---------------- c: adjustment robustness ----------------
    variants = [("rank_quadratic", "Rank\nquadratic*"), ("rank_cubic", "Rank\ncubic"),
                ("spline_raw", "B-spline\n(μm)"), ("shell_stratified", "Distance\nshells"),
                ("unadjusted", "Un-\nadjusted")]
    for i, (k, lab) in enumerate(variants):
        v = np.array([adj[s][k] for s in PRIMARY_ORDER], float)
        rng = np.random.default_rng(i)
        jit = (rng.random(len(v)) - 0.5) * 0.36
        axc.plot(i + jit, v, "o", ms=2.6, mfc="white", mec=S.MUTED, mew=0.6, ls="none", zorder=2)
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        axc.plot([i - 0.28, i + 0.28], [med, med], color=S.INK, lw=1.4, zorder=3)
        axc.plot([i, i], [q1, q3], color=S.INK, lw=0.8, zorder=3)
        axc.text(i, 0.475, f"{int((v > 0).sum())}/19 > 0", ha="center", va="bottom", fontsize=6.4, color=S.INK2)
    axc.axhline(0, color=S.INK2, lw=0.6, zorder=0)
    axc.set_xticks(range(len(variants)))
    axc.set_xticklabels([v[1] for v in variants], fontsize=6.6)
    axc.set_ylim(-0.06, 0.52)
    axc.set_ylabel("Partial Spearman ρ")
    S.panel(axc, "c", x=-0.30, y=1.02)
    S.save(fig, OUT, "Fig4_association")


if __name__ == "__main__":
    main()
