"""Fig. 4 -- field association within sections, strata and patient-level models."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E401, E402

import figstyle as S  # noqa: E402
from isdata import PATIENT_ORDER, PRIMARY_ORDER, REPLICATION_PATIENTS, j, ledger  # noqa: E402

OUT = Path(__file__).resolve().parents[2] / "results" / "figures" / "is"


def platform_for_primary(sid: str) -> str:
    return S.platform_key(sid)


def platform_for_extension(record: dict) -> str:
    if record["platform"] == "Slide-seqV2":
        return "slideseq"
    if record["cancer_type"] == "cscc":
        return "visium_cscc"
    return "visium_mel"


def main():
    S.apply()
    sn = j("spatial_null_check.json")["per_slide"]
    pl = j("patient_level_inference.json")
    adj = j("adjustment_robustness.json")["per_slide"]
    rep = j("replication_melanoma.json")
    rsn = rep["per_slide"]
    ext = j("extension_cohort.json")
    ep = ext["per_slide"]
    led = ledger()

    fig = plt.figure(figsize=(S.FULL_W, 165 * S.MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.18, 1.0], height_ratios=[1.45, 1.0],
                          wspace=0.48, hspace=0.48, left=0.16, right=0.90,
                          top=0.955, bottom=0.30)
    ax = fig.add_subplot(gs[:, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[1, 1])

    # ---------------- a: all sections by cohort / tumour stratum ----------------
    groups = [
        ("Primary cSCC", [
            {"sid": sid, "rho": sn[sid]["real_rho_partial"], "q": sn[sid]["bh_padj"],
             "platform": platform_for_primary(sid)}
            for sid in PRIMARY_ORDER if sid.startswith("CSCC")
        ]),
        ("Primary melanoma", [
            {"sid": sid, "rho": sn[sid]["real_rho_partial"], "q": sn[sid]["bh_padj"],
             "platform": platform_for_primary(sid)}
            for sid in PRIMARY_ORDER if sid.startswith("MEL")
        ]),
        ("Replication melanoma", [
            {"sid": sid, "rho": rsn[sid]["real_rho_partial"], "q": rsn[sid]["bh_padj"],
             "platform": "legacy_mel"}
            for sid in rsn
        ]),
        ("Extension cSCC", [
            {"sid": sid, "rho": r["real_rho_partial"], "q": r["bh_padj"],
             "platform": platform_for_extension(r)}
            for sid, r in ep.items() if r["cancer_type"] == "cscc"
        ]),
        ("Extension primary melanoma", [
            {"sid": sid, "rho": r["real_rho_partial"], "q": r["bh_padj"],
             "platform": platform_for_extension(r)}
            for sid, r in ep.items()
            if r["cancer_type"] == "melanoma" and r.get("site") in ("skin", "primary skin")
        ]),
        ("Extension metastatic melanoma", [
            {"sid": sid, "rho": r["real_rho_partial"], "q": r["bh_padj"],
             "platform": platform_for_extension(r)}
            for sid, r in ep.items()
            if r["cancer_type"] == "melanoma" and r.get("site") not in ("skin", "primary skin")
        ]),
    ]

    rng = np.random.default_rng(20261005)
    yticks, ylabels = [], []
    wrapped = {"Extension primary melanoma": "Extension primary\nmelanoma",
               "Extension metastatic melanoma": "Extension metastatic\nmelanoma"}
    for y, (label, records) in enumerate(groups):
        yticks.append(y)
        ylabels.append(f"{wrapped.get(label, label)}\n(n={len(records)})")
        jitter = rng.permutation(np.linspace(-0.32, 0.32, len(records)))
        for dy, record in zip(jitter, records):
            filled = record["q"] < 0.05
            ax.plot(record["rho"], y + dy, marker=S.SHAPE[record["platform"]],
                    ms=3.8, mec=S.INK, mew=0.65,
                    mfc=S.INK if filled else "white", ls="none", zorder=3)
        rhos = np.array([r["rho"] for r in records])
        q1, med, q3 = np.percentile(rhos, [25, 50, 75])
        ax.plot([q1, q3], [y, y], color=S.LIGHT, lw=3.2, solid_capstyle="butt", zorder=1)
        ax.plot([med, med], [y - 0.18, y + 0.18], color=S.INK, lw=1.0, zorder=2)

    ax.axvline(0, color=S.INK2, lw=0.6, zorder=0)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=6.8)
    ax.set_ylim(len(groups) - 0.55, -0.55)
    ax.set_xlim(-0.22, 0.68)
    ax.set_xlabel("Partial Spearman ρ (B$_{\\rm cell}$ field vs B$_{\\rm mAb}$ field)")
    ax.tick_params(axis="y", length=0)
    S.panel(ax, "a", x=-0.23, y=1.01)
    legend_keys = ["visium_cscc", "legacy_st", "visium_mel", "legacy_mel", "slideseq"]
    handles = [
        Line2D([0], [0], marker=S.SHAPE[k], color="none", mec=S.INK, mfc=S.INK, ms=4.0,
               label=S.SHAPE_LABEL[k]) for k in legend_keys
    ]
    handles += [
        Line2D([0], [0], marker="o", color="none", mec=S.INK, mfc="white", ms=4.0,
               label="Open: BH q ≥ 0.05"),
        Line2D([0], [0], color=S.LIGHT, lw=3.0, label="Interquartile range"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.46, -0.30),
              ncol=3, fontsize=6.1, columnspacing=0.7, handletextpad=0.3,
              borderaxespad=0.0)

    # ---------------- b: patient-level pooled estimates ----------------
    primary_model = pl["primary_association"]
    replication_model = rep["patient_level"]["association"]
    all12_model = rep["patient_level"]["pooled_primary_plus_replication"]
    extension_model = ext["patient_level"]["association"]
    primary_extension_model = ext["patient_level"]["pooled_primary_plus_extension"]
    pooled = [
        ("Primary", primary_model),
        ("Replication", replication_model),
        ("Extension", extension_model),
        ("Primary + replication", all12_model),
        ("Primary + extension", primary_extension_model),
    ]
    for y, (label, record) in enumerate(pooled):
        fit = record["nested_model"]
        lo, hi = fit["ci95"]
        axb.plot([lo, hi], [y, y], color=S.INK, lw=1.1)
        axb.plot(fit["mean"], y, marker="D", color=S.INK, ms=4.4)
        sf = record["sign_flip"]
        axb.text(0.52, y, f"{sf['n_positive']}/{sf['J']} > 0", va="center", ha="left",
                 fontsize=6.2, color=S.INK2)
    axb.axvline(0, color=S.INK2, lw=0.6)
    axb.set_yticks(range(len(pooled)))
    axb.set_yticklabels([item[0] for item in pooled], fontsize=6.8)
    axb.set_ylim(len(pooled) - 0.6, -0.6)
    axb.set_xlim(-0.03, 0.52)
    axb.set_xlabel("Pooled partial Spearman ρ")
    axb.tick_params(axis="y", length=0)
    S.panel(axb, "b", x=-0.42, y=1.03)

    # ---------------- c: adjustment robustness ----------------
    variants = [("rank_quadratic", "Rank\nquadratic*"), ("rank_cubic", "Rank\ncubic"),
                ("spline_raw", "B-spline\n(μm)"), ("shell_stratified", "Distance\nshells"),
                ("unadjusted", "Un-\nadjusted")]
    for i, (key, label) in enumerate(variants):
        values = np.array([adj[s][key] for s in PRIMARY_ORDER], float)
        jit = (rng.random(len(values)) - 0.5) * 0.36
        axc.plot(i + jit, values, "o", ms=2.6, mfc="white", mec=S.MUTED,
                 mew=0.6, ls="none", zorder=2)
        q1, med, q3 = np.percentile(values, [25, 50, 75])
        axc.plot([i - 0.28, i + 0.28], [med, med], color=S.INK, lw=1.3, zorder=3)
        axc.plot([i, i], [q1, q3], color=S.INK, lw=0.8, zorder=3)
        axc.text(i, 0.475, f"{int((values > 0).sum())}/19 > 0",
                 ha="center", va="bottom", fontsize=6.2, color=S.INK2)
    axc.axhline(0, color=S.INK2, lw=0.6)
    axc.set_xticks(range(len(variants)))
    axc.set_xticklabels([item[1] for item in variants], fontsize=6.4)
    axc.set_ylim(-0.06, 0.52)
    axc.set_ylabel("Partial Spearman ρ")
    S.panel(axc, "c", x=-0.28, y=1.03)

    S.save(fig, OUT, "Fig4_association")


if __name__ == "__main__":
    main()
