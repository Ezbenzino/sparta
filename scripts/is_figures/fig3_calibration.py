"""Fig. 3 (and Fig. S2) -- calibration of the graph-spectral surrogate test on the real section graphs."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy.stats import binom  # noqa: E402

import figstyle as S  # noqa: E402
from isdata import EXTERNAL_ORDER, PRIMARY_ORDER, j  # noqa: E402

OUT = Path(__file__).resolve().parents[2] / "results" / "figures" / "is"
MK = dict(spec=dict(marker="o", mfc=S.INK, mec=S.INK, ms=3.8),
          nscore=dict(marker="s", mfc="white", mec=S.INK, ms=3.6),
          naive=dict(marker="^", mfc="white", mec=S.MUTED, ms=3.8))
DODGE = 0.2   # vertical offset so coinciding markers stay distinguishable
LAB = dict(spec="Surrogates, raw spectrum", nscore="Surrogates, normal scores", naive="Point-level test")


def band(n, p=0.05):
    lo, hi = binom.ppf(0.025, n, p) / n, binom.ppf(0.975, n, p) / n
    return lo, hi


def per_section_panel(ax, nc, scen, sids, label_y=True, title=None):
    y = 0.0
    yt, yl = [], []
    n_sim = None
    for k, sid in enumerate(sids):
        if sid not in nc["per_slide"]:
            continue
        if sid == EXTERNAL_ORDER[0]:
            y += 0.6
        r = nc["per_slide"][sid]
        d = r[scen] if not scen.startswith("tau") else r["smoothness"][scen]
        n_sim = d["n_sim"]
        ax.plot([d["fpr_naive_05"]], [y], ls="none", mew=0.7, **MK["naive"])
        if "fpr_nscore_05" in d:
            ax.plot([d["fpr_nscore_05"]], [y + DODGE], ls="none", mew=0.7, **MK["nscore"])
        ax.plot([d["fpr_spectral_05"]], [y - DODGE], ls="none", mew=0.7, **MK["spec"])
        yt.append(y)
        yl.append(sid)
        y += 1.0
    if n_sim:
        lo, hi = band(n_sim)
        ax.axvspan(lo, hi, color="#ecebe6", zorder=0, lw=0)
    ax.axvline(0.05, color=S.INK2, lw=0.6, zorder=0)
    ax.set_yticks(yt)
    ax.set_yticklabels(yl if label_y else [], fontsize=6.2)
    ax.set_ylim(y - 0.4, -0.6)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(-0.012, 0.5)
    ax.set_xlabel("False-positive rate at α = 0.05")
    if title:
        ax.set_title(title, fontsize=7.2, loc="left")


def main():
    S.apply()
    nc = j("null_calibration_all.json")
    sm = nc.get("summary", {})
    sids = PRIMARY_ORDER + EXTERNAL_ORDER
    fig = plt.figure(figsize=(S.FULL_W, 104 * S.MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.05], height_ratios=[1.0, 1.0], wspace=0.40, hspace=0.75,
                          left=0.085, right=0.985, top=0.95, bottom=0.16)
    axa = fig.add_subplot(gs[:, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[1, 1])
    per_section_panel(axa, nc, "gaussian", sids)
    hs = [Line2D([0], [0], ls="none", mew=0.7, label=LAB[k], **MK[k]) for k in ("spec", "nscore", "naive")]
    hs.append(Line2D([0], [0], color="#ecebe6", lw=5, label="95% binomial range"))
    axa.legend(handles=hs, loc="upper left", bbox_to_anchor=(-0.02, -0.11), ncol=2, fontsize=6.4,
               handletextpad=0.3, columnspacing=0.8)
    S.panel(axa, "a", x=-0.25, y=1.0)

    # ---------------- b: pooled by scenario ----------------
    scen = [("gaussian", "Spectrum-matched"), ("marginal", "Observed marginal"), ("skewed", "Strongly skewed"),
            ("smoothness_tau_2", "Heat kernel, τ = 2"), ("smoothness_tau_10", "Heat kernel, τ = 10"),
            ("smoothness_tau_50", "Heat kernel, τ = 50")]
    scen = [s for s in scen if s[0] in sm]
    for i, (k, lab) in enumerate(scen):
        d = sm[k]
        axb.plot([d["pooled_fpr_naive_05"]], [i], ls="none", mew=0.7, **MK["naive"])
        if "pooled_fpr_nscore_05" in d:
            axb.plot([d["pooled_fpr_nscore_05"]], [i + DODGE], ls="none", mew=0.7, **MK["nscore"])
        axb.plot([d["pooled_fpr_spectral_05"]], [i - DODGE], ls="none", mew=0.7, **MK["spec"])
    axb.axvline(0.05, color=S.INK2, lw=0.6, zorder=0)
    axb.set_yticks(range(len(scen)))
    axb.set_yticklabels([s[1] for s in scen], fontsize=6.8)
    axb.set_ylim(len(scen) - 0.4, -0.6)
    axb.set_xlim(0, 0.45)
    axb.tick_params(axis="y", length=0)
    axb.spines["left"].set_visible(False)
    axb.set_xlabel("Pooled false-positive rate at α = 0.05")
    S.panel(axb, "b", x=-0.62, y=1.03)

    # ---------------- c: p-value distribution ----------------
    if "gaussian" in sm:
        dec = np.array(sm["gaussian"]["pooled_p_deciles"], float)
        dec = dec / dec.sum()
        edges = np.linspace(0, 1, 11)
        axc.bar(edges[:-1] + 0.05, dec, width=0.088, color=S.LIGHT, edgecolor="none", label=LAB["spec"])
        if "pooled_p_nscore_deciles" in sm["gaussian"]:
            dz = np.array(sm["gaussian"]["pooled_p_nscore_deciles"], float)
            dz = dz / dz.sum()
            axc.step(np.r_[edges[:-1], 1.0], np.r_[dz, dz[-1]], where="post", color=S.INK, lw=0.9,
                     label=LAB["nscore"])
        axc.axhline(0.1, color=S.INK2, lw=0.6)
        axc.set_xlim(0, 1)
        axc.set_ylim(0, 0.2)
        axc.set_xlabel("Surrogate-test p-value (spectrum-matched scenario)")
        axc.set_ylabel("Fraction of simulations")
        axc.legend(loc="upper right", fontsize=6.4)
    S.panel(axc, "c", x=-0.22, y=1.03)
    S.save(fig, OUT, "Fig3_calibration")

    # ---------------- Fig. S2: all scenarios per section ----------------
    keys = [("gaussian", "Spectrum-matched"), ("marginal", "Observed marginal"), ("skewed", "Strongly skewed"),
            ("tau_10", "Heat kernel, τ = 10")]
    keys = [k for k in keys if (k[0] in nc["per_slide"][sids[0]]) or k[0].startswith("tau")]
    fig2, axs = plt.subplots(1, len(keys), figsize=(S.FULL_W, 120 * S.MM))
    for i, (k, lab) in enumerate(keys):
        per_section_panel(axs[i], nc, k, sids, label_y=(i == 0), title=lab)
    axs[0].legend(handles=hs, loc="upper left", bbox_to_anchor=(0.0, -0.10), ncol=4, fontsize=6.4,
                  handletextpad=0.3, columnspacing=0.8)
    fig2.subplots_adjust(left=0.08, right=0.99, top=0.95, bottom=0.17, wspace=0.12)
    S.save(fig2, OUT / "supplement", "FigS2_calibration", eps=False)


if __name__ == "__main__":
    main()
