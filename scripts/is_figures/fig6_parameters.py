"""Fig. 6 -- the size-exclusion law is parameter-driven; size scan; selection-matched gap experiment."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import figstyle as S  # noqa: E402
from isdata import PATIENT_ORDER, PRIMARY_ORDER, VAL, j, ledger  # noqa: E402
from sparta.barrier import edge_pairs  # noqa: E402
from sparta.io_ import Paths, load_config, load_graph  # noqa: E402
from sparta.node_tables import load_nodes, scores_from_nodes  # noqa: E402

OUT = Path(__file__).resolve().parents[2] / "results" / "figures" / "is"
BETAS = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0])


def excluded_fraction(xl_edge, beta, xi0=20.0, r=5.5):
    xi = xi0 * np.exp(-beta * xl_edge)
    return float(np.mean(r >= xi))


def main():
    S.apply()
    cfg = load_config(None)
    P = Paths(cfg)
    led = ledger()
    fig = plt.figure(figsize=(S.FULL_W, 108 * S.MM))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.15], height_ratios=[1, 1], wspace=0.48, hspace=0.62,
                          left=0.075, right=0.985, top=0.95, bottom=0.10)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[1, 0:2])
    axd = fig.add_subplot(gs[:, 2])

    # ---------------- a: xi(x) for several beta ----------------
    x = np.linspace(0, 1, 300)
    xs = np.log(20 / 5.5) / 3.0
    axa.fill_between(x[x >= xs], 0, 20.0 * np.exp(-3.0 * x[x >= xs]), color="#fbe0d3", zorder=0, lw=0)
    for beta, shade, lw in [(1.0, S.LIGHT, 1.0), (2.0, S.MUTED, 1.0), (3.0, S.MAB, 1.6)]:
        axa.plot(x, 20.0 * np.exp(-beta * x), color=shade, lw=lw,
                 label=f"β = {beta:g}" + (" (default)" if beta == 3.0 else ""))
    axa.legend(loc="upper right", fontsize=6.4, handlelength=1.4, borderaxespad=0.1)
    axa.axhline(5.5, color=S.INK, lw=0.7)
    axa.text(0.015, 5.0, "r = 5.5 nm (IgG)", fontsize=6.2, ha="left", va="top", color=S.INK)
    axa.plot([xs, xs], [0, 5.5], color=S.INK, lw=0.7)
    axa.text(xs - 0.02, 0.5, f"X* = {xs:.2f}", fontsize=6.4, ha="right", va="bottom")
    axa.text(0.47, 0.4, "excluded", fontsize=6.0, ha="left", va="bottom", color=S.INK2)
    axa.set_xlim(0, 1)
    axa.set_ylim(0, 21)
    axa.set_xlabel("Edge crosslinking score (rank)")
    axa.set_ylabel("Effective pore radius ξ (nm)")
    S.panel(axa, "a", x=-0.30, y=1.02)

    # ---------------- b: excluded fraction vs beta, per section ----------------
    scan = j("size_exclusion_scan.json")
    BETAS_ = np.array(scan["betas"])
    curves = {sid: np.array(scan["per_slide"][sid]["excluded_fraction"]) for sid in PRIMARY_ORDER}
    for sid, c in curves.items():
        axb.plot(BETAS_, 100 * c, color=S.LIGHT, lw=0.6, zorder=1)
    med = np.median(np.vstack(list(curves.values())), axis=0)
    axb.plot(BETAS_, 100 * med, color=S.MAB, lw=1.5, zorder=2)
    at3 = np.array([curves[s][list(BETAS_).index(3.0)] for s in curves])
    axb.plot([3.0], [100 * np.median(at3)], "o", color=S.MAB, ms=4, zorder=3, mec="white", mew=0.6)
    axb.annotate(f"β = 3: {100 * at3.min():.1f}–{100 * at3.max():.1f}%\nin all 19 sections",
                 xy=(3.0, 100 * np.median(at3)), xytext=(3.9, 22), fontsize=6.4, color=S.INK,
                 arrowprops=dict(arrowstyle="-", lw=0.5, color=S.INK2))
    axb.set_xlabel("Mesh contraction β")
    axb.set_ylabel("Edges excluded (%)")
    axb.set_ylim(0, 100)
    axb.set_xlim(0, 8.3)
    S.panel(axb, "b", x=-0.30, y=1.02)

    # ---------------- c: size scan ----------------
    radii = None
    ys = []
    for sid in PRIMARY_ORDER:
        c = json.load(open(VAL.parent / "counterfactual" / f"{sid}.json", encoding="utf-8"))["s3"]
        radii = np.array(c["radii_nm"])
        yv = np.array(c["mean_core"])
        ys.append(yv)
        axc.plot(radii, yv, color=S.LIGHT, lw=0.6, zorder=1)
    ys = np.vstack(ys)
    axc.plot(radii, np.median(ys, axis=0), color=S.MAB, lw=1.5, marker="o", ms=3, mec="white", mew=0.5, zorder=2)
    axc.axvline(5.5, color=S.INK2, lw=0.6, zorder=0)
    axc.text(5.6, np.max(ys) * 0.96, "IgG (5.5 nm)", fontsize=6.6, ha="left", va="top", color=S.INK2)
    axc.text(0.5, np.max(ys) * 0.96, "thin lines: 19 sections; orange: median", fontsize=6.6, ha="left", va="top",
             color=S.INK2)
    axc.set_xlabel("Specified molecular radius r (nm)")
    axc.set_ylabel("Mean tumour-core B$_{\\rm mAb}$")
    axc.set_xlim(0, 10.3)
    S.panel(axc, "c", x=-0.12, y=1.02)

    # ---------------- d: selection-matched gap experiment ----------------
    s2 = j("s2_matched_selection.json")["per_slide"]
    from isdata import bh
    q = dict(zip(PRIMARY_ORDER, bh([s2[x]["p_vs_in_cut_matched"] for x in PRIMARY_ORDER])))
    y = 0.0
    yt, yl = [], []
    prev = None
    for sid in PRIMARY_ORDER:
        pat = led[sid]["patient"]
        if prev is not None and pat != prev:
            y += 0.5
        r = s2[sid]["ratio_vs_in_cut_matched"]
        pk = S.platform_key(sid)
        axd.plot([r], [y], marker=S.SHAPE[pk], ms=4.2, mec=S.INK, mew=0.8, mfc=S.INK if q[sid] < 0.05 else "white",
                 ls="none", zorder=3)
        yt.append(y)
        yl.append(sid)
        prev = pat
        y += 1.0
    XMIN, XMAX = 0.4, 2.0
    rng = np.random.default_rng(0)

    def sim_row(vals, y):
        vals = np.asarray(vals, float)
        jit = (rng.random(len(vals)) - 0.5) * 0.5
        inside = (vals >= XMIN) & (vals <= XMAX)
        axd.plot(vals[inside], y + jit[inside], "o", ms=3.2, mfc="white", mec=S.MUTED, mew=0.7, ls="none")
        hi = vals[vals > XMAX]
        if len(hi):
            axd.plot([XMAX - 0.015] * len(hi), y + jit[vals > XMAX], ">", ms=3.4, mfc=S.MUTED, mec=S.MUTED,
                     ls="none", clip_on=False)
            axd.text(XMAX + 0.06, y, f"{len(hi)} beyond axis", fontsize=5.6, ha="left", va="center",
                     color=S.INK2, clip_on=False)
        lo = vals[vals < XMIN]
        if len(lo):
            axd.plot([XMIN + 0.015] * len(lo), y + jit[vals < XMIN], "<", ms=3.4, mfc=S.MUTED, mec=S.MUTED,
                     ls="none", clip_on=False)

    p2 = VAL / "synthetic_s2_thickness.json"
    if p2.exists():
        th = j("synthetic_s2_thickness.json")["results"]
        y += 0.3
        for k in sorted(th, key=lambda s: int(s[:-2])):
            y += 1.6
            sim_row([r_["ratio_matched"] for r_ in th[k]["per_rep"]], y)
            yt.append(y)
            yl.append(f"Simulated capsule, {k[:-2]} μm")
    else:
        y += 1.4
        sim_row(j("synthetic_benchmark.json")["synthetic_s2_closed"]["ratios"], y)
        yt.append(y)
        yl.append("Simulated\nclosed capsules")
    axd.axvline(1.0, color=S.INK2, lw=0.6, zorder=0)
    axd.set_yticks(yt)
    axd.set_yticklabels(yl, fontsize=6.6)
    axd.set_ylim(y + 0.8, -0.8)
    axd.set_xlim(XMIN, XMAX + 0.85)
    axd.set_xticks([0.4, 0.8, 1.2, 1.6, 2.0])
    axd.tick_params(axis="y", length=0)
    axd.spines["left"].set_visible(False)
    axd.set_xlabel("Residual B$_{\\rm cell}$: scattered / contiguous")
    handles = [Line2D([0], [0], marker="o", color="none", mec=S.INK, mfc=S.INK, ms=4, label="BH q < 0.05"),
               Line2D([0], [0], marker="o", color="none", mec=S.INK, mfc="white", ms=4, label="BH q ≥ 0.05")]
    axd.legend(handles=handles, loc="upper right", fontsize=6.4, handletextpad=0.2, bbox_to_anchor=(1.02, 1.0))
    S.panel(axd, "d", x=-0.36, y=1.01)
    S.save(fig, OUT, "Fig6_parameters")


if __name__ == "__main__":
    main()
