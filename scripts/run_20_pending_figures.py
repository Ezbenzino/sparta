#!/usr/bin/env python
"""
run_20_pending_figures.py — 生成文稿待制图表（Fig 3/4/5 + Graphical abstract）
=============================================================================
从已有产物 JSON 出发，不重算任何数值。

Fig 3: 分子尺寸扫描（S3）+ 渗流/网格统计
Fig 4: 耦合与共享输入去除（逐切片配对图 + 散点）
Fig 5: S2 剂量反应 + S1 患者分层
Graphical abstract: 精简框架图
"""
from __future__ import annotations
import json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle
import matplotlib.patches as mpatches

from sparta.io_ import Paths, admitted_slides, load_config, patient_map

_cfg = load_config()
_P = Paths(_cfg)

ROOT = _P.root
FIG_DIR = _P.figure("")
VAL_DIR = _P.validation("")
CF_DIR = ROOT / "results" / "counterfactual"

# ── 色系（与 run_16 一致）──
CELL  = "#2a78d6"
MAB   = "#eb6834"
AQUA  = "#1baf7a"
INK   = "#0b0b0b"
INK2  = "#52514e"
MUTED = "#8f8e88"
GRID  = "#eeece6"
SURFACE = "#fcfcfb"
NEUTRAL = "#b9b7b0"

# ── 切片名单与患者映射：以台账为唯一事实来源 ──
ORDER = admitted_slides(_P)
PMAP = patient_map(_P)

# 患者色板（Fig 5 S1 分层用）——动态生成，适应队列变化
_PAT_PALETTE = [
    "#2a78d6", "#1baf7a", "#9b59b6", "#e74c3c", "#f39c12",
    "#1abc9c", "#eb6834", "#34495e", "#2ecc71", "#e67e22",
]
_ALL_PATIENTS = sorted(set(PMAP.values()))
PAT_COLORS = {p: _PAT_PALETTE[i % len(_PAT_PALETTE)]
              for i, p in enumerate(_ALL_PATIENTS)}

# ── 辅助 ──
def _load(p):
    with open(p, encoding="utf-8") as f:
        return json.load(f)

def _style(ax, ylab=None, xlab=None, grid_axis="y"):
    ax.set_facecolor(SURFACE)
    for sp in ("top","right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"):
        ax.spines[sp].set_color("#dedcd6")
        ax.spines[sp].set_linewidth(0.8)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=7.6, colors=INK2, length=3)
    if ylab: ax.set_ylabel(ylab, fontsize=8.6, color=INK2)
    if xlab: ax.set_xlabel(xlab, fontsize=8.6, color=INK2)

def _title(ax, letter, text, sub=None):
    ax.text(0, 1.12, f"{letter}  {text}", transform=ax.transAxes,
            fontsize=10.5, color=INK, va="bottom", ha="left")
    if sub:
        ax.text(0, 1.045, sub, transform=ax.transAxes, fontsize=8.2,
                color=INK2, va="bottom", ha="left")

def _cohort_ticks(ax, sids):
    ax.set_xticks(range(len(sids)))
    ax.set_xticklabels(sids, rotation=45, ha="right", fontsize=7.4)
    n_c = sum(1 for s in sids if s.startswith("CSCC"))
    if 0 < n_c < len(sids):
        ax.axvline(n_c - 0.5, color="#dedcd6", lw=1.0, ls=(0,(4,3)), zorder=0)

def _save(fig, name, dpi=400):
    for ext in ("png","pdf"):
        p = FIG_DIR / f"{name}.{ext}"
        fig.savefig(p, dpi=dpi, facecolor=SURFACE)
        print(f"  saved {p}")
    plt.close(fig)

def _mk(sid):
    if sid.startswith("CSCC"):
        return dict(marker="o", fc=CELL)
    return dict(marker="s", fc="none")


# ══════════════════════════════════════════════════════════════════════════
def fig3_size_scan():
    """Fig 3: 分子尺寸扫描（0.5–10 nm）+ 渗流/网格统计"""
    mesh_data = _load(VAL_DIR / "mesh_stats.json")
    mesh = mesh_data["per_slide"]
    cfs = {s: _load(CF_DIR / f"{s}.json") for s in ORDER}

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.6), facecolor=SURFACE)
    fig.subplots_adjust(left=0.075, right=0.97, top=0.74, bottom=0.22, wspace=0.28)

    # ---- (a) S3 分子尺寸扫描 ----
    a = ax[0]; _style(a, ylab=r"mean $B_{mAb}$ in tumour core",
                      xlab="hydrodynamic radius (nm)")
    for s in ORDER:
        s3 = cfs[s]["s3"]
        mk = _mk(s)
        a.plot(s3["radii_nm"], s3["mean_core"], color=MAB, lw=1.6, alpha=0.82,
               marker=mk["marker"], ms=4.5,
               markerfacecolor=mk["fc"] if mk["fc"]!="none" else SURFACE,
               markeredgecolor=MAB, markeredgewidth=1.2)
    a.axvline(5.5, color=MUTED, lw=1.0, ls=(0,(3,3)), zorder=1)
    a.text(5.65, a.get_ylim()[0]+0.3, "IgG\n5.5 nm", fontsize=7.6, color=INK2, va="bottom")
    a.text(0.02, 0.96, "same tissue, same graph —\nonly the molecule changes",
           transform=a.transAxes, fontsize=7.6, color=INK2, va="top")
    _title(a, "a", "Molecular size separates the two transport problems",
           "0.5 nm solute: nearly transparent   |   IgG: strongly obstructive")

    # ---- (b) 渗流/网格统计 ----
    b = ax[1]; _style(b, ylab="edges with mesh < IgG radius (%)")
    # Use the edge-level geometric exclusion fraction. `frac_crosslink` in
    # screen_decision.json is a variance-attribution share, a different metric.
    v = [mesh[s]["frac_size_excluded_pct"] for s in ORDER]
    colors = [MAB if s.startswith("CSCC") else "#f3a27f" for s in ORDER]
    b.bar(range(len(ORDER)), v, width=0.66, linewidth=1.0,
          edgecolor=SURFACE, color=colors)
    median_excluded = float(mesh_data["frac_excluded_median"]) * 100.0
    # Match the manuscript's half-up display at an exact 62.65% midpoint;
    # Python's built-in round() uses ties-to-even and serializes this as 62.6.
    median_label = f"{np.floor(median_excluded * 10 + 0.5 + 1e-9) / 10:.1f}"
    b.axhline(median_excluded, color=INK, lw=1.0, ls=(0,(4,3)), zorder=3)
    b.text(len(ORDER)-0.5, median_excluded + 1.2,
           f"median {median_label}%", ha="right", fontsize=7.4, color=INK)
    # 逐柱数字标签在 19 根柱子上必然重叠；median line 已表达主要信息，
    # 精确数字见 Supplementary Table S1。
    _cohort_ticks(b, ORDER)
    b.set_ylim(0, 75)
    _title(b, "b", "The antibody barrier is percolation-limited",
           "60–65% of edges exclude IgG; transport proceeds through the minority")
    # 图例
    b.legend(handles=[
        mpatches.Patch(color=MAB, label="primary cSCC"),
        mpatches.Patch(color="#f3a27f", label="metastatic melanoma"),
    ], loc="upper left", frameon=False, fontsize=7.0, ncol=1,
       borderaxespad=0.4, handlelength=1.0)

    fig.text(0.075, 0.955,
             "Figure 3  Molecular size scan and percolation statistics",
             fontsize=11.5, color=INK, ha="left", va="top")
    _save(fig, "fig3_size_scan")


# ══════════════════════════════════════════════════════════════════════════
def fig4_coupling_forest():
    """Fig 4: 耦合与共享输入去除（逐切片配对图 + 散点）"""
    dec = _load(VAL_DIR / "decoupling.json")
    chk = _load(VAL_DIR / "shared_ecm_check.json")["per_slide"]
    sids = [s for s in ORDER if s in dec and s in chk]
    chance = 0.0625

    fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.8), facecolor=SURFACE)
    fig.subplots_adjust(left=0.055, right=0.982, top=0.72, bottom=0.24, wspace=0.32)

    def _mk(s):
        return dict(marker="o", fc=CELL) if s.startswith("CSCC") else dict(marker="s", fc="none")

    # ---- (a) 偏相关散点 ----
    a = ax[0]; _style(a, ylab=r"partial $\rho$  ($B_{cell}$ vs $B_{mAb}$)")
    a.axhline(0, color=INK, lw=1.0, zorder=1)
    for i, s in enumerate(sids):
        r = dec[s]; m = _mk(s)
        a.scatter(i, r["rho_partial"], s=70, marker=m["marker"],
                  facecolors=m["fc"], edgecolors=CELL, linewidths=1.6, zorder=3)
        if r["p_partial"] >= 0.05:
            a.text(i, r["rho_partial"]+0.028, "n.s.", ha="center", fontsize=7, color=MUTED)
    _cohort_ticks(a, sids); a.set_ylim(-0.10, 0.46)
    n_pos = sum(1 for s in sids if dec[s]["rho_partial"] > 0)
    n_sig = sum(1 for s in sids if dec[s]["p_partial"] < 0.05)
    _title(a, "a", "Coupling, not dissociation",
           f"{n_pos}/{len(sids)} positive, {n_sig}/{len(sids)} significant")

    # ---- (b) 共享输入去除 — 配对相关系数图 ----
    b = ax[1]; _style(b, xlab=r"partial Spearman $\rho$", grid_axis="x")
    b.axvline(0, color=INK, lw=0.8, zorder=1)
    # 分层：Visium cSCC, 1st-gen ST cSCC, Melanoma. Each section is one row;
    # the x coordinate is the actual correlation, not its row index.
    strata = [
        ("cSCC (Visium)", [s for s in sids if s.startswith("CSCC") and int(s[4:])<=4]),
        ("cSCC (1st-gen ST)", [s for s in sids if s.startswith("CSCC") and int(s[4:])>=5]),
        ("Melanoma", [s for s in sids if s.startswith("MEL")]),
    ]
    y_pos = 0.0
    yticks_pos, yticks_lab = [], []
    xs_all = []
    for group_i, (sname, members) in enumerate(strata):
        for s in members:
            m0 = chk[s]["主配置"]
            m1 = chk[s]["切断共享ECM"]
            keep = m1["p_partial"] < 0.05 and m1["rho_partial"] > 0
            col = CELL if keep else MUTED
            rho_full = float(m0["rho_partial"])
            rho_removed = float(m1["rho_partial"])
            xs_all.extend((rho_full, rho_removed))
            b.plot([rho_full, rho_removed], [y_pos, y_pos],
                   color=col, lw=1.8 if keep else 1.0,
                   alpha=0.9 if keep else 0.55, zorder=2)
            b.scatter(rho_full, y_pos, s=38, marker="o", facecolors="none",
                      edgecolors=NEUTRAL, linewidths=1.2, zorder=3)
            b.scatter(rho_removed, y_pos, s=38, marker="o",
                      facecolors=col, edgecolors=col, linewidths=1.0, zorder=3)
            yticks_pos.append(y_pos); yticks_lab.append(s)
            y_pos += 1
        if members and group_i < len(strata) - 1:
            b.axhline(y_pos - 0.5, color="#dedcd6", lw=0.6, zorder=0)
            y_pos += 0.65
    b.set_yticks(yticks_pos); b.set_yticklabels(yticks_lab, fontsize=6.8)
    if xs_all:
        pad = 0.05
        b.set_xlim(min(-0.10, min(xs_all) - pad), max(0.45, max(xs_all) + pad))
    b.set_ylim(y_pos - 0.5, -0.5)
    b.tick_params(axis="y", length=0, pad=2)
    from matplotlib.lines import Line2D
    b.legend(handles=[
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="none",
               markeredgecolor=NEUTRAL, markersize=5, label="full model"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=CELL,
               markeredgecolor=CELL, markersize=5, label="shared matrix removed"),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.19), ncol=2,
       frameon=False, fontsize=6.6, columnspacing=0.8, handletextpad=0.35)
    # 统计标注
    cscc_all = [s for s in sids if s.startswith("CSCC")]
    cscc_surv = sum(1 for s in cscc_all
                    if chk[s]["切断共享ECM"]["p_partial"]<0.05
                    and chk[s]["切断共享ECM"]["rho_partial"]>0)
    mel_all = [s for s in sids if s.startswith("MEL")]
    mel_surv = sum(1 for s in mel_all
                   if chk[s]["切断共享ECM"]["p_partial"]<0.05
                   and chk[s]["切断共享ECM"]["rho_partial"]>0)
    _title(b, "b", "Does coupling persist after shared-input removal?",
           f"Significant positive: cSCC {cscc_surv}/{len(cscc_all)}; melanoma {mel_surv}/{len(mel_all)}")

    # ---- (c) 解离区 vs 随机期望 ----
    c = ax[2]; _style(c, ylab="dissociation-zone spots (%)")
    v = [dec[s]["frac_discordant_r"] * 100 for s in sids]
    c.bar(range(len(sids)), v, color=NEUTRAL, width=0.66,
          linewidth=1.6, edgecolor=SURFACE)
    c.axhline(chance*100, color=MAB, lw=2.0, zorder=3)
    c.text(len(sids)-0.4, chance*100+0.18,
           "expected under independence, 6.25%", fontsize=7.4, color=MAB, ha="right")
    for i, x in enumerate(v):
        c.text(i, x+0.12, f"{x:.1f}", ha="center", fontsize=7, color=INK)
    _cohort_ticks(c, sids); c.set_ylim(0, 7.6)
    n_below = sum(1 for x in v if x <= chance*100)
    _title(c, "c", "Fewer discordant spots than chance",
           f"{n_below}/{len(sids)} below the independence line")

    fig.text(0.055, 0.965,
             "Figure 4  The two barriers do not dissociate — "
             "coupling survives shared-input removal in the squamous cohort",
             fontsize=11.5, color=INK, ha="left", va="top")
    fig.text(0.055, 0.915,
             "● primary cSCC (6 patients, 15 sections)      "
             "■ metastatic melanoma (1 patient, 4 deposits)",
             fontsize=8, color=INK2, ha="left", va="top")
    _save(fig, "fig4_coupling_forest")


# ══════════════════════════════════════════════════════════════════════════
def fig5_counterfactuals():
    """Fig 5: S2 剂量反应 + S1 患者分层"""
    cfs = {s: _load(CF_DIR / f"{s}.json") for s in ORDER}

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.6), facecolor=SURFACE)
    fig.subplots_adjust(left=0.075, right=0.97, top=0.69, bottom=0.22, wspace=0.28)

    # ---- (a) S1 患者分层 ----
    a = ax[0]; _style(a, ylab="observed / permutation null")
    v = [cfs[s]["s1"]["fixed"].get("ratio_real_over_null", np.nan) for s in ORDER]
    sig = [cfs[s]["s1"]["fixed"]["p_emp"] < 0.05 for s in ORDER]
    colors = [PAT_COLORS[PMAP[s]] for s in ORDER]
    a.bar(range(len(ORDER)), v, width=0.66, linewidth=1.6, edgecolor=SURFACE,
          color=colors, alpha=0.85)
    a.axhline(1.0, color=INK, lw=1.2, zorder=3)
    vmax = float(np.nanmax(v)) if v else 1.0
    label_pad = max(vmax * 0.018, 0.04)
    for i, (x, k) in enumerate(zip(v, sig)):
        a.text(i, x + label_pad, f"{x:.2f}", ha="center", fontsize=6.8,
               color=INK if k else MUTED)
    _cohort_ticks(a, ORDER); a.set_ylim(0, max(vmax * 1.17, 1.2))
    # 患者图例
    seen = set()
    handles = []
    for s in ORDER:
        p = PMAP[s]
        if p not in seen:
            seen.add(p)
            handles.append(mpatches.Patch(color=PAT_COLORS[p], label=p))
    fig.legend(handles=handles, fontsize=6.4, frameon=False, ncol=min(7, len(handles)),
               loc="upper center", bbox_to_anchor=(0.51, 0.895),
               columnspacing=0.85, handletextpad=0.3)
    _title(a, "a", "S1  spatial rearrangement (by patient)",
           "composition held fixed; colours identify patients")

    # ---- (b) S2 剂量反应 ----
    b = ax[1]; _style(b, ylab="residual barrier ratio\n(scattered / contiguous)",
                      xlab="fraction of the blockade removed")
    for s in ORDER:
        pk = cfs[s]["s2"]["per_k"]
        xy = sorted((v_["k_over_cut"]*100, v_["ratio_vs_in_cut"]) for v_ in pk.values())
        xs_, ys_ = zip(*xy)
        b.plot(xs_, ys_, color=CELL, lw=1.6, alpha=0.82,
               marker="o" if s.startswith("CSCC") else "s", ms=4.5,
               markerfacecolor=CELL if s.startswith("CSCC") else SURFACE,
               markeredgecolor=CELL, markeredgewidth=1.2)
    b.axhline(1.0, color=INK, lw=1.2, zorder=3)
    b.axvline(20, color=MAB, lw=1.6, ls=(0,(4,3)), zorder=2)
    b.set_xlim(2.5, 32.5)
    b.set_xticks([5,10,20,30]); b.set_xticklabels(["5%","10%","20%","30%"])
    sig20 = sum(1 for s in ORDER
                if min(cfs[s]["s2"]["per_k"].values(),
                       key=lambda v_: abs(v_["k_over_cut"]*100-20)).get("p_vs_in_cut",1) < 0.05)
    mono = sum(1 for s in ORDER
               if (lambda pk:
                   (lambda s2: s2[-1]>s2[0])(
                       [v_["ratio_vs_in_cut"] for v_ in sorted(pk.values(),
                           key=lambda v_: v_["k_over_cut"])]))(
                       cfs[s]["s2"]["per_k"]))
    _title(b, "b", "S2  blockade continuity",
           f"{sig20}/{len(ORDER)} significant at pre-specified 20%; effect grows in {mono}/{len(ORDER)}")

    fig.text(0.075, 0.955,
             "Figure 5  Counterfactual experiments: "
             "spatial rearrangement (S1) and blockade continuity (S2)",
             fontsize=11.5, color=INK, ha="left", va="top")
    fig.text(0.97, 0.958, "● primary cSCC      ■ metastatic melanoma",
             fontsize=8, color=INK2, ha="right", va="top")
    _save(fig, "fig5_counterfactual")


# ══════════════════════════════════════════════════════════════════════════
def graphical_abstract():
    """Delegate to the canonical graphics-led builder.

    The landscape version that used to live here was superseded on 2026-10-06 by
    scripts/is_figures/graphical_abstract.py (13.0 x 5.2 cm, 2.5:1, graphics-led;
    writes results/figures/graphical_abstract.{pdf,png,tif}). Keeping a second
    implementation here would silently overwrite that file whenever run_20 is re-run.
    """
    import importlib.util

    path = Path(__file__).resolve().parent / 'is_figures' / 'graphical_abstract.py'
    spec = importlib.util.spec_from_file_location('_sparta_graphical_abstract', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.main()


# ══════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    print("[Fig 3] molecular size scan + percolation stats")
    fig3_size_scan()
    print("[Fig 4] paired coupling plot")
    fig4_coupling_forest()
    print("[Fig 5] counterfactuals S1+S2")
    fig5_counterfactuals()
    print("[Graphical abstract]")
    graphical_abstract()
    print("Done.")
