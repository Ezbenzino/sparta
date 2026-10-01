#!/usr/bin/env python
"""
run_20_pending_figures.py — 生成文稿待制图表（Fig 3/4/5 + Graphical abstract）
=============================================================================
从已有产物 JSON 出发，不重算任何数值。

Fig 3: 分子尺寸扫描（S3）+ 渗流/网格统计
Fig 4: 耦合与共享输入去除（森林图 + 散点）
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
    scr = _load(VAL_DIR / "screen_decision.json")["per_slide"]
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
    v = [scr[s]["frac_crosslink"] * 100 for s in ORDER]
    bars = b.bar(range(len(ORDER)), v, width=0.66, linewidth=1.6,
                 edgecolor=SURFACE, color=[MAB]*len(ORDER))
    b.axhline(62.7, color=INK, lw=1.0, ls=(0,(4,3)), zorder=3)
    b.text(len(ORDER)-0.5, 63.0, "median 62.7%", ha="right", fontsize=7.4, color=INK)
    for i, x in enumerate(v):
        b.text(i, x+0.6, f"{x:.1f}", ha="center", fontsize=6.8, color=INK)
    _cohort_ticks(b, ORDER)
    b.set_ylim(0, 75)
    _title(b, "b", "The antibody barrier is percolation-limited",
           "60–65% of edges exclude IgG; transport proceeds through the minority")
    # 图例
    b.text(0.02, 0.96, "● primary cSCC      ■ metastatic melanoma",
           transform=b.transAxes, fontsize=7.4, color=INK2, va="top")

    fig.text(0.075, 0.955,
             "Figure 3  Molecular size scan and percolation statistics",
             fontsize=11.5, color=INK, ha="left", va="top")
    _save(fig, "fig3_size_scan")


# ══════════════════════════════════════════════════════════════════════════
def fig4_coupling_forest():
    """Fig 4: 耦合与共享输入去除（分层森林图 + 散点）"""
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

    # ---- (b) 共享输入去除 — 分层森林图 ----
    b = ax[1]; _style(b, ylab=r"partial $\rho$", grid_axis="y")
    b.axhline(0, color=INK, lw=0.8, zorder=1)
    # 分层：Visium cSCC, 1st-gen ST cSCC, Melanoma
    strata = [
        ("cSCC (Visium)", [s for s in sids if s.startswith("CSCC") and int(s[4:])<=4]),
        ("cSCC (1st-gen ST)", [s for s in sids if s.startswith("CSCC") and int(s[4:])>=5]),
        ("Melanoma", [s for s in sids if s.startswith("MEL")]),
    ]
    y_pos = 0
    yticks_pos, yticks_lab = [], []
    for sname, members in strata:
        for s in members:
            m0 = chk[s]["主配置"]
            m1 = chk[s]["切断共享ECM"]
            keep = m1["p_partial"] < 0.05 and m1["rho_partial"] > 0
            col = CELL if keep else MUTED
            b.plot([0, 1], [m0["rho_partial"], m1["rho_partial"]],
                   color=col, lw=1.8 if keep else 1.0,
                   alpha=0.9 if keep else 0.5, zorder=2)
            mk = _mk(s)
            for x, v in ((0, m0["rho_partial"]), (1, m1["rho_partial"])):
                b.scatter(x, v, s=50, marker=mk["marker"],
                          facecolors=col if mk["fc"]!="none" else "none",
                          edgecolors=col, linewidths=1.3, zorder=3)
            yticks_pos.append(y_pos); yticks_lab.append(s)
            y_pos += 1
        if members:
            b.axhline(y_pos, color="#dedcd6", lw=0.6, zorder=0)
            y_pos += 0.5
    b.set_yticks(yticks_pos); b.set_yticklabels(yticks_lab, fontsize=6.8)
    b.set_xlim(-0.3, 1.4); b.set_xticks([0, 1])
    b.set_xticklabels(["full\nmodel", "shared matrix\nremoved"], fontsize=7.6)
    b.grid(axis="x", visible=False)
    b.invert_yaxis()
    # 统计标注
    cscc_all = [s for s in sids if s.startswith("CSCC")]
    cscc_surv = sum(1 for s in cscc_all
                    if chk[s]["切断共享ECM"]["p_partial"]<0.05
                    and chk[s]["切断共享ECM"]["rho_partial"]>0)
    mel_all = [s for s in sids if s.startswith("MEL")]
    mel_surv = sum(1 for s in mel_all
                   if chk[s]["切断共享ECM"]["p_partial"]<0.05
                   and chk[s]["切断共享ECM"]["rho_partial"]>0)
    _title(b, "b", "Is the coupling tissue-borne?",
           f"cSCC {cscc_surv}/{len(cscc_all)} survive   |   melanoma {mel_surv}/{len(mel_all)}")

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
    fig.subplots_adjust(left=0.075, right=0.97, top=0.74, bottom=0.22, wspace=0.28)

    # ---- (a) S1 患者分层 ----
    a = ax[0]; _style(a, ylab="observed / permutation null")
    v = [cfs[s]["s1"]["fixed"].get("ratio_real_over_null", np.nan) for s in ORDER]
    sig = [cfs[s]["s1"]["fixed"]["p_emp"] < 0.05 for s in ORDER]
    colors = [PAT_COLORS[PMAP[s]] for s in ORDER]
    a.bar(range(len(ORDER)), v, width=0.66, linewidth=1.6, edgecolor=SURFACE,
          color=colors, alpha=0.85)
    a.axhline(1.0, color=INK, lw=1.2, zorder=3)
    for i, (x, k) in enumerate(zip(v, sig)):
        a.text(i, x+0.06, f"{x:.2f}", ha="center", fontsize=6.8,
               color=INK if k else MUTED)
    _cohort_ticks(a, ORDER); a.set_ylim(0, 4.0)
    # 患者图例
    seen = set()
    handles = []
    for s in ORDER:
        p = PMAP[s]
        if p not in seen:
            seen.add(p)
            handles.append(mpatches.Patch(color=PAT_COLORS[p], label=p))
    a.legend(handles=handles, fontsize=6.5, frameon=False, ncol=2,
             loc="upper right", columnspacing=0.8, handletextpad=0.3)
    _title(a, "a", "S1  spatial rearrangement (by patient)",
           "composition held fixed, positions permuted")

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
    b.text(20.7, b.get_ylim()[0]+0.012, "pre-specified\nprimary test",
           fontsize=7.2, color=MAB, va="bottom")
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
    b.text(0.02, 0.96, f"{sig20}/{len(ORDER)} significant at 20%;\n"
           f"effect grows with removal in {mono}/{len(ORDER)}",
           transform=b.transAxes, fontsize=7.4, color=INK2, va="top")
    _title(b, "b", "S2  blockade continuity",
           "contiguous gap vs same material scattered")

    fig.text(0.075, 0.955,
             "Figure 5  Counterfactual experiments: "
             "spatial rearrangement (S1) and blockade continuity (S2)",
             fontsize=11.5, color=INK, ha="left", va="top")
    fig.text(0.97, 0.958, "● primary cSCC      ■ metastatic melanoma",
             fontsize=8, color=INK2, ha="right", va="top")
    _save(fig, "fig5_counterfactual")


# ══════════════════════════════════════════════════════════════════════════
def graphical_abstract():
    """Graphical abstract: 精简框架图"""
    fig, ax = plt.subplots(figsize=(10, 4.0), facecolor=SURFACE)
    ax.set_xlim(0, 10); ax.set_ylim(0, 4)
    ax.set_aspect("equal"); ax.axis("off")

    # ── 左：空间图 ──
    # 组织斑点
    rng = np.random.RandomState(42)
    pts = rng.randn(40, 2) * 0.6 + [1.8, 2.0]
    pts = pts[(pts[:,0]>0.8) & (pts[:,0]<2.8) & (pts[:,1]>0.8) & (pts[:,1]<3.2)]
    # 连线
    for i in range(len(pts)):
        for j in range(i+1, len(pts)):
            d = np.linalg.norm(pts[i]-pts[j])
            if d < 0.45:
                ax.plot([pts[i,0],pts[j,0]], [pts[i,1],pts[j,1]],
                        color=GRID, lw=0.5, zorder=1)
    # 斑点着色
    for i, (x,y) in enumerate(pts):
        if x < 1.5 and y > 2.0:
            c = CELL  # 源 (内皮)
        elif x > 2.2:
            c = MAB  # 汇 (肿瘤)
        else:
            c = NEUTRAL  # 基质
        ax.scatter(x, y, s=80, c=c, edgecolors=SURFACE, linewidths=0.8, zorder=2)
    ax.text(1.8, 3.45, "spatial graph", fontsize=9, color=INK, ha="center", weight="bold")
    ax.text(1.8, 0.45, "one substrate:\nextracellular matrix", fontsize=7.5, color=INK2, ha="center")

    # ── 中：两个算子 ──
    # B_cell 框
    box1 = FancyBboxPatch((4.2, 2.4), 2.0, 1.0, boxstyle="round,pad=0.1",
                          fc=SURFACE, ec=CELL, lw=2.0, zorder=2)
    ax.add_patch(box1)
    ax.text(5.2, 3.15, r"$B_{cell}$", fontsize=13, color=CELL, ha="center", weight="bold")
    ax.text(5.2, 2.75, "min-cut\nbarrier", fontsize=7.5, color=INK2, ha="center")
    # B_mAb 框
    box2 = FancyBboxPatch((4.2, 0.6), 2.0, 1.0, boxstyle="round,pad=0.1",
                          fc=SURFACE, ec=MAB, lw=2.0, zorder=2)
    ax.add_patch(box2)
    ax.text(5.2, 1.35, r"$B_{mAb}$", fontsize=13, color=MAB, ha="center", weight="bold")
    ax.text(5.2, 0.95, "screened Poisson\nsize exclusion", fontsize=7.5, color=INK2, ha="center")

    # 箭头
    ax.annotate("", xy=(4.2, 3.0), xytext=(2.9, 2.6),
                arrowprops=dict(arrowstyle="-|>", color=CELL, lw=1.5))
    ax.annotate("", xy=(4.2, 1.2), xytext=(2.9, 1.8),
                arrowprops=dict(arrowstyle="-|>", color=MAB, lw=1.5))

    # ── 右：结果 ──
    ax.text(7.8, 3.5, "findings", fontsize=9, color=INK, ha="center", weight="bold")
    findings = [
        ("percolation-limited", MAB),
        ("barriers co-located", CELL),
        ("deterministic, < 0.2 s", INK2),
    ]
    for i, (txt, c) in enumerate(findings):
        y = 3.0 - i * 0.55
        ax.scatter(7.2, y, s=30, c=c, edgecolors=SURFACE, linewidths=0.6, zorder=3)
        ax.text(7.4, y, txt, fontsize=7.8, color=c, va="center")

    ax.annotate("", xy=(7.0, 2.0), xytext=(6.2, 2.0),
                arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.2))

    fig.text(0.5, 0.97,
             "Two transport operators, one substrate: "
             "graph-transport modelling of cell and antibody delivery barriers",
             fontsize=10.5, color=INK, ha="center", va="top", style="italic")

    for ext in ("png","pdf"):
        p = FIG_DIR / f"graphical_abstract.{ext}"
        fig.savefig(p, dpi=400, facecolor=SURFACE, bbox_inches="tight")
        print(f"  saved {p}")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    print("[Fig 3] molecular size scan + percolation stats")
    fig3_size_scan()
    print("[Fig 4] coupling forest plot")
    fig4_coupling_forest()
    print("[Fig 5] counterfactuals S1+S2")
    fig5_counterfactuals()
    print("[Graphical abstract]")
    graphical_abstract()
    print("Done.")
