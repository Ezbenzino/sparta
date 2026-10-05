#!/usr/bin/env python
"""
run_16_figures.py —— Fig 2–5（全部从已有产物的汇总 JSON 出发，不重算）
==========================================================================

输入：results/validation/{screen_decision,bmab_sensitivity,decoupling,
      shared_ecm_check,benchmark_tools}.json + results/counterfactual/*.json
      + data/ledger.csv（取患者映射）
输出：results/figures/fig{2,3,4,5}_*.png / .pdf
上游模块：run_05 / run_06 / run_07 / run_12 / run_13 / run_14
下游模块：无（成文用）

设计约定（与 Fig 1 一致，全篇一套）
----------------------------------
· **色相承载"哪个屏障"**：cell = #2a78d6（蓝），antibody = #eb6834（橙）。
  这一对已过 CVD 校验（最差 ΔE 24.7，正常视觉 33.6）。
· **队列不用颜色区分**，用形状 + 分组位置 + 直接标注：
  cSCC = 实心圆，melanoma = 空心方。理由是颜色已经被"屏障模态"占用了，
  再拿它编码队列会让读者在两套语义之间来回猜。
· 第三个分类槽（抗原）用 aqua #1baf7a，仅出现在堆叠柱里，
  该色对浅底对比度 2.74:1 未达 3:1 —— 按 relief 规则一律带可见数值标注。
· 同一物理量在不同面板必须共用色标/轴范围，否则对比无意义。

⚠ 本脚本只画图，不做任何计算。任何数字对不上都应该回去查产物 JSON，
  不要在这里"修"。

用法
----
    python scripts/run_16_figures.py            # 出全部四张
    python scripts/run_16_figures.py --which 3  # 只出 Fig 3
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import (Paths, admitted_slides, load_config, load_json,  # noqa: E402
                        load_ledger, patient_map)

CELL, MAB, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8f8e88"
GRID, SURFACE = "#eeece6", "#fcfcfb"
NEUTRAL = "#b9b7b0"

# 出图顺序：cSCC 在前、黑色素瘤在后。名单在 main() 里从台账取
# （status == ingested），不再硬编码——队列一变就漏切片，而且不报错。
# 手改这个列表能对一时，对不了下一次；见 sparta/io_.py::admitted_slides
ORDER: list[str] = []


def _platform_map(P) -> dict:
    """切片 ID -> 平台。Fig 5 要按平台世代分开看，不能把两代混成一条中位线。"""
    return {r["slide_id"]: str(r.get("platform", "")) for r in load_ledger(P)}


def _bh(pvals):
    """Benjamini–Hochberg，与 run_12_paper_stats.py 同口径（每张切片一个检验）。

    图上写的显著数必须和正文一致。这里如果用未校正的 p，Fig 4b 会显示 17/19
    而 R4-S2 写的是 16/19——同一份数据两个数字，审稿人一对就出问题。
    """
    p = np.asarray(pvals, float); n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    adj = np.empty(n, float)
    adj[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


def _style(ax, ylab=None, xlab=None, grid_axis="y"):
    ax.set_facecolor(SURFACE)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color("#dedcd6"); ax.spines[sp].set_linewidth(0.8)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=7.6, colors=INK2, length=3)
    if ylab:
        ax.set_ylabel(ylab, fontsize=8.6, color=INK2)
    if xlab:
        ax.set_xlabel(xlab, fontsize=8.6, color=INK2)


def _title(ax, letter, text, sub=None):
    ax.text(0, 1.10, f"{letter}  {text}", transform=ax.transAxes,
            fontsize=10.5, color=INK, va="bottom", ha="left")
    if sub:
        ax.text(0, 1.035, sub, transform=ax.transAxes, fontsize=8.2,
                color=INK2, va="bottom", ha="left")


def _cohort_ticks(ax, sids):
    """x 轴刻度 + 用一条分隔线把两个队列隔开（队列靠位置区分，不靠颜色）。"""
    ax.set_xticks(range(len(sids)))
    ax.set_xticklabels(sids, rotation=45, ha="right", fontsize=7.4)
    n_c = sum(1 for s in sids if s.startswith("CSCC"))
    if 0 < n_c < len(sids):
        ax.axvline(n_c - 0.5, color="#dedcd6", lw=1.0, ls=(0, (4, 3)), zorder=0)


def _spread(vals, min_gap):
    """把一列会重叠的标签 y 坐标推开，保持原有上下顺序。

    直接把 8 个 slide 名标在线的末端一定会撞（CSCC02/03 差 0.006）。
    这里按值排序后逐个下推，只动标签位置，不动数据点。
    """
    order = sorted(range(len(vals)), key=lambda i: -vals[i])
    out = list(vals)
    for k, i in enumerate(order):
        if k == 0:
            continue
        prev = out[order[k - 1]]
        if prev - out[i] < min_gap:
            out[i] = prev - min_gap
    return out


def _save(fig, P, name, dpi):
    for ext in ("png", "pdf"):
        p = P.figure(f"{name}.{ext}")
        fig.savefig(p, dpi=dpi, facecolor=SURFACE)
        print(f"  已保存 {p}")


# ==========================================================================
def fig2(P, dpi):
    """驱动分解，以及它对 beta 的依赖。"""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    scr = load_json(P.validation("screen_decision.json"))["per_slide"]
    sens = load_json(P.validation("bmab_sensitivity.json"))
    sids = [s for s in ORDER if s in scr]

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.3), facecolor=SURFACE)
    fig.subplots_adjust(left=0.075, right=0.975, top=0.76, bottom=0.235, wspace=0.26)

    # ---- (a) 三个通道的份额（对数轴点图）----
    # 不用堆叠柱：三者跨越 0.5%–98%，两个数量级，堆叠柱里共享 ECM 与抗原
    # 会被压成看不见的一条线，等于把"抗原到底有多弱"这个信息丢掉。
    a = ax[0]; _style(a, xlab="share of B$_{mAb}$ variance (%)", grid_axis="x")
    ypos = np.arange(len(sids))[::-1]
    for key, lab, c, mk in (("frac_crosslink", "crosslinking / size exclusion", MAB, "o"),
                            ("frac_shared_ecm", "shared matrix (blocks both)", NEUTRAL, "D"),
                            ("frac_antigen", "target antigen", AQUA, "^")):
        v = np.array([scr[s_][key] * 100 for s_ in sids])
        a.scatter(v, ypos, s=62, marker=mk, color=c, edgecolors=SURFACE,
                  linewidths=1.2, label=lab, zorder=3)
    # CSCC10/14/15/16 的 var_antigen 恰好为 0——CD274 与 PDCD1LG2 都没过最小匹配
    # 基因数，签名被中性填充。**那是缺数据，不是测到接近零。** 对数轴画不了 0，
    # 若不显式标出来，这四张会安静地少一个绿三角，读者会当成"抗原份额极小"。
    X0 = 0.24                                    # 左边缘的"未检出"位置
    n_missing = 0
    for i, s_ in enumerate(sids):
        fa = scr[s_]["frac_antigen"] * 100
        missing = scr[s_].get("var_antigen", 1.0) == 0 or fa <= 0
        x_left = X0 if missing else fa
        a.plot([x_left, scr[s_]["frac_crosslink"] * 100],
               [ypos[i]] * 2, color="#e2e0da", lw=1.4, zorder=1)
        if missing:
            n_missing += 1
            a.scatter([X0], [ypos[i]], s=62, marker="^", facecolors="none",
                      edgecolors=AQUA, linewidths=1.3, zorder=3)
            a.plot([X0 * 0.86, X0 * 1.16], [ypos[i] - 0.34, ypos[i] + 0.34],
                   color=AQUA, lw=1.1, zorder=4)
        else:
            a.text(fa * 0.55, ypos[i], f"{fa:.1f}", ha="right", va="center",
                   fontsize=6.8, color=AQUA)
    if n_missing:
        a.text(0.5, -0.315, f"⌀ = antigen signature not detected "
                            f"(CD274 / PDCD1LG2 below the match threshold; "
                            f"{n_missing} sections) — missing, not near-zero",
               transform=a.transAxes, fontsize=6.9, color=AQUA, ha="center", va="top")
    a.set_xscale("log"); a.set_xlim(0.2, 260)
    a.axvline(0.31, color="#e2e0da", lw=0.9, ls=(0, (2, 3)), zorder=0)
    a.set_xticks([0.5, 1, 5, 10, 50, 100])
    a.set_xticklabels(["0.5", "1", "5", "10", "50", "100"])
    a.set_yticks(ypos); a.set_yticklabels(sids, fontsize=7.6)
    n_c = sum(1 for s_ in sids if s_.startswith("CSCC"))
    a.axhline(ypos[n_c - 1] - 0.5, color="#dedcd6", lw=1.0, ls=(0, (4, 3)), zorder=0)
    # 图例放到轴下方一行——放在轴内右下角会压住 MEL03/MEL04 的橙点
    a.legend(fontsize=7.2, frameon=False, ncol=3, handletextpad=0.4,
             columnspacing=1.1, loc="upper center", bbox_to_anchor=(0.5, -0.20))
    _title(a, "a", "Where B$_{mAb}$ variance comes from",
           "log axis — the three channels span two decades")

    # ---- (b) 对 beta 的依赖（多切片叠加）----
    b = ax[1]; _style(b, ylab="crosslinking share (%)",
                      xlab=r"$\beta$  (mesh-contraction steepness)")
    slides_sens = sens.get("slides", {})   # run_17 多切片格式
    sids_sens = [s for s in ORDER if s in slides_sens]
    if sids_sens:
        for sid in sids_sens:
            g = slides_sens[sid]["lam_x_beta"]["grid"]
            xs, ys = g["xs"], g["ys"]
            M = np.asarray(slides_sens[sid]["lam_x_beta"]["crosslink_pct"], float)
            j = xs.index(3.0)
            b.plot(ys, M[:, j], color=MAB, lw=1.8, alpha=0.9,
                   marker="o" if sid.startswith("CSCC") else "s", ms=4.5,
                   markerfacecolor=MAB if sid.startswith("CSCC") else SURFACE,
                   markeredgecolor=MAB, markeredgewidth=1.2)
        # 五条曲线在 beta=12 端点接近，图例映射切片身份，避免线尾标签互相覆盖。
        handles = [Line2D(
            [0], [0], color=MAB, lw=1.6,
            marker="o" if sid.startswith("CSCC") else "s",
            markerfacecolor=MAB if sid.startswith("CSCC") else SURFACE,
            markeredgecolor=MAB, markersize=4.5, label=sid,
        ) for sid in sids_sens]
        b.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.025, 0.98),
                 ncol=2, frameon=True, framealpha=0.88, facecolor="white",
                 edgecolor="none", fontsize=6.8, handlelength=1.1,
                 columnspacing=0.8, handletextpad=0.4, borderpad=0.35,
                 labelspacing=0.25)
        b.axvline(3.0, color=MUTED, lw=0.9, ls=(0, (3, 3)), zorder=1)
        b.set_xscale("log"); b.set_xticks(ys)
        b.set_xticklabels([f"{y:g}" for y in ys], fontsize=7.4)
        b.set_ylim(-4, 104)
        _title(b, "b", "…and that share is set by one scale parameter",
               r"$\beta$=3 marked; the share is a function of $\beta$ in every section")
    else:  # 旧单切片格式回退
        g = sens["grids"]["lam_x_beta"]
        xs, ys, M = g["grid"]["xs"], g["grid"]["ys"], np.asarray(g["crosslink_pct"], float)
        j = xs.index(3.0)
        b.plot(ys, M[:, j], color=MAB, lw=2.0, marker="o", ms=5.5,
               markeredgecolor=SURFACE, markeredgewidth=1.2, zorder=3)
        b.axvline(3.0, color=MUTED, lw=0.9, ls=(0, (3, 3)), zorder=1)
        b.annotate(f"default $\\beta$=3\n{M[ys.index(3.0), j]:.1f}%",
                   xy=(3.0, M[ys.index(3.0), j]), xytext=(10, -34),
                   textcoords="offset points", fontsize=8.2, color=INK,
                   arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
        for yv, lab in ((M[ys.index(1.0), j], r"$\beta$=1"),
                        (M[ys.index(0.25), j], r"$\beta$=0.25")):
            b.annotate(f"{lab}: {yv:.1f}%", xy=(ys[list(M[:, j]).index(yv)], yv),
                       xytext=(6, 8), textcoords="offset points", fontsize=7.6,
                       color=INK2)
        b.set_xscale("log"); b.set_xticks(ys)
        b.set_xticklabels([f"{y:g}" for y in ys], fontsize=7.4)
        b.set_ylim(-4, 104)
        _title(b, "b", "…and that share is set by one scale parameter",
               r"$\lambda$ = 3 held fixed; MEL01")

    fig.text(0.065, 0.955, "Figure 2  The antibody barrier decomposes onto one axis — "
             "whose dominance is a modelling choice, not a measurement",
             fontsize=11.5, color=INK, ha="left", va="top")
    _save(fig, P, "fig2_driver_decomposition", dpi); plt.close(fig)


# ==========================================================================
def fig3(P, dpi):
    """耦合的方向、强度，以及切断共享输入之后还剩多少。"""
    import matplotlib.pyplot as plt

    dec = load_json(P.validation("decoupling.json"))
    chk = load_json(P.validation("shared_ecm_check.json"))["per_slide"]
    pmap = patient_map(P)
    sids = [s for s in ORDER if s in dec]
    chance = 0.0625

    fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.4), facecolor=SURFACE)
    fig.subplots_adjust(left=0.055, right=0.982, top=0.74, bottom=0.235, wspace=0.30)

    def _mk(s):
        return dict(marker="o", facecolor=CELL) if s.startswith("CSCC") \
            else dict(marker="s", facecolor="none")

    # ---- (a) 偏相关 ----
    a = ax[0]; _style(a, ylab=r"partial $\rho$  (B$_{cell}$ vs B$_{mAb}$)")
    a.axhline(0, color=INK, lw=1.0, zorder=1)
    for i, s in enumerate(sids):
        r = dec[s]; m = _mk(s)
        a.scatter(i, r["rho_partial"], s=70, marker=m["marker"],
                  facecolors=m["facecolor"], edgecolors=CELL, linewidths=1.6, zorder=3)
        if r["p_partial"] >= 0.05:
            a.text(i, r["rho_partial"] + 0.028, "n.s.", ha="center", fontsize=7,
                   color=MUTED)
    _cohort_ticks(a, sids); a.set_ylim(-0.10, 0.46)
    a.text(0.98, 0.13, "positive = the two barriers rise together\nafter vessel distance is removed",
           transform=a.transAxes, fontsize=7.6, color=INK2, va="bottom", ha="right")
    n_pos = sum(1 for s in sids if dec[s]["rho_partial"] > 0)
    n_sig = sum(1 for s in sids if dec[s]["p_partial"] < 0.05)
    _title(a, "a", "Coupling, not dissociation",
           f"{n_pos}/{len(sids)} sections positive, {n_sig}/{len(sids)} significant")

    # ---- (b) 切断共享输入前后 ----
    b = ax[1]; _style(b, ylab=r"partial $\rho$", grid_axis="y")
    b.axhline(0, color=INK, lw=1.0, zorder=1)
    lab_y, lab_txt = [], []
    for i, s in enumerate(sids):
        m0, m1 = chk[s]["主配置"], chk[s]["切断共享ECM"]
        keep = m1["p_partial"] < 0.05 and m1["rho_partial"] > 0
        col = CELL if keep else MUTED
        b.plot([0, 1], [m0["rho_partial"], m1["rho_partial"]], color=col,
               lw=2.0 if keep else 1.2, alpha=1.0 if keep else 0.75, zorder=2)
        mk = _mk(s)
        for x, v in ((0, m0["rho_partial"]), (1, m1["rho_partial"])):
            b.scatter(x, v, s=58, marker=mk["marker"],
                      facecolors=col if mk["facecolor"] != "none" else "none",
                      edgecolors=col, linewidths=1.5, zorder=3)
        lab_y.append(m1["rho_partial"]); lab_txt.append((s, col, keep))
    # 选择性直标。19 条线的终点挤在 0.05–0.32 之间，逐条标名字必然叠成一团
    # （_spread 也铺不开）。只标真正需要按名字找的那几张：
    # 塌掉的（灰色，是例外，读者要知道是谁）+ 蓝色束的最高与最低。
    keep_y = [y for y, (_, _, k) in zip(lab_y, lab_txt) if k]
    hi, lo = (max(keep_y), min(keep_y)) if keep_y else (None, None)
    sel = [i for i, (y, (_, _, k)) in enumerate(zip(lab_y, lab_txt))
           if (not k) or y == hi or y == lo]
    gap = 0.034 if len(sel) <= 8 else 0.018
    ys = _spread([lab_y[i] for i in sel], gap)
    for y_, i in zip(ys, sel):
        txt, col, _k = lab_txt[i]
        if abs(y_ - lab_y[i]) > 1e-9:
            b.plot([1.02, 1.045], [lab_y[i], y_], color=col, lw=0.7, alpha=0.7, zorder=2)
        b.text(1.06, y_, txt, fontsize=6.8, va="center", color=col)
    b.text(1.06, (hi + lo) / 2 if keep_y else 0.2,
           f"+{len(sids) - len(sel)} more\n(all retained)", fontsize=6.4,
           va="center", color=MUTED)
    b.set_xlim(-0.18, 1.55); b.set_xticks([0, 1])
    b.set_xticklabels(["full model", "shared matrix\nterm removed"], fontsize=7.8)
    b.grid(axis="x", visible=False)
    b.text(0.5, -0.30, "solid = still significant and positive      grey = collapses",
           transform=b.transAxes, fontsize=7.4, color=INK2, ha="center")
    def _survive(prefix):
        sub = [s for s in sids if s.startswith(prefix)]
        if not sub:
            return 0, 0
        k = sum(1 for s in sub if chk[s]["切断共享ECM"]["p_partial"] < 0.05
                and chk[s]["切断共享ECM"]["rho_partial"] > 0)
        return k, len(sub)
    cscc_k, cscc_n = _survive("CSCC"); mel_k, mel_n = _survive("MEL")
    _title(b, "b", "Is the coupling tissue-borne?",
           f"cSCC {cscc_k}/{cscc_n} survive, melanoma {mel_k}/{mel_n}")

    # ---- (c) 解离区 vs 随机期望 ----
    c = ax[2]; _style(c, ylab="dissociation-zone spots (%)")
    v = np.array([dec[s]["frac_discordant_r"] * 100 for s in sids])
    c.bar(range(len(sids)), v, color=NEUTRAL, width=0.66,
          linewidth=1.6, edgecolor=SURFACE)
    c.axhline(chance * 100, color=MAB, lw=2.0, zorder=3)
    c.text(len(sids) - 0.4, chance * 100 + 0.18,
           "expected under independence, 6.25 %", fontsize=7.6, color=MAB, ha="right")
    for i, x in enumerate(v):
        c.text(i, x + 0.12, f"{x:.1f}", ha="center", fontsize=7, color=INK)
    _cohort_ticks(c, sids); c.set_ylim(0, 7.6)
    n_below = int(np.sum(v <= chance * 100))
    _title(c, "c", "Fewer discordant spots than chance",
           f"{n_below}/{len(sids)} below the independence line")

    fig.text(0.055, 0.968, "Figure 3  The two barriers do not dissociate, and in the squamous "
             "cohort the coupling survives removal of every shared input",
             fontsize=11.5, color=INK, ha="left", va="top")
    cscc_slides = [s for s in sids if s.startswith("CSCC")]
    mel_slides = [s for s in sids if s.startswith("MEL")]
    cscc_pats = len({pmap[s] for s in cscc_slides})
    mel_pats = len({pmap[s] for s in mel_slides})
    fig.text(0.055, 0.915,
             f"● primary cSCC ({cscc_pats} patients, {len(cscc_slides)} sections)      "
             f"■ metastatic melanoma ({mel_pats} patient, {len(mel_slides)} deposits)",
             fontsize=8, color=INK2, ha="left", va="top")
    _save(fig, P, "fig3_coupling", dpi); plt.close(fig)


# ==========================================================================
def fig4(P, dpi):
    """三个反事实实验。"""
    import matplotlib.pyplot as plt

    cf = {s: load_json(P.counterfactual(s)) for s in ORDER
          if P.counterfactual(s).exists()}
    sids = [s for s in ORDER if s in cf]

    fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.4), facecolor=SURFACE)
    fig.subplots_adjust(left=0.055, right=0.982, top=0.74, bottom=0.235, wspace=0.30)

    # ---- (a) S1 ----
    a = ax[0]; _style(a, ylab="observed / permutation null")
    v = [cf[s]["s1"]["fixed"].get("ratio_real_over_null", np.nan) for s in sids]
    sig = [cf[s]["s1"]["fixed"]["p_emp"] < 0.05 for s in sids]
    a.bar(range(len(sids)), v, width=0.66, linewidth=1.6, edgecolor=SURFACE,
          color=[CELL if k else NEUTRAL for k in sig])
    a.axhline(1.0, color=INK, lw=1.2, zorder=3)
    # 选择性直标：只标显著的那些。每根柱子都标数会在 1.0–1.2 那一带糊成一片，
    # 而不显著的柱子恰恰不需要精确读数——它们的信息就是"没过线"。
    for i, (x, k) in enumerate(zip(v, sig)):
        if k:
            a.text(i, x + 0.06, f"{x:.2f}", ha="center", fontsize=7, color=INK)
    # 上限跟着数据走：写死 3.9 会把 CSCC07（4.49）的柱子和数字一起切掉
    _cohort_ticks(a, sids); a.set_ylim(0, max([x for x in v if np.isfinite(x)] or [1]) * 1.18)
    a.text(0.98, 0.94, f"blue = significant after BH ({sum(sig)}/{len(sids)})",
           transform=a.transAxes, fontsize=7.6, color=INK2, va="top", ha="right")
    _title(a, "a", "S1  spatial rearrangement", "composition held fixed, positions permuted")

    # ---- (b) S2 剂量-反应 ----
    b = ax[1]; _style(b, ylab="residual barrier ratio\n(scattered / contiguous)",
                      xlab="fraction of the blockade removed")
    for s in sids:
        pk = cf[s]["s2"]["per_k"]
        xy = sorted((v_["k_over_cut"] * 100, v_["ratio_vs_in_cut"]) for v_ in pk.values())
        xs_, ys_ = zip(*xy)
        b.plot(xs_, ys_, color=CELL, lw=1.6, alpha=0.85,
               marker="o" if s.startswith("CSCC") else "s", ms=4.5,
               markerfacecolor=CELL if s.startswith("CSCC") else SURFACE,
               markeredgecolor=CELL, markeredgewidth=1.2)
    b.axhline(1.0, color=INK, lw=1.2, zorder=3)
    b.axvline(20, color=MAB, lw=1.6, ls=(0, (4, 3)), zorder=2)
    b.text(20.7, b.get_ylim()[0] + 0.012, "pre-specified\nprimary test", fontsize=7.4,
           color=MAB, va="bottom")
    b.set_xticks([5, 10, 20, 30]); b.set_xticklabels(["5 %", "10 %", "20 %", "30 %"])
    p20, mono = [], 0
    for s in sids:
        pk = cf[s]["s2"]["per_k"]
        k20 = min(pk.values(), key=lambda v_: abs(v_["k_over_cut"] * 100 - 20))
        p20.append(float(k20.get("p_vs_in_cut", 1.0)))
        series = [v_["ratio_vs_in_cut"] for v_ in sorted(pk.values(),
                                                         key=lambda v_: v_["k_over_cut"])]
        if len(series) >= 2 and series[-1] > series[0]:
            mono += 1
    # 用 BH 校正后的 p 计数，否则图上是 17/19 而正文是 16/19，审稿人一对就出问题
    sig20 = int((_bh(p20) < 0.05).sum())
    b.text(0.02, 0.96, f"{sig20}/{len(sids)} significant at 20 % (BH);\n"
           f"effect grows with removal fraction in {mono}/{len(sids)}",
           transform=b.transAxes, fontsize=7.6, color=INK2, va="top")
    _title(b, "b", "S2  blockade continuity", "contiguous gap vs same material scattered")

    # ---- (c) S3 ----
    c = ax[2]; _style(c, ylab="mean B$_{mAb}$ in tumour core",
                      xlab="hydrodynamic radius (nm)")
    for s in sids:
        s3 = cf[s]["s3"]
        c.plot(s3["radii_nm"], s3["mean_core"], color=MAB, lw=1.6, alpha=0.85,
               marker="o" if s.startswith("CSCC") else "s", ms=4.5,
               markerfacecolor=MAB if s.startswith("CSCC") else SURFACE,
               markeredgecolor=MAB, markeredgewidth=1.2)
    c.axvline(5.5, color=MUTED, lw=0.9, ls=(0, (3, 3)), zorder=1)
    c.text(5.7, c.get_ylim()[0] + 0.5, "IgG", fontsize=8, color=INK2)
    c.text(0.02, 0.96, "same tissue, same graph —\nonly the molecule changes",
           transform=c.transAxes, fontsize=7.6, color=INK2, va="top")
    _title(c, "c", "S3  molecular size", "monotone in every section, no ceiling")

    fig.text(0.055, 0.955, "Figure 4  Three counterfactuals: arrangement, continuity, and molecular size",
             fontsize=11.5, color=INK, ha="left", va="top")
    fig.text(0.982, 0.958, "● primary cSCC      ■ metastatic melanoma",
             fontsize=8, color=INK2, ha="right", va="top")
    _save(fig, P, "fig4_counterfactuals", dpi); plt.close(fig)


# ==========================================================================
def fig5(P, dpi):
    """与空间域方法的关系。"""
    import matplotlib.pyplot as plt

    bt = load_json(P.validation("benchmark_tools.json"))["per_slide"]
    sids = [s for s in ORDER if s in bt]

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.3), facecolor=SURFACE)
    fig.subplots_adjust(left=0.065, right=0.975, top=0.76, bottom=0.235, wspace=0.26)

    a = ax[0]; _style(a, ylab="fold enrichment of min-cut edges\non domain boundaries")
    v = [bt[s]["enrichment"] for s in sids]
    a.bar(range(len(sids)), v, color=CELL, width=0.66, linewidth=1.6, edgecolor=SURFACE)
    a.axhline(1.0, color=INK, lw=1.2, zorder=3)
    a.text(-0.42, 0.955, "no enrichment", fontsize=7.6, color=INK2, ha="left", va="top")
    # 选择性直标：分段中位线已经承担了主要信息，逐柱标数只会在 1.0–1.2 糊成一片。
    # 只标两个极值和所有低于 1（"域边界反而更少落在割上"）的切片——那几张才是要点。
    _fin = [x for x in v if np.isfinite(x)]
    _hi, _lo = (max(_fin), min(_fin)) if _fin else (None, None)
    for i, x in enumerate(v):
        if not np.isfinite(x):
            continue
        if x < 1.0 or x == _hi or x == _lo:
            a.text(i, x + 0.04, f"{x:.2f}", ha="center", fontsize=7,
                   color=MUTED if x < 1.0 else INK)
    # 按「队列 × 平台世代」的连续段各画一条中位线。
    # R6 的论点正是"富集随分辨率而变"，把两代混成一条中位线会把这件事抹掉；
    # 而 Visium 在 x 轴上被 cSCC/MEL 分成了两段，所以只能按**连续段**画，
    # 否则一条线会横跨整幅图、两个标签还会叠在一起。
    pf = _platform_map(P)
    blocks, cur = [], None
    for i, x in enumerate(sids):
        key = ("MEL" if x.upper().startswith("MEL") else "CSCC", pf.get(x, ""))
        if cur is None or key != cur[0]:
            cur = (key, [i]); blocks.append(cur)
        else:
            cur[1].append(i)
    for (_, _), idx in blocks[1:]:
        a.axvline(min(idx) - 0.5, color="#dedcd6", lw=1.0, ls=(0, (2, 3)), zorder=0)
    _cohort_ticks(a, sids); a.set_ylim(0, 2.25)

    def _med(idx):
        vv = [v[i] for i in idx if np.isfinite(v[i])]
        return float(np.median(vv)) if vv else float("nan")

    LBL = {"visium": "Visium", "legacy_st": "1st-gen ST"}
    for (coh, plat), idx in blocks:
        m = _med(idx)
        if not np.isfinite(m):
            continue
        a.plot([min(idx) - 0.42, max(idx) + 0.42], [m, m],
               color=MAB, lw=1.7, ls=(0, (5, 2)), zorder=4)
        a.text(np.mean(idx), 2.19, f"{LBL.get(plat, plat)}  {m:.2f}×",
               ha="center", va="top", fontsize=7.4, color=MAB)
    vis = [i for i, x in enumerate(sids) if pf.get(x) == "visium"]
    leg = [i for i, x in enumerate(sids) if pf.get(x) == "legacy_st"]
    m_vis, m_leg = _med(vis), _med(leg)
    _title(a, "a", "A standard segmentation sees part of it — at Visium resolution",
           f"orange dashes = median within each block · pooled {m_vis:.2f}× Visium, "
           f"{m_leg:.2f}× 1st-gen ST")

    b = ax[1]; _style(b, ylab="%")
    w, xs_ = 0.34, np.arange(len(sids))
    b.bar(xs_ - w / 2, [bt[s]["recall_of_cut_by_boundary"] * 100 for s in sids], w,
          color=CELL, linewidth=1.6, edgecolor=SURFACE,
          label="recall — cut edges on a domain boundary")
    b.bar(xs_ + w / 2, [bt[s]["precision_of_boundary_for_cut"] * 100 for s in sids], w,
          color=MAB, linewidth=1.6, edgecolor=SURFACE,
          label="precision — domain-boundary edges on the cut")
    for i, s in enumerate(sids):
        b.text(i + w / 2, bt[s]["precision_of_boundary_for_cut"] * 100 + 1.2,
               f"{bt[s]['precision_of_boundary_for_cut']*100:.0f}", ha="center",
               fontsize=7, color=INK)
    _cohort_ticks(b, sids); b.set_ylim(0, 66)
    b.legend(fontsize=7.4, frameon=False, loc="upper left")
    pr = [bt[s]["precision_of_boundary_for_cut"] * 100 for s in sids
          if np.isfinite(bt[s]["precision_of_boundary_for_cut"])]
    _title(b, "b", "…but it over-calls massively",
           f"only {min(pr):.0f}–{max(pr):.0f} % of domain boundaries limit transport")

    fig.text(0.065, 0.955, "Figure 5  Spatial-domain methods propose candidate boundaries; "
             "the min cut selects the one that limits flux and gives it a capacity",
             fontsize=11.5, color=INK, ha="left", va="top")
    _save(fig, P, "fig5_vs_domains", dpi); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Fig 2–5",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--which", nargs="+", type=int, default=[2, 3, 4, 5])
    ap.add_argument("--dpi", type=int, default=400)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")

    P = Paths(load_config(args.config))

    global ORDER
    ORDER = admitted_slides(P)
    if not ORDER:
        sys.exit("台账里没有 status=ingested 的切片，没东西可画。")
    print(f"[Fig] 切片 {len(ORDER)} 张：{', '.join(ORDER)}")

    for n in args.which:
        print(f"[Fig {n}]")
        {2: fig2, 3: fig3, 4: fig4, 5: fig5}[n](P, args.dpi)


if __name__ == "__main__":
    main()
