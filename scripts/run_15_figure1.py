#!/usr/bin/env python
"""
run_15_figure1.py —— Fig 1 框架示意图（一张图、两套边权语义、两个算子）
==========================================================================

输入：无（用 sparta.synthetic 现场生成的合成切片，图上明确标注 synthetic）
输出：results/figures/fig1_framework.png / .pdf
上游模块：sparta/synthetic.py、sparta/barrier.py
下游模块：无（成文用）

为什么用合成切片而不是真实切片
------------------------------
Fig 1 的任务是让读者在十秒内看懂"同一张图 + 两套边权语义 = 两个屏障"。
真实切片的噪声会淹没这个结构。合成切片把结构做干净，读者一眼就懂，
真实数据的结果放 Fig 2 起。按 CLAUDE.md 的要求，图上明确标注 synthetic。

配色（已过 CVD 校验）
---------------------
两个屏障用两个分类色承载身份：cell = #2a78d6（蓝），antibody = #eb6834（橙）。
这一对在 protan/deutan/tritan 下的最差 ΔE 为 24.7，正常视觉 33.6，
对比度均 >= 3:1（validate_palette.js, light mode 全通过）。
每个屏障场的**数值**用该色相的单色 light->dark 渐变承载，
所以"哪个屏障"看色相、"多强"看明暗，两条信息不混。

⚠ (d) 与 (e) 必须共用同一个色标，否则小分子/IgG 的对比没有意义。

用法
----
    python scripts/run_15_figure1.py
    python scripts/run_15_figure1.py --n 27 --dpi 600
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import Paths, load_config, set_seed  # noqa: E402

# ---- 分类色（身份）与由它派生的单色顺序渐变（数值）----
CELL = "#2a78d6"     # categorical slot 1
MAB = "#eb6834"      # categorical slot 2
INK = "#0b0b0b"      # text-primary
INK2 = "#52514e"     # text-secondary
MUTED = "#8f8e88"
SURFACE = "#fcfcfb"


def _ramp(name, hex_hi):
    """由一个分类色派生 light->dark 的单色顺序渐变（顺序色永远单色相）。"""
    from matplotlib.colors import LinearSegmentedColormap, to_rgb
    r, g, b = to_rgb(hex_hi)
    light = (1 - 0.92 * (1 - r), 1 - 0.92 * (1 - g), 1 - 0.92 * (1 - b))
    dark = (r * 0.42, g * 0.42, b * 0.42)
    return LinearSegmentedColormap.from_list(name, [light, (r, g, b), dark])


def _bare(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_aspect("equal")


def _panel_title(ax, letter, text, sub=None):
    """标题与副标题手工放置——set_title 与副标题会在同一条基线上打架。"""
    ax.text(0, 1.105, f"{letter}  {text}", transform=ax.transAxes,
            fontsize=10.5, color=INK, va="bottom", ha="left")
    if sub:
        ax.text(0, 1.035, sub, transform=ax.transAxes, fontsize=8.2,
                color=INK2, va="bottom", ha="left")


def main():
    ap = argparse.ArgumentParser(description="Fig 1 框架示意图",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--n", type=int, default=25, help="合成网格边长（spot 数）")
    ap.add_argument("--dpi", type=int, default=400)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    from sparta.barrier import (compute_b_cell, compute_b_cell_field,
                                compute_b_mab, edge_pairs)
    from sparta.synthetic import make_ring_grid

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    cc = cfg["barrier"]["b_cell"]
    cm = {k: v for k, v in cfg["barrier"]["b_mab"].items() if k != "r_nm"}

    sl = make_ring_grid(n=args.n, width_spots=2.0, seed=cfg["seed"])
    xy, A = sl.coords, sl.A
    pairs = edge_pairs(A, upper_only=True)
    seg = np.stack([xy[pairs[:, 0]], xy[pairs[:, 1]]], axis=1)

    bc = compute_b_cell(A, sl.ecm, sl.caf, sl.source, sl.sink, **cc)
    bcf = compute_b_cell_field(A, sl.ecm, sl.caf, sl.source, **cc)["b_cell_field"]
    b_igg = compute_b_mab(A, sl.ecm, sl.crosslink, sl.ag_target, sl.vessel,
                          r_nm=5.5, **cm)["b_mab"]
    b_small = compute_b_mab(A, sl.ecm, sl.crosslink, sl.ag_target, sl.vessel,
                            r_nm=0.5, **cm)["b_mab"]

    radii = [0.5, 1.0, 2.0, 3.0, 4.0, 5.5, 7.0, 10.0]
    core = [float(np.nanmean(compute_b_mab(A, sl.ecm, sl.crosslink, sl.ag_target,
                                           sl.vessel, r_nm=r, **cm)["b_mab"][sl.sink]))
            for r in radii]

    cmap_cell, cmap_mab = _ramp("cell", CELL), _ramp("mab", MAB)
    vmax_mab = float(max(np.nanmax(b_igg), np.nanmax(b_small)))   # (d)(e) 共用色标

    fig = plt.figure(figsize=(12.2, 8.4), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.30,
                          left=0.045, right=0.965, top=0.790, bottom=0.075)
    axes = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(3)]
    for ax in axes:
        ax.set_facecolor(SURFACE)

    # ---------------------------------------------------------------- (a)
    ax = axes[0]; _bare(ax)
    ax.add_collection(LineCollection(seg, colors="#dedcd6", linewidths=0.55, zorder=1))
    ax.scatter(*xy.T, s=13, c=sl.caf, cmap="Greys", vmin=-0.15, vmax=1.25,
               linewidths=0, zorder=2)
    ax.scatter(*xy[sl.source].T, s=26, facecolors="none", edgecolors=INK,
               linewidths=1.1, zorder=3)
    ax.scatter(*xy[sl.sink].T, s=17, c=INK, marker="s", linewidths=0, zorder=3)
    # 图例放在面板外——放进面板会被边缘的 source 节点压住
    ax.text(0.5, -0.075, "\u25cb immune entry (source)     \u25a0 tumour core (sink)",
            transform=ax.transAxes, fontsize=8.2, color=INK2, ha="center")
    ax.text(0.5, -0.145, "node shade = fibroblast / matrix score",
            transform=ax.transAxes, fontsize=8, color=MUTED, ha="center")
    _panel_title(ax, "a", "One spatial graph", "source, sink and vessel sets")

    # ---------------------------------------------------------------- (b)
    ax = axes[1]; _bare(ax)
    ax.add_collection(LineCollection(seg, colors="#e6e4de", linewidths=0.5, zorder=1))
    ax.scatter(*xy.T, s=8, c="#d8d6d0", linewidths=0, zorder=2)
    cut = np.array(bc["cut_edges"], int)
    if len(cut):
        ax.add_collection(LineCollection(np.stack([xy[cut[:, 0]], xy[cut[:, 1]]], axis=1),
                                         colors=CELL, linewidths=2.0, zorder=4))
    ax.scatter(*xy[sl.source].T, s=22, facecolors="none", edgecolors=INK,
               linewidths=1.0, zorder=3)
    ax.scatter(*xy[sl.sink].T, s=14, c=INK, marker="s", linewidths=0, zorder=3)
    ax.text(0.5, -0.075,
            f"max-flow = {bc['max_flow']:.4f}    →    "
            f"B$_{{cell}}$ = 1/max-flow = {bc['b_cell']:.1f}",
            transform=ax.transAxes, fontsize=8.4, color=INK, ha="center")
    _panel_title(ax, "b", "Operator 1 — source–sink min cut",
                 "the blockade is an object, not a score")

    # ---------------------------------------------------------------- (c)
    ax = axes[2]; _bare(ax)
    bcf_s = (bcf - np.nanmin(bcf)) / max(np.nanmax(bcf) - np.nanmin(bcf), 1e-12)
    sc = ax.scatter(*xy.T, s=15, c=bcf_s, cmap=cmap_cell, vmin=0, vmax=1,
                    linewidths=0, zorder=2)
    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02, ticks=[0, 0.5, 1])
    cb.set_label("B$_{cell}$ field (scaled)", fontsize=8.2, color=INK2)
    cb.ax.tick_params(labelsize=7, colors=INK2, length=2)
    cb.outline.set_visible(False)
    ax.text(0.5, -0.075, "min–max scaled within section",
            transform=ax.transAxes, fontsize=8, color=MUTED, ha="center")
    _panel_title(ax, "c", "Cellular barrier field",
                 "accumulated migration cost from source")

    # ---------------------------------------------------------------- (d)(e)
    for ax, field, lab, letter, sub in (
            (axes[3], b_small, "small molecule   r = 0.5 nm", "d", "same tissue, same graph"),
            (axes[4], b_igg, "IgG   r = 5.5 nm", "e", "only the molecule changed")):
        _bare(ax)
        sc = ax.scatter(*xy.T, s=15, c=field, cmap=cmap_mab, vmin=0, vmax=vmax_mab,
                        linewidths=0, zorder=2)
        cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
        cb.set_label("B$_{mAb}$", fontsize=8.2, color=INK2)
        cb.ax.tick_params(labelsize=7, colors=INK2, length=2)
        cb.outline.set_visible(False)
        ax.text(0.5, -0.075, lab, transform=ax.transAxes, fontsize=8.4,
                color=INK, ha="center")
        _panel_title(ax, letter, "Operator 2 — screened diffusion", sub)

    # ---------------------------------------------------------------- (f)
    ax = axes[5]
    ax.set_facecolor(SURFACE)
    ax.set_aspect("auto")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#dedcd6"); ax.spines[s].set_linewidth(0.8)
    ax.grid(axis="y", color="#eeece6", lw=0.8)
    ax.set_axisbelow(True)
    ax.plot(radii, core, color=MAB, lw=2.0, marker="o", ms=5.5,
            markeredgecolor=SURFACE, markeredgewidth=1.2, zorder=3)
    ax.margins(x=0.08, y=0.16)
    ax.annotate(f"{core[0]:.1f}", xy=(0.5, core[0]), xytext=(6, -3),
                textcoords="offset points", fontsize=8.4, color=INK, ha="left")
    i_igg = radii.index(5.5)
    ax.annotate(f"{core[i_igg]:.1f}", xy=(5.5, core[i_igg]), xytext=(0, -16),
                textcoords="offset points", fontsize=8.4, color=INK, ha="center")
    ax.axvline(5.5, color=MUTED, lw=0.9, ls=(0, (3, 3)), zorder=1)
    ax.text(5.5, ax.get_ylim()[0], " IgG", fontsize=8.2, color=INK2,
            va="bottom", ha="left")
    ax.set_xlabel("hydrodynamic radius (nm)", fontsize=8.6, color=INK2)
    ax.set_ylabel("mean B$_{mAb}$ in tumour core", fontsize=8.6, color=INK2)
    ax.tick_params(labelsize=7.6, colors=INK2, length=3)
    _panel_title(ax, "f", "Size selectivity", "one parameter sweep, same tissue")

    # ---------------------------------------------------------------- 顶部说明
    fig.text(0.045, 0.975,
             "SPARTA: one graph, two edge-weight semantics, two transport operators",
             fontsize=13.5, color=INK, ha="left", va="top")
    fig.text(0.045, 0.928,
             r"cell:  capacity$(u,v)=\sigma\,(a - b\,\overline{\mathrm{ECM}} - c\,\overline{\mathrm{CAF}})$"
             "        min cut",
             fontsize=9.4, color=CELL, ha="left", va="top")
    fig.text(0.045, 0.893,
             r"antibody:  conductance$(u,v)=g_0 e^{-\lambda \overline{\mathrm{ECM}}}\,\phi(r/\xi)$"
             r"        $(L+\mathrm{diag}(\kappa))\,\varphi=0,\ \ B_{mAb}=-\log\varphi$",
             fontsize=9.2, color=MAB, ha="left", va="top")
    fig.text(0.965, 0.975,
             "synthetic section — for illustration only\nno real data, no research claim",
             fontsize=8, color=MUTED, ha="right", va="top")

    P.figure("x")  # 确保 results/figures 存在
    for ext in ("png", "pdf"):
        p = P.figure(f"fig1_framework.{ext}")
        fig.savefig(p, dpi=args.dpi, facecolor=SURFACE)
        print(f"已保存 {p}")
    plt.close(fig)


if __name__ == "__main__":
    main()
