"""
可视化 —— 屏障线叠加图是本课题最有说服力的图件
=================================================

上游依赖：barrier.py（cut_edges、b_mab、b_meta）、graph.py（coords）
下游产出：results/figures/*

为什么单独做一个模块
--------------------
把最小割边集画在 H&E 上得到的那条"封锁线"，是本课题唯一能让审稿人一眼看懂
的图。它值得认真实现，而不是在 notebook 里临时拼几行 plt。

本模块只依赖 matplotlib + numpy，不依赖 scanpy/squidpy。
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "setup_cjk_font",
    "plot_barrier_line", "plot_barrier_field", "plot_permutation_null",
    "plot_ring_breaking", "plot_size_scan", "plot_decoupling",
]


def setup_cjk_font(verbose: bool = False) -> str | None:
    """让 matplotlib 能正确显示中文。

    注意：直接给 rcParams["font.sans-serif"] 赋值**不会**在字体不存在时报错，
    所以必须先查一遍系统里实际装了哪些字体。这是一个很容易写错的地方。

    返回实际选中的字体名；一个都没找到则返回 None，并给出安装提示。
    """
    import matplotlib
    from matplotlib import font_manager

    available = {f.name for f in font_manager.fontManager.ttflist}
    for fam in ("Noto Sans CJK SC", "Noto Sans CJK JP", "WenQuanYi Zen Hei",
                "WenQuanYi Micro Hei", "SimHei", "Microsoft YaHei",
                "PingFang SC", "Heiti SC", "Arial Unicode MS", "Source Han Sans SC"):
        if fam in available:
            matplotlib.rcParams["font.sans-serif"] = [fam, "DejaVu Sans"]
            matplotlib.rcParams["axes.unicode_minus"] = False
            if verbose:
                print(f"[viz] 中文字体：{fam}")
            return fam
    print("[viz] 未找到中文字体，图中的中文会显示为方块。")
    print("      Ubuntu/WSL: sudo apt install fonts-noto-cjk")
    print("      或改用英文标题（各绘图函数都接受 title 参数）。")
    return None

_CMAP_BARRIER = "magma_r"


def _new_ax(ax, figsize=(6, 6)):
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    return ax


# --------------------------------------------------------------------------
def plot_barrier_line(
    coords, cut_edges, *, he_image=None, extent=None, node_color=None,
    source=None, sink=None, ax=None, line_kw=None, title=None,
):
    """把最小割边集画成"封锁线"，可叠加在 H&E 图像上。

    参数
    ----
    coords    : (n, 2) 坐标。若要叠加 H&E，请传**图像像素坐标**而非微米坐标，
                并把 he_image 与 extent 一并传入。
    cut_edges : barrier.compute_b_cell 返回的 cut_edges
    he_image  : (H, W, 3) 的 H&E 缩略图数组；None 则画在白底上
    extent    : imshow 的 extent，用于把图像与坐标对齐
    node_color: (n,) 用于给 spot 上色的标量（如 CAF 分数）
    source/sink : 索引数组，会被标出来（这对读图很重要——
                  没有源汇的位置，一条封锁线是看不懂的）

    返回
    ----
    matplotlib Axes
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    ax = _new_ax(ax)
    coords = np.asarray(coords, float)

    if he_image is not None:
        ax.imshow(he_image, extent=extent, origin="upper")

    if node_color is not None:
        ax.scatter(coords[:, 0], coords[:, 1], c=node_color, s=6,
                   cmap="Greys", alpha=0.55, linewidths=0)
    else:
        ax.scatter(coords[:, 0], coords[:, 1], s=4, c="0.75", linewidths=0)

    if source is not None and len(source):
        ax.scatter(coords[source, 0], coords[source, 1], s=14, c="#2E86C1",
                   label="源：免疫入口", linewidths=0)
    if sink is not None and len(sink):
        ax.scatter(coords[sink, 0], coords[sink, 1], s=14, c="#C0392B",
                   label="汇：瘤巢核心", linewidths=0)

    if len(cut_edges):
        segs = [[coords[u], coords[v]] for u, v in cut_edges]
        kw = dict(colors="#F1C40F", linewidths=2.0, alpha=0.95, zorder=5)
        kw.update(line_kw or {})
        ax.add_collection(LineCollection(segs, **kw))
        ax.plot([], [], color=kw["colors"], lw=kw["linewidths"], label="最小割：封锁线")

    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    if title:
        ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.85)
    return ax


# --------------------------------------------------------------------------
def plot_barrier_field(coords, values, *, ax=None, title=None, cbar_label=None,
                       cmap=_CMAP_BARRIER, s=10, vmin=None, vmax=None):
    """画逐 spot 的屏障场热图（B_mAb / B_meta）。"""
    import matplotlib.pyplot as plt

    ax = _new_ax(ax)
    coords = np.asarray(coords, float)
    finite = np.isfinite(values)
    sc = ax.scatter(coords[finite, 0], coords[finite, 1], c=np.asarray(values)[finite],
                    cmap=cmap, s=s, linewidths=0, vmin=vmin, vmax=vmax)
    if (~finite).any():
        ax.scatter(coords[~finite, 0], coords[~finite, 1], s=s, c="0.85",
                   linewidths=0, label="不可达")
        ax.legend(loc="upper right", fontsize=8)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    if title:
        ax.set_title(title)
    cb = plt.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
    if cbar_label:
        cb.set_label(cbar_label)
    return ax


# --------------------------------------------------------------------------
def plot_permutation_null(res, *, ax=None, title="S1 空间重排对照"):
    """画 S1 的零分布与真实值。"""
    ax = _new_ax(ax, figsize=(5.5, 3.6))
    null, real = np.asarray(res["null"]), res["b_real"]
    ax.hist(null, bins=40, color="0.75", edgecolor="white", label="空间重排零分布")
    ax.axvline(real, color="#C0392B", lw=2.2, label=f"真实值 = {real:.3f}")
    ax.set_xlabel("B_cell"); ax.set_ylabel("频数")
    ax.set_title(f"{title}（mode={res['mode']}, z={res['z']:.1f}, p={res['p_emp']:.4f}）")
    ax.legend(fontsize=8)
    return ax


def plot_ring_breaking(res, *, ax=None):
    """画 S2 的三条曲线：连续缺口 vs 全局随机 vs 屏障内分散。"""
    ax = _new_ax(ax, figsize=(5.5, 3.8))
    ks = sorted(res["per_k"])
    b0 = res["b0"]
    tgt = [res["per_k"][k]["targeted"] / b0 for k in ks]
    rg = [res["per_k"][k]["rand_global_mean"] / b0 for k in ks]
    rc = [res["per_k"][k]["rand_in_cut_mean"] / b0 for k in ks]
    ax.plot(ks, tgt, "o-", color="#C0392B", label="连续缺口（定向）")
    ax.plot(ks, rc, "s--", color="#E67E22", label="屏障内分散（关键对照）")
    ax.plot(ks, rg, "^:", color="0.55", label="全局随机")
    ax.axhline(1.0, color="0.8", lw=1, ls="-")
    ax.set_xlabel("移除的 spot 数 k"); ax.set_ylabel("剩余屏障 / 基线屏障")
    ax.set_title("S2 环带断裂")
    ax.legend(fontsize=8)
    return ax


def plot_size_scan(res, *, ax=None):
    """画 S3 的尺寸-屏障曲线，并标出 IgG 与小分子的位置。"""
    ax = _new_ax(ax, figsize=(5.5, 3.8))
    r = res["radii_nm"]
    ax.plot(r, res["mean_all"], "o-", color="0.45", label="全片均值")
    if np.isfinite(res["mean_core"]).any():
        ax.plot(r, res["mean_core"], "s-", color="#C0392B", label="瘤巢核心均值")
    ax.axvline(5.5, color="#2E86C1", ls="--", lw=1.2)
    ax.text(5.6, ax.get_ylim()[1] * 0.92, "IgG / ADC\n(5.5 nm)", fontsize=8, color="#2E86C1")
    ax.axvline(0.5, color="#27AE60", ls="--", lw=1.2)
    ax.text(0.6, ax.get_ylim()[1] * 0.78, "小分子", fontsize=8, color="#27AE60")
    ax.set_xlabel("分子流体力学半径 (nm)"); ax.set_ylabel("B_mAb")
    ax.set_title("S3 分子尺寸扫描")
    ax.legend(fontsize=8)
    return ax


def plot_decoupling(b_cell_like, b_mab, *, labels=None, ax=None,
                    xlabel="B_cell（T 细胞迁移屏障）", ylabel="B_mAb（抗体传质屏障）"):
    """解耦散点图 —— 本课题的中心图。

    横轴 T 细胞屏障、纵轴抗体屏障。右下象限（T 细胞进得去、抗体进不去）
    就是本课题主张存在、而现有框架完全看不到的那类区域。
    """
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr

    ax = _new_ax(ax, figsize=(5.2, 5.0))
    x, y = np.asarray(b_cell_like, float), np.asarray(b_mab, float)
    ok = np.isfinite(x) & np.isfinite(y)
    ax.scatter(x[ok], y[ok], s=8, c="#34495E", alpha=0.5, linewidths=0)
    rho, p = spearmanr(x[ok], y[ok]) if ok.sum() > 3 else (np.nan, np.nan)

    mx, my = np.nanmedian(x[ok]), np.nanmedian(y[ok])
    ax.axvline(mx, color="0.8", lw=1); ax.axhline(my, color="0.8", lw=1)
    ax.text(0.98, 0.03, f"Spearman ρ = {rho:.3f}\np = {p:.2e}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9)
    # 标注关键象限
    ax.text(0.02, 0.97, "抗体挡住\nT 细胞进得去", transform=ax.transAxes,
            ha="left", va="top", fontsize=8, color="#C0392B")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_title("双屏障解耦")
    return ax
