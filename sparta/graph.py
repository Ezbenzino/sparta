"""
M3 空间图构建与源汇定义
========================

上游依赖：M2 签名打分（{slide_id}.scored.h5ad，分数在 adata.obs 中）
下游产出：邻接矩阵 A + 源汇索引，供 barrier.py 使用

本模块同样遵循「核心函数只吃数组」的原则：build_graph 与 define_source_sink
的核心版本接受坐标与分数数组；*_from_adata 是薄适配器。

两种建图方式
------------
  build_graph_radius : 半径建图（cKDTree）。通用，第一代 ST 的非规则坐标也能用。
  build_graph_knn    : k 近邻建图。当 spot 密度不均时更稳。
squidpy 的 spatial_neighbors 也可以用，但它不是必需依赖——这里给出的实现
只用 scipy，装不上 squidpy 也不影响主流程。

关于 radius_um 的取值
---------------------
  Visium      : spot 中心间距 100 μm，取 150 μm 可覆盖六个一阶邻居
  第一代 ST   : 间距约 200 μm，取 250–300 μm
这个参数必须做敏感性分析——它直接决定图的连通性，进而决定所有屏障值。
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree

__all__ = [
    "build_graph_radius",
    "build_graph_knn",
    "graph_report",
    "define_source_sink",
    "rank_normalize_array",
]


# --------------------------------------------------------------------------
def build_graph_radius(coords_um: np.ndarray, radius_um: float = 150.0):
    """半径建图：距离小于 radius_um 的 spot 互为邻居。

    参数
    ----
    coords_um : (n, 2) 微米坐标。务必先把像素坐标换算成微米，
                否则 radius_um 这个有物理意义的参数就失去了意义。
    radius_um : 邻接半径（微米）

    返回
    ----
    A : (n, n) 稀疏对称二值邻接矩阵（scipy CSR）
    D : (n, n) 稀疏对称距离矩阵（微米），供需要真实距离的算子使用
    """
    coords_um = np.asarray(coords_um, dtype=float)
    n = len(coords_um)
    tree = cKDTree(coords_um)
    pairs = tree.query_pairs(r=radius_um, output_type="ndarray")
    if len(pairs) == 0:
        raise ValueError(
            f"radius_um={radius_um} 下没有找到任何邻居对。"
            f"请检查坐标单位是否为微米（当前坐标范围 "
            f"{coords_um.min():.1f}–{coords_um.max():.1f}）"
        )
    d = np.linalg.norm(coords_um[pairs[:, 0]] - coords_um[pairs[:, 1]], axis=1)

    rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
    A = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n)).tocsr()
    A.data[:] = 1.0
    D = sp.coo_matrix((np.concatenate([d, d]), (rows, cols)), shape=(n, n)).tocsr()
    return A, D


def build_graph_knn(coords_um: np.ndarray, k: int = 6):
    """k 近邻建图。spot 密度不均时比半径建图更稳（每个点保证有 k 个邻居）。"""
    coords_um = np.asarray(coords_um, dtype=float)
    n = len(coords_um)
    tree = cKDTree(coords_um)
    dist, idx = tree.query(coords_um, k=k + 1)  # 第 0 列是自己
    rows = np.repeat(np.arange(n), k)
    cols = idx[:, 1:].ravel()
    d = dist[:, 1:].ravel()

    A = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    A = ((A + A.T) > 0).astype(float).tocsr()   # 对称化
    D = sp.coo_matrix((d, (rows, cols)), shape=(n, n))
    D = D.maximum(D.T).tocsr()
    return A, D


def graph_report(A: sp.spmatrix) -> dict:
    """图的健康检查。建图之后必须看一眼这个报告。

    最重要的一项是连通分量数：如果图碎成很多块，说明 radius_um 太小，
    有效阻抗与图距离都会失去意义（跨分量是无穷大）。
    """
    n = A.shape[0]
    deg = np.asarray(A.sum(axis=1)).ravel()
    n_comp, labels = connected_components(A, directed=False)
    sizes = np.bincount(labels)
    return dict(
        n_nodes=n,
        n_edges=int(A.nnz // 2),
        mean_degree=float(deg.mean()),
        min_degree=int(deg.min()),
        n_isolated=int((deg == 0).sum()),
        n_components=int(n_comp),
        largest_component_frac=float(sizes.max() / n),
    )


# --------------------------------------------------------------------------
def rank_normalize_array(x: np.ndarray) -> np.ndarray:
    """把任意分数秩标准化到 [0, 1]。

    为什么必须做这一步：不同签名的原始分数量纲差异很大（score_genes 的输出
    可能是 -0.3 到 2.1，也可能是 0.01 到 0.05）。直接把它们代进边权公式，
    等于给某些项施加了隐性的高权重。秩转换消除了这个问题，也顺带抹平了
    离群值的影响。
    """
    from scipy.stats import rankdata

    x = np.asarray(x, dtype=float)
    if len(x) <= 1:
        return np.full_like(x, 0.5)
    return (rankdata(x) - 1) / (len(x) - 1)


# --------------------------------------------------------------------------
def define_source_sink(
    A: sp.spmatrix,
    endothelial: np.ndarray,
    t_nk: np.ndarray,
    malignant: np.ndarray,
    *,
    q_vessel: float = 0.80,
    q_immune_nbr: float = 0.60,
    q_malig: float = 0.70,
    q_core: float = 0.50,
    fallback_border_coords: np.ndarray | None = None,
) -> dict:
    """定义血管节点、免疫入口（源）与瘤巢核心（汇）。

    这是全流程中最"人为"的一步，也是审稿人一定会追问的一步。
    四个分位数阈值必须做网格扫描，证明结论不依赖于某一组特定取值。

    参数
    ----
    A            : (n, n) 邻接矩阵
    endothelial  : (n,) 内皮签名分数，已秩标准化
    t_nk         : (n,) T/NK 签名分数，已秩标准化
    malignant    : (n,) 恶性签名分数，已秩标准化
    q_vessel     : 内皮分数分位阈值，超过者视为血管节点
    q_immune_nbr : 血管节点中，其邻域 T 细胞分数超过此分位者视为免疫入口
    q_malig      : 恶性分数分位阈值
    q_core       : 在恶性节点中，距非肿瘤区图距离超过此分位者视为瘤巢核心
    fallback_border_coords :
        可选。若切片内皮信号过弱（准入标准 C6 不合格的降级路径），
        传入坐标数组，则改用「组织最外圈 spot」作为血管代理源，
        并在返回值中把 used_fallback 标记为 True。
        这一替代必须在论文方法部分明确写出。

    返回
    ----
    dict，含 source / sink / vessel 三组索引，以及诊断信息
    """
    n = A.shape[0]
    used_fallback = False

    vessel = np.where(endothelial >= np.quantile(endothelial, q_vessel))[0]

    # 降级路径：内皮信号不可用时改用组织边缘
    if len(vessel) < 5 and fallback_border_coords is not None:
        used_fallback = True
        c = np.asarray(fallback_border_coords, float)
        center = c.mean(axis=0)
        d = np.linalg.norm(c - center, axis=1)
        vessel = np.where(d >= np.quantile(d, 0.90))[0]

    # 免疫入口：血管节点中，邻域 T 细胞浸润也高的
    deg = np.maximum(np.asarray(A.sum(axis=1)).ravel(), 1.0)
    t_nbr = (A @ t_nk) / deg
    if len(vessel):
        thr_imm = np.quantile(t_nbr[vessel], q_immune_nbr)
        source = vessel[t_nbr[vessel] >= thr_imm]
    else:
        source = np.array([], int)

    # 瘤巢核心：恶性分数高，且距非肿瘤区的图距离大
    malig = np.where(malignant >= np.quantile(malignant, q_malig))[0]
    non_tumor = np.setdiff1d(np.arange(n), malig)
    if len(non_tumor) and len(malig):
        d2edge = dijkstra(A, directed=False, indices=non_tumor, min_only=True)
        d2edge = np.where(np.isfinite(d2edge), d2edge, 0.0)
        thr_core = np.quantile(d2edge[malig], q_core)
        sink = malig[d2edge[malig] >= thr_core]
    else:
        sink = malig

    # 源汇不能重叠，否则最大流无解
    overlap = np.intersect1d(source, sink)
    if len(overlap):
        sink = np.setdiff1d(sink, overlap)

    return dict(
        source=source,
        sink=sink,
        vessel=vessel,
        used_fallback=used_fallback,
        n_source=len(source),
        n_sink=len(sink),
        n_vessel=len(vessel),
    )


# --------------------------------------------------------------------------
# AnnData 适配器
# --------------------------------------------------------------------------
def build_graph_from_adata(adata, radius_um: float = 150.0, key: str = "spatial_um"):
    """从 adata.obsm[key] 取微米坐标建图，结果写回 adata.obsp。"""
    if key not in adata.obsm:
        raise KeyError(
            f"adata.obsm 中没有 '{key}'。请先在 M1 质控阶段把像素坐标换算为微米。"
        )
    A, D = build_graph_radius(adata.obsm[key], radius_um=radius_um)
    adata.obsp["sparta_adj"] = A
    adata.obsp["sparta_dist"] = D
    adata.uns["sparta_graph_report"] = graph_report(A)
    return A, D
