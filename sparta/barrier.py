"""
M4 三分量屏障算子 —— 本项目的核心模块
==========================================

上游依赖：graph.py（提供邻接矩阵 A 与源汇索引）、signatures.py（提供各项签名分数）
下游产出：{slide_id}.barrier.npz（B_cell / B_mAb / B_meta 逐 spot 值）
          {slide_id}.mincut.json（最小割边集，供 viz.py 绘制"封锁线"）

设计说明（重要）
----------------
本模块中的三个核心函数**只接受 numpy / scipy 数组，不接受 AnnData 对象**。
这样做有三个理由：
  1. 核心数学与数据格式解耦，换成 Xenium、CosMx 或任何其他平台都不用改这里；
  2. 单元测试可以在合成图上运行，不需要安装 scanpy（tests/test_barrier.py 即如此）；
  3. 零基础读者更容易看懂——输入就是几个数组，没有隐藏状态。

如果你手上是 AnnData，用本模块末尾的 *_from_adata 适配器，它只负责把
adata.obs 里的列取出来变成数组，然后调用核心函数。

三个屏障分量的物理含义
----------------------
B_cell  : T 细胞（直径约 10 μm 的有形细胞）从血管迁移到瘤巢的阻力。
          用【源汇最小割】——最窄封锁截面的通过容量的倒数。
B_mAb   : 抗 PD-1/PD-L1 抗体（流体力学半径约 5.5 nm 的大分子）扩散到某点的阻力。
          用【尺寸依赖有效阻抗】+ 结合位点消耗修正。
B_meta  : 小分子药物到得了但失效的程度（缺氧静止 + 外排）。
          用【自血管出发的图距离】× 局部代谢状态。
"""

from __future__ import annotations

import json
from typing import Iterable, Sequence

import networkx as nx
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components, dijkstra, laplacian
from scipy.sparse.linalg import cg, spsolve

__all__ = [
    "edge_pairs",
    "exact_min_cut",
    "compute_b_cell",
    "compute_b_cell_field",
    "compute_b_mab",
    "compute_b_meta",
    "compute_all_barriers",
    "save_barriers",
]

_EPS = 1e-9


# --------------------------------------------------------------------------
# 工具
# --------------------------------------------------------------------------
def edge_pairs(A: sp.spmatrix, upper_only: bool = True) -> np.ndarray:
    """从稀疏邻接矩阵中提取边列表。

    参数
    ----
    A          : (n, n) 稀疏对称邻接矩阵。只用其非零结构，不用其数值。
    upper_only : True 时只返回 u < v 的边（每条无向边出现一次）；
                 False 时返回全部有向对（每条无向边出现两次）。

    返回
    ----
    (m, 2) 的 int 数组，每行是一条边的两个端点索引。
    """
    Ac = A.tocoo()
    u, v = Ac.row, Ac.col
    keep = u < v if upper_only else u != v
    return np.column_stack([u[keep], v[keep]]).astype(int)


def _pair_mean(x: np.ndarray, pairs: np.ndarray) -> np.ndarray:
    """取每条边两端节点分数的均值。返回 (m,) 数组。"""
    return 0.5 * (x[pairs[:, 0]] + x[pairs[:, 1]])


def _sigmoid(z: np.ndarray | float) -> np.ndarray | float:
    return 1.0 / (1.0 + np.exp(-np.clip(z, -60, 60)))


def _as_idx(x: Iterable[int] | np.ndarray) -> np.ndarray:
    return np.asarray(list(x), dtype=int) if not isinstance(x, np.ndarray) else x.astype(int)


# --------------------------------------------------------------------------
# B_cell —— T 细胞迁移屏障（源汇最小割）
# --------------------------------------------------------------------------
# 2026-10-05：最小割改用整数容量求解。
# networkx 的 preflow-push 在浮点容量上给出的最大流数值是对的（与精确值相对误差 < 1e-12），
# 但它用残量图可达性导出的"源侧集合"在浮点舍入下并不可靠：实测 22 张切片里有 12 张
# 返回的割边总容量与最大流相差 0.03%–5.4%（有的割边缺失，有的多出），也就是说
# 返回的 cut_edges 并不是真正的最小割。networkx 文档也提示浮点容量可能出问题，
# 建议乘一个常数换成整数。这里把容量取到 2^-40 的整数倍后用 Python 整数精确求解：
# 最大流与原实现的差 < 1e-12（相对），而割边集的容量严格等于最大流。
_CUT_SCALE = 2 ** 40


def exact_min_cut(nodes, pairs: np.ndarray, cap: np.ndarray, source, sink):
    """在无向图（每条边拆成两个方向）上求多源多汇最大流与一个精确的最小割。

    参数
    ----
    nodes  : 参与计算的节点编号（可迭代）
    pairs  : (m, 2) 无向边端点
    cap    : (m,) 边容量（浮点，>= 0）
    source, sink : 源 / 汇节点编号

    返回
    ----
    (max_flow, reach, unreach)：max_flow 为浮点；reach / unreach 是最小割两侧的节点集合
    （含超级源 "s" / 超级汇 "t"）。割边 = reach -> unreach 的有限容量边，其容量和等于 max_flow
    （在 2^-40 的取整精度内）。
    """
    ci = np.rint(np.asarray(cap, float) * _CUT_SCALE).astype(np.int64).tolist()
    G = nx.DiGraph()
    G.add_nodes_from(int(u) for u in nodes)
    for (u, v), c in zip(np.asarray(pairs, int).tolist(), ci):
        G.add_edge(u, v, capacity=c)
        G.add_edge(v, u, capacity=c)
    big = int(sum(ci)) * 10 + 1
    for i in source:
        G.add_edge("s", int(i), capacity=big)
    for i in sink:
        G.add_edge(int(i), "t", capacity=big)
    f, (reach, unreach) = nx.minimum_cut(G, "s", "t", capacity="capacity")
    return float(f) / _CUT_SCALE, reach, unreach, G


def compute_b_cell(
    A: sp.spmatrix,
    ecm: np.ndarray,
    caf: np.ndarray,
    source: Sequence[int],
    sink: Sequence[int],
    *,
    a: float = 3.0,
    b_ecm: float = 8.0,
    c_caf: float = 4.0,
    return_capacity: bool = False,
) -> dict:
    """T 细胞迁移屏障 = 1 / 最大流。

    直觉
    ----
    把组织想成一张水管网。每条边是一段管子，管子的"通径"（capacity）由该处
    ECM 有多致密、CAF 有多密集决定——越致密越细。然后问：从血管（源）到瘤巢
    中心（汇），最多能通过多少水？这个最大流量的倒数就是屏障强度。

    最大流—最小割定理保证：最大流 == 最窄封锁截面上所有管子通径之和。
    所以我们同时拿到了一个数值（屏障强度）和一条几何（封锁线的位置）。

    参数
    ----
    A       : (n, n) 稀疏对称邻接矩阵
    ecm     : (n,) core matrisome 签名分数，已秩标准化到 [0, 1]
    caf     : (n,) 成纤维细胞签名分数，已秩标准化到 [0, 1]
    source  : 源节点索引（免疫入口 spot）
    sink    : 汇节点索引（瘤巢核心 spot）
    a       : 基线通径（截距）。a 越大，所有管子都越粗
    b_ecm   : ECM 对通径的抑制系数
    c_caf   : CAF 对通径的抑制系数
    return_capacity : 是否在返回值中附带每条边的容量（调试与可视化用）

    返回
    ----
    dict，含：
      b_cell    : float，屏障强度 = 1 / (最大流 + eps)
      max_flow  : float，最大流量本身
      cut_edges : list[(u, v)]，最小割边集 —— 画在 H&E 上就是"封锁线"
      cut_nodes : np.ndarray，最小割边涉及的全部节点（S2 环带断裂实验要用）

    注意
    ----
    源集或汇集为空时返回 b_cell = nan，并给出 max_flow = nan。调用方需要处理
    这种情况（通常意味着该切片没通过准入标准 C6）。
    """
    source, sink = _as_idx(source), _as_idx(sink)
    n = A.shape[0]

    if len(source) == 0 or len(sink) == 0:
        return dict(b_cell=np.nan, max_flow=np.nan, cut_edges=[], cut_nodes=np.array([], int))

    # The manuscript defines the section-level barrier on the largest
    # connected component.  Keep source/sink nodes and capacities on that
    # component only; otherwise disconnected fragments with both endpoints
    # contribute additional, independent flows to the section total.
    n_comp, labels = connected_components(A, directed=False)
    component_sizes = np.bincount(labels, minlength=n_comp)
    largest_label = int(np.argmax(component_sizes))
    component_nodes = np.flatnonzero(labels == largest_label)
    source = source[labels[source] == largest_label]
    sink = sink[labels[sink] == largest_label]
    if len(source) == 0 or len(sink) == 0:
        return dict(b_cell=np.nan, max_flow=np.nan, cut_edges=[], cut_nodes=np.array([], int))

    pairs = edge_pairs(A, upper_only=True)
    pairs = pairs[labels[pairs[:, 0]] == largest_label]
    if len(pairs) == 0:
        return dict(b_cell=np.nan, max_flow=np.nan, cut_edges=[], cut_nodes=np.array([], int))

    # 边容量：ECM 与 CAF 越强 -> 通径越小 -> 细胞越难通过
    z = a - b_ecm * _pair_mean(ecm, pairs) - c_caf * _pair_mean(caf, pairs)
    cap = _sigmoid(z)  # 值域 (0, 1)

    # 源汇重叠时无解（会得到无穷大流），这是源汇定义写错的信号
    if set(source.tolist()) & set(sink.tolist()):
        raise ValueError("源集与汇集有重叠节点，请检查 graph.define_source_sink 的分位数阈值")

    # 有向化（每条无向边两个方向）+ 超级源汇，用整数容量精确求最大流 / 最小割（见 exact_min_cut）
    max_flow, reach, unreach, G = exact_min_cut(component_nodes, pairs, cap, source, sink)

    # 最小割边集：从可达侧指向不可达侧的边，排除超级源汇
    cut_edges = [
        (u, v)
        for u in reach
        if u != "s"
        for v in G.successors(u)
        if v in unreach and v != "t"
    ]
    cut_nodes = np.unique(np.array(cut_edges, dtype=int).ravel()) if cut_edges else np.array([], int)

    out = dict(
        b_cell=float(1.0 / (max_flow + _EPS)),
        max_flow=float(max_flow),
        cut_edges=cut_edges,
        cut_nodes=cut_nodes,
    )
    if return_capacity:
        out["edge_pairs"] = pairs
        out["edge_capacity"] = cap
    return out


def compute_b_cell_field(
    A: sp.spmatrix,
    ecm: np.ndarray,
    caf: np.ndarray,
    source: Sequence[int],
    *,
    a: float = 3.0,
    b_ecm: float = 8.0,
    c_caf: float = 4.0,
) -> dict:
    """T 细胞迁移屏障的**逐 spot 版本**（compute_b_cell 的场版本）。

    为什么需要这个函数
    ------------------
    compute_b_cell 返回的是切片级标量（整张切片的最小割），它回答
    "这张切片整体挡不挡"。但解耦分析需要回答的是"**这个位置**挡不挡"——
    横轴 B_cell、纵轴 B_mAb 的散点图，每个点是一个 spot，所以必须有逐点值。

    定义
    ----
    B_cell_field(v) = 从免疫入口集合出发，到达 v 的最小累积迁移代价，
                      其中每条边的代价 = 1 / capacity(u,v)。
    capacity 与 compute_b_cell 完全一致，所以两者物理语义统一：
    一个是"整张切片最窄截面的总容量"，一个是"到达某点最省力的那条路要多费劲"。

    与 compute_b_cell 的关系
    ------------------------
    二者互补而非替代：
      · 论文里报告"这张切片的屏障强度"用 compute_b_cell（并附最小割封锁线）；
      · 论文里做逐 spot 的解耦散点图、与 B_mAb 比较，用本函数。
    不要用局部 ECM/CAF 分数当逐点代理——那几乎是二值的，做分位数切分时
    会退化（这个坑在合成数据演示中出现过）。
    ⚠ 2026-08-26 审查发现 run_06_validate.py --dim decoupling 正是踩了这个坑
      （用 0.5*(ECM+CAF) 当代理且未传 control），产生的负相关是伪结果。已修。

    这个口径的已知局限（论文必须写）
    --------------------------------
    本函数是**最短路**（cost = 1/capacity 的 Dijkstra），不是拓扑量：
    只要存在一条低代价的绕行路线，累积代价就不增加——即"可以免费绕开封锁带"。
    这与本项目在 B_mAb 上修掉过的那个 bug 同型（见 CLAUDE.md「历史 bug」第 2 条）。
    因此：
      · 切片级的拓扑主张一律用 compute_b_cell（最小割），本函数不承担；
      · 本函数只用于需要逐点值的场合（解耦散点、与 B_mAb 的相关分析），
        并且在结果解释时要意识到它偏"深度"而非"封锁"。

    返回
    ----
    dict: b_cell_field (n,)、reachable (n,)
    """
    source = _as_idx(source)
    n = A.shape[0]
    if len(source) == 0:
        return dict(b_cell_field=np.full(n, np.nan), reachable=np.zeros(n, bool))

    pairs = edge_pairs(A, upper_only=True)
    z = a - b_ecm * _pair_mean(ecm, pairs) - c_caf * _pair_mean(caf, pairs)
    cap = _sigmoid(z)
    cost = 1.0 / np.maximum(cap, 1e-9)     # 通径越小，走过去越费劲

    C = sp.coo_matrix(
        (np.concatenate([cost, cost]),
         (np.concatenate([pairs[:, 0], pairs[:, 1]]),
          np.concatenate([pairs[:, 1], pairs[:, 0]]))),
        shape=(n, n),
    ).tocsr()

    d = dijkstra(C, directed=False, indices=source, min_only=True)
    reachable = np.isfinite(d)
    if reachable.any():
        d[~reachable] = d[reachable].max()
    else:
        d[:] = np.nan
    return dict(b_cell_field=d, reachable=reachable)


# --------------------------------------------------------------------------
# B_mAb —— 抗体传质屏障（尺寸依赖有效阻抗 + 结合位点消耗）
# --------------------------------------------------------------------------
def compute_b_mab(
    A: sp.spmatrix,
    ecm: np.ndarray,
    crosslink: np.ndarray,
    ag_target: np.ndarray,
    vessel: Sequence[int],
    *,
    r_nm: float = 5.5,
    g0: float = 1.0,
    lam: float = 3.0,
    xi0_nm: float = 20.0,
    beta: float = 3.0,
    kd_eff: float = 0.5,
    kappa_w: float = 1.0,
    g_floor: float = 1e-6,
) -> dict:
    """抗体传质屏障，逐 spot 输出。

    直觉
    ----
    把组织想成一张电阻网络。每条边的"电导"表示该处对一个特定尺寸分子的通透性：
      · ECM 越致密，电导越低（指数衰减）；
      · 局部网孔尺寸 xi 越小、分子半径 r 越大，尺寸排阻越严重，电导越低。
    在这张网络上求解【扩散-吸收方程】（屏蔽泊松）：
          (L + diag(k)) phi = 0 ，血管节点处 phi = 1
    L 是电导拉普拉斯矩阵（扩散），k 是逐节点吸收率（结合位点屏障——
    高抗原细胞把抗体结合内吞掉）。解出的 phi 就是稳态浓度场，
    B_mAb = -log(phi)：浓度越低，屏障越高。

    当 k=0 时该方程退化为调和场，正是标准的电阻网络类比，图论中的
    有效阻抗是它的一个特例。所以"图论传输阻力"这个叙事完全成立，
    而且这里用的是更完整的那个版本。

    关键参数 r_nm 就是"你在算哪种分子"：
      r_nm = 5.5  -> IgG / ADC（约 150 kDa 抗体，文献锚定值，不参与拟合）
      r_nm = 0.5  -> 小分子药物
    把 r_nm 从 0.5 扫到 10 就是 S3 分子尺寸扫描实验。

    参数
    ----
    A         : (n, n) 稀疏对称邻接矩阵
    ecm       : (n,) core matrisome 分数，[0, 1]
    crosslink : (n,) 交联酶（LOX 家族等）分数，[0, 1]，用于估计局部网孔尺寸
    ag_target : (n,) 模型吸收项的表达代理，[0, 1]。当前取 CD274/PDCD1LG2
                （PD-L1/PD-L2 配体），并非 PD-1 受体或抗 PD-1 药物特异性靶点。
                此代理及 kd_eff 未经蛋白浓度/结合动力学标定。
    vessel    : 血管节点索引（源）
    r_nm      : 分子流体力学半径（纳米）
    g0        : 基线电导
    lam       : ECM 对电导的指数衰减系数
    xi0_nm    : 无交联时的基线网孔尺寸（纳米）
    beta      : 交联对网孔尺寸的收缩系数
    kd_eff    : 有效解离常数，控制结合位点消耗的强度。
                kappa(v) = Ag(v) / (Ag(v) + kd_eff)，即该节点截留抗体的比例
    kappa_w   : 消耗项的权重。沿路径的累积吸收按 -log(1-kappa_w*kappa) 计
    g_floor   : 电导下限，避免图断开导致数值爆炸

    返回
    ----
    dict，含：
      b_mab      : (n,) 逐 spot 屏障值（越大越难到达）
      reachable  : (n,) bool，该节点是否与任一血管节点连通
      r_nm       : 回显所用的分子半径，便于结果归档

    注意
    ----
    与血管不连通的节点（孤立组织碎片）其 b_mab 在数学上是无穷大。本实现把它们
    赋为有限节点中的最大值，并在 reachable 中标记为 False，由调用方决定是否
    纳入统计。不要静默地把它们当成正常值。
    """
    vessel = _as_idx(vessel)
    n = A.shape[0]
    if len(vessel) == 0:
        return dict(b_mab=np.full(n, np.nan), reachable=np.zeros(n, bool), r_nm=r_nm)

    pairs = edge_pairs(A, upper_only=True)
    ecm_e = _pair_mean(ecm, pairs)
    xl_e = _pair_mean(crosslink, pairs)

    # 局部有效网孔尺寸：交联越强，网孔越小
    xi = xi0_nm * np.exp(-beta * xl_e)
    s = r_nm / np.maximum(xi, 1e-6)
    # 尺寸排阻因子：分子半径接近或超过网孔尺寸时通透性归零
    phi = np.where(s < 1.0, (1.0 - s) ** 2, 0.0)

    g = g0 * np.exp(-lam * ecm_e) * phi + g_floor

    # 组装电导矩阵（对称）
    Gw = sp.coo_matrix(
        (np.concatenate([g, g]),
         (np.concatenate([pairs[:, 0], pairs[:, 1]]),
          np.concatenate([pairs[:, 1], pairs[:, 0]]))),
        shape=(n, n),
    ).tocsr()

    # 连通性检查：只在包含血管的连通分量内求解
    n_comp, labels = connected_components(Gw, directed=False)
    vessel_comps = set(labels[vessel].tolist())
    reachable = np.isin(labels, list(vessel_comps))

    # ---- 扩散-吸收方程（屏蔽泊松 / screened Poisson）----
    # 抗体在组织中的稳态浓度场 phi 满足
    #       (L + diag(k)) phi = 0，血管节点处 phi = 1（Dirichlet 边界）
    # 其中 L 是电导拉普拉斯矩阵（描述扩散），k 是逐节点的吸收率
    # （描述结合位点屏障：高抗原细胞把抗体结合内吞掉）。
    # B_mAb = -log(phi)：浓度越低，屏障越高。
    #
    # 为什么用这个而不是"沿最短路径累积吸收"
    # ------------------------------------
    # 早期实现用 Dijkstra 求"最少被吸收的那条路线"。问题是：抗体可以
    # **免费绕开**高抗原区域——只要绕路上的吸收低，路径再长代价也不增加。
    # 结果是一个扇区的高抗原完全挡不住深部，与真实的扩散物理不符。
    # 扩散-吸收方程没有这个毛病：绕路本身要付出扩散阻力，两者自动权衡。
    # 而且当 k=0 时它退化为调和场，正是标准的电阻网络类比，
    # 有效阻抗是它的一个特例——所以论文里"图论传输阻力"的叙事完全保留。
    #
    # 计算量：一次稀疏线性求解，n<=3000 时秒级。
    kappa = ag_target / (ag_target + kd_eff + _EPS)
    k_abs = np.clip(kappa_w * kappa, 0.0, None)

    b_mab = np.full(n, np.nan)
    for comp in vessel_comps:
        idx = np.where(labels == comp)[0]
        if len(idx) < 2:
            continue
        sub = Gw[idx][:, idx]
        L = laplacian(sub).tocsr()
        M = (L + sp.diags(k_abs[idx] + 1e-9)).tocsr()

        is_v = np.zeros(len(idx), bool)
        loc = {g_: i for i, g_ in enumerate(idx)}
        for x in vessel:
            if labels[x] == comp:
                is_v[loc[x]] = True
        if not is_v.any():
            continue

        free = ~is_v
        if not free.any():
            b_mab[idx] = 0.0
            continue

        # (L+diag(k)) phi = 0 ，phi[vessel]=1  ->  M_ff phi_f = -M_fv * 1
        M_ff = M[free][:, free]
        rhs = -np.asarray(M[free][:, is_v].sum(axis=1)).ravel()
        try:
            phi_f = spsolve(M_ff.tocsc(), rhs)
        except Exception:  # noqa: BLE001 —— 退化情形改用迭代解
            phi_f, _ = cg(M_ff.tocsr(), rhs, rtol=1e-8, maxiter=5000)

        phi = np.ones(len(idx))
        phi[free] = np.clip(np.nan_to_num(phi_f, nan=1e-12), 1e-12, 1.0)
        b_mab[idx] = -np.log(phi)

    # 不可达节点赋为最大值并标记
    if np.any(np.isfinite(b_mab)):
        b_mab[~np.isfinite(b_mab)] = np.nanmax(b_mab[np.isfinite(b_mab)])
    return dict(b_mab=b_mab, reachable=reachable, r_nm=r_nm)


# --------------------------------------------------------------------------
# B_meta —— 小分子代谢庇护
# --------------------------------------------------------------------------
def compute_b_meta(
    A: sp.spmatrix,
    hypoxia: np.ndarray,
    proliferation: np.ndarray,
    efflux: np.ndarray,
    vessel: Sequence[int],
    *,
    D: sp.spmatrix | None = None,
    spacing_um: float = 100.0,
    d0_um: float = 130.0,
    tau_um: float = 40.0,
    w_hyp: float = 1.0,
    w_quiesc: float = 1.0,
    w_efflux: float = 1.0,
) -> dict:
    """小分子代谢庇护，逐 spot 输出。

    直觉
    ----
    小分子扩散能力强，基本"到得了"，问题是"到了没用"：距最近功能性血管
    100–150 μm 之外形成缺氧区，细胞进入静止期从而对抗有丝分裂药物不敏感；
    加上外排泵把药物泵出去。所以这一项 = 距血管有多远（几何）× 局部状态有多差（生物）。

    参数
    ----
    A             : (n, n) 稀疏对称邻接矩阵
    hypoxia       : (n,) 缺氧签名分数，[0, 1]
    proliferation : (n,) 增殖签名分数，[0, 1]
    efflux        : (n,) 外排泵签名分数，[0, 1]
    vessel        : 血管节点索引
    D             : (n, n) 稀疏距离矩阵，边权以**微米**计（build_graph_* 的第二个
                    返回值，load_graph 也直接给）。传了就用加权图距，这是正式口径。
    spacing_um    : 相邻 spot 中心间距（Visium 约 100 μm，第一代 ST 约 200 μm）。
                    **仅在 D 为 None 时**用来把跳数换算成微米（合成图没有 D）
    d0_um         : sigmoid 拐点，对应氧扩散极限（文献值 100–150 μm）
    tau_um        : sigmoid 陡峭度
    w_*           : 三个状态项的权重，默认等权

    返回
    ----
    dict，含 b_meta (n,)、d_vessel_um (n,)、d_vessel_mode（"weighted" 或 "hops"）

    为什么不再用 hops × spacing_um（2026-08-27 改）
    -----------------------------------------------
    跳数近似的偏差**随平台阵列几何反号**，实测（hops×spacing ÷ 直线距离，中位数）：
      · Visium 六方阵列 r=150   -> 1.109（高估 11%）
      · 第一代 ST 交错阵列 r=300 -> 0.894（低估 11%）
    而 d0_um=130 锚的是直线氧扩散极限，于是两个队列的 B_meta 之间凭空多了
    约 22% 的系统偏移——R3b 比的正是队列差异，这个偏移会被读成生物学差异。
    加权图距把这个平台相关常数整个去掉：Visium 侧 d 只动 0.4–1.6%
    （psi 均值 0.4217->0.4297 与 0.4629->0.4640，逐点 |Δpsi| 中位 0.0009–0.0089），
    第一代 ST 侧则一次到位。
    """
    vessel = _as_idx(vessel)
    n = A.shape[0]
    if len(vessel) == 0:
        return dict(b_meta=np.full(n, np.nan), d_vessel_um=np.full(n, np.nan))

    # 距最近血管的距离。正式口径是加权图距（D 的边权已是微米）；
    # 只有拿不到 D 的场合（合成图、老的单元测试）才退回 hops × spacing_um。
    if D is not None:
        d_um = dijkstra(D, directed=False, indices=vessel, min_only=True)
        mode = "weighted"
    else:
        Ab = A.copy().astype(float)
        Ab.data[:] = 1.0
        d_um = dijkstra(Ab, directed=False, indices=vessel, min_only=True) * spacing_um
        mode = "hops"
    finite = np.isfinite(d_um)
    if finite.any():
        d_um[~finite] = d_um[finite].max() + spacing_um   # 不连通的点放到最远处之外
    else:
        d_um[:] = 0.0

    psi = _sigmoid((d_um - d0_um) / max(tau_um, 1e-6))
    state = (
        w_hyp * hypoxia
        + w_quiesc * (1.0 - proliferation)
        + w_efflux * efflux
    ) / max(w_hyp + w_quiesc + w_efflux, _EPS)

    return dict(b_meta=psi * state, d_vessel_um=d_um, d_vessel_mode=mode)


# --------------------------------------------------------------------------
# 一站式入口
# --------------------------------------------------------------------------
def compute_all_barriers(A, scores: dict, source, sink, vessel, cfg: dict, D=None) -> dict:
    """按配置一次算出三个分量。

    参数
    ----
    scores : dict，键需包含 ecm / caf / crosslink / ag_target /
             hypoxia / proliferation / efflux，值均为 (n,) 且已秩标准化到 [0,1]
    cfg    : 配置字典，结构见 configs/default.yaml 的 barrier 段
    D      : (n, n) 微米距离矩阵。**真实切片一律要传**——不传的话 B_meta 会退回
             跳数近似，其偏差随平台阵列几何反号（详见 compute_b_meta 的说明）
    """
    bc = compute_b_cell(A, scores["ecm"], scores["caf"], source, sink, **cfg["b_cell"])
    bm = compute_b_mab(
        A, scores["ecm"], scores["crosslink"], scores["ag_target"], vessel, **cfg["b_mab"]
    )
    bt = compute_b_meta(
        A, scores["hypoxia"], scores["proliferation"], scores["efflux"], vessel,
        D=D, **cfg["b_meta"]
    )
    return {**bc, **bm, **bt}


def save_barriers(path_npz: str, res: dict, path_cut_json: str | None = None) -> None:
    """把屏障结果落盘。cut_edges 是变长列表，单独存 JSON。"""
    np.savez_compressed(
        path_npz,
        b_cell=np.array([res.get("b_cell", np.nan)]),
        max_flow=np.array([res.get("max_flow", np.nan)]),
        b_mab=res.get("b_mab", np.array([])),
        b_meta=res.get("b_meta", np.array([])),
        d_vessel_um=res.get("d_vessel_um", np.array([])),
        reachable=res.get("reachable", np.array([])),
        cut_nodes=res.get("cut_nodes", np.array([], int)),
    )
    if path_cut_json:
        with open(path_cut_json, "w", encoding="utf-8") as f:
            json.dump(
                {"cut_edges": [[int(u), int(v)] for u, v in res.get("cut_edges", [])]},
                f, ensure_ascii=False, indent=1,
            )


# --------------------------------------------------------------------------
# AnnData 适配器（薄封装，不含数学）
# --------------------------------------------------------------------------
_SCORE_KEYS = {
    "ecm": "ECM_core_n",
    "caf": "CAF_n",
    "crosslink": "ECM_crosslink_n",
    "ag_target": "Ag_target_n",
    "hypoxia": "Hypoxia_n",
    "proliferation": "Proliferation_n",
    "efflux": "Efflux_n",
}


def scores_from_adata(adata, keys: dict | None = None,
                      missing_out: list | None = None) -> dict:
    """从 adata.obs 中取出秩标准化后的签名分数，转成核心函数需要的数组字典。

    缺失的键会用 0.5（中性值）填充并打印告警——不要静默失败，那会让你在
    分析阶段拿到一堆看似合理其实无意义的数字。

    参数
    ----
    missing_out : 可选 list。若传入，被中性填充的键名（如 ``"ag_target"``）会被
                  append 进去，供下游脚本写进结果 JSON，使"这张切片的某个屏障
                  分量是常数 0.5"这件事在产物里可追溯。仅打印 warning 不够——
                  事后没人会回查日志，而 0.5 常数填充在结果里长得跟生物学阴性一模一样。
    """
    keys = keys or _SCORE_KEYS
    out = {}
    for k, col in keys.items():
        if col in adata.obs:
            out[k] = np.asarray(adata.obs[col].values, dtype=float)
        else:
            print(f"[warn] adata.obs 缺少 '{col}'，用中性值 0.5 填充。请检查 M2 打分是否完成。")
            out[k] = np.full(adata.n_obs, 0.5)
            if missing_out is not None:
                missing_out.append(k)
    return out
