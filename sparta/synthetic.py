"""
合成图生成器 —— 用于单元测试与教学
=====================================

上游依赖：无
下游用途：tests/test_barrier.py；也可在 notebooks 里用来直观理解三个算子在做什么

为什么要有这个模块
------------------
屏障算子写错了，在真实数据上几乎不可能被察觉——你会得到一堆看起来合理的数字，
然后基于它们写完整篇论文。合成图是唯一能在动真实数据之前发现这类错误的手段：
这里的每张图我们都**事先知道正确答案**。

四张关键的合成图
----------------
  make_ring_grid(gap=0)      闭合的高阻抗环带      -> 屏障应该很强
  make_ring_grid(gap=5)      有缺口的环带          -> 屏障应该明显下降（拓扑敏感性）
  make_scattered_grid()      同等数量的散在阻抗点  -> 屏障应该很弱
  make_dissociated_grid()    两道屏障解离的场景    -> 阳性对照（见该函数文档）

如果算子分不出前三者，它就退化成了局部统计量，整个方法学论证落空。

关于交联度
----------
各生成器里的交联度设为"与 ECM 相关但不等同"（约 r=0.7），因为真实组织里
LOX 家族与胶原相关但不是同一个量。写成 crosslink = ecm.copy() 会让
validate.dissociation_drivers 的共线性检查判为完全共线，分解结果退化。
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

__all__ = ["make_grid", "make_ring_grid", "make_scattered_grid",
           "make_dissociated_grid", "SyntheticSlide"]


class SyntheticSlide:
    """一张合成切片。字段名与真实数据的用法保持一致，便于直接喂给 barrier.py。"""

    def __init__(self, coords, A, ecm, caf, crosslink, ag, hyp, prol, eff, source, sink, vessel,
                 blocked_idx=None, spacing_um=100.0):
        self.coords = coords            # (n, 2) 微米坐标
        self.A = A                      # (n, n) 稀疏邻接
        self.ecm = ecm
        self.caf = caf
        self.crosslink = crosslink
        self.ag_target = ag
        self.hypoxia = hyp
        self.proliferation = prol
        self.efflux = eff
        self.source = source            # 免疫入口（网格外圈）
        self.sink = sink                # 瘤巢核心（网格中心）
        self.vessel = vessel            # 血管节点
        self.blocked_idx = blocked_idx if blocked_idx is not None else np.array([], int)
        self.spacing_um = spacing_um
        self.n = len(coords)

    @property
    def scores(self) -> dict:
        """打包成 barrier.compute_all_barriers 需要的分数字典。"""
        return dict(
            ecm=self.ecm, caf=self.caf, crosslink=self.crosslink,
            ag_target=self.ag_target, hypoxia=self.hypoxia,
            proliferation=self.proliferation, efflux=self.efflux,
        )


def make_grid(n: int = 25, spacing_um: float = 100.0, connectivity: int = 4):
    """构造 n×n 规则网格。

    参数
    ----
    n            : 每边的节点数，总节点数为 n²
    spacing_um   : 相邻节点间距（微米）
    connectivity : 4 = 上下左右；8 = 含对角线

    返回
    ----
    coords (n², 2)、A 稀疏邻接矩阵
    """
    ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    coords = np.column_stack([ii.ravel(), jj.ravel()]).astype(float) * spacing_um

    if connectivity == 4:
        offs = [(0, 1), (1, 0)]
    elif connectivity == 8:
        offs = [(0, 1), (1, 0), (1, 1), (1, -1)]
    else:
        raise ValueError("connectivity 只支持 4 或 8")

    rows, cols = [], []
    idx = lambda i, j: i * n + j
    for i in range(n):
        for j in range(n):
            for di, dj in offs:
                a, b = i + di, j + dj
                if 0 <= a < n and 0 <= b < n:
                    rows += [idx(i, j), idx(a, b)]
                    cols += [idx(a, b), idx(i, j)]
    A = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n * n, n * n)).tocsr()
    A.data[:] = 1.0
    return coords, A


def _shell(coords, center, radius_um, width_um):
    """返回位于 [radius-width/2, radius+width/2] 环形壳层内的节点索引。"""
    d = np.linalg.norm(coords - center, axis=1)
    return np.where(np.abs(d - radius_um) <= width_um / 2.0)[0]


def _base_slide(n, spacing_um, connectivity):
    coords, A = make_grid(n, spacing_um, connectivity)
    N = len(coords)
    center = coords.mean(axis=0)
    d = np.linalg.norm(coords - center, axis=1)

    # 源：网格最外一圈（模拟组织边缘的血管入口）
    on_border = (
        (coords[:, 0] == coords[:, 0].min()) | (coords[:, 0] == coords[:, 0].max())
        | (coords[:, 1] == coords[:, 1].min()) | (coords[:, 1] == coords[:, 1].max())
    )
    source = np.where(on_border)[0]
    # 汇：最中心的一小块（模拟瘤巢核心）
    sink = np.where(d <= 1.6 * spacing_um)[0]
    # 血管：与源相同（合成场景下简化）
    vessel = source.copy()
    return coords, A, N, center, source, sink, vessel


def make_ring_grid(
    n: int = 25,
    spacing_um: float = 100.0,
    radius_spots: float = 7.0,
    width_spots: float = 1.2,
    gap_spots: int = 0,
    high: float = 1.0,
    low: float = 0.0,
    connectivity: int = 4,
    seed: int = 0,
) -> SyntheticSlide:
    """构造一张"中心瘤巢被高阻抗环带包围"的合成切片。

    参数
    ----
    radius_spots : 环带半径，以 spot 为单位
    width_spots  : 环带厚度
    gap_spots    : 在环带上开一个缺口，缺口跨越的角度对应约这么多个 spot。
                   gap_spots=0 表示闭合环带
    high / low   : 环带内外的 ECM/CAF 分数

    返回
    ----
    SyntheticSlide
    """
    rng = np.random.default_rng(seed)
    coords, A, N, center, source, sink, vessel = _base_slide(n, spacing_um, connectivity)

    ring = _shell(coords, center, radius_spots * spacing_um, width_spots * spacing_um)

    if gap_spots > 0 and len(ring) > 0:
        # 按极角开一个连续缺口
        ang = np.arctan2(coords[ring, 1] - center[1], coords[ring, 0] - center[0])
        ang = np.mod(ang, 2 * np.pi)
        gap_frac = gap_spots / max(len(ring), 1)
        start = 0.0
        width = 2 * np.pi * gap_frac
        in_gap = (ang >= start) & (ang < start + width)
        ring = ring[~in_gap]

    ecm = np.full(N, low, dtype=float)
    ecm[ring] = high
    caf = ecm.copy()
    # 交联度与 ECM 相关但不等同（真实组织里 LOX 家族与胶原相关约 r=0.6–0.8）。
    # 早期版本写成 crosslink = ecm.copy()，会让 dissociation_drivers 的
    # 共线性检查判定为完全共线，分解结果退化——那不是真实情形。
    crosslink = np.clip(0.7 * ecm + 0.3 * rng.uniform(0, 1, N), 0, 1)

    ag = np.full(N, 0.2)
    ag[sink] = 0.8                       # 瘤巢抗原表达高
    hyp = np.clip(np.linalg.norm(coords - center, axis=1) / (n * spacing_um / 2), 0, 1)
    hyp = 1.0 - hyp                      # 越靠中心越缺氧
    prol = np.clip(rng.uniform(0.3, 0.7, N), 0, 1)
    eff = np.clip(rng.uniform(0.2, 0.5, N), 0, 1)

    return SyntheticSlide(coords, A, ecm, caf, crosslink, ag, hyp, prol, eff,
                          source, sink, vessel, blocked_idx=ring, spacing_um=spacing_um)


def make_scattered_grid(
    n: int = 25,
    spacing_um: float = 100.0,
    n_blocked: int | None = None,
    high: float = 1.0,
    low: float = 0.0,
    connectivity: int = 4,
    seed: int = 0,
    radius_spots: float = 7.0,
    width_spots: float = 1.2,
) -> SyntheticSlide:
    """构造对照切片：与环带**相同数量**的高阻抗节点，但随机散布。

    这是最关键的对照——它与环带图的细胞组成完全一致（同样多的高 ECM/CAF 节点），
    只有空间排布不同。如果算子分不出这两张图，说明它测的是组成而不是拓扑。
    """
    rng = np.random.default_rng(seed)
    coords, A, N, center, source, sink, vessel = _base_slide(n, spacing_um, connectivity)

    if n_blocked is None:
        ref = make_ring_grid(n=n, spacing_um=spacing_um, radius_spots=radius_spots,
                             width_spots=width_spots, connectivity=connectivity, seed=seed)
        n_blocked = len(ref.blocked_idx)

    # 不在源/汇上放阻抗点，避免平凡地把汇整个封死
    candidates = np.setdiff1d(np.arange(N), np.union1d(source, sink))
    blocked = rng.choice(candidates, size=min(n_blocked, len(candidates)), replace=False)

    ecm = np.full(N, low, dtype=float)
    ecm[blocked] = high
    caf = ecm.copy()
    crosslink = np.clip(0.7 * ecm + 0.3 * rng.uniform(0, 1, N), 0, 1)

    ag = np.full(N, 0.2)
    ag[sink] = 0.8
    hyp = 1.0 - np.clip(np.linalg.norm(coords - center, axis=1) / (n * spacing_um / 2), 0, 1)
    prol = np.clip(rng.uniform(0.3, 0.7, N), 0, 1)
    eff = np.clip(rng.uniform(0.2, 0.5, N), 0, 1)

    return SyntheticSlide(coords, A, ecm, caf, crosslink, ag, hyp, prol, eff,
                          source, sink, vessel, blocked_idx=blocked, spacing_um=spacing_um)


def make_dissociated_grid(
    n: int = 25,
    spacing_um: float = 100.0,
    radius_spots: float = 7.0,
    width_spots: float = 2.0,
    gap_frac: float = 0.30,
    ag_sector_frac: float = 0.30,
    rim_ag: float = 0.98,
    connectivity: int = 4,
    seed: int = 0,
) -> SyntheticSlide:
    """构造一张**两道屏障解离**的合成切片 —— 本课题中心主张的阳性对照。

    场景（这正是论文要描述的那类"热而无效"的肿瘤）
    ------------------------------------------------
    在同一个**扇区**内：
      · ECM/CAF 环带在这里开着缺口 -> T 细胞顺着缺口进入瘤巢，B_cell 低，
        病理上看起来是一个"热"区域；
      · 但这里的瘤缘一段**抗原表达极高** -> 到达的抗体在瘤缘就被结合内吞殆尽
        （结合位点屏障），深部几乎接触不到药物，B_mAb 高。
    其余扇区则相反：ECM 环带完整（T 细胞进不去），抗原低（抗体本可以进去）。

    为什么解离因素必须是**扇区**而不是整圈
    --------------------------------------
    这是在跑合成实验时发现的：如果把高抗原做成一整圈径向对称的壳层，它与
    "距血管的径向距离"完全共线，而解耦分析必须校正这个几何混杂
    （见 validate.decoupling_stats 的说明），于是解离信号会被校正过程一并
    抹掉，阳性对照失效。
    真实的 PD-L1 表达本来就是斑片状而非径向均匀的，所以扇区设计也更贴近现实。

    为什么必须有这张图
    ------------------
    在真实数据上声称"发现了解离"之前，必须先证明框架在**已知存在解离**的
    场景下能把它检出来。否则无法区分"真的没有解离"与"我的方法测不出解离"。
    这是一个必要的阳性对照，tests/test_barrier.py 中有对应断言。

    参数
    ----
    gap_frac       : ECM 环带缺口所占的角度比例
    ag_sector_frac : 高抗原瘤缘所占的角度比例（与缺口扇区重叠）
    rim_ag         : 该扇区瘤缘的抗原水平（模拟 PD-L1 高表达的血管周肿瘤细胞）
    """
    rng = np.random.default_rng(seed)
    coords, A, N, center, source, sink, vessel = _base_slide(n, spacing_um, connectivity)
    d = np.linalg.norm(coords - center, axis=1)
    ang = np.mod(np.arctan2(coords[:, 1] - center[1], coords[:, 0] - center[0]), 2 * np.pi)

    # 细胞屏障：环带在 [0, gap) 扇区开缺口 -> 该扇区 B_cell 低
    ring = _shell(coords, center, radius_spots * spacing_um, width_spots * spacing_um)
    gap_w = 2 * np.pi * gap_frac
    ring = ring[~((ang[ring] >= 0.0) & (ang[ring] < gap_w))]

    ecm = np.zeros(N); ecm[ring] = 1.0
    caf = ecm.copy()
    # 这张图里交联度**独立于** ECM（低且随机）——瓶颈是抗原而不是网孔。
    # 这正是"只影响抗体不影响细胞"的因素能否被检出的关键场景。
    crosslink = np.clip(rng.uniform(0.05, 0.25, N), 0, 1)

    # 抗体屏障：**同一个扇区**的瘤缘高抗原 -> 抗体在此被消耗
    rim = _shell(coords, center, (radius_spots - 2.0) * spacing_um, 1.6 * spacing_um)
    ag_w = 2 * np.pi * ag_sector_frac
    rim_hot = rim[(ang[rim] >= 0.0) & (ang[rim] < ag_w)]
    ag = np.full(N, 0.05)
    ag[rim_hot] = rim_ag

    hyp = 1.0 - np.clip(d / (n * spacing_um / 2), 0, 1)
    prol = np.clip(rng.uniform(0.3, 0.7, N), 0, 1)
    eff = np.clip(rng.uniform(0.2, 0.5, N), 0, 1)

    sl = SyntheticSlide(coords, A, ecm, caf, crosslink, ag, hyp, prol, eff,
                        source, sink, vessel, blocked_idx=ring, spacing_um=spacing_um)
    sl.hot_sector = np.where((ang >= 0.0) & (ang < min(gap_w, ag_w)))[0]
    return sl
