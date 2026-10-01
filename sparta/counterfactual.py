"""
M5 反事实实验 —— 本项目性价比最高的模块
=========================================

上游依赖：graph.py（A、源汇）、barrier.py（三个算子）
下游产出：results/counterfactual/*.json 与配套图件

这三个实验成本极低（CPU 分钟级），但分别封堵了审稿人最可能提出的三个致命质疑：

  S1 空间重排对照   -> "你的分数是不是只是细胞类型比例的复杂重写？"
  S2 环带断裂       -> "为什么非要用空间数据？"
  S3 分子尺寸扫描   -> "两个屏障凭什么说是解耦的？"

三者都不依赖任何外部数据，是自洽的内部证据。务必在主线分析跑通后立即完成，
不要拖到最后。
"""

from __future__ import annotations

import numpy as np

from .barrier import compute_b_cell, compute_b_mab

__all__ = ["s1_spatial_permutation", "s2_ring_breaking", "s3_size_scan"]


# --------------------------------------------------------------------------
# S1 空间重排对照
# --------------------------------------------------------------------------
def s1_spatial_permutation(
    A,
    ecm: np.ndarray,
    caf: np.ndarray,
    source,
    sink,
    cfg_cell: dict,
    *,
    n_perm: int = 500,
    seed: int = 0,
    mode: str = "fixed",
    endothelial: np.ndarray | None = None,
    t_nk: np.ndarray | None = None,
    malignant: np.ndarray | None = None,
    cfg_source_sink: dict | None = None,
) -> dict:
    """S1：保持细胞组成不变、只打乱空间位置，检验屏障是否来自空间排布。

    做法
    ----
    把每个 spot 的分数向量整体做随机置换（即"把这些细胞随机搬家"），
    组成完全不变，只有位置变了。重算屏障，构建零分布。
    若真实屏障显著高于零分布，说明屏障确实来自空间排布而非组成。

    两种模式（论文中两种都要报告）
    ------------------------------
    mode="fixed"  : 源汇节点位置保持不变，只置换 ECM/CAF。
                    直接检验"阻力物质的空间排布是否重要"，是最贴合本课题
                    核心主张的检验。
    mode="follow" : 置换全部分数（含内皮/T细胞/恶性），并从置换后的分数
                    重新推导源汇。这是更彻底的随机化，但源汇位置也随之变动，
                    因此零分布的方差更大。需要额外传入
                    endothelial / t_nk / malignant 与 cfg_source_sink。

    返回
    ----
    dict: b_real, null(数组), p_emp, z, mode, n_perm
    """
    rng = np.random.default_rng(seed)
    real = compute_b_cell(A, ecm, caf, source, sink, **cfg_cell)
    b_real = real["b_cell"]

    n = A.shape[0]
    null = []

    if mode == "follow":
        from .graph import define_source_sink

        if any(x is None for x in (endothelial, t_nk, malignant)):
            raise ValueError("mode='follow' 需要额外传入 endothelial / t_nk / malignant")
        cfg_ss = cfg_source_sink or {}

    for _ in range(n_perm):
        p = rng.permutation(n)
        if mode == "fixed":
            b = compute_b_cell(A, ecm[p], caf[p], source, sink, **cfg_cell)["b_cell"]
        elif mode == "follow":
            ss = define_source_sink(A, endothelial[p], t_nk[p], malignant[p], **cfg_ss)
            if ss["n_source"] == 0 or ss["n_sink"] == 0:
                continue
            b = compute_b_cell(A, ecm[p], caf[p], ss["source"], ss["sink"], **cfg_cell)["b_cell"]
        else:
            raise ValueError("mode 只能是 'fixed' 或 'follow'")
        if np.isfinite(b):
            null.append(b)

    null = np.asarray(null, dtype=float)
    if len(null) == 0:
        return dict(b_real=b_real, null=null, p_emp=np.nan, z=np.nan,
                    mode=mode, n_perm=0)

    # 单侧经验 p 值：真实屏障有多罕见地高
    p_emp = (np.sum(null >= b_real) + 1) / (len(null) + 1)
    z = (b_real - null.mean()) / (null.std() + 1e-12)
    return dict(b_real=float(b_real), null=null, p_emp=float(p_emp), z=float(z),
                mode=mode, n_perm=len(null))


# --------------------------------------------------------------------------
# S2 环带断裂
# --------------------------------------------------------------------------
def _ablate(scores: np.ndarray, idx, low_q: float = 0.05) -> np.ndarray:
    """把指定节点的阻力分数降到切片的低分位水平。

    注意这里是"把这些位置的屏障物质去掉"，不是"把组织挖掉"。
    前者是有治疗学意义的反事实（比如用药把 CAF 正常化），后者没有生物学对应物。
    """
    out = scores.copy()
    if len(idx):
        out[np.asarray(idx, int)] = np.quantile(scores, low_q)
    return out


def _nearest_within(coords, pool, seed_node, k):
    """在 pool 内取距离 seed_node 最近的 k 个节点（含 seed 自身），构成一段连续弧。"""
    d = np.linalg.norm(coords[pool] - coords[seed_node], axis=1)
    return np.asarray(pool)[np.argsort(d)[:k]]


def s2_ring_breaking(
    A,
    ecm: np.ndarray,
    caf: np.ndarray,
    source,
    sink,
    cfg_cell: dict,
    *,
    coords: np.ndarray | None = None,
    k_list=(3, 5, 8, 12),
    k_frac=None,
    n_rand: int = 200,
    seed: int = 0,
    low_q: float = 0.05,
) -> dict:
    """S2：在最小割上开一个**连续缺口** vs 移除同样多的、但分散的屏障物质。

    为什么必须是"连续缺口"
    ----------------------
    本课题的主张是"屏障是连续封锁线"。检验它的正确方式不是随机拿掉几个
    最小割上的节点——最小割集通常有几十个节点，随机拿掉 3–5 个大概率
    在环带各处打几个小洞，环带整体仍然闭合，屏障不会崩溃。
    正确的做法是拿掉**空间上连续的一段弧**，制造一个真正的缺口。
    （这一点是在合成图上跑实验时才暴露出来的——直接按"随机移除最小割节点"
      实现会得到与随机对照无异的结果，从而错误地否定拓扑假设。）

    两个对照
    --------
    对照 A（全局随机）：从全片随机取 k 个 spot 移除屏障物质。
        回答"随便去掉几个点会怎样"。
    对照 B（屏障内分散，关键对照）：从最小割集内随机取 k 个**不连续**的节点移除。
        材料完全相同、数量完全相同，只有"是否连续"不同。
        这个对照直接把"拓扑"从"组成"里分离出来，是 S2 最有说服力的部分。

    参数
    ----
    coords : (n, 2) 空间坐标。连续弧策略需要它；不给则退化为对照 B 的做法
             并在返回值中标记 contiguous=False。
    k_list : 移除的节点数（绝对值）。**只建议用于合成图与复现旧结果。**
    k_frac : 移除的节点数占割集的比例，如 (0.05, 0.10, 0.20, 0.30)。
             给了它就忽略 k_list。**真实切片一律用这个。**

    ⚠ k 必须随割集规模缩放（2026-08-26 审查发现）
    ---------------------------------------------
    合成图的割集只有 40–76 个节点，k=3–12 已经能开出一段真缺口（效应量 7.45x）。
    真实切片的割集有 150–750 个节点，同样的 k 只动了割集的 0.5%–8%，
    在 614 个节点的封锁带上拿掉连续的 3 个节点当然什么都不会发生——
    实测效应量随割集变大系统性衰减（割集 156 -> 1.07–1.16x；247 -> 1.09–1.30x；
    614 -> 1.01–1.03x 且全不显著）。
    这**不是拓扑假设不成立**，是实验规模没跟着割集缩放。
    详见 docs/review_for_journal.md §3 B2。

    效应量的读法
    ------------
    ratio_vs_in_cut = 分散移除后剩余屏障 / 连续缺口后剩余屏障。
    大于 1 表示"同样的材料，连起来比散开更能挡"。
    p_vs_in_cut 是更稳健的统计量：连续缺口比多大比例的等量分散移除更有效。

    一个必须知道的依赖关系
    ----------------------
    连续性的效应强度依赖于屏障带的**厚度**。一层薄的屏障带随便打几个洞就能
    穿过，连续性的优势不明显；多层厚的屏障带则必须开出一段连续缺口才能突破。
    真实的促纤维增生带在 Visium 分辨率下通常有数个 spot 厚，属于后者。
    分析真实数据时应报告屏障带的厚度分布，否则这个实验的效应大小无法解释。

    返回
    ----
    dict，per_k[k] 含 targeted / rand_global / rand_in_cut 三组结果与效应比
    """
    rng = np.random.default_rng(seed)
    base = compute_b_cell(A, ecm, caf, source, sink, **cfg_cell)
    b0 = base["b_cell"]
    cut_nodes = np.asarray(base["cut_nodes"], int)

    n = A.shape[0]
    protected = np.union1d(np.asarray(source, int), np.asarray(sink, int))
    pool_all = np.setdiff1d(np.arange(n), protected)
    cut_pool = np.setdiff1d(cut_nodes, protected)

    def _b(idx):
        return compute_b_cell(A, _ablate(ecm, idx, low_q), _ablate(caf, idx, low_q),
                              source, sink, **cfg_cell)["b_cell"]

    contiguous = coords is not None and len(cut_pool) > 0

    if k_frac is not None:
        ks = sorted({max(3, int(round(float(f) * len(cut_pool)))) for f in k_frac})
        k_mode = "frac"
    else:
        ks = sorted({int(k) for k in k_list})
        k_mode = "absolute"

    out = {"b0": float(b0), "n_cut_nodes": int(len(cut_nodes)),
           "n_cut_pool": int(len(cut_pool)), "k_mode": k_mode,
           "k_frac": [float(f) for f in k_frac] if k_frac is not None else None,
           "contiguous": bool(contiguous), "per_k": {}}

    for k in ks:
        if len(cut_pool) < k:
            continue

        # 定向：在最小割上开一段连续弧
        if contiguous:
            # 取多个种子各试一次，报告最有效的一段（即"最脆弱的弧"）
            seeds = rng.choice(cut_pool, size=min(8, len(cut_pool)), replace=False)
            cands = [_nearest_within(coords, cut_pool, s_, k) for s_ in seeds]
            vals = [_b(c) for c in cands]
            j = int(np.argmin(vals))
            tgt, b_t = cands[j], vals[j]
        else:
            tgt = rng.choice(cut_pool, size=k, replace=False)
            b_t = _b(tgt)

        # 对照 A：全局随机
        b_ra = np.array([_b(rng.choice(pool_all, size=k, replace=False))
                         for _ in range(n_rand)])
        # 对照 B：屏障内分散（关键对照）
        b_rc = np.array([_b(rng.choice(cut_pool, size=k, replace=False))
                         for _ in range(max(n_rand // 2, 20))])

        drop = lambda b: float((b0 - b) / (b0 + 1e-12))
        out["per_k"][int(k)] = dict(
            k_over_cut=float(k / max(len(cut_pool), 1)),
            targeted=float(b_t),
            drop_targeted=drop(b_t),
            rand_global_mean=float(b_ra.mean()),
            drop_rand_global=drop(b_ra.mean()),
            rand_in_cut_mean=float(b_rc.mean()),
            drop_rand_in_cut=drop(b_rc.mean()),
            # 效应量用"剩余屏障之比"而非"降幅之比"。
            # 降幅接近 100% 时比值会饱和：降 78% vs 降 54% 的降幅比只有 1.4，
            # 但剩余屏障 0.16 vs 0.34 实际相差 2.1 倍。后者才是可解释的量。
            ratio_vs_global=float(b_ra.mean() / (b_t + 1e-12)),
            ratio_vs_in_cut=float(b_rc.mean() / (b_t + 1e-12)),
            p_vs_in_cut=float((np.sum(b_rc <= b_t) + 1) / (len(b_rc) + 1)),
            nodes=[int(x) for x in tgt],
        )
    return out


# --------------------------------------------------------------------------
# S3 分子尺寸扫描
# --------------------------------------------------------------------------
def s3_size_scan(
    A,
    ecm: np.ndarray,
    crosslink: np.ndarray,
    ag_target: np.ndarray,
    vessel,
    cfg_mab: dict,
    *,
    radii_nm=(0.5, 1.0, 2.0, 4.0, 5.5, 8.0, 10.0),
    core_idx=None,
) -> dict:
    """S3：把分子半径从小分子扫到大分子，观察屏障场如何变化。

    这是"不同尺寸的分子面对不同的屏障"这一核心概念最直观的证据，
    实现成本几乎为零（改一个参数重跑），但它把整篇文章的物理主张
    可视化成了一条曲线。

    预期结果：小分子（r≈0.5 nm）的屏障场几乎不受空间结构影响；
    大分子（r≈5.5 nm，IgG）在瘤巢深部出现明显的低暴露区。

    参数
    ----
    core_idx : 瘤巢核心节点索引。若给出，额外报告核心区的屏障均值——
               这是最有临床意义的数字（"抗体到底能不能到达瘤巢深部"）。

    返回
    ----
    dict: radii, mean_all, mean_core, fields(每个半径的完整屏障场)
    """
    radii = list(radii_nm)
    fields, mean_all, mean_core = [], [], []
    for r in radii:
        res = compute_b_mab(A, ecm, crosslink, ag_target, vessel, r_nm=r, **cfg_mab)
        b = res["b_mab"]
        fields.append(b)
        mean_all.append(float(np.nanmean(b)))
        mean_core.append(float(np.nanmean(b[np.asarray(core_idx, int)]))
                         if core_idx is not None and len(core_idx) else np.nan)

    # 屏障对尺寸的敏感度：核心区屏障从最小半径到 IgG 半径的增幅
    idx_igg = int(np.argmin([abs(r - 5.5) for r in radii]))
    sensitivity = (mean_core[idx_igg] - mean_core[0]) if core_idx is not None else np.nan

    return dict(
        radii_nm=radii,
        mean_all=mean_all,
        mean_core=mean_core,
        fields=np.array(fields),
        igg_vs_small_core_delta=float(sensitivity) if np.isfinite(sensitivity) else np.nan,
    )
