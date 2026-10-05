"""
逐 spot 反事实干预排序 —— 把"测量屏障"升级为"定位屏障"
==========================================================

上游依赖：barrier.py（两个算子）、graph.py（A、源汇、血管）
下游产出：results/validation/intervention_ranking.json（run_44）

动机
----
compute_b_cell / compute_b_mAb 回答的是"这张切片挡不挡"。本模块回答的是
"**哪个位置**挡"：对每个 spot，把它周围的边"完全打通"（把该处的阻力物质
降到切片低分位水平，与 S2 的 _ablate 同一口径），重算两个屏障，看各降多少。
对全部 spot 循环一遍，就得到一张按"干预价值"排序的表——反事实意义上的
靶点排序。

设计口径（与 S2 保持一致，避免两套干预语义）
--------------------------------------------
· 干预 = 把该 spot 的 ecm / caf / crosslink 降到全片 low_q 分位
  （"把这些位置的屏障物质去掉"，不是"把组织挖掉"；前者有治疗学对应物
  ——基质正常化、LOX 抑制——后者没有）。
· ag_target（吸收项）**不动**：干预的是递送屏障，不是靶点表达。
· B_cell 的干预读数 = 切片级最小割标量的相对下降；
  B_mAb 的干预读数 = reachable 集合上 b_mab 均值的相对下降
  （与 run_23 stromal_intervention 同一口径，两者可直接对比）。

两个自带的理论锚点（不是对照，是刻度）
--------------------------------------
· full_cut_breach：把整条最小割封锁带上的物质全部去掉后的 B_cell 下降。
  top-k 单点干预的下降可以表示成它的百分比——"打通最有效的 k 个点，
  相当于把封锁线撕开多大比例的口子"。
· all_material_removed：把全片阻力物质都去掉后的 B_mAb 下降
  （无屏障极限）。单点干预读数同样可以表示成它的百分比。

复杂度
------
每张切片 n 次 min-cut + n 次稀疏线性求解。实测（2026-10，CPU）：
n=664 -> 0.03 s/次；n=3798 -> 0.15 s/次。22 张切片全量约 0.5–1 小时，
无需任何候选池裁剪——全 spot 测试本身就构成经验零分布
（"随机选一个点干预"的效果分布），不需要另设随机对照。

输出语义警示（论文里必须写）
----------------------------
排序结果**是模型内的反事实**，不是实测通透性；它说的是"在这个传输模型里，
把这些位置的基质参数设为低值会降低屏障分数多少"。Δ 的绝对量随参数
（b_ecm / lam / beta 等）缩放，跨切片比较用相对值（fraction of
full_cut_breach / of all_material_removed）。

一个必须知道的方向性细节（B_mAb 会出现少量负 Δ）
-------------------------------------------------
B_cell 的 Δ 严格非负：打通只抬高边容量，最大流不减。
B_mAb 则可能出现极小的负 Δ（合成图上约 1e-4 × max）：打通高吸收区附近的
电导会让抗体更多地在吸收区被截留，使**别处**浓度略降、屏障略微上升。
这是屏蔽泊松（扩散-吸收）方程的真实性质而非数值噪声，因此不要 clip 掉，
报告时说明即可（tests/test_intervention.py 测试 4 锁住其量级）。
"""

from __future__ import annotations

import numpy as np

from .barrier import compute_b_cell, compute_b_mab

__all__ = ["spot_intervention_ranking", "joint_intervention_curve", "intervention_anchors"]

_EPS = 1e-12


def _ablate_one(x: np.ndarray, i: int, low_val: float) -> np.ndarray:
    """单点干预：把第 i 个值换成 low_val（与 S2 的 _ablate 同一口径）。返回副本。"""
    out = x.copy()
    out[i] = low_val
    return out


def spot_intervention_ranking(
    A,
    ecm: np.ndarray,
    caf: np.ndarray,
    crosslink: np.ndarray,
    ag_target: np.ndarray,
    source,
    sink,
    vessel,
    cfg_cell: dict,
    cfg_mab: dict,
    *,
    low_q: float = 0.05,
) -> dict:
    """逐 spot 反事实干预排序（B_cell 与 B_mAb 各一列）。

    参数
    ----
    其余同 compute_b_cell / compute_b_mab。ecm / caf / crosslink / ag_target
    均为 (n,) 且已秩标准化到 [0, 1]。
    low_q : 干预时把阻力物质降到这个分位（默认 0.05，与 S2 一致）。

    返回
    ----
    dict：
      b_cell_base / b_mab_base     : 基线屏障
      delta_b_cell (n,)            : 干预后 B_cell 的绝对下降（>= 0 期望）
      delta_b_mab (n,)             : 干预后 B_mAb 均值的绝对下降
      full_cut_breach              : 整条割集物质全去掉后的 B_cell 下降
      all_material_removed_mab     : 全片阻力物质去掉后的 B_mAb 下降
      cut_nodes                    : 基线最小割节点（用于富集分析）
    """
    n = A.shape[0]

    # ---- 基线 ----
    base_cell = compute_b_cell(A, ecm, caf, source, sink, **cfg_cell)
    b_cell0 = base_cell["b_cell"]
    cut_nodes = np.asarray(base_cell.get("cut_nodes", []), int)

    base_mab = compute_b_mab(A, ecm, crosslink, ag_target, vessel, **cfg_mab)
    reachable = np.asarray(base_mab["reachable"], bool)
    b_mab0 = float(np.nanmean(base_mab["b_mab"][reachable])) if reachable.any() else np.nan

    # 干预目标值：全片 low_q 分位（对 ecm / caf / crosslink 各自取）
    q_ecm = float(np.quantile(ecm, low_q))
    q_caf = float(np.quantile(caf, low_q))
    q_xl = float(np.quantile(crosslink, low_q))

    delta_b_cell = np.zeros(n)
    delta_b_mab = np.zeros(n)

    if np.isfinite(b_cell0) or np.isfinite(b_mab0):
        for i in range(n):
            if np.isfinite(b_cell0):
                bc = compute_b_cell(A, _ablate_one(ecm, i, q_ecm),
                                    _ablate_one(caf, i, q_caf),
                                    source, sink, **cfg_cell)
                if np.isfinite(bc["b_cell"]):
                    delta_b_cell[i] = b_cell0 - bc["b_cell"]
            if np.isfinite(b_mab0):
                bm = compute_b_mab(A, _ablate_one(ecm, i, q_ecm),
                                   _ablate_one(crosslink, i, q_xl),
                                   ag_target, vessel, **cfg_mab)
                m = float(np.nanmean(bm["b_mab"][reachable]))
                if np.isfinite(m):
                    delta_b_mab[i] = b_mab0 - m

    anc = intervention_anchors(A, ecm, caf, crosslink, ag_target, source, sink, vessel,
                               cfg_cell, cfg_mab, low_q=low_q, cut_nodes=cut_nodes,
                               b_cell0=b_cell0, b_mab0=b_mab0, reachable=reachable)
    full_cut_breach = anc["full_cut_breach"]
    all_removed_cell = anc["all_removed_cell"]
    all_removed_mab = anc["all_removed_mab"]

    return dict(
        n=int(n),
        b_cell_base=float(b_cell0) if np.isfinite(b_cell0) else float("nan"),
        b_mab_base=b_mab0,
        delta_b_cell=delta_b_cell,
        delta_b_mab=delta_b_mab,
        full_cut_breach=float(full_cut_breach),
        all_removed_cell=float(all_removed_cell),
        all_removed_mab=float(all_removed_mab),
        cut_nodes=cut_nodes,
        low_q=float(low_q),
    )


def joint_intervention_curve(
    A,
    ecm: np.ndarray,
    caf: np.ndarray,
    crosslink: np.ndarray,
    ag_target: np.ndarray,
    source,
    sink,
    vessel,
    cfg_cell: dict,
    cfg_mab: dict,
    ranking_idx: np.ndarray,
    *,
    ks=(1, 5, 20, 50),
    low_q: float = 0.05,
    operator: str = "cell",
) -> dict:
    """按排序**联合**打通 top-k 个 spot，测联合效应（不是把独立效应累加）。

    为什么需要它
    ------------
    ``spot_intervention_ranking`` 的逐 spot Δ 是**各自单独**干预得到的，
    直接相加会高估（各点效果彼此重叠、且屏障是次可加的）。把 top-k 作为一个
    集合同时打通，才回答"按这个排序打到第 k 个，总共能撕开多大口子"。

    参数
    ----
    ranking_idx : (n,) 排序后的 spot 索引（从最有价值到最无价值），
                  即 np.argsort(-delta)。
    operator    : "cell" 或 "mab"，决定用哪个算子与哪组锚点。
    ks          : 要评估的 k 值。

    返回
    ----
    dict: ks, joint_drop（各 k 的绝对下降）, joint_frac_of_anchor,
          anchor（cell=full_cut_breach, mab=all_removed_mab）, baseline
    """
    n = A.shape[0]
    q_ecm = float(np.quantile(ecm, low_q))
    q_caf = float(np.quantile(caf, low_q))
    q_xl = float(np.quantile(crosslink, low_q))

    if operator == "cell":
        base = compute_b_cell(A, ecm, caf, source, sink, **cfg_cell)
        b0 = base["b_cell"]
        cut_nodes = np.asarray(base.get("cut_nodes", []), int)
        e2, c2 = ecm.copy(), caf.copy()
        if len(cut_nodes):
            e2[cut_nodes] = q_ecm
            c2[cut_nodes] = q_caf
        anchor = float(b0 - compute_b_cell(A, e2, c2, source, sink, **cfg_cell)["b_cell"])
    elif operator == "mab":
        base = compute_b_mab(A, ecm, crosslink, ag_target, vessel, **cfg_mab)
        reach = np.asarray(base["reachable"], bool)
        b0 = float(np.nanmean(base["b_mab"][reach]))
        e2 = np.full(n, q_ecm)
        x2 = np.full(n, q_xl)
        bm = compute_b_mab(A, e2, x2, ag_target, vessel, **cfg_mab)
        anchor = float(b0 - float(np.nanmean(bm["b_mab"][reach])))
    else:
        raise ValueError("operator 只能是 'cell' 或 'mab'")

    out = []
    for k in ks:
        if k > n:
            continue
        idx = np.asarray(ranking_idx[:k], int)
        if operator == "cell":
            e2, c2 = ecm.copy(), caf.copy()
            e2[idx] = q_ecm
            c2[idx] = q_caf
            b = compute_b_cell(A, e2, c2, source, sink, **cfg_cell)["b_cell"]
            drop = float(b0 - b) if np.isfinite(b0) and np.isfinite(b) else float("nan")
        else:
            e2, x2 = ecm.copy(), crosslink.copy()
            e2[idx] = q_ecm
            x2[idx] = q_xl
            bm = compute_b_mab(A, e2, x2, ag_target, vessel, **cfg_mab)
            m = float(np.nanmean(bm["b_mab"][reach]))
            drop = float(b0 - m)
        out.append(dict(
            k=int(k),
            k_over_n=float(k / n),
            joint_drop=drop,
            joint_frac_of_anchor=(float(drop / anchor)
                                  if np.isfinite(anchor) and anchor > 0 else float("nan")),
        ))
    return dict(operator=operator, baseline=float(b0), anchor=float(anchor),
                low_q=float(low_q), points=out)


def intervention_anchors(A, ecm, caf, crosslink, ag_target, source, sink, vessel, cfg_cell, cfg_mab,
                         *, low_q: float = 0.05, cut_nodes=None, b_cell0=None, b_mab0=None,
                         reachable=None) -> dict:
    """两个刻度锚点（不需要逐 spot 循环，可单独刷新）。

    full_cut_breach  : 最小割封锁带上全部节点的物质降到 low_q 分位后，B_cell 的下降
    all_removed_cell : 全片物质降到 low_q 分位后，B_cell 的下降
    all_removed_mab  : 全片物质降到 low_q 分位后，reachable 上 B_mAb 均值的下降
    未传入的基线量在这里现算。cut_nodes 必须来自精确最小割（barrier.exact_min_cut）。
    """
    n = A.shape[0]
    if b_cell0 is None or cut_nodes is None:
        base = compute_b_cell(A, ecm, caf, source, sink, **cfg_cell)
        b_cell0 = base["b_cell"]
        cut_nodes = np.asarray(base.get("cut_nodes", []), int)
    if b_mab0 is None or reachable is None:
        bm0 = compute_b_mab(A, ecm, crosslink, ag_target, vessel, **cfg_mab)
        reachable = np.asarray(bm0["reachable"], bool)
        b_mab0 = float(np.nanmean(bm0["b_mab"][reachable])) if reachable.any() else np.nan
    q_ecm = float(np.quantile(ecm, low_q))
    q_caf = float(np.quantile(caf, low_q))
    q_xl = float(np.quantile(crosslink, low_q))
    full_cut_breach = np.nan
    if np.isfinite(b_cell0) and len(cut_nodes):
        ecm2, caf2 = ecm.copy(), caf.copy()
        ecm2[cut_nodes] = q_ecm
        caf2[cut_nodes] = q_caf
        bc = compute_b_cell(A, ecm2, caf2, source, sink, **cfg_cell)
        if np.isfinite(bc["b_cell"]):
            full_cut_breach = float(b_cell0 - bc["b_cell"])
    ecm_all = np.full(n, q_ecm)
    caf_all = np.full(n, q_caf)
    xl_all = np.full(n, q_xl)
    all_removed_cell = np.nan
    if np.isfinite(b_cell0):
        bc = compute_b_cell(A, ecm_all, caf_all, source, sink, **cfg_cell)
        if np.isfinite(bc["b_cell"]):
            all_removed_cell = float(b_cell0 - bc["b_cell"])
    all_removed_mab = np.nan
    if np.isfinite(b_mab0):
        bm = compute_b_mab(A, ecm_all, xl_all, ag_target, vessel, **cfg_mab)
        m = float(np.nanmean(bm["b_mab"][reachable]))
        if np.isfinite(m):
            all_removed_mab = float(b_mab0 - m)
    return dict(full_cut_breach=float(full_cut_breach), all_removed_cell=float(all_removed_cell),
                all_removed_mab=float(all_removed_mab), cut_nodes=np.asarray(cut_nodes, int),
                b_cell0=float(b_cell0), b_mab0=float(b_mab0))
