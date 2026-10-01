"""
barrier.py 的合成图单元测试 —— 全项目唯一不可省略的一环
=========================================================

运行方式（两种都行）：
    pytest tests/test_barrier.py -v          # 有 pytest 时
    python tests/test_barrier.py             # 没装 pytest 时也能跑

这三个测试必须全绿，才允许进入真实数据分析。

它们分别检验：
  测试 1  能否区分「连续环带」与「同等数量的散在阻抗点」
          -> 检验算子测的是拓扑而非组成
  测试 2  能否感知环带上的缺口（最核心）
          -> 检验拓扑敏感性。若不敏感，算子已退化为局部统计量
  测试 3  B_mAb 是否随分子尺寸单调变化
          -> 检验尺寸依赖项实现正确，这是 B_cell 与 B_mAb 解耦的物理基础
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from sparta.barrier import compute_b_cell, compute_b_mab, compute_b_meta  # noqa: E402
from sparta.synthetic import (make_ring_grid, make_scattered_grid,        # noqa: E402
                              make_dissociated_grid)

CFG_CELL = dict(a=3.0, b_ecm=8.0, c_caf=4.0)
CFG_MAB = dict(g0=1.0, lam=3.0, xi0_nm=20.0, beta=3.0, kd_eff=0.5, kappa_w=1.0)


# ---------------------------------------------------------------------------
def test_ring_vs_scattered():
    """测试 1：连续环带的屏障必须显著强于同等数量的散在阻抗节点。"""
    ring = make_ring_grid(n=25, gap_spots=0, seed=0)
    scat = make_scattered_grid(n=25, seed=0)

    # 前提：两张图的高阻抗节点数量必须一致，否则这个对照不成立
    assert len(ring.blocked_idx) == len(scat.blocked_idx), (
        f"对照不成立：环带 {len(ring.blocked_idx)} 个阻抗节点，"
        f"散在 {len(scat.blocked_idx)} 个"
    )

    r = compute_b_cell(ring.A, ring.ecm, ring.caf, ring.source, ring.sink, **CFG_CELL)
    s = compute_b_cell(scat.A, scat.ecm, scat.caf, scat.source, scat.sink, **CFG_CELL)

    ratio = r["b_cell"] / s["b_cell"]
    print(f"  环带 b_cell = {r['b_cell']:.4f} (max_flow={r['max_flow']:.3f})")
    print(f"  散在 b_cell = {s['b_cell']:.4f} (max_flow={s['max_flow']:.3f})")
    print(f"  比值 = {ratio:.1f}x   阻抗节点数 = {len(ring.blocked_idx)}")

    assert ratio > 5.0, f"环带屏障仅为散在的 {ratio:.1f} 倍，算子对空间排布不敏感"
    assert len(r["cut_edges"]) > 0, "闭合环带应产生非空的最小割边集"


# ---------------------------------------------------------------------------
def test_ring_with_gap():
    """测试 2（最核心）：环带开缺口后屏障必须显著下降。

    这检验的正是本课题的核心主张——屏障是拓扑性质。一圈连续的阻抗细胞构成
    封锁线；一旦开一个口子，即便阻抗细胞总数几乎不变，通路就恢复了。
    """
    closed = make_ring_grid(n=25, gap_spots=0, seed=0)
    gapped = make_ring_grid(n=25, gap_spots=5, seed=0)

    c = compute_b_cell(closed.A, closed.ecm, closed.caf, closed.source, closed.sink, **CFG_CELL)
    g = compute_b_cell(gapped.A, gapped.ecm, gapped.caf, gapped.source, gapped.sink, **CFG_CELL)

    drop = 1.0 - g["b_cell"] / c["b_cell"]
    print(f"  闭合环带 b_cell = {c['b_cell']:.4f}  (阻抗节点 {len(closed.blocked_idx)})")
    print(f"  有缺口   b_cell = {g['b_cell']:.4f}  (阻抗节点 {len(gapped.blocked_idx)})")
    print(f"  屏障下降 = {drop * 100:.1f}%")

    assert g["b_cell"] < 0.5 * c["b_cell"], (
        f"开缺口后屏障只降了 {drop * 100:.1f}%，算子对拓扑缺口不敏感——"
        "它实际上退化成了局部统计量"
    )


# ---------------------------------------------------------------------------
def test_size_dependence():
    """测试 3：同一张图上，大分子的传质屏障必须高于小分子。"""
    sl = make_ring_grid(n=25, gap_spots=0, seed=0)

    big = compute_b_mab(sl.A, sl.ecm, sl.crosslink, sl.ag_target, sl.vessel,
                        r_nm=5.5, **CFG_MAB)      # IgG / ADC
    small = compute_b_mab(sl.A, sl.ecm, sl.crosslink, sl.ag_target, sl.vessel,
                          r_nm=0.5, **CFG_MAB)    # 小分子

    core = sl.sink  # 瘤巢核心，位于环带内侧
    mb, ms = np.nanmean(big["b_mab"][core]), np.nanmean(small["b_mab"][core])
    print(f"  瘤巢核心 B_mAb(IgG  r=5.5nm) = {mb:.3f}")
    print(f"  瘤巢核心 B_mAb(小分子 r=0.5nm) = {ms:.3f}")
    print(f"  差值 = {mb - ms:.3f}")

    assert mb > ms, "大分子在瘤巢核心的屏障未高于小分子，尺寸排阻项实现有误"

    # 单调性：半径越大，核心屏障越高（允许在完全阻断后进入平台）
    radii = [0.5, 2.0, 4.0, 5.5, 8.0]
    vals = [np.nanmean(compute_b_mab(sl.A, sl.ecm, sl.crosslink, sl.ag_target,
                                     sl.vessel, r_nm=r, **CFG_MAB)["b_mab"][core])
            for r in radii]
    print("  尺寸扫描: " + "  ".join(f"r={r}->{v:.2f}" for r, v in zip(radii, vals)))
    diffs = np.diff(vals)
    assert np.all(diffs >= -1e-6), f"B_mAb 随分子半径非单调：{vals}"


# ---------------------------------------------------------------------------
def test_b_meta_distance_monotonic():
    """补充测试：B_meta 应随距血管距离增大而增大（在状态项相同时）。"""
    sl = make_ring_grid(n=25, gap_spots=0, seed=0)
    n = sl.n
    flat = np.full(n, 0.5)
    res = compute_b_meta(sl.A, flat, flat, flat, sl.vessel,
                         spacing_um=sl.spacing_um, d0_um=400.0, tau_um=150.0)
    d, b = res["d_vessel_um"], res["b_meta"]
    order = np.argsort(d)
    # 距离排序后，b_meta 应基本单调不减（状态项恒定，只剩几何项）
    corr = np.corrcoef(d, b)[0, 1]
    print(f"  d_vessel 与 b_meta 的相关系数 = {corr:.4f}")
    assert corr > 0.95, f"B_meta 与距血管距离的相关性仅 {corr:.3f}，几何项实现有误"


# ---------------------------------------------------------------------------
def test_b_meta_weighted_distance():
    """d_vessel_um 必须用加权图距，不能再拿"跳数 × 间距"顶替。

    为什么单独测这一条（2026-08-27）
    --------------------------------
    hops × spacing_um 只在"所有边等长"时才对。真实阵列不是这样：
      · Visium 六方 r=150   —— 边基本都是 100 μm，近似还算准
      · 第一代 ST 交错阵列 r=300 —— 边混了 200 与 282.8 μm，却都按 200 算
    偏差因此**随平台反号**（实测 +11% 与 −11%），而 d0_um=130 锚的是直线
    氧扩散极限，于是两个队列的 B_meta 之间凭空多出约 22% 的系统偏移。
    这里用等边长与混边长两种图钉死行为：等边长时两种口径必须一致，
    混边长时加权口径必须给出更大的距离。
    """
    from sparta.graph import build_graph_radius

    sp_um = 100.0
    xy = np.array([[i * sp_um, j * sp_um] for i in range(12) for j in range(12)], float)
    flat = np.full(len(xy), 0.5)
    ves = [0]
    kw = dict(d0_um=400.0, tau_um=150.0)

    # (1) 只连正交邻居 -> 所有边都是 100 μm -> 两种口径必须一致
    A4, D4 = build_graph_radius(xy, radius_um=1.2 * sp_um)
    d_hop = compute_b_meta(A4, flat, flat, flat, ves, spacing_um=sp_um, **kw)
    d_wgt = compute_b_meta(A4, flat, flat, flat, ves, D=D4, spacing_um=sp_um, **kw)
    assert d_hop["d_vessel_mode"] == "hops" and d_wgt["d_vessel_mode"] == "weighted"
    dev = np.abs(d_wgt["d_vessel_um"] - d_hop["d_vessel_um"]).max()
    print(f"  等边长图：两种口径最大偏差 {dev:.3g} μm")
    assert dev < 1e-6, f"边全等长时加权口径不该有差异，却差了 {dev} μm"

    # (2) 加上对角邻居（141.4 μm）-> 跳数口径把对角边也当 100 μm，必然低估
    A8, D8 = build_graph_radius(xy, radius_um=1.5 * sp_um)
    h8 = compute_b_meta(A8, flat, flat, flat, ves, spacing_um=sp_um, **kw)["d_vessel_um"]
    w8 = compute_b_meta(A8, flat, flat, flat, ves, D=D8, spacing_um=sp_um, **kw)["d_vessel_um"]
    far = h8 > 0
    ratio = float(np.median(h8[far] / w8[far]))
    print(f"  混边长图：hops×spacing ÷ 加权图距 中位数 = {ratio:.3f}（应 < 1）")
    assert ratio < 0.98, f"混边长图上跳数口径应明显低估，实测比值 {ratio:.3f}"
    assert np.all(w8 >= h8 - 1e-9), "加权图距不应小于跳数口径"


# ---------------------------------------------------------------------------
def test_dissociation_is_detectable():
    """阳性对照（最重要的一个）：在**已知存在解离**的场景下，框架必须能检出解离。

    在真实数据上声称"发现了两道屏障解离"之前，必须先证明方法在已知解离的
    合成场景下检得出来。否则无法区分"真的没解离"与"我的方法测不出解离"。

    合成场景：ECM 环带开大缺口（T 细胞进得去，B_cell 低）+ 瘤巢外缘高抗原
    （抗体被结合位点屏障消耗，B_mAb 高）。即临床上的"热而无效"。
    """
    diss = make_dissociated_grid(n=25, seed=0)
    ring = make_ring_grid(n=25, gap_spots=0, width_spots=2.0, seed=0)

    bc_diss = compute_b_cell(diss.A, diss.ecm, diss.caf, diss.source, diss.sink, **CFG_CELL)
    bc_ring = compute_b_cell(ring.A, ring.ecm, ring.caf, ring.source, ring.sink, **CFG_CELL)

    core = diss.sink
    bm_diss = compute_b_mab(diss.A, diss.ecm, diss.crosslink, diss.ag_target,
                            diss.vessel, r_nm=5.5, **CFG_MAB)["b_mab"]
    # 同一张解离图，但把抗原拉平（去掉结合位点屏障）作为对照
    flat_ag = np.full(diss.n, 0.05)
    bm_flat = compute_b_mab(diss.A, diss.ecm, diss.crosslink, flat_ag,
                            diss.vessel, r_nm=5.5, **CFG_MAB)["b_mab"]

    print(f"  解离图 B_cell = {bc_diss['b_cell']:.4f}  (闭合环带对照 {bc_ring['b_cell']:.4f})")
    print(f"  解离图 瘤巢核心 B_mAb = {np.nanmean(bm_diss[core]):.3f}")
    print(f"  抗原拉平后 核心 B_mAb = {np.nanmean(bm_flat[core]):.3f}")

    # 细胞屏障必须明显低于闭合环带（T 细胞确实进得去）
    assert bc_diss["b_cell"] < 0.5 * bc_ring["b_cell"], (
        "解离场景下细胞屏障未明显降低，缺口没起作用"
    )
    # 抗体屏障必须因结合位点屏障而升高（抗体确实进不去）
    assert np.nanmean(bm_diss[core]) > np.nanmean(bm_flat[core]), (
        "结合位点屏障未使核心区抗体屏障升高，kappa 项实现有误"
    )

    # ------------------------------------------------------------------
    # 下面是**诊断输出而非断言**。原因见 README「已知的核心风险」一节：
    # B_cell 与 B_mAb 之间除了共享几何（距血管距离）之外，还共享同一个 ECM 项
    # ——ECM 既进了最小割的边容量，也进了扩散的边电导。因此一个 ECM 缺口
    # 会同时降低两个屏障。要在全片尺度上观察到解离，抗原/交联这类
    # **只影响一个屏障**的因素必须强到能压过这条共享通道。
    # 这是一个关于真实组织的经验问题，不应该靠调合成数据的参数来"证明"。
    # ------------------------------------------------------------------
    from sparta.validate import dissociation_drivers

    dd = dissociation_drivers(
        diss.A, diss.ecm, diss.caf, diss.crosslink, diss.ag_target,
        diss.source, diss.vessel, cfg_cell=CFG_CELL, cfg_mab=CFG_MAB,
    )
    print("  [驱动分解] B_mAb 方差中：")
    print(f"    共享 ECM 通道解释 {dd['frac_shared_ecm'] * 100:5.1f}%")
    print(f"    抗原(BSB)特异贡献 {dd['frac_antigen'] * 100:5.1f}%")
    print(f"    交联/尺寸特异贡献 {dd['frac_crosslink'] * 100:5.1f}%")
    print(f"    -> 解离潜力比 = {dd['dissociation_potential']:.3f}"
          f"（特异因素 / 共享因素；>1 才有望在全片尺度观察到解离）")


# ---------------------------------------------------------------------------
def test_empty_source_returns_nan():
    """健壮性：源集为空时应返回 nan 而不是崩溃或静默返回 0。"""
    sl = make_ring_grid(n=15, seed=0)
    r = compute_b_cell(sl.A, sl.ecm, sl.caf, [], sl.sink, **CFG_CELL)
    assert np.isnan(r["b_cell"]), "源集为空时应返回 nan"
    m = compute_b_mab(sl.A, sl.ecm, sl.crosslink, sl.ag_target, [], **CFG_MAB)
    assert np.all(np.isnan(m["b_mab"])), "血管集为空时 b_mab 应全为 nan"


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        print(f"\n[RUN ] {t.__name__}")
        print(f"       {(t.__doc__ or '').strip().splitlines()[0]}")
        try:
            t()
            print(f"[PASS] {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"[FAIL] {t.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"[ERROR] {t.__name__}: {type(e).__name__}: {e}")
    print("\n" + "=" * 60)
    print(f"{len(tests) - failed} / {len(tests)} passed")
    sys.exit(1 if failed else 0)
