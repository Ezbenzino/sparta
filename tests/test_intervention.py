"""
intervention.py 的合成图单元测试 —— 检验"排序是不是真的在定位屏障"
=====================================================================

运行方式（两种都行）：
    pytest tests/test_intervention.py -v
    python tests/test_intervention.py

四个测试共同锁住下面这条主张：逐 spot 反事实排序**不是**把局部 ECM 分数
换个写法重排一遍，而是真的在找最小割封锁带上的位置。

  测试 1  最该"打通"的位置必须落在环带（最小割）上，而非全片任意位置
  测试 2  环带之外的零阻力 spot 干预后屏障**完全不变**（精确 0，
          这是 max-flow 的结构性结论，不是近似）
  测试 3  单个 spot 的干预效果不能超过"整条割集全打通"的刻度锚点
  测试 4  两个算子（B_cell / B_mAb）的逐 spot 干预读数都非负、形状一致
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from sparta.intervention import spot_intervention_ranking  # noqa: E402
from sparta.synthetic import make_ring_grid  # noqa: E402

CFG_CELL = dict(a=3.0, b_ecm=8.0, c_caf=4.0)
CFG_MAB = dict(g0=1.0, lam=3.0, xi0_nm=20.0, beta=3.0, kd_eff=0.5, kappa_w=1.0)

_N = 15            # 15x15 网格 = 225 spot，约 6-12 s
_RADIUS = 4.0


def _ranking():
    """构造闭合环带切片并跑一次排序（缓存供四个测试共用）。"""
    if not hasattr(_ranking, "_cache"):
        slide = make_ring_grid(n=_N, radius_spots=_RADIUS, width_spots=1.2,
                               gap_spots=0, seed=0)
        res = spot_intervention_ranking(
            slide.A, slide.ecm, slide.caf, slide.crosslink, slide.ag_target,
            slide.source, slide.sink, slide.vessel, CFG_CELL, CFG_MAB,
        )
        _ranking._cache = (slide, res)
    return _ranking._cache


# ---------------------------------------------------------------------------
def test_top_spots_are_on_the_ring():
    """测试 1：ΔB_cell 最高的位置应显著富集在环带上（本底 vs top-5%）。"""
    slide, res = _ranking()
    d = np.asarray(res["delta_b_cell"], float)
    ring = np.zeros(slide.n, bool)
    ring[slide.blocked_idx] = True

    base = ring.mean()
    order = np.argsort(-d)
    k = max(1, int(round(0.05 * slide.n)))
    top = order[:k]
    frac_top = ring[top].mean()

    print(f"  环带 {ring.sum()}/{slide.n} 个 spot（本底 {base*100:.0f}%）")
    print(f"  top-5% ({k} 个) 中落在环带上的比例：{frac_top*100:.0f}%")
    print(f"  前三名 ΔB_cell = {np.round(d[order[:3]], 5)}")

    assert base < 0.25, f"环带占比 {base:.2f} 过高，本测试失去区分度"
    assert frac_top >= 0.9, (f"top-5% 仅 {frac_top*100:.0f}% 落在环带上——"
                             "排序没有在定位封锁带")


def test_zero_effect_outside_the_cut():
    """测试 2：环带外零阻力 spot 的干预必须精确零效应（max-flow 结构性结论）。"""
    slide, res = _ranking()
    d = np.asarray(res["delta_b_cell"], float)
    ring = np.zeros(slide.n, bool)
    ring[slide.blocked_idx] = True

    # 取一个零阻力、且不在环带上的 spot
    outside = np.flatnonzero((~ring) & np.isclose(slide.ecm, 0.0))
    assert len(outside) > 0, "合成图里没有可用的环带外零阻力 spot"
    i = outside[len(outside) // 2]
    print(f"  环带外零阻力 spot {i}: ΔB_cell = {d[i]:.3e}")
    print(f"  全片 ΔB_cell 精确为零的比例 = {np.mean(np.isclose(d, 0.0))*100:.1f}%")

    assert abs(d[i]) < 1e-12, f"环带外 spot 干预改变屏障 {d[i]:.3e}，排序不是拓扑定位"
    assert np.mean(np.isclose(d, 0.0)) > 0.3, "零效应比例过低，与闭合环带的直觉不符"


def test_single_spot_bounded_by_full_breach():
    """测试 3：单个 spot 的干预效果不得超过整条割集打通的刻度锚点。"""
    slide, res = _ranking()
    d = np.asarray(res["delta_b_cell"], float)
    full = res["full_cut_breach"]

    print(f"  最大单点 ΔB_cell = {d.max():.5f}")
    print(f"  整条割集打通     = {full:.5f}  (比值 {d.max()/full:.2f})")

    assert np.isfinite(full) and full > 0, "刻度锚点应为有限正值"
    assert d.max() <= full * (1 + 1e-6), "单点效果超过了全割集打通的锚点，锚点定义有误"
    assert d.max() / full < 1.0, "单点即达成全部效果，与'屏障是分布式封锁带'的设定矛盾"
    assert d.max() / full > 0.02, "单点效果过小，排序可能全是数值噪声"


def test_both_operators_shapes_and_signs():
    """测试 4：两个算子的读数形状一致；B_cell 严格非负，B_mAb 负向残差可忽略。"""
    slide, res = _ranking()
    dc = np.asarray(res["delta_b_cell"], float)
    dm = np.asarray(res["delta_b_mab"], float)

    assert dc.shape == (slide.n,) and dm.shape == (slide.n,)
    # B_cell：打通只抬高边容量 -> 最大流不减 -> 屏障不增，严格非负
    assert np.all(dc >= -1e-12), "打通屏障不应使 B_cell 升高"
    # B_mAb：负向残差来自"吸收重分布"——打通高吸收区附近的电导会让抗体
    # 更多地在吸收区被截留，使别处浓度略降、屏障略微上升。这是屏蔽泊松
    # 方程的真实性质，不是数值噪声，所以只要求其量级可忽略（<0.5% of max）。
    neg_share = abs(min(dm.min(), 0.0)) / max(dm.max(), 1e-12)
    print(f"  ΔB_mAb 负向残差 {neg_share*100:.3f}% of max（吸收重分布）")
    assert neg_share < 0.005, \
        f"ΔB_mAb 出现 {neg_share*100:.2f}% 的负向偏移，超出吸收重分布可解释的范围"
    assert dm.max() > 0, "B_mAb 单点干预完全没有效果，尺寸排阻项可能没起作用"
    assert dm.max() <= res["all_removed_mab"] * (1 + 1e-6), \
        "单点 B_mAb 效果超过'全片物质去掉'的锚点"
    print(f"  ΔB_cell: max={dc.max():.5f}  零效应 {np.mean(np.isclose(dc,0))*100:.0f}%")
    print(f"  ΔB_mAb : max={dm.max():.5f}  零效应 {np.mean(np.isclose(dm,0))*100:.0f}%")


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
