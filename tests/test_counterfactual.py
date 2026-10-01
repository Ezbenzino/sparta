"""
counterfactual.py 的合成图测试
================================

运行方式：
    pytest tests/test_counterfactual.py -v
    python tests/test_counterfactual.py

这里检验的是三个反事实实验在"已知正确答案"的合成图上是否给出正确结论。
合成图上环带是人为放进去的，所以：
  S1 必须显著（屏障确实来自空间排布）
  S2 必须显示定向移除远强于随机移除
  S3 必须显示大分子在核心区的屏障远高于小分子
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from sparta.counterfactual import s1_spatial_permutation, s2_ring_breaking, s3_size_scan  # noqa: E402
from sparta.synthetic import make_ring_grid                                               # noqa: E402

CFG_CELL = dict(a=3.0, b_ecm=8.0, c_caf=4.0)
CFG_MAB = dict(g0=1.0, lam=3.0, xi0_nm=20.0, beta=3.0, kd_eff=0.5, kappa_w=1.0)


def test_s1_permutation_detects_structure():
    """S1：合成图上有真实环带，空间重排必须使屏障显著下降。"""
    sl = make_ring_grid(n=21, gap_spots=0, seed=0)
    res = s1_spatial_permutation(sl.A, sl.ecm, sl.caf, sl.source, sl.sink,
                                 CFG_CELL, n_perm=120, seed=1, mode="fixed")
    print(f"  真实 b_cell = {res['b_real']:.4f}")
    print(f"  零分布均值  = {res['null'].mean():.4f} ± {res['null'].std():.4f}")
    print(f"  z = {res['z']:.2f}   经验 p = {res['p_emp']:.4f}")

    assert res["z"] > 3.0, f"z 值仅 {res['z']:.2f}，S1 未能检出人为放入的环带结构"
    assert res["p_emp"] < 0.05, f"经验 p = {res['p_emp']:.4f}，未达显著"


def _run_s2(width_spots, k_list=(5, 10)):
    sl = make_ring_grid(n=21, gap_spots=0, width_spots=width_spots, seed=0)
    res = s2_ring_breaking(sl.A, sl.ecm, sl.caf, sl.source, sl.sink, CFG_CELL,
                           coords=sl.coords, k_list=k_list, n_rand=60, seed=1)
    return sl, res


def test_s2_contiguous_gap_beats_scattered():
    """S2：连续缺口必须比等量的分散移除更有效（关键对照：材料相同，只差是否连续）。"""
    for width in (1.2, 2.5):
        sl, res = _run_s2(width)
        print(f"\n  --- 屏障带厚度 = {width} spot，最小割节点 {res['n_cut_nodes']} 个 ---")
        sig = []
        for k, r in res["per_k"].items():
            print(f"  k={k:>2}: 连续缺口 {r['targeted']:.4f} | "
                  f"全局随机 {r['rand_global_mean']:.4f} | "
                  f"屏障内分散 {r['rand_in_cut_mean']:.4f}")
            print(f"        剩余屏障比 vs 全局 {r['ratio_vs_global']:.2f}x | "
                  f"vs 屏障内分散 {r['ratio_vs_in_cut']:.2f}x | p={r['p_vs_in_cut']:.3f}")
            sig.append(r["p_vs_in_cut"] < 0.05 and r["ratio_vs_in_cut"] > 1.0)
        assert res["per_k"], "没有产生任何 k 的结果"
        assert any(sig), (
            f"厚度 {width}：连续缺口未显著优于等量分散移除，拓扑主张在此设定下不成立"
        )


def test_s2_thickness_amplifies_contiguity():
    """S2 补充：屏障带越厚，连续性的优势应该越大（这是拓扑主张的一个可证伪推论）。"""
    _, thin = _run_s2(1.2, k_list=(5,))
    _, thick = _run_s2(3.0, k_list=(5,))
    rt = thin["per_k"][5]["ratio_vs_in_cut"]
    rk = thick["per_k"][5]["ratio_vs_in_cut"]
    print(f"  薄带(1.2 spot) 剩余屏障比 = {rt:.2f}x")
    print(f"  厚带(3.0 spot) 剩余屏障比 = {rk:.2f}x")
    assert rk >= rt * 0.9, (
        f"厚带的连续性优势({rk:.2f}x)未达到薄带({rt:.2f}x)的水平，"
        "与拓扑主张的推论不符"
    )


def test_s3_size_scan_monotone():
    """S3：核心区屏障必须随分子半径单调不减，且 IgG 远高于小分子。"""
    sl = make_ring_grid(n=21, gap_spots=0, seed=0)
    res = s3_size_scan(sl.A, sl.ecm, sl.crosslink, sl.ag_target, sl.vessel,
                       CFG_MAB, core_idx=sl.sink)
    for r, m in zip(res["radii_nm"], res["mean_core"]):
        print(f"  r={r:>4.1f} nm  ->  核心区 B_mAb = {m:.3f}")
    print(f"  IgG 相对小分子的核心区屏障增幅 = {res['igg_vs_small_core_delta']:.3f}")

    mc = np.asarray(res["mean_core"])
    assert np.all(np.diff(mc) >= -1e-6), f"核心区屏障随半径非单调：{mc}"
    assert res["igg_vs_small_core_delta"] > 1.0, "IgG 与小分子的核心区屏障差异过小"


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
