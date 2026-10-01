"""
统计口径的回归测试（2026-08-26 审查后新增）
==============================================

运行方式：
    pytest tests/test_statistics.py -v
    python tests/test_statistics.py

为什么需要这一组测试
--------------------
原有的 14 个测试全部检验**算子的定性行为**（环带 vs 散在、大分子 vs 小分子），
它们跑在合成图上，跑得很好，也确实抓出过三个真 bug。
但它们里没有任何一个会因为下面这些错误而变红：

  · `decoupling_stats` 忘了传 control，产出 rho_partial=NaN，
    而下游却把未校正的 rho 当成"偏相关"写进论文；
  · "解离区占比 6%"被当成发现报告，而两个屏障相互独立时的期望本来就是 6.25%；
  · S2 的 k 是绝对值，割集从 76 个节点长到 614 个节点后，
    同一个 k 只动了割集的 0.5%，效应量被稀释到看不见。

这三件事在 2026-08-26 的投稿前审查里全部实际发生了
（见 docs/review_for_journal.md §3）。本文件把它们变成会报警的测试。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from sparta.counterfactual import s2_ring_breaking          # noqa: E402
from sparta.synthetic import make_ring_grid                 # noqa: E402
from sparta.validate import decoupling_stats                # noqa: E402

CFG_CELL = dict(a=3.0, b_ecm=8.0, c_caf=4.0)
Q = 0.75


def test_discordant_fraction_has_a_chance_baseline():
    """两个**独立**的屏障场，解离区占比必须≈(1-q)^2，富集倍数≈1。

    失败意味着：可以把随机水平的解离区占比当成"发现"报出去。
    """
    rng = np.random.default_rng(0)
    x, y = rng.normal(size=4000), rng.normal(size=4000)
    d = decoupling_stats(x, y, q=Q)

    print(f"  随机期望 chance_discordant = {d['chance_discordant']:.4f}")
    print(f"  实测 frac_discordant       = {d['frac_discordant']:.4f}")
    print(f"  富集倍数 enrichment        = {d['enrichment_vs_chance']:.3f}x")

    assert abs(d["chance_discordant"] - (1 - Q) ** 2) < 1e-12, "随机期望的定义写错了"
    assert 0.8 < d["enrichment_vs_chance"] < 1.25, (
        f"独立场的富集倍数应≈1，实测 {d['enrichment_vs_chance']:.3f}x —— "
        "解离区占比的基线定义有问题"
    )


def test_enrichment_detects_true_dissociation():
    """构造一对真正反相关的场，富集倍数必须明显 > 1。

    这一条是上一条的配对：证明这个指标不是对什么都给 1x，
    否则"没测到解离"就分不清是真没有还是指标是瞎的。
    """
    rng = np.random.default_rng(1)
    x = rng.normal(size=4000)
    y = -x + 0.35 * rng.normal(size=4000)          # 强反相关
    d = decoupling_stats(x, y, q=Q)

    print(f"  反相关场 frac_discordant = {d['frac_discordant']:.4f}"
          f"（随机期望 {d['chance_discordant']:.4f}）")
    print(f"  富集倍数 = {d['enrichment_vs_chance']:.2f}x， rho = {d['rho']:+.3f}")

    assert d["enrichment_vs_chance"] > 2.0, (
        f"已知强反相关的场只给出 {d['enrichment_vs_chance']:.2f}x 富集，指标不灵敏"
    )


def test_decoupling_marks_whether_control_was_applied():
    """没传 control 时必须 controlled=False 且 rho_partial=NaN；传了才有偏相关。

    失败意味着：可以把未校正的相关伪装成"控制距血管距离后的偏相关"。
    这正是 2026-08-26 之前 run_06_validate.py 的实际情况。
    """
    rng = np.random.default_rng(2)
    depth = rng.uniform(size=3000)                  # 共同混杂：距血管深度
    x = depth + 0.3 * rng.normal(size=3000)
    y = depth + 0.3 * rng.normal(size=3000)

    d0 = decoupling_stats(x, y, q=Q)                       # 不传 control
    d1 = decoupling_stats(x, y, q=Q, control=depth)        # 传 control

    print(f"  不传 control: controlled={d0['controlled']}, "
          f"rho={d0['rho']:+.3f}, rho_partial={d0['rho_partial']}")
    print(f"  传 control  : controlled={d1['controlled']}, "
          f"rho={d1['rho']:+.3f}, rho_partial={d1['rho_partial']:+.3f}")

    assert d0["controlled"] is False and np.isnan(d0["rho_partial"]), (
        "没传 control 却给出了偏相关，下游会误以为已经校正过几何"
    )
    assert d1["controlled"] is True and np.isfinite(d1["rho_partial"]), "传了 control 却没算偏相关"
    assert abs(d1["rho_partial"]) < abs(d1["rho"]), (
        "校正共同混杂后相关性没有下降，_residualize 可能没起作用"
    )


def test_s2_k_scales_with_cut_set():
    """S2 的移除规模必须相对割集定义，且必须把占比记进结果。

    背景：`k_list=[3,5,8,12]` 是给合成图（割集 40–76 节点）设计的。
    真实切片的割集有 150–750 个节点，同一个 k 只动了割集的 0.5%–8%，
    效应量被稀释到测不出来——这不是拓扑假设不成立，是实验规模没缩放。

    本测试用两张割集规模相差一倍以上的合成图检查：
      · 用 k_frac 时，两张图上实际移除的割集比例必须一致；
      · 用绝对 k_list 时，比例必然不一致——所以结果里**必须**记录 k_over_cut，
        否则这种稀释在下游完全不可见。
    """
    frac = 0.30
    small = make_ring_grid(n=15, radius_spots=4.0, width_spots=2.0, seed=0)
    big = make_ring_grid(n=31, radius_spots=10.0, width_spots=2.0, seed=0)

    got = {}
    for name, sl in (("small", small), ("big", big)):
        r_frac = s2_ring_breaking(sl.A, sl.ecm, sl.caf, sl.source, sl.sink, CFG_CELL,
                                  coords=sl.coords, k_frac=(frac,), n_rand=20, seed=1)
        r_abs = s2_ring_breaking(sl.A, sl.ecm, sl.caf, sl.source, sl.sink, CFG_CELL,
                                 coords=sl.coords, k_list=(5,), n_rand=20, seed=1)
        kf = next(iter(r_frac["per_k"].values()))
        ka = next(iter(r_abs["per_k"].values()))
        got[name] = dict(pool=r_frac["n_cut_pool"],
                         frac_over_cut=kf["k_over_cut"], frac_ratio=kf["ratio_vs_in_cut"],
                         abs_over_cut=ka["k_over_cut"], abs_ratio=ka["ratio_vs_in_cut"])
        print(f"  {name:<6} 割集可动 {got[name]['pool']:>4} 个｜"
              f"k_frac={frac} -> 实际移除 {kf['k_over_cut']*100:>5.1f}%（{kf['ratio_vs_in_cut']:.2f}x）｜"
              f"绝对 k=5 -> 实际移除 {ka['k_over_cut']*100:>5.1f}%（{ka['ratio_vs_in_cut']:.2f}x）")

    assert got["big"]["pool"] > 1.5 * got["small"]["pool"], "两张图的割集规模差得不够，测不出问题"

    for name in ("small", "big"):
        assert abs(got[name][ "frac_over_cut"] - frac) < 0.06, (
            f"{name}: k_frac={frac} 但实际移除了 {got[name]['frac_over_cut']:.3f} 的割集"
        )

    abs_lo, abs_hi = sorted([got["small"]["abs_over_cut"], got["big"]["abs_over_cut"]])
    print(f"  -> 绝对 k 在两张图上移除的割集比例相差 {abs_hi / abs_lo:.1f} 倍"
          f"（{abs_lo*100:.1f}% vs {abs_hi*100:.1f}%）")
    assert abs_hi / abs_lo > 1.5, (
        "绝对 k 在两种割集规模上移除的比例竟然一致，说明这张合成图没有复现真实的规模差异"
    )
    for name in ("small", "big"):
        assert got[name]["abs_over_cut"] > 0 and got[name]["frac_over_cut"] > 0, (
            "结果里没有记录 k_over_cut —— 移除规模的稀释在下游将完全不可见"
        )


TESTS = [
    test_discordant_fraction_has_a_chance_baseline,
    test_enrichment_detects_true_dissociation,
    test_decoupling_marks_whether_control_was_applied,
    test_s2_k_scales_with_cut_set,
]

if __name__ == "__main__":
    ok = 0
    for fn in TESTS:
        print(f"\n[RUN ] {fn.__name__}")
        print(f"       {fn.__doc__.strip().splitlines()[0]}")
        try:
            fn()
            print(f"[PASS] {fn.__name__}")
            ok += 1
        except AssertionError as e:
            print(f"[FAIL] {fn.__name__}: {e}")
    print("\n" + "=" * 60)
    print(f"{ok} / {len(TESTS)} passed")
    sys.exit(0 if ok == len(TESTS) else 1)
