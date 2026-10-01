#!/usr/bin/env python
"""
run_14_shared_ecm_check.py —— 决策岔口：耦合是组织的性质，还是模型自带的？
==============================================================================

输入：{sid}.scored.h5ad + {sid}.graph.npz（不需要 run_04 的产物）
输出：results/validation/shared_ecm_check.json（**新文件，不覆盖任何已有产物**）
上游模块：run_03_graph.py
下游模块：无 —— 这是一个**决策节点**，它的结果决定这篇文章投哪个刊

它回答的问题
------------
PIVOT 结论是"两道屏障高度耦合"。审稿人一定会问：

    "它们正相关，会不会只是因为你把同一个 core matrisome 分数
     同时放进了最小割的边容量和扩散的边电导？"

这个质疑是合理的——按当前公式，ECM 确实是两个算子的共同输入，
所以"两屏障正相关"**有一部分是结构上必然的**。

本脚本把共享通道切断后再问一次：
  · 主配置    B_cell = f(ECM, CAF)          B_mAb = g(ECM, 交联, 抗原)
  · 变体配置  B_cell = f(     CAF)          B_mAb = g(     交联, 抗原)
变体里两个算子**没有任何共同输入**。距血管距离（control）两种口径下都用
主配置算，保证校正的是同一个几何量。

判读
----
  · 变体下偏相关**仍显著为正** -> 耦合是**组织的性质**：在真实皮肤肿瘤里，
    挡细胞的结构与挡大分子的结构本来就长在一起。PIVOT 升级为强结论，
    值得冲 Briefings in Bioinformatics。
  · 变体下偏相关**接近 0 或转负** -> 耦合主要由共享输入造成。
    这不是失败，但措辞必须收窄为"在本框架的边权定义下二者不可分离"，
    并如实报告本检查的结果。投 Bioinformatics。

**两种结果都要写进论文。** 主动报告远好过被审稿人问出来。

用法
----
    python scripts/run_14_shared_ecm_check.py --slides MEL01 MEL02 MEL03 MEL04 \
                                              CSCC01 CSCC02 CSCC03 CSCC04
    python scripts/run_14_shared_ecm_check.py --synthetic     # 先看输出长什么样
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import (Paths, load_config, load_graph, patient_map,  # noqa: E402
                        save_json, set_seed, stamp_run)

VARIANT_DEFAULT = "configs/no_shared_ecm.yaml"


def _cohort(sid: str) -> str:
    """从切片 ID 取队列前缀（MEL01 -> MEL）。队列间的差异必须单独看，
    中位数会把'一个队列全中、另一个队列全不中'这种结构平掉。"""
    return "".join(c for c in sid if not c.isdigit()) or sid


def _summarise(rows: dict, hdr_width: int = 76, pmap: dict | None = None) -> dict:
    """按显著性计数 + 分队列给判定。

    为什么不能只看符号和中位数
    --------------------------
    2026-08-27 第一次真实运行就撞上了：全体中位 +0.153、7/8 为正，
    按"中位>0.05 且 75% 为正"会判 TISSUE；但拆开看是
    CSCC 4/4 显著（保留 72–93%）而 MEL 只有 1/4 显著（另外三张塌到 ±0.03）。
    把这两个队列平均成一个数字，等于把结论里最有信息量的那部分丢掉。
    """
    main_rp = np.array([r["主配置"]["rho_partial"] for r in rows.values()], float)
    var_rp = np.array([r["切断共享ECM"]["rho_partial"] for r in rows.values()], float)
    var_p = np.array([r["切断共享ECM"]["p_partial"] for r in rows.values()], float)
    sig_pos = (var_rp > 0) & (var_p < 0.05)
    med_main, med_var = float(np.median(main_rp)), float(np.median(var_rp))

    print("-" * hdr_width)
    print(f"主配置    偏相关：中位 {med_main:+.3f}，"
          f"{int((main_rp > 0).sum())}/{len(main_rp)} 为正")
    print(f"切断共享  偏相关：中位 {med_var:+.3f}，"
          f"{int((var_rp > 0).sum())}/{len(var_rp)} 为正，"
          f"**{int(sig_pos.sum())}/{len(var_rp)} 显著为正（p<0.05）**")

    sids = list(rows)
    if pmap:
        pats = {}
        for i, sid in enumerate(sids):
            pats.setdefault(pmap.get(sid, sid), []).append(i)
        print(f"\n  样本量：{len(sids)} 张切片 / **{len(pats)} 位患者**")
        for k, idx in pats.items():
            n_ok = int(sig_pos[idx].sum())
            print(f"    {k:<12} {len(idx)} 张，切断后显著为正 {n_ok}/{len(idx)}"
                  f"　中位 {np.median(var_rp[idx]):+.3f}")
        if len(pats) < len(sids):
            print("    ⚠ 同患者多切片不是独立样本；正文的 n 必须按患者报。")

    by_cohort = {}
    for c in sorted({_cohort(s) for s in sids}):
        idx = [i for i, s in enumerate(sids) if _cohort(s) == c]
        by_cohort[c] = dict(
            n=len(idx),
            median_rho_partial_variant=float(np.median(var_rp[idx])),
            n_significant_positive=int(sig_pos[idx].sum()),
            retention=[float(var_rp[i] / main_rp[i]) if main_rp[i] else None for i in idx],
            slides=[sids[i] for i in idx])
    print("\n  分队列（中位数会把队列间的差异平掉，所以必须拆开看）：")
    for c, v in by_cohort.items():
        print(f"    {c:<6} n={v['n']}  切断后中位 {v['median_rho_partial_variant']:+.3f}"
              f"  显著为正 {v['n_significant_positive']}/{v['n']}")

    all_in = [v["n_significant_positive"] == v["n"] for v in by_cohort.values()]
    frac_sig = sig_pos.mean()

    print("\n" + "=" * hdr_width)
    print("旧的单一判定（**已降级为参考，不要写进正文**，理由见下）：")
    if frac_sig >= 0.75 and all(all_in):
        verdict = "TISSUE"
        print("判定：TISSUE —— 耦合是组织的性质，不是模型自带的")
        print("  切断两个算子的全部共同输入后，每个队列都仍然显著正相关。")
        print("  -> 值得冲 Briefings in Bioinformatics；本检查作为主图之一。")
    elif any(all_in) and frac_sig >= 0.4:
        verdict = "SPLIT"
        keep = [c for c, v in by_cohort.items() if v["n_significant_positive"] == v["n"]]
        drop = [c for c, v in by_cohort.items() if v["n_significant_positive"] < v["n"]]
        print("判定：SPLIT —— 队列之间不一致，这本身就是结果")
        print(f"  {'/'.join(keep)}：切断共同输入后仍然全部显著正相关"
              f" -> 在该队列里耦合是组织的性质。")
        print(f"  {'/'.join(drop)}：切断后大部分塌掉"
              f" -> 在该队列里耦合主要由共享的基质项带来。")
        print("  正文必须**分队列陈述**，不要合并成一句'皮肤肿瘤中两屏障耦合'。")
        print("  这个对比是可解释、可证伪的（组织结构组织度不同），比一个均匀结论更有信息量，")
        print("  但队列 = 瘤种 + 原发/转移 + 平台世代，三者完全共变，无法在本设计下分离——必须写明。")
    elif med_var < -0.05:
        verdict = "REVERSED"
        print("判定：REVERSED —— 切断共享输入后转为负相关")
        print("  共享项掩盖了一个真实的反向关系。在改写任何结论之前先停下来讨论。")
    else:
        verdict = "MODEL"
        print("判定：MODEL —— 耦合主要由共享输入造成")
        print("  措辞必须收窄为：'在本框架的边权定义下，两道屏障不可分离；")
        print("  其耦合的主要来源是两个算子共享的细胞外基质项，而非几何深度。'")
        print("  -> 投 Bioinformatics，不要冲 BIB。本检查如实写进 Results。")
    print("=" * hdr_width)

    by_patient = {}
    for sid in sids:
        by_patient.setdefault((pmap or {}).get(sid, sid), []).append(sid)

    # ------------------------------------------------------------------
    # 为什么不再给单一判定（2026-08-27 决定）
    # ------------------------------------------------------------------
    # 上面那条 TISSUE/SPLIT/MODEL 规则要求「某个队列的**每一张**切片都显著」。
    # 这个门槛的严苛程度随 n 单调上升：切片越多越不可能满贯，
    # **证据越多反而越走不到 TISSUE**。这是规则构造上的缺陷，跟数字是多少无关。
    # 实际后果：8 张时它给 SPLIT（靠「cSCC 4/4 满贯」），19 张时给 MODEL
    # （cSCC 变成 11–12/15），两次都没有描述数据——中间没有任何证据反转，
    # 变的只有 n。
    #
    # 换一条阈值当然能让标签变好看，但那正是事后调参。所以这里不换规则，
    # 而是取消「用一个标签概括」这件事本身：按队列 × 按患者陈述，
    # 效应量、保留率、显著性计数全部摆出来。这样没有可调的旋钮。
    #
    # 另一个理由是样本结构：黑色素瘤只有 1 位患者，任何整体标签都会被
    # 这一位患者的 4 张切片牵着走。分层陈述才诚实。
    # ------------------------------------------------------------------
    strat = {}
    for c, v in by_cohort.items():
        idx = [i for i, x in enumerate(sids) if _cohort(x) == c]
        keep = [var_rp[i] / main_rp[i] for i in idx if main_rp[i]]
        pats_c = {}
        for i in idx:
            pats_c.setdefault((pmap or {}).get(sids[i], sids[i]), []).append(i)
        strat[c] = dict(
            n_slides=len(idx), n_patients=len(pats_c),
            median_rho_partial_main=float(np.median(main_rp[idx])),
            median_rho_partial_variant=float(np.median(var_rp[idx])),
            n_significant_positive=int(sig_pos[idx].sum()),
            median_retention=float(np.median(keep)) if keep else None,
            per_patient={k: dict(
                n_slides=len(ii),
                n_significant_positive=int(sig_pos[ii].sum()),
                median_rho_partial_variant=float(np.median(var_rp[ii])),
                majority_significant=bool(int(sig_pos[ii].sum()) * 2 > len(ii)),
                slides=[sids[i] for i in ii]) for k, ii in pats_c.items()})

    print("\n分层陈述（**正文用这一份**）：")
    for c, v in strat.items():
        maj = sum(1 for q in v["per_patient"].values() if q["majority_significant"])
        ret = "—" if v["median_retention"] is None else f"{v['median_retention']:.2f}"
        print(f"  {c}：{v['n_slides']} 张 / {v['n_patients']} 位患者　"
              f"切断后中位 {v['median_rho_partial_variant']:+.3f}　"
              f"显著为正 {v['n_significant_positive']}/{v['n_slides']} 张、"
              f"{maj}/{v['n_patients']} 位患者过半　保留率中位 {ret}")
        for k, q in sorted(v["per_patient"].items()):
            print(f"      {k:<12} {q['n_significant_positive']}/{q['n_slides']}"
                  f"　中位 {q['median_rho_partial_variant']:+.3f}"
                  f"　{'过半' if q['majority_significant'] else '未过半'}")
    print("  单一判定已降级：它随 n 变严，8 张给 SPLIT、19 张给 MODEL，两次都没描述数据。")
    print("=" * hdr_width)

    return dict(verdict=verdict, verdict_status="deprecated_2026-08-27",
                verdict_note="判定规则要求队列内每张切片都显著，严苛程度随 n 单调上升，"
                             "证据越多越走不到 TISSUE。已改为按队列×按患者分层陈述，"
                             "正文不使用这个字段。",
                stratified=strat,
                n_slides=len(rows),
                n_patients=len(by_patient), patients=by_patient,
                median_rho_partial_main=med_main,
                median_rho_partial_variant=med_var,
                n_positive_variant=int((var_rp > 0).sum()),
                n_significant_positive_variant=int(sig_pos.sum()),
                by_cohort=by_cohort)


def _reverdict(cfg, tag=None):
    """只读回已有结果重新判定，不重算。改了判定规则之后刷新用。

    tag 对应 --tag：两种距离口径各存一份，刷新时也要分别指定，
    否则 hops 那一份永远停在旧的汇总结构上。
    """
    from sparta.io_ import load_json
    P = Paths(cfg)
    p = P.validation(f"shared_ecm_check_{tag}.json" if tag else "shared_ecm_check.json")
    if not p.exists():
        sys.exit(f"找不到 {p}，先跑一次完整的 run_14。")
    d = load_json(p)
    rows = d["per_slide"]
    print(f"重新判定 {p}（{len(rows)} 张切片，不重算）\n")
    d.update(_summarise(rows, pmap=patient_map(P)))
    save_json(p, d)
    print(f"\n已回写 {p}")


def main():
    ap = argparse.ArgumentParser(
        description="决策岔口：切断共享 ECM 通道后两屏障还相关吗",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slides", nargs="+", default=None)
    ap.add_argument("--synthetic", action="store_true",
                    help="用合成切片演示（不需要真实数据），先看看输出长什么样")
    ap.add_argument("--config", default=None, help="主配置")
    ap.add_argument("--variant", default=VARIANT_DEFAULT, help="切断共享通道的变体配置")
    ap.add_argument("--d-vessel", default="weighted", choices=["weighted", "hops"],
                    help="偏相关里当 control 的距血管距离用哪个口径。"
                         "weighted=沿图的加权最短路（正式口径，2026-08-27 起）；"
                         "hops=旧的 跳数×spacing_um。加这个开关是为了能证明本检查的"
                         "结论不是那次几何修正带来的假象——两种口径都要报。")
    ap.add_argument("--tag", default=None,
                    help="输出文件后缀，例如 --tag hops 会写 shared_ecm_check_hops.json，"
                         "免得两种口径互相覆盖")
    ap.add_argument("--reverdict", action="store_true",
                    help="不重算，只读回已有的 shared_ecm_check.json 重新判定并回写。"
                         "改了判定规则之后用它刷新，一秒钟的事")
    args = ap.parse_args()

    if args.reverdict:
        _reverdict(load_config(args.config), args.tag)
        return

    if not (args.slides or args.synthetic):
        ap.error("需要 --slides 或 --synthetic 之一")

    if not args.synthetic:
        try:
            import scanpy as sc
        except ImportError:
            sys.exit("需要 scanpy 读取签名分数。想先看输出长什么样："
                     "python scripts/run_14_shared_ecm_check.py --synthetic")

    from sparta.barrier import (compute_b_cell_field, compute_b_mab,
                                compute_b_meta, scores_from_adata)
    from sparta.validate import decoupling_stats

    cfg = load_config(args.config)
    var = load_config(args.variant)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    q = cfg["validate"]["decouple_q"]
    chance = (1.0 - q) ** 2

    def _cells(c):
        return (c["barrier"]["b_cell"],
                {k: v for k, v in c["barrier"]["b_mab"].items() if k != "r_nm"},
                c["barrier"]["b_mab"]["r_nm"])

    main_cell, main_mab, r_nm = _cells(cfg)
    var_cell, var_mab, _ = _cells(var)

    print("共享 ECM 通道检查")
    print(f"  主配置  ：b_cell.b_ecm={main_cell['b_ecm']}  c_caf={main_cell['c_caf']}"
          f"  |  b_mab.lam={main_mab['lam']}  beta={main_mab['beta']}")
    print(f"  变体    ：b_cell.b_ecm={var_cell['b_ecm']}  c_caf={var_cell['c_caf']}"
          f"  |  b_mab.lam={var_mab['lam']}  beta={var_mab['beta']}")
    if var_cell["b_ecm"] != 0 or var_mab["lam"] != 0:
        print("  [warn] 变体配置里 b_ecm 或 lam 不为 0，共享通道并未真正切断，"
              "结论不可用作稳健性证据。")
    print(f"  解离区随机期望 = (1-{q})^2 = {chance*100:.2f}%\n")

    if args.synthetic:
        from sparta.synthetic import make_dissociated_grid, make_ring_grid
        print("【合成演示模式】切片由 sparta.synthetic 现场生成，"
              "不含真实数据，也不构成任何研究结论。\n")
        synth = {"SYN_ring": make_ring_grid(n=25, width_spots=2.0, seed=0),
                 "SYN_diss": make_dissociated_grid(n=25, seed=0)}
        args.slides = list(synth)

    hdr = (f"{'slide':<10}{'口径':<12}{'rho':>9}{'rho_partial':>13}"
           f"{'解离区%':>9}{'富集':>8}")
    print(hdr); print("-" * len(hdr))

    rows = {}
    for sid in args.slides:
        if args.synthetic:
            sl = synth[sid]
            A, source, vessel = sl.A, sl.source, sl.vessel
            S = sl.scores
        else:
            try:
                adata = sc.read_h5ad(P.scored(sid))
                A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
            except FileNotFoundError as e:
                print(f"{sid:<10} 跳过：{e}")
                continue
            S = scores_from_adata(adata)

        # control 两种口径共用，且一律用主配置算——校正的必须是同一个几何量
        dv = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"], vessel,
                            D=(D if args.d_vessel == "weighted" else None),
                            **cfg["barrier"]["b_meta"])["d_vessel_um"]

        rec = {}
        for label, cc, cm in (("主配置", main_cell, main_mab),
                              ("切断共享ECM", var_cell, var_mab)):
            bcf = compute_b_cell_field(A, S["ecm"], S["caf"], source, **cc)["b_cell_field"]
            bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel,
                               r_nm=r_nm, **cm)["b_mab"]
            d = decoupling_stats(bcf, bm, q=q, control=dv)
            rec[label] = dict(rho=d["rho"], rho_partial=d["rho_partial"],
                              p_partial=d["p_partial"],
                              frac_discordant_r=d["frac_discordant_r"],
                              enrichment_vs_chance_r=d["enrichment_vs_chance_r"],
                              n=d["n"])
            print(f"{sid if label == '主配置' else '':<10}{label:<12}"
                  f"{d['rho']:>+9.3f}{d['rho_partial']:>+13.3f}"
                  f"{d['frac_discordant_r']*100:>9.1f}"
                  f"{d['enrichment_vs_chance_r']:>7.2f}x")
        rows[sid] = rec

    if not rows:
        sys.exit("没有任何切片可用。请先跑 run_01 -> run_02 -> run_03。")

    summary = _summarise(rows, hdr_width=len(hdr), pmap=patient_map(P))

    out = dict(q=q, chance_discordant=chance,
               config_main=str(args.config or "configs/default.yaml"),
               config_variant=str(args.variant),
               d_vessel_mode=args.d_vessel,
               per_slide=rows, **summary,
               meta=stamp_run(cfg, {"module": "M14-shared-ecm-check",
                                    "d_vessel_mode": args.d_vessel}))
    name = f"shared_ecm_check_{args.tag}.json" if args.tag else "shared_ecm_check.json"
    p = P.validation(name)
    save_json(p, out)
    print(f"\n已写出 {p}（新文件，未覆盖任何已有产物）")


if __name__ == "__main__":
    main()
