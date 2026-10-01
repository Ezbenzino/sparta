#!/usr/bin/env python
"""
run_12_paper_stats.py —— 论文统计汇总（把散落的临时脚本收编进管线）
=====================================================================

输入：results/counterfactual/*.json + results/validation/{screen_decision,decoupling,
      benchmark_ext}.json
输出：results/validation/paper_stats.json + 终端汇总表
上游模块：run_05_counterfactual.py / run_06_validate.py / run_07_screen.py / run_10_benchmark_ext.py
下游模块：无（成文用）

为什么有这个脚本
----------------
2026-08-26 审查发现：论文里的 FDR 校正、S2 汇总、跨队列汇总表，
全部只能由 scripts/_*.py 这批一次性脚本产生，而它们硬编码了 `D:\\sparta\\...`
绝对路径——README 里那条 run_01 -> run_10 的管线**复现不出论文的数字**。
本脚本把这些统计收编进正式管线，路径一律走 io_.Paths。

它不做任何新计算，只汇总已有产物，因此可以随时重跑，秒级完成。

三个口径上的修正（相对 _fdr_correction.py）
-------------------------------------------
1. **检验族按切片划分。** 旧版把同一张切片的 fixed/follow 当 2 个独立检验、
   把 4 个 k 当 4 个独立检验，family 被人为放大到 16 / 32。
   这些检验共享同一张图、同一个割集，强相关，BH 的前提不成立。
   新版：每张切片先取该切片内最小的 p 作为代表（S1 分模式分别成族；
   S2 每张切片一个代表），family = 切片数。
2. **报告 p 值的分辨率下限。** 经验 p 的下限是 1/(n_perm+1)。
   MEL 跑了 500 次（下限 0.002），CSCC 只跑了 150 次（下限 0.0066），
   两者不可直接比较——本脚本把下限一并打印出来。
3. **效应量优先于 z。** z 大可能只是零分布方差小。
   S1 报 b_real/null_mean 的比值，S2 报 ratio_vs_in_cut 与移除的割集比例。

用法
----
    python scripts/run_12_paper_stats.py
    python scripts/run_12_paper_stats.py --slides MEL01 MEL02 CSCC01
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import (Paths, admitted_slides, load_config, load_json,  # noqa: E402
                        patient_map, save_json, stamp_run)

# 切片名单以台账为准（status == ingested），不再硬编码。
# 硬编码的表现是队列一变就漏切片，而且**安静地漏**：图和统计少几张，不报错。
# 见 sparta/io_.py::admitted_slides


def bh(pvals) -> np.ndarray:
    """Benjamini–Hochberg 校正。返回与输入同序的 adjusted p。"""
    p = np.asarray(pvals, float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    adj = np.minimum(1.0, p[order] * n / np.arange(1, n + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = adj
    return out


def main():
    ap = argparse.ArgumentParser(description="论文统计汇总",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slides", nargs="+", default=None,
                    help="不传则取台账里 status=ingested 的全部切片")
    ap.add_argument("--primary-frac", type=float, default=0.20,
                    help="S2 的主检验用哪个移除比例（占割集）。**预先指定，不要按结果挑**。"
                         "每张切片取离它最近的那个 k，其余 k 作为剂量-反应报告")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    P = Paths(cfg)
    slides = args.slides or admitted_slides(P)
    if not slides:
        sys.exit("台账里没有 status=ingested 的切片。先跑 run_00b_ingest.py。")
    print(f"[stats] 切片 {len(slides)} 张：{', '.join(slides)}")
    out = {"slides": [], "s1": [], "s2": [], "decoupling": [], "screen": [],
           "meta": stamp_run(cfg, {"module": "M12-paper-stats"})}

    cf = {}
    for sid in slides:
        p = P.counterfactual(sid)
        if p.exists():
            cf[sid] = load_json(p)
            out["slides"].append(sid)
    if not cf:
        sys.exit("没有找到任何 results/counterfactual/*.json。请先跑 run_05。")

    # ---- 患者层 ----------------------------------------------------------
    # 同一患者的多张切片不是独立样本。所有"n=8"的说法都要同时给出患者数，
    # 否则就是把伪重复当样本量（2026-08-27 查 GEO 才发现本项目 8 张 = 3 位患者）。
    pmap = patient_map(P)
    pats = {}
    for sid in out["slides"]:
        pats.setdefault(pmap.get(sid, sid), []).append(sid)
    out["patients"] = {k: v for k, v in pats.items()}
    print("=" * 92)
    print(f"【样本量】{len(out['slides'])} 张切片，来自 **{len(pats)} 位患者**")
    print("=" * 92)
    for k, v in pats.items():
        print(f"  {k:<12} {len(v)} 张：{', '.join(v)}")
    if len(pats) < len(out["slides"]):
        print("\n  ⚠ 存在同患者多切片。下面所有按切片成族的 FDR 都是**切片层**的，")
        print("    正文报告时必须同时给出患者层的结论，不能只说 'n/8 significant'。")

    # ------------------------------------------------------------------ S1
    print("=" * 92)
    print("【S1】空间重排对照 —— 屏障是否来自阻力物质的空间排布")
    print("=" * 92)
    print(f"{'slide':<8}{'mode':<8}{'b_real':>9}{'null均值':>10}{'效应量':>9}"
          f"{'z':>8}{'p':>9}{'p下限':>8}{'p_fdr':>9}")
    print("-" * 92)
    for mode in ("fixed", "follow"):
        rows = []
        for sid, d in cf.items():
            r = (d.get("s1") or {}).get(mode)
            if not r or not np.isfinite(r.get("p_emp", np.nan)):
                continue
            ns = r.get("null_summary", {})
            mu = ns.get("mean") or np.nan
            n_perm = int(r.get("n_perm") or ns.get("n") or 0)
            rows.append(dict(slide=sid, mode=mode, b_real=r["b_real"], null_mean=mu,
                             ratio=float(r["b_real"] / mu) if mu else np.nan,
                             z=r.get("z"), p=r["p_emp"],
                             p_floor=1.0 / (n_perm + 1) if n_perm else np.nan,
                             n_perm=n_perm))
        # 检验族 = 该模式下的切片数（不与另一个模式混族：它们检验的零假设不同）
        for r, a in zip(rows, bh([r["p"] for r in rows])):
            r["p_fdr"] = float(a)
            sig = "*" if r["p_fdr"] < 0.05 else " "
            print(f"{r['slide']:<8}{r['mode']:<8}{r['b_real']:>9.4f}{r['null_mean']:>10.4f}"
                  f"{r['ratio']:>8.2f}x{r['z']:>8.1f}{r['p']:>9.4f}"
                  f"{r['p_floor']:>8.4f}{r['p_fdr']:>8.4f}{sig}")
        out["s1"] += rows
    floors = sorted({round(r["p_floor"], 4) for r in out["s1"] if np.isfinite(r["p_floor"])})
    if len(floors) > 1:
        print(f"\n⚠ 各切片的 p 值分辨率下限不一致：{floors}。"
              f"\n  说明置换次数在队列之间不同（config 里写的是 "
              f"{cfg['counterfactual']['s1']['n_perm']}），两队列的 p 不可直接比较。"
              f"\n  投稿前必须用同一个 n_perm 重跑 run_05。")
    print("\n注：follow 模式把源汇也一起重新推导，其零分布测的是"
          "\"完全没有空间结构\"，")
    print("    与 fixed 模式（只打乱阻力物质排布）不是同一个零假设，不能合并计数，")
    print("    也不能用 follow 的显著去补 fixed 的不显著。")

    # ------------------------------------------------------------------ S2
    print("\n" + "=" * 92)
    print("【S2】环带断裂 —— 连续缺口 vs 等量分散移除")
    print("=" * 92)
    print(f"{'slide':<8}{'割集':>6}  {'k模式':<10}{'k':>5}{'占割集':>8}"
          f"{'ratio':>8}{'p':>8}")
    print("-" * 92)
    per_slide = []
    for sid, d in cf.items():
        s2 = d.get("s2") or {}
        pk = s2.get("per_k") or {}
        if not pk:
            continue
        # 旧结果（2026-08-26 之前）里没有 k_over_cut，这里按割集规模补算，
        # 好让"到底动了割集的百分之几"在汇总表里一眼可见。
        denom = max(int(s2.get("n_cut_pool") or s2.get("n_cut_nodes") or 1), 1)

        def _kov(k, r):
            v = r.get("k_over_cut")
            return float(v) if v is not None else int(k) / denom

        # 每张切片只出一个代表检验（4 个 k 是同一个割集的嵌套子集，不是 4 个独立检验）。
        # 但**不能取 p 最小的那个**——那是从 4 个里挑出来的，等于隐性多重比较，
        # 审稿人一句 "you picked the best k" 就能把显著性打掉。
        # 正确做法：预先指定一个移除比例（--primary-frac，默认 20%），
        # 所有切片都用同一个比例，其余 k 作为剂量-反应曲线报告。
        target = args.primary_frac
        best_k, best_r = min(pk.items(),
                             key=lambda kv: abs(_kov(kv[0], kv[1]) - target))
        first = sorted(pk, key=int)[0]
        for k, r in sorted(pk.items(), key=lambda kv: int(kv[0])):
            head = (f"{sid:<8}{s2.get('n_cut_nodes', 0):>6}  " if k == first else " " * 16)
            mark = f"  <- 代表（预设 {args.primary_frac:.0%}）" if k == best_k else ""
            print(f"{head}{s2.get('k_mode', 'absolute'):<10}{int(k):>5}"
                  f"{_kov(k, r) * 100:>7.1f}%"
                  f"{r['ratio_vs_in_cut']:>8.3f}{r['p_vs_in_cut']:>8.3f}{mark}")
        per_slide.append(dict(slide=sid, k_mode=s2.get("k_mode", "absolute"),
                              n_cut_nodes=s2.get("n_cut_nodes"),
                              n_cut_pool=s2.get("n_cut_pool"),
                              best_k=int(best_k),
                              k_over_cut=_kov(best_k, best_r),
                              ratio=best_r["ratio_vs_in_cut"],
                              p=best_r["p_vs_in_cut"],
                              all_k={int(k): dict(k_over_cut=_kov(k, v),
                                                  ratio=v["ratio_vs_in_cut"],
                                                  p=v["p_vs_in_cut"])
                                     for k, v in pk.items()}))
    for r, a in zip(per_slide, bh([r["p"] for r in per_slide])):
        r["p_fdr"] = float(a)
    out["s2"] = per_slide
    print(f"\n每片代表检验（预设移除比例 {args.primary_frac:.0%}）的 BH 校正：")
    for r in per_slide:
        print(f"  {r['slide']:<8} k={r['best_k']:>3}（{r['k_over_cut']*100:>4.1f}% 割集）"
              f"  ratio={r['ratio']:.3f}x  p={r['p']:.3f}  p_fdr={r['p_fdr']:.3f}"
              f"{'  *' if r['p_fdr'] < 0.05 else ''}")

    # 剂量-反应：效应量是否随移除比例单调上升。这是拓扑主张的一个可证伪推论，
    # 比任何单个 k 的显著性都更有说服力——随机噪声不会给出单调的剂量曲线。
    # 列用 config 里预设的 k_frac，而不是从实测比例里去重——
    # 实测比例会因为 round() 把 20.6% 单独拆成一列，表就花了
    fracs = [float(f) for f in (cfg["counterfactual"]["s2"].get("k_frac")
                                or [0.05, 0.10, 0.20, 0.30])]
    print("\n剂量-反应（同一张切片、不同移除比例）：")
    hdr2 = f"  {'slide':<8}" + "".join(f"{f:>9.0%}" for f in fracs)
    print(hdr2); print("  " + "-" * (len(hdr2) - 2))
    mono = 0
    for r in per_slide:
        cells, series = "", []
        for f in fracs:
            hit = sorted((v for v in r["all_k"].values()
                          if abs(v["k_over_cut"] - f) < f * 0.35),
                         key=lambda v: abs(v["k_over_cut"] - f))
            if hit:
                cells += f"{hit[0]['ratio']:>9.3f}"; series.append(hit[0]["ratio"])
            else:
                cells += f"{'—':>9}"
        if len(series) >= 2 and series[-1] > series[0]:
            mono += 1
        print(f"  {r['slide']:<8}{cells}")
    print(f"  -> 最大移除比例下的效应量高于最小比例的切片：{mono}/{len(per_slide)}")

    p_floors = sorted({round(min(v["p"] for v in r["all_k"].values()), 4) for r in per_slide})
    print(f"\n注：S2 经验 p 的分辨率下限约 {min(p_floors):.3f}"
          f"（= 1/(n_rand//2 + 1)）。落在下限上的 p 只说明'比全部随机对照都强'，"
          f"不代表更小的 p 值。")
    n_sig = sum(1 for r in per_slide if r["p_fdr"] < 0.05)
    ratios = [r["ratio"] for r in per_slide]
    print("-" * 92)
    print(f"按切片成族（n={len(per_slide)}）：FDR 后显著 {n_sig}/{len(per_slide)}；"
          f"代表效应量中位 {np.median(ratios):.3f}x，范围 {min(ratios):.3f}–{max(ratios):.3f}x")
    modes = {r["k_mode"] for r in per_slide}
    if "absolute" in modes:
        covs = [(r["slide"], (r["k_over_cut"] or 0) * 100) for r in per_slide]
        print("\n⚠ 有切片仍在用绝对 k。各片实际移除的割集比例："
              + "、".join(f"{s} {c:.1f}%" for s, c in covs))
        print("  比例差异越大，跨切片比较效应量就越没有意义。"
              "请在 config 里启用 counterfactual.s2.k_frac 后重跑 run_05。")

    # ------------------------------------------------------------------ 解耦
    dec = load_json(P.validation("decoupling.json")) if P.validation("decoupling.json").exists() else {}
    scr = load_json(P.validation("screen_decision.json")) if P.validation("screen_decision.json").exists() else {}
    q = cfg["validate"]["decouple_q"]
    chance = (1 - q) ** 2
    print("\n" + "=" * 92)
    print(f"【解耦】B_cell 与 B_mAb（解离区的随机期望 = (1-{q})^2 = {chance*100:.2f}%）")
    print("=" * 92)
    print(f"{'slide':<8}{'来源':<16}{'rho':>9}{'rho_partial':>13}{'解离区%':>9}{'富集':>8}")
    print("-" * 92)
    for sid in out["slides"]:
        got = None
        if isinstance(dec, dict) and sid in dec and dec[sid].get("controlled"):
            got, src = dec[sid], "run_06(已校正)"
        elif isinstance(dec, dict) and sid in dec:
            got, src = dec[sid], "run_06(未校正!)"
        elif scr.get("per_slide", {}).get(sid):
            got, src = scr["per_slide"][sid], "run_07"
        if not got:
            continue
        fdr_ = got.get("frac_discordant_r")
        fdr_ = fdr_ if fdr_ is not None and np.isfinite(fdr_) else got.get("frac_discordant", np.nan)
        rp = got.get("rho_partial", np.nan)
        rp = float(rp) if rp is not None and np.isfinite(rp) else np.nan
        enr = fdr_ / chance if np.isfinite(fdr_) else np.nan
        rp_s = f"{rp:+.3f}" if np.isfinite(rp) else "未校正"
        print(f"{sid:<8}{src:<16}{got.get('rho', np.nan):>+9.3f}{rp_s:>13}"
              f"{fdr_*100:>9.1f}{enr:>7.2f}x")
        out["decoupling"].append(dict(slide=sid, source=src, rho=got.get("rho"),
                                      rho_partial=rp, frac_discordant=fdr_,
                                      chance=chance, enrichment=enr))
    if out["decoupling"]:
        rps = [r["rho_partial"] for r in out["decoupling"] if r["rho_partial"] is not None
               and np.isfinite(r["rho_partial"])]
        enrs = [r["enrichment"] for r in out["decoupling"] if np.isfinite(r["enrichment"])]
        if rps:
            print(f"\n偏相关：{sum(1 for r in rps if r > 0)}/{len(rps)} 为正，"
                  f"范围 {min(rps):+.3f} ~ {max(rps):+.3f}")
        if enrs:
            print(f"解离区富集：中位 {np.median(enrs):.2f}x"
                  f"（<=1.0x 表示解离区比随机还少，不能称为解离）")

    # ------------------------------------------------------------------ 决策点
    if scr.get("per_slide"):
        out["screen"] = [dict(slide=s, **{k: v for k, v in r.items()
                                          if k in ("dissociation_potential", "frac_shared_ecm",
                                                   "frac_antigen", "frac_crosslink",
                                                   "ecm_crosslink_corr", "crosslink_collinear")})
                         for s, r in scr["per_slide"].items()]
        pots = [r["dissociation_potential"] for r in out["screen"]]
        print("\n" + "=" * 92)
        print("【驱动分解】B_mAb 变异来源（提醒：交联占比几乎只是 beta 的单调函数，"
              "见 review §3 B3）")
        print("=" * 92)
        print(f"  中位解离潜力比 {np.median(pots):.2f}（n={len(pots)}）；"
              f"交联占比 {min(r['frac_crosslink'] for r in out['screen'])*100:.1f}–"
              f"{max(r['frac_crosslink'] for r in out['screen'])*100:.1f}%")
        print(f"  当前 beta = {cfg['barrier']['b_mab']['beta']}，"
              f"身份已下调为[定性尺度参数]——报告时必须写明该占比对 beta 的依赖。")

    # ------------------------------------------------------------------ 患者层结论
    if len(pats) < len(out["slides"]):
        print("\n" + "=" * 92)
        print("【患者层】把切片层结果折叠到患者（同患者取全部切片一致才算该患者阳性）")
        print("=" * 92)
        for label, rows_, key, thr in (
                ("S1 fixed 显著", [r for r in out["s1"] if r["mode"] == "fixed"], "p_fdr", 0.05),
                ("S2 显著", out["s2"], "p_fdr", 0.05)):
            byp = {}
            for r in rows_:
                byp.setdefault(pmap.get(r["slide"], r["slide"]), []).append(r[key] < thr)
            allpos = sum(1 for v in byp.values() if all(v))
            anypos = sum(1 for v in byp.values() if any(v))
            detail = "；".join(f"{k} {sum(v)}/{len(v)}" for k, v in byp.items())
            print(f"  {label:<14} 全部切片阳性的患者 {allpos}/{len(byp)}"
                  f"（至少一张阳性 {anypos}/{len(byp)}）　{detail}")
        if out["decoupling"]:
            byp = {}
            for r in out["decoupling"]:
                rp = r.get("rho_partial")
                if rp is not None and np.isfinite(rp):
                    byp.setdefault(pmap.get(r["slide"], r["slide"]), []).append(rp > 0)
            allpos = sum(1 for v in byp.values() if all(v))
            print(f"  {'偏相关为正':<14} 全部切片为正的患者 {allpos}/{len(byp)}")

    p = P.validation("paper_stats.json")
    save_json(p, out)
    print(f"\n已写出 {p}")


if __name__ == "__main__":
    main()
