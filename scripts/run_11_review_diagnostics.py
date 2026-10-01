#!/usr/bin/env python
"""
run_11_review_diagnostics.py —— 投稿前审查诊断（只读，不改任何产物）
=====================================================================

输入：{sid}.scored.h5ad + {sid}.graph.npz + results/counterfactual/{sid}.json
输出：results/validation/review_diagnostics.json（新文件，不覆盖任何已有产物）
上游模块：run_03_graph.py / run_05_counterfactual.py
下游模块：无（供 docs/review_for_journal.md 的三个论断做数值核验）

它核验三件事，全部是审稿人会自己动手算的：

  D1 尺寸排阻的"构造性阈值"
      crosslink 是片内秩标准化到 [0,1] 的量，而 xi = xi0*exp(-beta*x) 有纳米单位。
      本诊断报告：每张切片的 xi 实际范围、被完全排阻（phi=0）的边占比，
      以及 b_mab 落在数值上限 -log(1e-12)=27.631 的 spot 占比。
      若各切片的排阻边占比高度一致，说明它由秩分布而非生物学决定。

  D2 解耦统计的正确口径
      用 compute_b_cell_field + control=d_vessel 重算偏相关，
      并把"解离区占比"与随机期望 (1-q)*(1-q) 对照。
      q=0.75 时随机期望是 6.25% —— 低于它就不能称为"解离富集"。

  D3 S1 的效应量
      从 counterfactual JSON 读出 b_real / null_mean 的比值（比 z 更可解释），
      同时打印 n_perm 与经验 p 的分辨率下限 1/(n_perm+1)。

用法
----
    python scripts/run_11_review_diagnostics.py --slides MEL01 MEL02 MEL03 MEL04 \
                                                CSCC01 CSCC02 CSCC03 CSCC04
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import Paths, load_config, load_graph, save_json, set_seed  # noqa: E402

CAP = float(-np.log(1e-12))          # compute_b_mab 的数值上限
CHANCE_Q = 0.75                      # 与 validate.decouple_q 对齐


def main():
    ap = argparse.ArgumentParser(description="投稿前审查诊断（只读）",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slides", nargs="+", required=True)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    try:
        import scanpy as sc
    except ImportError:
        sys.exit("需要 scanpy 读取 .scored.h5ad。")

    from sparta.barrier import (compute_b_cell_field, compute_b_mab, compute_b_meta,
                                edge_pairs, scores_from_adata)
    from sparta.validate import decoupling_stats

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    cfg_cell = cfg["barrier"]["b_cell"]
    mab = cfg["barrier"]["b_mab"]
    cfg_mab = {k: v for k, v in mab.items() if k != "r_nm"}
    r_nm, xi0, beta = mab["r_nm"], mab["xi0_nm"], mab["beta"]

    # 构造性阈值：xi(x) = xi0*exp(-beta*x) 何时降到 r_nm 以下
    x_star = float(np.log(xi0 / r_nm) / beta) if xi0 > r_nm else 0.0
    chance = (1.0 - CHANCE_Q) ** 2

    print(f"[D1] 参数：r={r_nm} nm, xi0={xi0} nm, beta={beta}")
    print(f"     xi 的构造性范围 = [{xi0*np.exp(-beta):.2f}, {xi0:.1f}] nm（与切片无关）")
    print(f"     完全排阻的秩阈值 x* = ln(xi0/r)/beta = {x_star:.3f}\n")

    out = {"params": dict(r_nm=r_nm, xi0_nm=xi0, beta=beta,
                          xi_min_nm=float(xi0*np.exp(-beta)), x_star=x_star,
                          chance_discordant=chance), "per_slide": {}}

    hdr = (f"{'slide':<8}{'排阻边%':>9}{'xi中位':>9}{'b_mab封顶%':>11}"
           f"{'rho':>8}{'rho_part':>10}{'解离区%':>9}{'随机期望%':>10}")
    print(hdr); print("-" * len(hdr))

    for sid in args.slides:
        try:
            adata = sc.read_h5ad(P.scored(sid))
            A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        except FileNotFoundError as e:
            print(f"{sid:<8} 跳过：{e}")
            continue
        S = scores_from_adata(adata)

        # ---- D1 边级排阻 ----
        pairs = edge_pairs(A, upper_only=True)
        xl_e = 0.5 * (S["crosslink"][pairs[:, 0]] + S["crosslink"][pairs[:, 1]])
        xi = xi0 * np.exp(-beta * xl_e)
        excluded = float(np.mean(r_nm / np.maximum(xi, 1e-6) >= 1.0))

        bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel,
                           r_nm=r_nm, **cfg_mab)["b_mab"]
        at_cap = float(np.mean(np.isclose(bm, CAP, rtol=1e-6)))

        # ---- D2 正确口径的解耦 ----
        bcf = compute_b_cell_field(A, S["ecm"], S["caf"], source, **cfg_cell)["b_cell_field"]
        dv = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"], vessel,
                            D=D, **cfg["barrier"]["b_meta"])["d_vessel_um"]
        dec = decoupling_stats(bcf, bm, q=CHANCE_Q, control=dv)

        # ---- D3 S1 效应量 ----
        s1 = {}
        try:
            cf = json.load(open(P.counterfactual(sid), encoding="utf-8"))
            for mode, r in cf.get("s1", {}).items():
                ns = r.get("null_summary", {})
                mu = ns.get("mean", np.nan)
                s1[mode] = dict(z=r.get("z"), n_perm=r.get("n_perm"),
                                p_emp=r.get("p_emp"),
                                p_floor=1.0 / (int(r.get("n_perm", 0)) + 1) if r.get("n_perm") else None,
                                ratio_real_over_null=float(r["b_real"] / mu) if mu else None)
        except FileNotFoundError:
            pass

        out["per_slide"][sid] = dict(
            frac_edges_size_excluded=excluded,
            xi_median_nm=float(np.median(xi)),
            frac_bmab_at_cap=at_cap,
            rho=dec["rho"], rho_partial=dec["rho_partial"],
            p_partial=dec["p_partial"],
            frac_discordant_r=dec["frac_discordant_r"],
            s1=s1, n_spots=int(A.shape[0]), n_edges=int(len(pairs)),
        )
        print(f"{sid:<8}{excluded*100:>9.1f}{np.median(xi):>9.2f}{at_cap*100:>11.1f}"
              f"{dec['rho']:>8.3f}{dec['rho_partial']:>10.3f}"
              f"{dec['frac_discordant_r']*100:>9.1f}{chance*100:>10.2f}")

    print("-" * len(hdr))
    print("读法：")
    print("  · 排阻边% 若各片高度一致 -> 由秩分布决定，不是生物学信号（见 review §3 B3）")
    print("  · rho_part 为正 -> 控制几何后两屏障仍正相关，即耦合而非解离（见 §3 B1）")
    print("  · 解离区% 若 <= 随机期望 -> 不能称为解离富集")
    for mode in ("fixed", "follow"):
        vals = [v["s1"].get(mode, {}).get("ratio_real_over_null")
                for v in out["per_slide"].values() if v["s1"].get(mode)]
        vals = [v for v in vals if v]
        if vals:
            print(f"  · S1 {mode} 效应量（真实/零分布均值）中位 {np.median(vals):.2f}x "
                  f"范围 {min(vals):.2f}–{max(vals):.2f}x")

    p = P.validation("review_diagnostics.json")
    save_json(p, out)
    print(f"\n已写出 {p}（新文件，未覆盖任何已有产物）")


if __name__ == "__main__":
    main()
