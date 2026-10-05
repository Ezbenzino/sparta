#!/usr/bin/env python
"""
run_22_null_crosslink.py —— Size-exclusion channel 的 null model 检验
=====================================================================
回答审稿人质疑："97.5% size-exclusion 占比是不是 rank-normalisation 构造出来的？"

做法：
  对每张切片，把 crosslink signature 随机打乱（permutation）N 次，
  每次重跑 dissociation_drivers，记录 frac_crosslink_raw 的 null 分布。
  如果真实值落在 null 分布的极端尾部（empirical p < 0.05），
  说明 size-exclusion 贡献确实来自真实 crosslink 空间分布，
  不是单纯 rank-normalisation 的数学产物；
  如果真实值和 null 分布重叠，说明这个数字确实是 by construction。

输出：results/validation/null_crosslink_check.json
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import (Paths, load_config, load_graph, patient_map,  # noqa: E402
                        save_json, set_seed, stamp_run)
from sparta.barrier import scores_from_adata  # noqa: E402
from sparta.validate import dissociation_drivers  # noqa: E402

N_PERM = 50
SEED = 20261001


def main():
    import scanpy as sc

    cfg = load_config(None)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    rng = np.random.default_rng(SEED)

    cfg_cell = cfg["barrier"]["b_cell"]
    cfg_mab = {k: v for k, v in cfg["barrier"]["b_mab"].items() if k != "r_nm"}
    r_nm = cfg["barrier"]["b_mab"]["r_nm"]

    # all ingested slides (mirrors run_25)
    ledger = Path(__file__).resolve().parents[1] / "data" / "ledger.csv"
    with open(ledger, encoding="utf-8") as f:
        SLIDES = [r["slide_id"] for r in csv.DictReader(f)
                  if r.get("status") == "ingested"]

    print(f"Null crosslink permutation test  (N={N_PERM} per slide, seed={SEED})")
    print(f"Slides: {SLIDES}\n")

    rows = {}
    for sid in SLIDES:
        try:
            adata = sc.read_h5ad(P.scored(sid))
            A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        except FileNotFoundError as e:
            print(f"{sid}: SKIP ({e})")
            continue
        S = scores_from_adata(adata)
        ecm = S["ecm"]; caf = S["caf"]; xl = S["crosslink"]; ag = S["ag_target"]

        # 真实值
        real = dissociation_drivers(A, ecm, caf, xl, ag, source, vessel,
                                    cfg_cell=cfg_cell, cfg_mab=cfg_mab, r_nm=r_nm)
        real_xl = real["frac_crosslink_raw"]

        # null 分布
        null_xl = []
        for _ in range(N_PERM):
            xl_perm = rng.permutation(xl)
            r = dissociation_drivers(A, ecm, caf, xl_perm, ag, source, vessel,
                                     cfg_cell=cfg_cell, cfg_mab=cfg_mab, r_nm=r_nm)
            null_xl.append(r["frac_crosslink_raw"])
        null_xl = np.array(null_xl)

        # empirical p：真实值比 null 中多少比例的 permutation 大
        emp_p = float((null_xl >= real_xl).mean())
        rows[sid] = dict(
            real_frac_crosslink=float(real_xl),
            null_mean=float(null_xl.mean()),
            null_std=float(null_xl.std()),
            null_ci95=[float(np.percentile(null_xl, 2.5)),
                       float(np.percentile(null_xl, 97.5))],
            empirical_p_one_sided=emp_p,
            n_perm=N_PERM,
            ecm_crosslink_corr=real["ecm_crosslink_corr"],
            crosslink_collinear=real["crosslink_collinear"],
        )
        print(f"{sid:<8} real={real_xl*100:5.1f}%  "
              f"null={null_xl.mean()*100:5.1f}±{null_xl.std()*100:4.1f}%  "
              f"95%CI=[{np.percentile(null_xl,2.5)*100:4.1f}, {np.percentile(null_xl,97.5)*100:4.1f}]%  "
              f"p={emp_p:.3f}")

    # 汇总
    reals = np.array([r["real_frac_crosslink"] for r in rows.values()])
    nulls = np.array([r["null_mean"] for r in rows.values()])
    ps = np.array([r["empirical_p_one_sided"] for r in rows.values()])
    summary = dict(
        median_real=float(np.median(reals)),
        median_null=float(np.median(nulls)),
        median_empirical_p=float(np.median(ps)),
        n_slides=len(rows),
        n_per_significant=int((ps < 0.05).sum()),
    )
    print("\n=== Summary ===")
    print(f"Median real size-exclusion share: {summary['median_real']*100:.1f}%")
    print(f"Median null share (permuted crosslink): {summary['median_null']*100:.1f}%")
    print(f"Median one-sided empirical p: {summary['median_empirical_p']:.3f}")
    print(f"Slides with p<0.05: {summary['n_per_significant']}/{summary['n_slides']}")

    out = dict(per_slide=rows, summary=summary,
               seed=SEED, n_perm=N_PERM,
               meta=stamp_run(cfg, {"module": "M22-null-crosslink-permutation"}))
    p = P.validation("null_crosslink_check.json")
    save_json(p, out)
    print(f"\nWrote {p}")


if __name__ == "__main__":
    main()
