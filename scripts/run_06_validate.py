#!/usr/bin/env python
"""
run_06_validate.py —— M6 四维验证
======================================

输入：全部切片的 .barrier.npz + H&E 图像 + 外部队列
输出：results/validation/*
上游模块：run_04_barrier.py / run_05_counterfactual.py
下游模块：无（成文用）

维度① 跨切片一致性  ② 同片 H&E 形态学  ③ 临床队列与头对头  ④ 生物学一致性

维度③ 的三件事缺一不可：多变量校正后仍显著、与已有签名头对头、
增量价值检验（似然比）。本课题在这一维度的优势在于 B_mAb 是已有方法
完全不覆盖的维度——重点报告它的独立贡献，而非三分量合并分数。

用法
----
    python scripts/run_06_validate.py --slides MEL01 MEL02 SCC01 --dim consistency
    python scripts/run_06_validate.py --slides MEL01 --dim morphology --he data/raw/MEL01/he.png
    python scripts/run_06_validate.py --dim sensitivity --slides MEL01
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import (Paths, load_config, load_graph, load_json, save_json,  # noqa: E402
                        set_seed, stamp_run)


def main():
    ap = argparse.ArgumentParser(description="M6 四维验证",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slides", nargs="+", required=True, help="切片 ID 列表")
    ap.add_argument("--dim", required=True,
                    choices=["consistency", "decoupling", "morphology", "sensitivity"],
                    help="要跑的验证维度。临床队列（维度③）请见 notebooks/，"
                         "因为不同队列的数据格式差异太大，不适合做成统一脚本")
    ap.add_argument("--he", default=None, help="morphology 维度：H&E 图像路径")
    ap.add_argument("--tag", default=None,
                    help="给输出文件加后缀，例如 --tag MEL 会写 consistency_MEL.json。"
                         "跨队列一致性没有意义，必须分队列跑并分别落盘，"
                         "否则第二次运行会把第一次的结果冲掉")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    from sparta.validate import (correlate_with_morphology, cross_slide_consistency,
                                 decoupling_stats, he_patch_features,
                                 parameter_sensitivity)

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)

    def load_b(sid):
        z = np.load(P.barrier(sid), allow_pickle=True)
        return {k: z[k] for k in z.files}

    def merge_into(name: str, rows: dict) -> dict:
        """读回已有结果再合并写出。

        2026-08-26 修：旧实现每次都整文件覆盖，导致跑完 CSCC 就把 MEL 的结果冲掉
        （decoupling.json 只剩 MEL、consistency.json 只剩 CSCC）。
        写论文汇总表时会以为另一个队列没跑过。
        """
        p = P.validation(name)
        old = load_json(p) if p.exists() else {}
        if not isinstance(old, dict):
            old = {}
        old.update(rows)
        return old

    # ---------------------------------------------------------- 维度①
    if args.dim == "consistency":
        fields = {s: load_b(s)["b_mab"] for s in args.slides}
        res = cross_slide_consistency(fields)
        print(f"[①] {len(res['slides'])} 张切片的 B_mAb 分布一致性")
        print(f"    平均相似度 {res['mean_similarity']:.3f}，最低 {res['min_similarity']:.3f}")
        print("    注意：一致性只在**同一队列内**有意义。跨队列（MEL vs CSCC）"
              "的分布差异同时包含平台与瘤种差异，不要合并成一个数字报告。")
        name = f"consistency_{args.tag}.json" if args.tag else "consistency.json"
        if not args.tag:
            print("    [warn] 没给 --tag：本次结果会写进 consistency.json 并覆盖上一次。"
                  "分队列跑时请用 --tag MEL / --tag CSCC。")
        save_json(P.validation(name),
                  {k: v for k, v in res.items() if k != "similarity"} |
                  {"similarity": res["similarity"].tolist(),
                   "slides_in_this_run": list(args.slides)})
        print(f"    已写出 {P.validation(name)}")

    # ---------------------------------------------------------- 解耦性
    elif args.dim == "decoupling":
        try:
            import scanpy as sc
        except ImportError:
            sys.exit("decoupling 维度需要 scanpy 读取签名分数。")
        from sparta.barrier import compute_b_cell_field, scores_from_adata

        # ⚠ 2026-08-26 修正的口径错误（docs/review_for_journal.md §3 B1）
        # 旧实现：proxy = 0.5*(ECM_core_n + CAF_n)，且**没有传 control**。
        #   1) ECM/CAF 局部代理是 validate.decoupling_stats 的 docstring 明文禁止的；
        #   2) 没传 control 时返回的 rho_partial 是 NaN，controlled=False，
        #      但产物却被当成"控制距血管距离后的偏相关"写进了论文骨架。
        # 正确口径：逐点 B_cell 用 compute_b_cell_field，control 用距血管图距离。
        rows = {}
        for s in args.slides:
            adata = sc.read_h5ad(P.scored(s))
            A, D, source, sink, vessel, _ = load_graph(P.graph(s))
            b = load_b(s)
            S = scores_from_adata(adata)
            bcf = compute_b_cell_field(A, S["ecm"], S["caf"], source,
                                       **cfg["barrier"]["b_cell"])["b_cell_field"]
            d = decoupling_stats(bcf, b["b_mab"], q=cfg["validate"]["decouple_q"],
                                 control=np.asarray(b["d_vessel_um"], float))
            print(f"[解耦] {s}: ρ={d['rho']:+.3f}｜偏相关 {d['rho_partial']:+.3f} "
                  f"(p={d['p_partial']:.2e})｜解离区 {d['frac_discordant_r']*100:.1f}%"
                  f"（随机期望 {d['chance_discordant']*100:.2f}%，"
                  f"富集 {d['enrichment_vs_chance_r']:.2f}x）")
            rows[s] = {k: v for k, v in d.items()
                       if k not in ("idx_discordant", "idx_discordant_r")}
        save_json(P.validation("decoupling.json"), merge_into("decoupling.json", rows))
        print("    读法：偏相关为正 = 控制几何后两屏障仍同向，即耦合而非解离；")
        print("          解离区富集 <= 1.0x = 解离区比随机还少，不能称为解离。")

    # ---------------------------------------------------------- 维度②
    elif args.dim == "morphology":
        if not args.he:
            sys.exit("morphology 维度需要 --he 指定 H&E 图像路径")
        try:
            import scanpy as sc
            from skimage.io import imread
        except ImportError:
            sys.exit("需要 scanpy 与 scikit-image。pip install scanpy scikit-image")
        sid = args.slides[0]
        adata = sc.read_h5ad(P.scored(sid))
        img = imread(args.he)
        coords_px = np.asarray(adata.obsm["spatial"], int)  # 注意：这里用像素坐标
        print(f"[②] {sid}：图像 {img.shape}，{len(coords_px)} 个 spot，"
              f"patch {cfg['validate']['he_patch_px']} px")
        morph = he_patch_features(img, coords_px, patch_px=cfg["validate"]["he_patch_px"])
        b = load_b(sid)
        corr = correlate_with_morphology(b["b_mab"], morph)
        for k, v in corr.items():
            print(f"    B_mAb vs {k:22s} ρ={v['rho']:+.3f}  p={v['p']:.2e}  n={v['n']}")
        print("    预期方向：高屏障区应对应更高的基质占比、更低的核密度。")
        print("    方向不符时先检查 patch 是否裁在了正确位置，不要硬解释。")
        save_json(P.validation(f"morphology_{sid}.json"), corr)

    # ---------------------------------------------------------- 敏感性
    elif args.dim == "sensitivity":
        try:
            import scanpy as sc
        except ImportError:
            sys.exit("sensitivity 维度需要 scanpy。")
        from sparta.barrier import compute_b_cell, scores_from_adata
        sid = args.slides[0]
        adata = sc.read_h5ad(P.scored(sid))
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        S = scores_from_adata(adata)
        grid = {k: v for k, v in cfg["validate"]["sensitivity_grid"].items()
                if k in ("b_ecm", "c_caf")}
        rows = parameter_sensitivity(
            lambda **kw: compute_b_cell(A, S["ecm"], S["caf"], source, sink,
                                        a=cfg["barrier"]["b_cell"]["a"], **kw)["b_cell"],
            grid=grid,
        )
        for r in rows:
            print(f"    {r}")
        save_json(P.validation(f"sensitivity_{sid}.json"), rows)
        print("    把这张表画成热图，作为补充材料的必备一图。")


if __name__ == "__main__":
    main()
