#!/usr/bin/env python
"""
run_02_score.py —— M2 签名打分
================================

输入：data/interim/{slide_id}.qc.h5ad
输出：data/interim/{slide_id}.scored.h5ad
上游模块：run_01_qc.py
下游模块：run_03_graph.py

用签名打分刻画每个 spot 的细胞/基质/代谢状态，替代反卷积。
打分后必须秩标准化（写入 *_n 列），barrier.py 消费的是这些 _n 列。

做完这一步必须做的检查：把每个分数画在切片上目视确认空间分布合理
（恶性集中在瘤巢、内皮呈条索/点状、CAF 在间质带）。
这一步能发现基因名版本不一致、样本标签错配等问题，花十分钟能省一个月。

用法
----
    python scripts/run_02_score.py --slide MEL01 --tumor-type melanoma
    python scripts/run_02_score.py --slide SCC01 --tumor-type cscc --plot
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import Paths, load_config, save_json, set_seed, stamp_run  # noqa: E402


def _require_scanpy():
    """真实数据流程需要 scanpy；合成流程（run_00_demo）不需要。"""
    try:
        import scanpy  # noqa: F401
    except ImportError:
        sys.exit(
            "需要 scanpy 才能处理真实空间数据。\n"
            "  pip install scanpy squidpy\n"
            "如果只想先验证环境与算子，请跑 python scripts/run_00_demo.py（不需要 scanpy）。"
        )


def main():
    ap = argparse.ArgumentParser(description="M2 签名打分",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slide", required=True)
    ap.add_argument("--tumor-type", default=None, choices=["melanoma", "cscc", "bcc", "brca"],
                    help="不传则用 config 里的 signatures.tumor_type")
    ap.add_argument("--hypoxia-gmt", default=None,
                    help="MSigDB HALLMARK_HYPOXIA 的 gmt 文件。不传则用代码内的占位集合")
    ap.add_argument("--config", default=None)
    ap.add_argument("--plot", action="store_true", help="输出各签名的空间分布图（强烈建议）")
    args = ap.parse_args()

    _require_scanpy()
    import scanpy as sc
    from sparta.signatures import SCORE_COLUMNS, rank_normalize, score_all

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    tumor = args.tumor_type or cfg["signatures"]["tumor_type"]

    adata = sc.read_h5ad(P.qc(args.slide))
    print(f"[M2] {args.slide}：{adata.n_obs} spot，瘤种 {tumor}")

    # 缺氧基因集的来源：命令行 > config 的 signatures.hypoxia_gmt > 代码内占位集
    #
    # 2026-08-27 修：原先只认命令行，而 run_batch.py 与 run_ingest_legacy_cscc.py
    # 都不传这个参数，于是 config 里明明写着 HALLMARK_HYPOXIA，实际用的却是代码里
    # 那 10 个基因的 placeholder，只在日志里留一行警告。更糟的是事后无法追查——
    # 从前 run_02 不把「用了哪个集合」写进 uns。
    # 现在：config 兜底、来源写进 stamp；config 指了文件却读不到就直接报错，
    # 绝不静默退回 placeholder。两个队列用上不同的缺氧集会直接污染 B_meta 与 R3b，
    # 而这种污染在结果里长得跟生物学差异一模一样。
    gmt = args.hypoxia_gmt or cfg["signatures"].get("hypoxia_gmt")
    hyp, hyp_src = None, "placeholder"
    if gmt:
        gp = Path(gmt)
        if not gp.is_absolute():
            gp = Path(cfg["paths"]["root"]) / gp
        if not gp.exists():
            sys.exit(f"[M2] config 里指定的缺氧基因集不存在：{gp}\n"
                     f"     不会退回占位集——那会让不同批次的切片用上不同的基因集。\n"
                     f"     要么把文件放回去，要么显式把 signatures.hypoxia_gmt 置空。")
        with open(gp, encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if "HYPOXIA" in parts[0].upper():
                    hyp = parts[2:]
                    break
        if not hyp:
            sys.exit(f"[M2] {gp} 里找不到名字含 HYPOXIA 的通路行，请检查这个 gmt。")
        hyp_src = str(gp)
        print(f"[M2] 缺氧基因集：{len(hyp)} 个基因（来自 {gp}）")
    else:
        print("[M2] 警告：使用代码内的缺氧占位集合（placeholder，10 个基因）。"
              "正式分析请把 config 的 signatures.hypoxia_gmt 指向 MSigDB HALLMARK_HYPOXIA。")

    score_all(adata, tumor, hypoxia_genes=hyp,
              min_genes=cfg["signatures"]["min_genes_matched"],
              ctrl_size=cfg["signatures"]["ctrl_size"], seed=cfg["seed"])
    rank_normalize(adata)

    got = [c for c in SCORE_COLUMNS if c + "_n" in adata.obs]
    print(f"[M2] 已生成 {len(got)} 个秩标准化签名：{', '.join(got)}")
    missing = [c for c in ("Malignant", "CAF", "Endothelial", "ECM_core", "Ag_target")
               if c + "_n" not in adata.obs]
    if missing:
        print(f"[M2] 严重警告：核心签名缺失 {missing}，下游屏障算子将使用中性填充值，"
              f"结果不可用。请检查基因名版本。")

    adata.uns["sparta_run"] = stamp_run(cfg, {"module": "M2", "slide": args.slide,
                                              "tumor_type": tumor,
                                              "hypoxia_source": hyp_src,
                                              "hypoxia_n_genes": len(hyp) if hyp else 10})
    adata.write(P.scored(args.slide))
    print(f"[M2] 已写出 {P.scored(args.slide)}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        keys = [c for c in got]
        ncol = 4
        nrow = int(np.ceil(len(keys) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.2 * nrow))
        for ax, k in zip(np.ravel(axes), keys):
            xy = adata.obsm.get("spatial_um", adata.obsm.get("spatial"))
            sc_ = ax.scatter(xy[:, 0], xy[:, 1], c=adata.obs[k + "_n"], s=4,
                             cmap="viridis", linewidths=0)
            ax.set_title(k, fontsize=9); ax.set_aspect("equal")
            ax.set_xticks([]); ax.set_yticks([])
        for ax in np.ravel(axes)[len(keys):]:
            ax.axis("off")
        fig.tight_layout()
        p = P.figure(f"scores_{args.slide}.png")
        fig.savefig(p, dpi=130); plt.close(fig)
        print(f"[M2] 已保存 {p} —— 请目视检查空间分布是否合理，这一步不能跳过")


if __name__ == "__main__":
    main()
