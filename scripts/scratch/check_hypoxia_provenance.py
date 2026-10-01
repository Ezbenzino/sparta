#!/usr/bin/env python
"""
check_hypoxia_provenance.py —— 查清每张切片的 Hypoxia 签名到底用了哪个基因集

为什么需要这个脚本
------------------
run_02_score.py 从前只认命令行的 --hypoxia-gmt，而 run_batch.py 与
run_ingest_legacy_cscc.py 都不传它。于是 config 里明明写着
signatures.hypoxia_gmt: HALLMARK_HYPOXIA（200 基因），实际很可能用的是
signatures.py 里那 10 个基因的 placeholder ——而且从前不写进 uns，事后查不出来。

两个队列若用了不同的缺氧集，B_meta 会整体偏移，在 R3b 的队列比较里
长得跟生物学差异一模一样。所以必须先查清事实，再决定要不要重跑。

做法：拿 {slide}.qc.h5ad 分别用两个集合重算一遍 score_genes，
      再跟 {slide}.scored.h5ad 里存着的 Hypoxia 比 Spearman。
      哪个 r≈1，当初用的就是哪个。不改任何文件。

用法
----
    python scripts/scratch/check_hypoxia_provenance.py
    python scripts/scratch/check_hypoxia_provenance.py --slides MEL01 CSCC01
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402

from sparta.io_ import Paths, load_config  # noqa: E402

PLACEHOLDER = ["VEGFA", "SLC2A1", "CA9", "LDHA", "PGK1", "ADM",
               "NDRG1", "P4HA1", "BNIP3", "ANKRD37"]


def read_gmt(path):
    with open(path, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if "HYPOXIA" in parts[0].upper():
                return parts[2:]
    return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slides", nargs="+", default=None)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    import scanpy as sc
    from scipy.stats import spearmanr

    cfg = load_config(args.config)
    P = Paths(cfg)
    gmt = Path(cfg["paths"]["root"]) / cfg["signatures"]["hypoxia_gmt"]
    hall = read_gmt(gmt) if gmt.exists() else []
    print(f"HALLMARK 集：{len(hall)} 个基因（{gmt}）")
    print(f"placeholder：{len(PLACEHOLDER)} 个基因\n")

    slides = args.slides or sorted(
        p.name.split(".")[0] for p in P.interim.glob("*.scored.h5ad"))

    out = []
    print(f"{'切片':10} {'n_spot':>7} {'r(placeholder)':>15} {'r(HALLMARK)':>13}   判定")
    for sid in slides:
        sp, qp = P.scored(sid), P.qc(sid)
        if not (Path(sp).exists() and Path(qp).exists()):
            print(f"{sid:10} —— 缺 scored 或 qc，跳过")
            continue
        a = sc.read_h5ad(sp)
        if "Hypoxia" not in a.obs:
            print(f"{sid:10} —— scored 里没有 Hypoxia 列，跳过")
            continue
        stored = np.asarray(a.obs["Hypoxia"].values, float)

        rs = {}
        for tag, genes in (("placeholder", PLACEHOLDER), ("HALLMARK", hall)):
            if not genes:
                rs[tag] = float("nan"); continue
            q = sc.read_h5ad(qp)
            present = [g for g in genes if g in q.var_names]
            if len(present) < cfg["signatures"]["min_genes_matched"]:
                rs[tag] = float("nan"); continue
            sc.tl.score_genes(q, gene_list=present, score_name="_H",
                              ctrl_size=cfg["signatures"]["ctrl_size"],
                              random_state=cfg["seed"])
            rs[tag] = float(spearmanr(stored, np.asarray(q.obs["_H"].values, float)).statistic)

        rp, rh = rs["placeholder"], rs["HALLMARK"]
        if np.isnan(rp) and np.isnan(rh):
            verdict = "无法判定"
        elif np.nanmax([rp, rh]) < 0.99:
            verdict = "两个都不像 —— 可能用了第三个集合，需人工查"
        else:
            verdict = "placeholder" if (rp > rh) else "HALLMARK"
        print(f"{sid:10} {a.n_obs:7d} {rp:15.4f} {rh:13.4f}   {verdict}")
        out.append({"slide": sid, "n_spot": int(a.n_obs),
                    "r_placeholder": None if np.isnan(rp) else round(rp, 4),
                    "r_hallmark": None if np.isnan(rh) else round(rh, 4),
                    "verdict": verdict})

    print("\n注：r≈1.000 表示当初就是用的那个集合（同样的 seed 与 ctrl_size 能精确复现）。")

    # 结果自己落盘，不依赖 PowerShell 重定向——`*>` 会按控制台代码页
    # 把 UTF-8 字节再解一遍，中文全变乱码，两次踩过了。
    res = Path(cfg["paths"]["root"]) / "results" / "validation"
    res.mkdir(parents=True, exist_ok=True)
    fp = res / "hypoxia_provenance.json"
    fp.write_text(json.dumps(
        {"hallmark_gmt": str(gmt), "n_hallmark": len(hall),
         "n_placeholder": len(PLACEHOLDER), "slides": out},
        ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"已写出 {fp}")


if __name__ == "__main__":
    main()
