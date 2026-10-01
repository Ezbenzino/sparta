#!/usr/bin/env python
"""
run_00b_ingest.py —— 原始数据摄入与数据台账（第 1 步的主力工具）
==================================================================

输入：data/raw/{任意格式的下载}
输出：data/interim/{slide_id}.raw.h5ad  +  data/ledger.csv 中的一行
上游模块：无
下游模块：run_01_qc.py

这个脚本解决什么问题
--------------------
公开空间数据的格式是混乱的：第一代 ST 给 TSV（坐标编码在 spot 名里）、
Visium 给 spaceranger 目录、GEO 补充文件经常是 genes×spots 需要转置。
本脚本自动识别格式、纠正方向、提取坐标，统一转成标准 h5ad，
并把结果记进**数据台账**（data/ledger.csv）。

数据台账是 P1 阶段最重要的产出。它让你随时能回答：
  · 我一共下了几个数据集？其中几个三要素齐全？
  · 哪几个因为缺坐标/缺图像被剔除了？
  · 每张切片是什么瘤种、什么平台、治疗过没有？

用法
----
    # 只看一眼，不转换（筛选阶段用这个，很快）
    python scripts/run_00b_ingest.py --scan data/raw

    # 转换单个数据集并登记
    python scripts/run_00b_ingest.py --input data/raw/GSE144240_P2 \
        --slide SCC01 --cancer cscc --platform legacy_st \
        --source "GEO GSE144240" --treatment naive --site primary

    # 方向判断错了（n_obs 是基因数而不是 spot 数）时手工纠正
    python scripts/run_00b_ingest.py --input ... --slide SCC01 --transpose

    python scripts/run_00b_ingest.py --help
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402

from sparta.io_ import Paths, load_config, set_seed  # noqa: E402
from sparta.loaders import inspect_raw, load_raw  # noqa: E402

LEDGER_COLS = [
    # patient 与 replicate 是 2026-08-27 补的：同一患者的多张切片不是独立样本，
    # 台账里没有这一层的话，下游统计会把伪重复当样本量。详见 io_.patient_map()
    "slide_id", "patient", "replicate",
    "cancer_type", "platform", "source", "accession", "raw_path",
    "format", "coord_source", "has_counts", "has_coords", "has_image",
    "n_spots", "n_genes", "median_umi",
    "treatment", "site", "status", "notes",
]


def _ledger_path(P):
    return P.root / "data" / "ledger.csv"


def _load_ledger(P):
    p = _ledger_path(P)
    if p.exists():
        df = pd.read_csv(p)
        # 老台账没有新加的列，补空列而不是崩掉
        for c in LEDGER_COLS:
            if c not in df.columns:
                df[c] = ""
        return df
    return pd.DataFrame(columns=LEDGER_COLS)


def _save_ledger(P, df):
    p = _ledger_path(P)
    p.parent.mkdir(parents=True, exist_ok=True)
    df[LEDGER_COLS].to_csv(p, index=False)
    return p


def do_scan(root: Path):
    """扫描一个目录下的所有候选数据集，只检查不加载。"""
    subs = [d for d in sorted(Path(root).iterdir())
            if d.is_dir() or d.suffix == ".h5ad"]
    if not subs:
        print(f"{root} 下没有找到任何数据集目录或 .h5ad 文件")
        return

    print(f"扫描 {root}：{len(subs)} 个候选\n")
    print(f"{'数据集':<28}{'格式':<15}{'计数':<6}{'坐标':<6}{'图像':<6}{'判定'}")
    print("-" * 78)
    usable = 0
    for d in subs:
        info = inspect_raw(d)
        t = lambda b: "  ✓  " if b else "  ✗  "
        ok = info["has_counts"] and info["has_coords"]
        verdict = "可用" if ok and info["has_image"] else \
                  ("缺图像(C3)" if ok else "不可用")
        usable += bool(ok and info["has_image"])
        print(f"{d.name[:27]:<28}{info['format']:<15}"
              f"{t(info['has_counts'])}{t(info['has_coords'])}{t(info['has_image'])}{verdict}")
        for n in info["notes"]:
            if n.startswith("⚠"):
                print(f"    {n}")
    print("-" * 78)
    print(f"三要素齐全：{usable} / {len(subs)}")
    print("\n缺坐标(C2)的直接剔除；缺图像(C3)的可以先留着做数值分析，")
    print("但它们参与不了形态学验证，也画不出封锁线叠加图——统计时要单列。")


def main():
    ap = argparse.ArgumentParser(
        description="原始数据摄入与数据台账",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--scan", metavar="DIR",
                    help="只扫描该目录下所有候选数据集，不做转换（筛选阶段用）")
    ap.add_argument("--input", help="单个数据集目录或 .h5ad")
    ap.add_argument("--slide", help="切片 ID，用作后续所有文件的前缀")
    ap.add_argument("--cancer", default="melanoma",
                    choices=["melanoma", "cscc", "bcc", "other"])
    ap.add_argument("--platform", default="visium",
                    choices=["visium", "legacy_st", "visium_hd", "xenium", "other"])
    ap.add_argument("--source", default="", help='数据出处，如 "GEO GSE144240"')
    ap.add_argument("--accession", default="", help="具体的样本编号")
    ap.add_argument("--patient", default="",
                    help="患者 ID。**同一患者的多张切片必须填同一个值。**"
                         "统计以患者为独立单位，不填会把技术重复当独立样本")
    ap.add_argument("--replicate", default="",
                    help="同一患者内的重复编号，如 rep1 / rep2")
    ap.add_argument("--treatment", default="unknown",
                    choices=["naive", "post_treatment", "unknown"],
                    help="准入标准 C7：治疗后切片的屏障结构已改变，须单列或剔除")
    ap.add_argument("--site", default="unknown",
                    choices=["primary", "metastasis", "unknown"])
    ap.add_argument("--transpose", action="store_true",
                    help="强制转置（当自动判断把 spot 数和基因数弄反时用）")
    ap.add_argument("--no-transpose", action="store_true", help="强制不转置")
    ap.add_argument("--coords-file", default=None, help="显式指定坐标文件")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)

    if args.scan:
        do_scan(Path(args.scan))
        return

    if not (args.input and args.slide):
        ap.error("除 --scan 外，--input 与 --slide 都是必需的")

    src = Path(args.input)
    info = inspect_raw(src)
    print(f"[摄入] {src}")
    print(f"       格式={info['format']}  坐标来源={info['coord_source']}")
    for n in info["notes"]:
        print(f"       · {n}")

    if not info["has_counts"]:
        sys.exit("[摄入] 找不到计数矩阵，无法继续。")
    if not info["has_coords"]:
        sys.exit("[摄入] 没有空间坐标（准入标准 C2 不合格）。该切片不可用，"
                 "请在台账中记为 rejected。")

    tr = True if args.transpose else (False if args.no_transpose else None)
    try:
        adata = load_raw(src, coords_file=args.coords_file, transpose=tr)
    except ImportError:
        sys.exit("[摄入] 需要 anndata / scanpy。pip install anndata scanpy")
    except Exception as e:  # noqa: BLE001
        sys.exit(f"[摄入] 加载失败：{type(e).__name__}: {e}")

    import numpy as np
    tot = np.asarray(adata.X.sum(axis=1)).ravel()
    med = float(np.median(tot))
    print(f"[摄入] 结果：{adata.n_obs} spot × {adata.n_vars} 基因，中位 UMI {med:.0f}")

    # 方向自检 —— 这是最容易出错的一步，主动提醒
    if adata.n_obs > adata.n_vars:
        print(f"[摄入] ⚠ spot 数({adata.n_obs}) 大于基因数({adata.n_vars})。")
        print("       真实空间数据通常基因数远多于 spot 数。若这里反了，")
        print("       请加 --transpose 重跑，否则后续全部结果都是错的。")

    out = P.interim / f"{args.slide}.raw.h5ad"
    adata.write(out)
    print(f"[摄入] 已写出 {out}")

    # ---- 更新数据台账 ----
    led = _load_ledger(P)
    row = {
        "slide_id": args.slide, "patient": args.patient, "replicate": args.replicate,
        "cancer_type": args.cancer, "platform": args.platform,
        "source": args.source, "accession": args.accession, "raw_path": str(src),
        "format": info["format"], "coord_source": info["coord_source"],
        "has_counts": info["has_counts"], "has_coords": info["has_coords"],
        "has_image": info["has_image"],
        "n_spots": adata.n_obs, "n_genes": adata.n_vars, "median_umi": round(med, 1),
        "treatment": args.treatment, "site": args.site,
        "status": "ingested",
        "notes": " | ".join(info["notes"])[:300],
    }
    led = led[led["slide_id"] != args.slide]
    led = pd.concat([led, pd.DataFrame([row])], ignore_index=True)
    p = _save_ledger(P, led)
    print(f"[摄入] 台账已更新：{p}（当前 {len(led)} 条）")
    if not args.patient:
        print("[warn] 没给 --patient。同一患者的多张切片若不标注，"
              "下游统计会把它们当独立样本（伪重复）。强烈建议补上。")
    else:
        npat = led["patient"].replace("", pd.NA).dropna().nunique()
        print(f"[摄入] 台账现覆盖 {npat} 位患者 / {len(led)} 张切片")
    print(f"[摄入] 下一步：python scripts/run_01_qc.py --slide {args.slide} "
          f"--input {out} --platform {args.platform}"
          + (" --has-he" if info["has_image"] else "")
          + (" --treatment-known" if args.treatment != "unknown" else ""))


if __name__ == "__main__":
    main()
