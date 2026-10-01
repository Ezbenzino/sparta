#!/usr/bin/env python
"""
run_01_qc.py —— M1 数据质控与准入核查
============================================

输入：data/raw/{slide_id}/  （Visium 目录或 h5ad）
输出：data/interim/{slide_id}.qc.h5ad + results/figures/qc_{slide_id}.png
上游模块：无
下游模块：run_02_score.py

这一步做两件事：标准质控，以及**逐条核查准入标准 C1–C7**。

准入核查比质控本身更重要。相当比例的公开切片缺少空间坐标或配对 H&E，
或者内皮信号弱到无法定义血管源集（C6）。这些问题必须在建库阶段发现，
而不是在第五个月建图时才发现——那时数据工作已经白做了两个月。

用法
----
    python scripts/run_01_qc.py --slide MEL01 --input data/raw/MEL01
    python scripts/run_01_qc.py --slide MEL01 --input data/raw/MEL01 --platform legacy_st
    python scripts/run_01_qc.py --help
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


ADMISSION_DESC = {
    "C1": "有 count 矩阵",
    "C2": "有空间坐标",
    "C3": "有配对 H&E 图像",
    "C4": "spot 数达标",
    "C5": "中位 UMI 达标",
    "C6": "存在可识别内皮/血管信号（最易忽略、最致命）",
    "C7": "已知治疗状态与原发/转移",
}

ENDOTHELIAL = ["PECAM1", "VWF", "CDH5", "CLDN5", "ENG", "EGFL7"]


def check_admission(adata, cfg, platform: str, has_he: bool, treatment_known: bool,
                    cohort: str | None = None) -> tuple[dict, dict]:
    """逐条核查 C1–C7，返回 ({条目: (通过与否, 说明)}, 实际生效的阈值)。

    cohort : 队列名。若 configs/default.yaml 的 admission_overrides 里有同名段，
             则用它覆盖默认阈值。**放宽阈值必须走这条路**——写在 config 里、
             注明日期与理由，而不是靠命令行 --force 悄悄放行。
             生效阈值会一并写进 {slide}.admission.json，使准入过程可审计。
    """
    import numpy as np
    adm = dict(cfg["admission"])
    ov = (cfg.get("admission_overrides") or {}).get(cohort or "", {})
    _THRESHOLD_KEYS = ("min_spots_visium", "min_spots_legacy_st",
                       "min_median_umi", "min_endothelial_spots")
    applied = []
    for k in _THRESHOLD_KEYS:
        if k in ov:
            adm[k] = ov[k]
            applied.append(k)
    res = {}

    res["C1"] = (adata.n_vars > 0 and adata.n_obs > 0, f"{adata.n_obs} spot × {adata.n_vars} 基因")

    has_coord = ("spatial" in adata.obsm) or ("spatial_um" in adata.obsm)
    res["C2"] = (has_coord, "有坐标" if has_coord else "缺少 obsm['spatial']，无法建图")

    res["C3"] = (has_he, "有配对 H&E" if has_he else "缺少 H&E，无法做形态学验证与屏障线叠加图")

    min_spots = adm["min_spots_legacy_st"] if platform == "legacy_st" else adm["min_spots_visium"]
    res["C4"] = (adata.n_obs >= min_spots, f"{adata.n_obs} spot（阈值 {min_spots}）")

    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    tot = np.asarray(counts.sum(axis=1)).ravel()
    med = float(np.median(tot))
    res["C5"] = (med >= adm["min_median_umi"], f"中位 UMI {med:.0f}（阈值 {adm['min_median_umi']}）")

    present = [g for g in ENDOTHELIAL if g in adata.var_names]
    if present:
        X = adata[:, present].X
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        n_pos = int((X.sum(axis=1) > 0).sum())
    else:
        n_pos = 0
    res["C6"] = (n_pos >= adm["min_endothelial_spots"],
                 f"{n_pos} 个 spot 有内皮信号（阈值 {adm['min_endothelial_spots']}，"
                 f"匹配到 {len(present)}/{len(ENDOTHELIAL)} 个标记基因）")

    res["C7"] = (treatment_known, "已知" if treatment_known else "未知——治疗后切片的屏障结构已改变")
    # 只报真正改过的阈值键——reason / decided_on 是给人看的元信息，不是阈值
    return res, dict(adm, _cohort=cohort, _overridden=applied,
                     _override_reason=ov.get("reason"),
                     _override_decided_on=ov.get("decided_on"))


def main():
    ap = argparse.ArgumentParser(description="M1 数据质控与准入核查",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slide", required=True, help="切片 ID，用作所有中间文件的前缀")
    ap.add_argument("--input", required=True, help="Visium 目录或 .h5ad 文件路径")
    ap.add_argument("--platform", default="visium", choices=["visium", "legacy_st"],
                    help="平台。影响 spot 数阈值与坐标间距换算")
    ap.add_argument("--config", default=None, help="配置文件路径")
    ap.add_argument("--has-he", action="store_true", help="声明该切片有配对 H&E 图像")
    ap.add_argument("--treatment-known", action="store_true", help="声明已知治疗状态")
    ap.add_argument("--cohort", default=None,
                    help="队列名，用于查 configs/default.yaml 的 admission_overrides。"
                         "例：cscc_gse144239。放宽阈值请走这条路，不要用 --force")
    ap.add_argument("--force", action="store_true",
                    help="即使准入不合格也继续（会在 admission.json 里记 forced=true，"
                         "并要求用 --force-reason 写明理由）")
    ap.add_argument("--force-reason", default=None,
                    help="--force 的书面理由，会原样写进 admission.json 供审稿核查")
    args = ap.parse_args()

    _require_scanpy()
    import scanpy as sc

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)

    src = Path(args.input)
    print(f"[M1] 读入 {src}")
    adata = sc.read_h5ad(src) if src.suffix == ".h5ad" else sc.read_visium(src)
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()

    # --- 准入核查（在过滤之前做，反映数据的原始状态）---
    print("\n[M1] 准入核查 C1–C7")
    adm, thresholds = check_admission(adata, cfg, args.platform, args.has_he,
                                      args.treatment_known, cohort=args.cohort)
    for k in sorted(adm):
        ok, msg = adm[k]
        print(f"   {'✓' if ok else '✗'} {k} {ADMISSION_DESC[k]}：{msg}")
    failed = [k for k, (ok, _) in adm.items() if not ok]
    if thresholds.get("_overridden"):
        print(f"   ! 已按队列 '{args.cohort}' 放宽阈值：{thresholds['_overridden']}"
              f"（理由见 configs/default.yaml 的 admission_overrides）")

    # 2026-08-26 修：**无论通过与否都写 admission.json**。
    # 旧实现只在失败分支写，而失败分支又被 --force 跳过，
    # 结果是 8 张切片里只有 1 张留下了准入记录，Methods 里"全部通过准入"
    # 这句话在磁盘上没有任何证据。这是审稿人与编辑最容易致命一击的地方。
    record = {
        "slide": args.slide,
        "checks": {k: {"pass": bool(v[0]), "detail": v[1]} for k, v in adm.items()},
        "failed": failed,
        "admitted": bool(not failed or args.force),
        "forced": bool(failed and args.force),
        "force_reason": args.force_reason,
        "cohort": args.cohort,
        "thresholds_applied": {k: v for k, v in thresholds.items()
                               if not k.startswith("_")},
        "thresholds_overridden": thresholds.get("_overridden", []),
        "meta": stamp_run(cfg, {"module": "M1-admission", "slide": args.slide}),
    }
    save_json(P.interim / f"{args.slide}.admission.json", record)

    if failed and not args.force:
        sys.exit(f"\n[M1] 该切片未通过 {', '.join(failed)}，已记录但不进入主分析。"
                 f"\n     队列级放宽请用 --cohort（写进 config，可审计）；"
                 f"\n     确需单张放行请加 --force --force-reason \"...\"。")
    if failed and args.force:
        if not args.force_reason:
            sys.exit("\n[M1] --force 必须同时给 --force-reason，"
                     "否则准入记录里会留下一条没有理由的放行。")
        print(f"\n[M1] ⚠ 强制放行（未通过 {', '.join(failed)}）。"
              f"理由已记入 admission.json：{args.force_reason}")
        print("     Methods 必须如实写明这一条，不能写成'全部通过准入'。")

    # --- 质控 ---
    q = cfg["qc"]
    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_counts=q["min_counts"])
    sc.pp.filter_genes(adata, min_cells=q["min_cells_per_gene"])
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=False)
    print(f"\n[M1] spot {n0} -> {adata.n_obs}，基因 {adata.n_vars}")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 坐标统一为微米。这一步做错会让所有"微米"参数失去意义
    if "spatial" in adata.obsm:
        xy = np.asarray(adata.obsm["spatial"], float)
        # 用最近邻间距标定：把观测到的中位最近邻距离对齐到平台标称间距
        from scipy.spatial import cKDTree
        d = cKDTree(xy).query(xy, k=2)[0][:, 1]
        scale = q["spacing_um"] / max(float(np.median(d)), 1e-9)
        adata.obsm["spatial_um"] = xy * scale
        print(f"[M1] 坐标换算：中位最近邻间距 {np.median(d):.2f} -> {q['spacing_um']} μm "
              f"(scale={scale:.4f})")

    adata.uns["sparta_admission"] = {k: bool(v[0]) for k, v in adm.items()}
    adata.uns["sparta_run"] = stamp_run(cfg, {"module": "M1", "slide": args.slide,
                                              "platform": args.platform})
    out = P.qc(args.slide)
    adata.write(out)
    print(f"[M1] 已写出 {out}")


if __name__ == "__main__":
    main()
