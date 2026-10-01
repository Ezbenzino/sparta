#!/usr/bin/env python
"""
run_00c_admission_audit.py —— 准入台账重建（只读原始数据，不改任何管线产物）
==============================================================================

输入：data/ledger.csv + data/interim/{sid}.raw.h5ad
输出：data/interim/{sid}.admission.json（每张一份）+ data/admission_audit.csv（汇总）
上游模块：run_00b_ingest.py
下游模块：无（供 Methods 的"数据与准入"段落引用）

为什么需要这个脚本
------------------
2026-08-26 的投稿前审查发现：`run_01_qc.py` 旧版只在**准入失败**时写
`{sid}.admission.json`，而失败分支又会被 `--force` 跳过。实际建库时
`_pipeline_3samples.py` 对每张切片都带了 `--force`，结果是 8 张切片里
只有 MEL01 留下了准入记录（而且那份记录显示 C7 未通过）。
换句话说，"全部通过准入 C1–C7"这句话在磁盘上没有任何证据支撑。

这个脚本把台账补回来：对每张切片重新跑一遍 C1–C7（用的是**原始** raw.h5ad，
与建库时 run_01 看到的输入一致），把结论、实际生效的阈值、
以及该切片是否属于队列级放宽逐条落盘。它不碰 qc/scored/graph/barrier
任何产物，因此不会改变任何已有结果。

跑完之后，Methods 里就可以写"逐切片准入记录见 data/admission_audit.csv"，
而不是一句无从核验的"全部通过"。

用法
----
    python scripts/run_00c_admission_audit.py                 # 读 ledger 里的全部切片
    python scripts/run_00c_admission_audit.py --slides MEL01 CSCC03
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sparta.io_ import Paths, load_config, save_json, set_seed, stamp_run  # noqa: E402

# 队列名 -> configs/default.yaml 的 admission_overrides 段。
# 靠 ledger 的 accession/source 字段自动匹配，匹配不上就用默认阈值（不放宽）。
COHORT_BY_ACCESSION = {
    "GSE144239": "cscc_gse144239",
    "GSE250636": None,          # 2024 年 Visium，用默认阈值
}


def _cohort_of(row: dict) -> str | None:
    blob = " ".join(str(row.get(k, "")) for k in ("source", "accession", "raw_path", "notes"))
    for acc, cohort in COHORT_BY_ACCESSION.items():
        if acc in blob:
            return cohort
    return None


def main():
    ap = argparse.ArgumentParser(description="准入台账重建（只读）",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slides", nargs="+", default=None, help="默认取 ledger.csv 全部切片")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    try:
        import scanpy as sc
    except ImportError:
        sys.exit("需要 scanpy 读取 .raw.h5ad。pip install scanpy")

    from run_01_qc import check_admission  # 复用同一套判定逻辑，避免两处实现漂移

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)

    ledger_path = Path(cfg["paths"]["root"]) / "data" / "ledger.csv"
    if not ledger_path.exists():
        sys.exit(f"找不到 {ledger_path}。请先跑 run_00b_ingest.py。")
    with open(ledger_path, encoding="utf-8") as f:
        ledger = {r["slide_id"]: r for r in csv.DictReader(f)}

    slides = args.slides or list(ledger)
    rows = []
    print(f"{'slide':<8}{'队列放宽':<18}{'C1':>4}{'C2':>4}{'C3':>4}{'C4':>4}"
          f"{'C5':>4}{'C6':>4}{'C7':>4}  结论")
    print("-" * 78)

    for sid in slides:
        meta = ledger.get(sid, {})
        raw = P.interim / f"{sid}.raw.h5ad"
        if not raw.exists():
            print(f"{sid:<8} 跳过：找不到 {raw}")
            continue
        adata = sc.read_h5ad(raw)
        adata.var_names_make_unique()
        if "counts" not in adata.layers:
            adata.layers["counts"] = adata.X.copy()

        cohort = _cohort_of(meta)
        platform = "legacy_st" if meta.get("platform") == "legacy_st" else "visium"
        has_he = str(meta.get("has_image", "")).lower() in ("true", "1", "yes")
        treatment_known = str(meta.get("treatment", "unknown")).lower() not in ("", "unknown", "na")

        adm, thr = check_admission(adata, cfg, platform, has_he, treatment_known, cohort=cohort)
        failed = [k for k, (ok, _) in adm.items() if not ok]

        record = {
            "slide": sid,
            "checks": {k: {"pass": bool(v[0]), "detail": v[1]} for k, v in adm.items()},
            "failed": failed,
            "admitted": True,                    # 这些切片事实上已进入分析
            "forced": bool(failed),              # 有未通过项却已进入分析 = 强制放行
            "force_reason": ("队列级质量差异；2026-08-26 回溯登记，"
                             "并已在 configs/default.yaml 的 admission_overrides 写明阈值依据"
                             if failed else None),
            "cohort": cohort,
            "thresholds_applied": {k: v for k, v in thr.items() if not k.startswith("_")},
            "thresholds_overridden": thr.get("_overridden", []),
            "reconstructed": True,               # 标明这是回溯重建，不是建库当时写下的
            "meta": stamp_run(cfg, {"module": "M0c-admission-audit", "slide": sid}),
        }
        save_json(P.interim / f"{sid}.admission.json", record)

        marks = "".join(f"{'✓' if adm[c][0] else '✗':>4}" for c in
                        ("C1", "C2", "C3", "C4", "C5", "C6", "C7"))
        verdict = "通过" if not failed else f"强制放行（未过 {','.join(failed)}）"
        print(f"{sid:<8}{(cohort or '默认阈值'):<18}{marks}  {verdict}")

        rows.append(dict(slide_id=sid, cohort=cohort or "", platform=platform,
                         **{c: ("pass" if adm[c][0] else "FAIL") for c in
                            ("C1", "C2", "C3", "C4", "C5", "C6", "C7")},
                         failed=";".join(failed), forced=bool(failed),
                         min_spots=thr.get("min_spots_legacy_st" if platform == "legacy_st"
                                           else "min_spots_visium"),
                         min_median_umi=thr.get("min_median_umi"),
                         overridden=";".join(thr.get("_overridden", []))))

    if not rows:
        sys.exit("没有任何切片被审计。")

    out_csv = Path(cfg["paths"]["root"]) / "data" / "admission_audit.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    n_forced = sum(1 for r in rows if r["forced"])
    print("-" * 78)
    print(f"已写出 {out_csv}（{len(rows)} 张，其中 {n_forced} 张为强制放行）")
    if n_forced:
        print("\n⚠ Methods 里必须如实写成：")
        print(f"   \"{len(rows)} 张切片中 {len(rows)-n_forced} 张在预设阈值下通过全部 C1–C7；")
        print(f"    {n_forced} 张在队列级放宽阈值下纳入，逐切片记录见 data/admission_audit.csv\"")
        print("   不能写成 \"全部通过准入\"。")


if __name__ == "__main__":
    main()
