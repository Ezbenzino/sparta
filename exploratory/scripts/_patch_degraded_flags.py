#!/usr/bin/env python
"""一次性补丁：把"哪些切片用了中性 0.5 填充"补进已有验证 JSON。
=====================================================================
背景：scores_from_adata 缺列时只打 warning，不写进结果 JSON。事后看
spatial_null_check.json 时，CSCC14/15/16 的 p=0.156/0.497 会被误读成
生物学阴性，实际是 Ag_target/Efflux 基因在第一代 ST 里没测到、被填成
常数 0.5 导致 B_mAb/B_meta 退化。

本脚本只读各切片 scored.h5ad 的 obs.columns（秒级），不重跑任何计算，
把 filled_keys / degraded 字段补进：
  - results/validation/spatial_null_check.json
  - results/validation/s2_matched_selection.json
未来重跑 run_27/run_28 时代码已自动写这些字段，本脚本不需要再跑。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import warnings  # noqa: E402
warnings.filterwarnings("ignore")
import scanpy as sc  # noqa: E402

# 与 sparta.barrier._SCORE_KEYS 一致
_N_COL = ["ECM_core_n", "CAF_n", "ECM_crosslink_n", "Ag_target_n",
          "Hypoxia_n", "Proliferation_n", "Efflux_n"]
_N_TO_KEY = {c: k for k, c in {
    "ecm": "ECM_core_n", "caf": "CAF_n", "crosslink": "ECM_crosslink_n",
    "ag_target": "Ag_target_n", "hypoxia": "Hypoxia_n",
    "proliferation": "Proliferation_n", "efflux": "Efflux_n",
}.items()}


def filled_for(sid: str) -> list[str]:
    p = ROOT / "data" / "interim" / f"{sid}.scored.h5ad"
    if not p.exists():
        return []
    a = sc.read_h5ad(p)
    cols = set(a.obs.columns)
    missing = [_N_TO_KEY[c] for c in _N_COL if c not in cols]
    del a
    return missing


def patch(path: Path, *, slide_field: str = "slide") -> None:
    if not path.exists():
        print(f"SKIP {path} (not found)")
        return
    obj = json.loads(path.read_text(encoding="utf-8"))
    # decoupling.json 的结构是顶层直接是 slide id（无 per_slide 包装）
    if "per_slide" in obj:
        per = obj["per_slide"]
    else:
        per = obj
    n_deg = 0
    deg_map = {}
    for sid, rec in per.items():
        if not isinstance(rec, dict):
            continue
        filled = filled_for(sid)
        rec["filled_keys"] = filled
        rec["degraded"] = bool(filled)
        if filled:
            n_deg += 1
            deg_map[sid] = filled
    if "per_slide" in obj:
        s = obj.setdefault("summary", {})
        s["n_degraded"] = n_deg
        s["degraded_slides"] = deg_map
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Patched {path.name}: {n_deg}/{len(per)} degraded")
    for sid, keys in deg_map.items():
        print(f"    {sid}: {keys}")


if __name__ == "__main__":
    v = ROOT / "results" / "validation"
    for name in ["spatial_null_check.json", "s2_matched_selection.json",
                  "decoupling.json", "biological_validation.json",
                  "null_crosslink_check.json"]:
        patch(v / name)
