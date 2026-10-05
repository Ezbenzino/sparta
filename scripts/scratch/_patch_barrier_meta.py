"""一次性补丁：把 degraded 标记写进各切片的 barrier_meta.json。
run_04 旧版没写这个字段，现在基于 scored.h5ad 的 obs.columns 静态补。"""
from __future__ import annotations
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import warnings; warnings.filterwarnings("ignore")
import scanpy as sc

_N_COL = ["ECM_core_n", "CAF_n", "ECM_crosslink_n", "Ag_target_n",
          "Hypoxia_n", "Proliferation_n", "Efflux_n"]
_N_TO_KEY = {c: k for k, c in {
    "ecm": "ECM_core_n", "caf": "CAF_n", "crosslink": "ECM_crosslink_n",
    "ag_target": "Ag_target_n", "hypoxia": "Hypoxia_n",
    "proliferation": "Proliferation_n", "efflux": "Efflux_n",
}.items()}

interim = ROOT / "data" / "interim"
for meta_path in sorted(interim.glob("*.barrier_meta.json")):
    sid = meta_path.name.replace(".barrier_meta.json", "")
    scored = interim / f"{sid}.scored.h5ad"
    if not scored.exists():
        continue
    a = sc.read_h5ad(scored)
    cols = set(a.obs.columns)
    filled = [_N_TO_KEY[c] for c in _N_COL if c not in cols]
    del a
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    meta["filled_keys"] = filled
    meta["degraded"] = bool(filled)
    meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"{sid:<8} degraded={bool(filled)}  filled={filled}")
