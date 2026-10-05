"""把 BH 校正写进 spatial_null_check.json，并在 ext_validation 里加 sink 局限说明。"""
import json
from pathlib import Path
import numpy as np

V = Path(r"D:\sparta\results\validation")

# --- spatial_null_check.json: 加 BH 校正 ---
p = V / "spatial_null_check.json"
o = json.load(open(p, encoding="utf-8"))
slides = sorted(o["per_slide"].keys())
ps = np.array([o["per_slide"][s]["empirical_p_one_sided"] for s in slides])
m = len(ps)
order = np.argsort(ps)
adj = np.empty(m)
prev = 1.0
for rank, idx in enumerate(order[::-1], start=1):
    i = m - rank
    val = min(prev, ps[idx] * m / (i + 1))
    adj[idx] = min(val, 1.0)
    prev = val
for s, a in zip(slides, adj):
    o["per_slide"][s]["bh_padj"] = float(a)
o["summary"]["n_sig_raw_p05"] = int((ps < 0.05).sum())
o["summary"]["n_sig_bh_p05"] = int((adj < 0.05).sum())
o["summary"]["bh_note"] = (
    f"Raw p<0.05: {int((ps<0.05).sum())}/{m}; BH-FDR adjusted p<0.05: {int((adj<0.05).sum())}/{m}. "
    f"CSCC10 (raw p=0.048) is not BH-significant (adj=0.061); it is also ag_target-degraded."
)
json.dump(o, open(p, "w", encoding="utf-8"), indent=2, ensure_ascii=False)
print(f"spatial_null_check.json: raw {int((ps<0.05).sum())}/{m} -> BH {int((adj<0.05).sum())}/{m}")

# --- ext_validation.json: 加 sink 局限 note ---
p2 = V / "ext_validation.json"
e = json.load(open(p2, encoding="utf-8"))
e["summary"]["caveat"] = (
    "External slides were scored with tumor_type=melanoma (run_batch fallback bug, "
    "since fixed). R3 spot-level fields (b_cell_field from source only, b_mab from vessel) "
    "do NOT depend on sink/Malignant, so rho_partial values are valid. However, the "
    "slice-level B_cell min-cut and S1/S2 counterfactuals DO depend on sink, which was "
    "defined using melanoma markers (MLANA/PMEL/TYR...) that are not expressed in breast "
    "tissue or lymph node; those slice-level external results should be treated as "
    "exploratory pending re-scoring with tumor_type=brca."
)
# 外部 3 个 p 值的 BH（小样本）
ext_ps = [r["spatial_null"]["empirical_p_one_sided"] for r in e["per_slide"].values()]
order = np.argsort(ext_ps)
ext_adj = np.empty(len(ext_ps))
prev = 1.0
for rank, idx in enumerate(order[::-1], start=1):
    i = len(ext_ps) - rank
    ext_adj[idx] = min(prev, ext_ps[idx] * len(ext_ps) / (i + 1))
    prev = ext_adj[idx]
for (s, r), a in zip(e["per_slide"].items(), ext_adj):
    r["spatial_null"]["bh_padj"] = float(a)
e["summary"]["n_null_sig_bh"] = int((ext_adj < 0.05).sum())
json.dump(e, open(p2, "w", encoding="utf-8"), indent=2, ensure_ascii=False)
print(f"ext_validation.json: BH {int((ext_adj<0.05).sum())}/3 (BRCA02 adj={ext_adj[1]:.4f})")
