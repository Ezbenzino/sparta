"""更新 ext_validation.json 的 caveat——sink 已用 brca 标记重建。"""
import json
from pathlib import Path
import numpy as np

V = Path(r"D:\sparta\results\validation")
p = V / "ext_validation.json"
e = json.load(open(p, encoding="utf-8"))

# BH 校正（外部 3 个 p 值）
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

e["summary"]["caveat"] = (
    "External slides re-scored with tumor_type=brca (EPCAM/KRT8/KRT18/KRT19/MUC1) on "
    "2026-10-03, after fixing run_batch's silent fallback to melanoma. Sink (tumor core) "
    "is now defined on epithelial-malignant markers. Limitations: (1) BRCA01 and BRCA02 are "
    "two adjacent sections from the SAME patient (BRCA_P1), not independent tumor replicates; "
    "(2) LN01 is a non-tumor lymph node, included as a negative control — its S2 ratio is "
    "near 1.0 (1.016x, p=0.12), consistent with no tumor-core barrier; (3) BRCA02 null "
    "p=0.050 is borderline after BH adjustment (adj=0.050), treat as sensitive; "
    "(4) epithelial-malignant markers have not been histopathologically confirmed, so "
    "'tumor core' is a model-defined epithelial region, not a pathologist-verified region. "
    "These data support cross-dataset reproduction of model-field coupling, not independent "
    "biological validation of a tumor barrier."
)
e["summary"]["n_null_sig_bh"] = int((ext_adj < 0.05).sum())
e["summary"]["median_rho_partial"] = float(np.median([r["decoupling"]["rho_partial"] for r in e["per_slide"].values()]))

json.dump(e, open(p, "w", encoding="utf-8"), indent=2, ensure_ascii=False)
print("ext_validation.json caveat updated.")
for s, r in e["per_slide"].items():
    print(f"  {s}: rho_partial={r['decoupling']['rho_partial']:.3f} "
          f"null_p={r['spatial_null']['empirical_p_one_sided']:.3f} "
          f"bh={r['spatial_null']['bh_padj']:.3f}")
