"""对 spatial_null_check.json 的 19 个逐切片 p 值做 BH 校正。"""
import json
import numpy as np
from pathlib import Path

V = Path(r"D:\sparta\results\validation")
o = json.load(open(V / "spatial_null_check.json", encoding="utf-8"))
slides = sorted(o["per_slide"].keys())
ps = np.array([o["per_slide"][s]["empirical_p_one_sided"] for s in slides])

# Benjamini-Hochberg
m = len(ps)
order = np.argsort(ps)
adj = np.empty(m)
prev = 1.0
for rank, idx in enumerate(order[::-1], start=1):
    i = m - rank
    val = min(prev, ps[idx] * m / (i + 1))
    adj[idx] = min(val, 1.0)
    prev = val

print(f"{'slide':<8}{'raw_p':>10}{'BH_padj':>12}{'sig_raw':>9}{'sig_bh':>8}")
print("-" * 50)
for s, p, a in zip(slides, ps, adj):
    print(f"{s:<8}{p:>10.4f}{a:>12.4f}{'*' if p<0.05 else ' ':>9}{'*' if a<0.05 else ' ':>8}")

print(f"\nraw p<0.05: {(ps<0.05).sum()}/{m}")
print(f"BH padj<0.05: {(adj<0.05).sum()}/{m}")

# 外部队列
ext = json.load(open(V / "ext_validation.json", encoding="utf-8"))
print("\n--- 外部队列 ---")
for s, r in ext["per_slide"].items():
    p = r["spatial_null"]["empirical_p_one_sided"]
    print(f"{s:<8} raw p={p:.4f}")
ext_ps = [r["spatial_null"]["empirical_p_one_sided"] for r in ext["per_slide"].values()]
ext_adj = np.empty(len(ext_ps))
order = np.argsort(ext_ps)
prev = 1.0
for rank, idx in enumerate(order[::-1], start=1):
    i = len(ext_ps) - rank
    ext_adj[idx] = min(prev, ext_ps[idx] * len(ext_ps) / (i + 1))
print(f"外部 BH padj<0.05: {(ext_adj<0.05).sum()}/{len(ext_ps)}")
for s, a in zip(ext["per_slide"].keys(), ext_adj):
    print(f"  {s}: BH padj={a:.4f}")
