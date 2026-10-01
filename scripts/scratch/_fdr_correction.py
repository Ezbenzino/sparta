# 统计严谨性：S1/S2 置换检验 p 值 BH-FDR 校正（8 张双队列）
import json, sys
import numpy as np
sys.path.insert(0, r"D:\sparta")
from sparta.io_ import Paths, load_config

P = Paths(load_config(None))
slides = ["MEL01", "MEL02", "MEL03", "MEL04", "CSCC01", "CSCC02", "CSCC03", "CSCC04"]


def bh(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    # BH: adj_i = min(1, p_i * n / i)，再从大到小后缀累计最小
    adj_ranked = np.minimum(1.0, ranked * n / np.arange(1, n + 1))
    adj_ranked = np.minimum.accumulate(adj_ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = adj_ranked
    return adj


# ---- 收集 S1 ----
s1_rows = []
for sid in slides:
    cf = json.load(open(P.counterfactual(sid), encoding="utf-8"))
    for mode, r in cf["s1"].items():
        if "p_emp" in r and np.isfinite(r["p_emp"]):
            s1_rows.append(dict(slide=sid, mode=mode, z=r["z"], p=r["p_emp"], b=r["b_real"]))
s1_p = [r["p"] for r in s1_rows]
s1_adj = bh(s1_p)
for r, a in zip(s1_rows, s1_adj):
    r["p_fdr"] = float(a)
print(f"S1: {len(s1_rows)} 个检验")
for r in s1_rows:
    sig = "***" if r["p_fdr"] < 0.001 else ("**" if r["p_fdr"] < 0.01 else ("*" if r["p_fdr"] < 0.05 else ""))
    print(f"  {r['slide']:<7} mode={r['mode']:<7} z={r['z']:+6.1f}  p={r['p']:.4f}  p_fdr={r['p_fdr']:.4f}  {sig}")

# ---- 收集 S2 ----
s2_rows = []
for sid in slides:
    cf = json.load(open(P.counterfactual(sid), encoding="utf-8"))
    for k, r in cf["s2"]["per_k"].items():
        if "p_vs_in_cut" in r and np.isfinite(r["p_vs_in_cut"]):
            s2_rows.append(dict(slide=sid, k=int(k), ratio=r["ratio_vs_in_cut"],
                                p=r["p_vs_in_cut"]))
s2_p = [r["p"] for r in s2_rows]
s2_adj = bh(s2_p)
for r, a in zip(s2_rows, s2_adj):
    r["p_fdr"] = float(a)
print(f"\nS2: {len(s2_rows)} 个检验 (continuous gap vs in-barrier diffuse)")
n_sig = sum(1 for r in s2_rows if r["p_fdr"] < 0.05)
print(f"  校正后显著 (p_fdr<0.05): {n_sig}/{len(s2_rows)}")
print(f"  中位 ratio={np.median([r['ratio'] for r in s2_rows]):.3f}")
for r in s2_rows:
    sig = "*" if r["p_fdr"] < 0.05 else ""
    if r["p_fdr"] >= 0.05:
        print(f"  {r['slide']:<7} k={r['k']:>2}  ratio={r['ratio']:.3f}  p={r['p']:.3f}  p_fdr={r['p_fdr']:.3f}  {sig}")

out = dict(s1=s1_rows, s2=s2_rows,
           summary=dict(
               s1_n=len(s1_rows),
               s1_sig_fdr05=sum(1 for r in s1_rows if r["p_fdr"] < 0.05),
               s2_n=len(s2_rows),
               s2_sig_fdr05=n_sig,
               s2_median_ratio=float(np.median([r["ratio"] for r in s2_rows]))))
json.dump(out, open(r"D:\sparta\results\validation\fdr_correction.json", "w"),
          ensure_ascii=False, indent=1, default=float)
print(f"\nsaved fdr_correction.json")
