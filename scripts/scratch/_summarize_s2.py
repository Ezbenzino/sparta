"""S2 环带断裂：跨队列汇总 + 屏障带厚度分析，决定其论文定位。"""
import json, glob, numpy as np

slides = ["MEL01","MEL02","MEL03","MEL04","CSCC01","CSCC02","CSCC03","CSCC04"]
rows = []
for sid in slides:
    p = rf"D:\sparta\results\counterfactual\{sid}.json"
    try:
        d = json.load(open(p, encoding="utf-8"))
    except FileNotFoundError:
        print(f"{sid}: 缺产物"); continue
    s2 = d.get("s2") or d.get("S2") or {}
    if not s2:
        print(f"{sid}: 无 S2 段"); continue
    for k, v in s2.get("per_k", {}).items():
        rows.append((sid, int(k), v.get("ratio_vs_in_cut"), v.get("p_vs_in_cut"),
                     v.get("drop_targeted"), v.get("n_cut_nodes") if False else None))
        # 打印 n_cut_nodes 一次
    nc = s2.get("n_cut_nodes")
    rows[-1] = (sid, int(k), v.get("ratio_vs_in_cut"), v.get("p_vs_in_cut"), nc)

print(f"{'切片':<7}{'k':>4}{'连续缺口后降幅':>14}{'ratio(分散/连续)':>18}{'p':>8}")
print("-" * 60)
for sid in slides:
    d = json.load(open(rf"D:\sparta\results\counterfactual\{sid}.json", encoding="utf-8"))
    s2 = d.get("s2") or d.get("S2") or {}
    nc = s2.get("n_cut_nodes")
    print(f"{sid:<7}  (最小割 {nc} 节点)")
    for k, v in s2.get("per_k", {}).items():
        ratio = v.get("ratio_vs_in_cut")
        p = v.get("p_vs_in_cut")
        drop_t = v.get("drop_targeted", 0)
        sig = "***" if p < 0.01 else ("**" if p < 0.05 else ("*" if p < 0.1 else ""))
        print(f"{'':7}{k:>4}{drop_t*100:>13.1f}%{ratio:>18.2f}{p:>8.3f} {sig}")

# 汇总
print("\n==== 汇总 ====")
all_r = [r[2] for r in rows if r[2] is not None]
sig_r = [r for r in rows if r[2] is not None and r[3] is not None and r[3] < 0.05]
print(f"有效组合 {len(all_r)} 个，其中 ratio_vs_in_cut>1: {sum(1 for r in all_r if r>1)}")
print(f"显著(p<0.05): {len(sig_r)} 个，中位 ratio = {np.median(all_r):.2f}")
print(f"ratio 范围: {min(all_r):.2f}–{max(all_r):.2f}")
