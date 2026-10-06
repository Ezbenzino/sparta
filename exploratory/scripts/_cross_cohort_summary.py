"""跨队列综合对比：MEL（转移黑色素瘤）vs CSCC（原发皮肤鳞癌）核心证据汇总。"""
import json, numpy as np

def load(p):
    return json.load(open(rf"D:\sparta\results\validation\{p}", encoding="utf-8"))

# 1) 决策点
scr = load("screen_decision.json")
print("="*78)
print("【1】决策点：解离潜力筛查")
print("="*78)
for sid, r in scr["per_slide"].items():
    print(f"  {sid}: potential={r['dissociation_potential']:6.2f}  "
          f"ECM={r['frac_shared_ecm']*100:5.1f}%  抗原={r['frac_antigen']*100:5.1f}%  "
          f"交联={r['frac_crosslink']*100:5.1f}%  col={r['crosslink_collinear']}")
print(f"  → 判定 {scr['verdict']}，中位潜力 {scr['median_potential']:.2f}")

# 2) 解耦
dec = load("decoupling.json")
print("\n" + "="*78)
print("【2】解耦验证（B_cell 与 B_mAb 偏相关，控制距血管距离）")
print("="*78)
per = dec.get("per_slide", dec)
for sid, r in (per.items() if isinstance(per, dict) else []):
    print(f"  {sid}: ρ={r.get('rho', r.get('rho_partial','?'))}  p={r.get('p','?')}")

# 3) S1 排布
print("\n" + "="*78)
print("【3】S1 空间重排（B_cell 对排布的敏感性，fixed 模式）")
print("="*78)
for sid in ["MEL01","MEL02","MEL03","MEL04","CSCC01","CSCC02","CSCC03","CSCC04"]:
    try:
        cf = json.load(open(rf"D:\sparta\results\counterfactual\{sid}.json", encoding="utf-8"))
    except FileNotFoundError:
        continue
    s1 = cf.get("s1") or {}
    fixed = s1.get("fixed", {}) if isinstance(s1, dict) else {}
    print(f"  {sid}: z={fixed.get('z','?')!s:>7}  p={fixed.get('p_emp','?')}")

# 4) S2
print("\n" + "="*78)
print("【4】S2 环带断裂（连续缺口 vs 分散移除，32 组合）")
print("="*78)
all_r = []
for sid in ["MEL01","MEL02","MEL03","MEL04","CSCC01","CSCC02","CSCC03","CSCC04"]:
    try:
        cf = json.load(open(rf"D:\sparta\results\counterfactual\{sid}.json", encoding="utf-8"))
    except FileNotFoundError:
        continue
    s2 = cf.get("s2") or {}
    for k, v in s2.get("per_k", {}).items():
        if v.get("ratio_vs_in_cut"): all_r.append(v["ratio_vs_in_cut"])
print(f"  ratio>1 比例: {sum(1 for r in all_r if r>1)}/{len(all_r)}  "
      f"中位 {np.median(all_r):.2f}  范围 {min(all_r):.2f}-{max(all_r):.2f}")

# 5) B_mAb 敏感性
sens = load("bmab_sensitivity.json")
print("\n" + "="*78)
print("【5】B_mAb 参数敏感性（合理参数邻域）")
print("="*78)
for k, v in sens["summary"].items():
    print(f"  {k}: 交联% {v['crosslink_pct_min']:.1f}-{v['crosslink_pct_max']:.1f}  "
          f"潜力 {v['potential_min']:.2f}-{v['potential_max']:.2f}  "
          f"GO {v['n_above_GO']}/{v['n_cells']}")

# 6) 真实 benchmark
bm = load("benchmark_real.json")
print("\n" + "="*78)
print("【6】真实数据 benchmark（B_cell vs 局部统计量对排布敏感性）")
print("="*78)
for sid, r in bm.items():
    print(f"  {sid}: B_cell z={r['z_b_cell']:+.1f}  "
          f"邻域z={r['z_nbr']:+.1f}  Ripley z={r['z_ripley']:+.1f}  "
          f"ρ(B_cell~CAF)={r['rho_Bcell_vs_CAF']:+.2f}")
