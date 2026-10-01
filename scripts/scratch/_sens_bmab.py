# B_mAb 参数敏感性网格：判定"交联主导"是真实发现还是参数假象
import sys, json
sys.path.insert(0, r"D:\sparta")
import numpy as np
import scanpy as sc

from sparta.io_ import Paths, load_config, load_graph
from sparta.barrier import scores_from_adata
from sparta.validate import dissociation_drivers

cfg = load_config(None)
P = Paths(cfg)
SID = "MEL01"

adata = sc.read_h5ad(P.scored(SID))
A, D, source, sink, vessel, _ = load_graph(P.graph(SID))
S = scores_from_adata(adata)

cfg_cell = cfg["barrier"]["b_cell"]
base_mab = {k: v for k, v in cfg["barrier"]["b_mab"].items() if k != "r_nm"}
r_nm = cfg["barrier"]["b_mab"]["r_nm"]

print(f"基线 B_mAb 参数: {base_mab}")
print("=" * 100)

def run(mab_cfg, tag):
    dd = dissociation_drivers(A, S["ecm"], S["caf"], S["crosslink"], S["ag_target"],
                              source, vessel, cfg_cell=cfg_cell, cfg_mab=mab_cfg, r_nm=r_nm)
    pot = dd["dissociation_potential"]
    return pot, dd["frac_shared_ecm"], dd["frac_antigen"], dd["frac_crosslink"], dd["crosslink_collinear"]

def report(tag, mab_cfg):
    pot, fe, fa, fc, col = run(mab_cfg, tag)
    flag = "共线!" if col else ""
    print(f"{tag:<34} 潜力={pot:>7.2f}  ECM%={fe*100:>6.1f}  "
          f"抗原%={fa*100:>6.1f}  交联%={fc*100:>6.1f}  {flag}")

print("\n[1] 单参数扫描（其余取默认）\n" + "-" * 100)
print(f"{'参数组合':<34}{'潜力':>9}{'ECM%':>8}{'抗原%':>8}{'交联%':>8}")
print("-" * 100)

report("默认 (lam=3, beta=3, xi0=20)", base_mab)

# lam：ECM 衰减
for v in [1.0, 8.0, 20.0]:
    m = dict(base_mab); m["lam"] = v
    report(f"lam={v}", m)
# beta：交联收缩
for v in [0.5, 1.5, 6.0, 12.0]:
    m = dict(base_mab); m["beta"] = v
    report(f"beta={v}", m)
# xi0：基线网孔
for v in [10.0, 40.0, 80.0]:
    m = dict(base_mab); m["xi0_nm"] = v
    report(f"xi0={v}", m)
# kappa_w：结合位点权重
for v in [0.1, 0.5, 3.0]:
    m = dict(base_mab); m["kappa_w"] = v
    report(f"kappa_w={v}", m)
# kd_eff：解离常数
for v in [0.1, 2.0]:
    m = dict(base_mab); m["kd_eff"] = v
    report(f"kd_eff={v}", m)

print("\n[2] 联合扫描：把交联推到最不利/最有利的极端\n" + "-" * 100)
# 最不利：lam 大（ECM 也强）、beta 小（交联弱）
m = dict(base_mab); m["lam"] = 20.0; m["beta"] = 0.5
report("lam=20, beta=0.5 (ECM强,交联弱)", m)
# 最有利：lam 小、beta 大
m = dict(base_mab); m["lam"] = 1.0; m["beta"] = 12.0
report("lam=1, beta=12 (ECM弱,交联强)", m)
# 网孔超大：交联几乎不影响网孔
m = dict(base_mab); m["xi0_nm"] = 80.0; m["beta"] = 0.5
report("xi0=80, beta=0.5 (交联几乎失效)", m)
# 网孔很小：尺寸排阻普遍强
m = dict(base_mab); m["xi0_nm"] = 10.0; m["beta"] = 3.0
report("xi0=10 (网孔普遍小)", m)
# 抗原权重拉到最大
m = dict(base_mab); m["kappa_w"] = 10.0; m["kd_eff"] = 0.05
report("kappa_w=10, kd=0.05 (抗原主导)", m)

print("\n[3] 结论判据：")
print("  若在所有合理参数下交联%始终 >50% 且潜力>1 -> 交联主导是稳健发现")
print("  若在某些参数下 ECM%/抗原%反超 -> 交联主导是参数产物，需谨慎")
