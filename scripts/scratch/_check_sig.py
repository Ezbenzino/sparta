# 验证 ECM_crosslink 与 ECM_core 签名在 MEL01 中的基因匹配与空间分布
import sys
sys.path.insert(0, r"D:\sparta")
import numpy as np
import anndata as ad

a = ad.read_h5ad(r"D:\sparta\data\interim\MEL01.scored.h5ad")
genes = set(a.var_names)

crosslink = ["LOX", "LOXL1", "LOXL2", "LOXL3", "PLOD1", "PLOD2", "TGM2"]
ecm_core = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2",
            "FN1", "LAMB1", "TNC", "THBS2", "FBN1", "VCAN", "BGN", "LUM", "DCN"]

print("ECM_crosslink 基因匹配:")
for g in crosslink:
    print(f"  {g}: {'√' if g in genes else '×'}")
print("ECM_core 基因匹配:")
for g in ecm_core:
    print(f"  {g}: {'√' if g in genes else '×'}")

# 检查秩标准化后的签名分布（应有空间异质性）
for col in ["ECM_core_n", "ECM_crosslink_n", "Ag_target_n"]:
    if col in a.obs:
        v = a.obs[col].values
        print(f"{col}: min={v.min():.3f} max={v.max():.3f} std={v.std():.3f} "
              f"IQR={np.percentile(v,75)-np.percentile(v,25):.3f}")
    else:
        print(f"{col}: 缺失")

# crosslink 与 ecm 的空间相关
if "ECM_core_n" in a.obs and "ECM_crosslink_n" in a.obs:
    from scipy.stats import spearmanr
    r, p = spearmanr(a.obs["ECM_core_n"], a.obs["ECM_crosslink_n"])
    print(f"ECM_core_n vs ECM_crosslink_n: Spearman ρ={r:.3f} (p={p:.2e})")
