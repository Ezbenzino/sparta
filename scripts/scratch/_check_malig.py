import warnings; warnings.filterwarnings("ignore")
import scanpy as sc
import numpy as np
for s in ["BRCA01", "BRCA02", "LN01"]:
    a = sc.read_h5ad(rf"D:\sparta\data\interim\{s}.scored.h5ad")
    m = a.obs["Malignant_n"]
    # 看 Malignant_n 和 EPCAM/KRT8 的关系——应该高相关
    print(f"{s}: Malignant_n mean={m.mean():.3f} std={m.std():.3f} "
          f"q90={m.quantile(0.9):.3f} q10={m.quantile(0.1):.3f}")
    del a
