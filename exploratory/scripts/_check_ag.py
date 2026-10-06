# 临时诊断脚本：检查靶抗原与外排泵基因在数据中的存在性
import sys
sys.path.insert(0, r"D:\sparta")
import anndata as ad

a = ad.read_h5ad(r"D:\sparta\data\interim\MEL01.scored.h5ad")
genes = list(a.var_names)
print("n_genes:", len(genes))

for g in ["CD274", "PDCD1LG2", "CD80", "CD86", "ABCB1", "ABCG2",
          "ABCC1", "MDR1", "ERBB2", "TACSTD2", "PECAM1", "VWF", "CDH5"]:
    hit = [x for x in genes if x == g]
    print(g, "存在" if hit else "缺失")

print("var_names 前 25 个:", genes[:25])
sig_cols = [c for c in a.obs.columns if c.endswith("_n")]
print("obs 中秩标准化签名列:", sig_cols)
print("obs 中 Ag_target 列:", [c for c in a.obs.columns if "Ag" in c])
print("obs 中 Efflux 列:", [c for c in a.obs.columns if "fflux" in c])
