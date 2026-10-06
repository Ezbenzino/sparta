"""维度③概念验证（2/2）：GSE91061 屏障/免疫签名 vs ICB 响应。

用与空间分析完全相同的签名基因集（ECM_crosslink/CAF/Hypoxia = B_mAb 组成；
CD8T/T_NK/Ag_target = 免疫激活对照），在治疗前 bulk RNA-seq 上打分，
检验"屏障成分高 → 非响应"假说。
"""
import gzip, json, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
sys.path.insert(0, r"D:\sparta")
from sparta.signatures import GENE_SETS

META = json.load(open(r"D:\sparta\data\external\GSE91061_meta.json", encoding="utf-8"))
S2E = json.load(open(r"D:\sparta\data\external\GSE91061_sym2entrez.json", encoding="utf-8"))

# ---- 读 rld ----
print("reading rld ...")
df = pd.read_csv(r"D:\sparta\data\external\GSE91061_rld.csv.gz", index_col=0)
df.columns = [c.strip().strip('"') for c in df.columns]
df.index = df.index.astype(str).str.strip().str.strip('"')
print("matrix:", df.shape)

# 只保留 Pre-treatment 且 response 已知
pre = {t: m for t, m in META.items() if "_Pre" in t and m["response"] != "UNK"}
print(f"pre-treatment with known response: {len(pre)}")
print("resp:", {r: sum(1 for m in pre.values() if m['response']==r) for r in set(m['response'] for m in pre.values())})

# ---- 签名分数（基因级 z-score 均值）----
def sig_score(sample_expr: pd.Series, syms) -> float:
    ezs = [S2E[s] for s in syms if s in S2E and str(S2E[s]) in df.index]
    if not ezs:
        return np.nan
    sub = df.loc[[str(e) for e in ezs], sample_expr.name]  # 列=样本
    # 跨样本 z（在全部样本上标准化该基因）
    gz = (sub - df.loc[[str(e) for e in ezs]].mean(axis=1)) / df.loc[[str(e) for e in ezs]].std(axis=1)
    return gz.mean()

# 构建基因集符号（用 melanoma 的 CAF/ECM_crosslink；hypoxia 从 gmt）
S = GENE_SETS["melanoma"]
hypoxia = []
for line in open(r"D:\sparta\data\external\h.all.v7.1.HALLMARK_HYPOXIA.symbols.gmt", encoding="utf-8"):
    hypoxia = line.rstrip("\n").split("\t")[2:]

sigdefs = {
    "B_mAb_crosslink (ECM_crosslink)": S["ECM_crosslink"],
    "B_mAb_CAF": S["CAF"],
    "B_mAb_hypoxia (HALLMARK)": hypoxia,
    "B_cell_ECM_core": S["ECM_core"],
    "immune_CD8T": S["CD8T"],
    "immune_T_NK": S["T_NK"],
    "immune_Ag_target (PD-L1/2)": S["Ag_target"],
    "Efflux (drug resistance)": S["Efflux"],
}

# 全部样本的基因级 z 矩阵（显式 axis=0 按行标准化）
zdf = df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1), axis=0)
scores = {}
for name, syms in sigdefs.items():
    ezs = [S2E[s] for s in syms if s in S2E and str(S2E[s]) in zdf.index]
    present = len(ezs) / max(len(syms), 1)
    sub = zdf.loc[[str(e) for e in ezs]]
    scores[name] = sub.mean(axis=0)  # 每样本分数
    print(f"  {name}: {len(ezs)}/{len(syms)} genes ({present*100:.0f}%)")

S_df = pd.DataFrame(scores)

# ---- 组别 ----
pre_cols = [t for t in pre]
resp = np.array([1 if pre[t]["response"] == "PRCR" else 0 for t in pre_cols])
# NR 定义：PD 严格；含 SD 版本
resp_strict = np.array([1 if pre[t]["response"] == "PRCR" else (0 if pre[t]["response"] == "PD" else np.nan) for t in pre_cols])
print("R:", int(resp.sum()), "NR(PD+SD):", int((resp == 0).sum()))
print("R(strict):", int((resp_strict == 1).sum()), "NR(strict PD):", int((resp_strict == 0).sum()))

# ---- 检验：屏障复合 vs 免疫复合 ----
barrier_sig = ["B_mAb_crosslink (ECM_crosslink)", "B_mAb_CAF", "B_mAb_hypoxia (HALLMARK)"]
immune_sig = ["immune_CD8T", "immune_T_NK", "immune_Ag_target (PD-L1/2)"]
S_df["Barrier composite"] = S_df[barrier_sig].mean(axis=1)
S_df["Immune composite"] = S_df[immune_sig].mean(axis=1)

def auc_report(col, y, label):
    mask = ~np.isnan(y)
    x = S_df.loc[pre_cols, col].values[mask]
    yy = y[mask]
    r = np.array(x)[yy == 1]; nr = np.array(x)[yy == 0]
    u = stats.mannwhitneyu(r, nr, alternative="two-sided")
    auc = u[0] / (len(r) * len(nr))
    # AUC 定义：分数高 → R（AUC>0.5 表示 R 分数更高）；对 NR 预测用 1-auc
    print(f"  {label}: R_median={np.median(r):.3f} NR_median={np.median(nr):.3f}  "
          f"AUC(R vs NR)={auc:.3f}  p={u[1]:.2e}  (barrier 高→NR 的 AUC={1-auc:.3f})")
    return auc, u[1]

print("\n=== 全部 Pre 样本 (R=PRCR vs NR=PD+SD) ===")
for c, lab in [("Barrier composite", "Barrier (crosslink+CAF+hypoxia)"),
               ("Immune composite", "Immune (CD8T+TNK+PDL1)")]:
    auc_report(c, resp, lab)
print("\n=== 严格 NR (R=PRCR vs NR=PD, 剔除 SD) ===")
for c, lab in [("Barrier composite", "Barrier (crosslink+CAF+hypoxia)"),
               ("Immune composite", "Immune (CD8T+TNK+PDL1)")]:
    auc_report(c, resp_strict, lab)

# 单项也报告
print("\n=== 单项签名 (全部样本) ===")
for c in sigdefs:
    mask = ~np.isnan(resp)
    x = S_df.loc[pre_cols, c].values[mask]
    yy = resp[mask]
    r = np.array(x)[yy == 1]; nr = np.array(x)[yy == 0]
    u = stats.mannwhitneyu(r, nr, alternative="two-sided")
    auc = u[0] / (len(r) * len(nr))
    print(f"  {c}: AUC(R vs NR)={auc:.3f} p={u[1]:.2e}")

# 保存
S_df.to_csv(r"D:\sparta\results\validation\gse91061_scores.csv")
json.dump({"pre_columns": pre_cols, "response": [int(v) for v in resp],
           "response_strict": [None if np.isnan(v) else int(v) for v in resp_strict]},
          open(r"D:\sparta\results\validation\gse91061_labels.json", "w"))
print("\nsaved gse91061_scores.csv + labels")
