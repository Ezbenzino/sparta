"""维度③概念验证：跨队列 ICB 响应（GSE78220 + GSE91061）。

每个队列独立用相同的签名基因集打分（基因级 z 均值），
比较 Responder vs Non-responder 的 Barrier/Immune 复合分数。
"""
import json, sys
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, r"D:\sparta")
from sparta.signatures import GENE_SETS

S = GENE_SETS["melanoma"]
hypoxia = []
for line in open(r"D:\sparta\data\external\h.all.v7.1.HALLMARK_HYPOXIA.symbols.gmt", encoding="utf-8"):
    hypoxia = line.rstrip("\n").split("\t")[2:]
sigdefs = {
    "crosslink": S["ECM_crosslink"],
    "CAF": S["CAF"],
    "hypoxia": hypoxia,
    "CD8T": S["CD8T"],
    "T_NK": S["T_NK"],
    "PDL1": S["Ag_target"],
}
barrier_sig = ["crosslink", "CAF", "hypoxia"]
immune_sig = ["CD8T", "T_NK", "PDL1"]


def zscore_rows(mat: pd.DataFrame):
    return mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1), axis=0)


def score_matrix(zmat: pd.DataFrame, sigdefs, sym_lookup):
    out = {}
    for name, syms in sigdefs.items():
        idx = [sym_lookup[s] for s in syms if s in sym_lookup and sym_lookup[s] in zmat.index]
        out[name] = zmat.loc[idx].mean(axis=0)
    df = pd.DataFrame(out)
    df["Barrier"] = df[barrier_sig].mean(axis=1)
    df["Immune"] = df[immune_sig].mean(axis=1)
    return df


def test(df_scores, y, label):
    r = df_scores.loc[y == 1]; nr = df_scores.loc[y == 0]
    res = {}
    for col in ["Barrier", "Immune"]:
        xr = r[col].values; xn = nr[col].values
        u = stats.mannwhitneyu(xr, xn, alternative="two-sided")
        auc = u[0] / (len(xr) * len(xn))
        res[col] = dict(R_med=float(np.median(xr)), NR_med=float(np.median(xn)),
                        AUC_R_vs_NR=float(auc), p=float(u[1]))
        print(f"  [{label}] {col}: R={np.median(xr):+.3f} NR={np.median(xn):+.3f} "
              f"AUC(R>NR)={auc:.3f} p={u[1]:.3f}")
    return res


summary = {}

# ---------------- 队列 1: GSE78220 (symbol) ----------------
print("=== GSE78220 (Hugo 2016, anti-PD-1, n=28) ===")
meta = json.load(open(r"D:\sparta\data\external\GSE78220_meta.json", encoding="utf-8"))
df1 = pd.read_excel(r"D:\sparta\data\external\GSE78220_PatientFPKM.xlsx")
df1.columns = [str(c) for c in df1.columns]
df1 = df1.set_index("Gene")
# 列 = PtN.baseline → PtN
df1.columns = [c.split(".")[0] for c in df1.columns]
# 去重列（若有）
df1 = df1.loc[:, ~df1.columns.duplicated()]
# 行去重（symbol）
df1 = df1[~df1.index.duplicated()]
# 只保留有响应的患者
kept = [t for t in df1.columns if t in meta and meta[t]["response"] != "UNK"]
df1 = df1[kept]
y1 = np.array([1 if meta[t]["response"] in ("Complete Response", "Partial Response")
               else 0 for t in kept])
print(f"  n={len(kept)}, R={int(y1.sum())}, NR={int((y1==0).sum())}")
sym1 = {s: s for s in df1.index}  # symbol 直接匹配
z1 = zscore_rows(df1)
s1 = score_matrix(z1, sigdefs, sym1)
summary["GSE78220"] = test(s1, y1, "GSE78220")

# ---------------- 队列 2: GSE91061 (entrez) ----------------
print("=== GSE91061 (Riaz 2017, anti-PD-1, pre-treatment) ===")
meta2 = json.load(open(r"D:\sparta\data\external\GSE91061_meta.json", encoding="utf-8"))
S2E = json.load(open(r"D:\sparta\data\external\GSE91061_sym2entrez.json", encoding="utf-8"))
df2 = pd.read_csv(r"D:\sparta\data\external\GSE91061_rld.csv.gz", index_col=0)
df2.columns = [c.strip().strip('"') for c in df2.columns]
df2.index = df2.index.astype(str).str.strip().str.strip('"')
pre = {t: m for t, m in meta2.items() if "_Pre" in t and m["response"] != "UNK"}
df2 = df2[[t for t in pre]]
# y: R=PRCR, NR=PD+SD
y2 = np.array([1 if pre[t]["response"] == "PRCR" else 0 for t in pre])
print(f"  n={len(pre)}, R={int(y2.sum())}, NR(PD+SD)={int((y2==0).sum())}")
# 对每个 signature 直接构建 symbol→row 映射
def sym_lookup_from(syms):
    return {s: str(S2E[s]) for s in syms if s in S2E and str(S2E[s]) in df2.index}
all_syms = set()
for v in sigdefs.values(): all_syms |= set(v)
sym2 = sym_lookup_from(all_syms)
z2 = zscore_rows(df2)
s2 = score_matrix(z2, sigdefs, sym2)
summary["GSE91061"] = test(s2, y2, "GSE91061(PD+SD)")

# 严格版
y2s = np.array([1 if pre[t]["response"] == "PRCR" else (0 if pre[t]["response"] == "PD" else np.nan)
                for t in pre])
mask = ~np.isnan(y2s)
s2s = s2[mask]; y2ss = y2s[mask].astype(int)
print(f"  [严格] n={int(mask.sum())}, R={int(y2ss.sum())}, NR(PD)={int((y2ss==0).sum())}")
summary["GSE91061_strict"] = test(s2s, y2ss, "GSE91061(PD only)")

json.dump(summary, open(r"D:\sparta\results\validation\icb_crosscohort_summary.json", "w"), indent=1, default=float)
print("\nsaved icb_crosscohort_summary.json")
