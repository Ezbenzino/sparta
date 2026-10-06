# 维度③跨队列图：屏障/免疫分数 vs ICB 响应（两队列分面板）
import json, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
sys.path.insert(0, r"D:\sparta")
from sparta.signatures import GENE_SETS

S = GENE_SETS["melanoma"]
hypoxia = []
for line in open(r"D:\sparta\data\external\h.all.v7.1.HALLMARK_HYPOXIA.symbols.gmt", encoding="utf-8"):
    hypoxia = line.rstrip("\n").split("\t")[2:]
sigdefs = {"crosslink": S["ECM_crosslink"], "CAF": S["CAF"], "hypoxia": hypoxia,
           "CD8T": S["CD8T"], "T_NK": S["T_NK"], "PDL1": S["Ag_target"]}
barrier_sig = ["crosslink", "CAF", "hypoxia"]
immune_sig = ["CD8T", "T_NK", "PDL1"]

def zscore_rows(mat): return mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1), axis=0)
def score_matrix(zmat, lookup):
    out = {}
    for name, syms in sigdefs.items():
        idx = [lookup[s] for s in syms if s in lookup and lookup[s] in zmat.index]
        out[name] = zmat.loc[idx].mean(axis=0)
    d = pd.DataFrame(out); d["Barrier"] = d[barrier_sig].mean(axis=1); d["Immune"] = d[immune_sig].mean(axis=1)
    return d

# 队列1 GSE78220
meta = json.load(open(r"D:\sparta\data\external\GSE78220_meta.json", encoding="utf-8"))
df1 = pd.read_excel(r"D:\sparta\data\external\GSE78220_PatientFPKM.xlsx").set_index("Gene")
df1.columns = [str(c).split(".")[0] for c in df1.columns]
df1 = df1.loc[:, ~df1.columns.duplicated()]; df1 = df1[~df1.index.duplicated()]
kept1 = [t for t in df1.columns if t in meta]
df1 = df1[kept1]
y1 = np.array([1 if meta[t]["response"] in ("Complete Response", "Partial Response") else 0 for t in kept1])
s1 = score_matrix(zscore_rows(df1), {s: s for s in df1.index})

# 队列2 GSE91061
meta2 = json.load(open(r"D:\sparta\data\external\GSE91061_meta.json", encoding="utf-8"))
S2E = json.load(open(r"D:\sparta\data\external\GSE91061_sym2entrez.json", encoding="utf-8"))
df2 = pd.read_csv(r"D:\sparta\data\external\GSE91061_rld.csv.gz", index_col=0)
df2.columns = [c.strip().strip('"') for c in df2.columns]; df2.index = df2.index.astype(str).str.strip().str.strip('"')
pre = {t: m for t, m in meta2.items() if "_Pre" in t and m["response"] != "UNK"}
df2 = df2[[t for t in pre]]
y2 = np.array([1 if pre[t]["response"] == "PRCR" else 0 for t in pre])
look2 = {s: str(S2E[s]) for s in S2E if s in set().union(*[set(v) for v in sigdefs.values()]) and str(S2E[s]) in df2.index}
s2 = score_matrix(zscore_rows(df2), look2)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), sharey=True)
def panel(ax, sdf, y, title):
    for col, off, c in [("Barrier", -0.12, "#E67E22"), ("Immune", 0.12, "#2980B9")]:
        r = sdf.loc[y == 1, col].values; nr = sdf.loc[y == 0, col].values
        bp = ax.boxplot([r, nr], positions=[1 + off, 2 + off], widths=0.22, patch_artist=True,
                        boxprops=dict(facecolor=c, alpha=0.6))
        ax.scatter(np.full(len(r), 1 + off) + np.random.normal(0, 0.03, len(r)), r, color=c, s=14, alpha=0.7)
        ax.scatter(np.full(len(nr), 2 + off) + np.random.normal(0, 0.03, len(nr)), nr, color=c, s=14, alpha=0.7, marker="x")
        u = stats.mannwhitneyu(r, nr, alternative="two-sided")
        auc = u[0] / (len(r) * len(nr))
        ax.text(1.5 + off, ax.get_ylim()[1] * 0.92, f"{col}\nAUC={auc:.2f}, p={u[1]:.2f}",
                ha="center", fontsize=8, color=c)
    ax.set_xticks([1, 2]); ax.set_xticklabels(["Responder", "Non-responder"])
    ax.set_title(title); ax.axhline(0, color="gray", lw=0.8, ls="--")

panel(axes[0], s1, y1, f"GSE78220 (Hugo 2016, anti-PD-1)\nn={len(y1)} (R={int(y1.sum())}, NR={int((y1==0).sum())})")
panel(axes[1], s2, y2, f"GSE91061 (Riaz 2017, anti-PD-1 pre)\nn={len(y2)} (R={int(y2.sum())}, NR={int((y2==0).sum())})")
fig.suptitle("Barrier (crosslink+CAF+hypoxia) and Immune (CD8T+TNK+PDL1) composite scores vs ICB response", fontsize=10)
fig.tight_layout()
p = r"D:\sparta\results\figures\icb_bulk_crosscohort.png"
fig.savefig(p, dpi=150); plt.close(fig)
print("saved", p)
