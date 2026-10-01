# benchmark_ext 英文汇总图：B_cell vs 替代度量的敏感性 + 信息冗余
import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, r"D:\sparta")
from sparta.io_ import Paths, load_config
P = Paths(load_config(None))
d = json.load(open(r"D:\sparta\results\validation\benchmark_ext.json", encoding="utf-8"))
sids = list(d)
cohort = {s: ("Melanoma" if s.startswith("MEL") else "cSCC") for s in sids}
cmap = {"Melanoma": "#C0392B", "cSCC": "#2E8B57"}

# 指标列表（z 得分）
metrics = [("z_b_cell", "B_cell (min-cut)"),
           ("z_dens", "CAF density"),
           ("z_nbr", "CAF neigh. enrich."),
           ("z_ripley", "Ripley's L (CAF)"),
           ("z_tdist", "T-to-core dist."),
           ("z_tinf", "T infiltr. fraction")]
keys = [m[0] for m in metrics]
labels = [m[1] for m in metrics]
zmat = np.array([[d[s][k] for k in keys] for s in sids])

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), gridspec_kw={"width_ratios": [1.4, 1]})

# 左：z 得分热图
im = axes[0].imshow(zmat, cmap="RdBu_r", vmin=-4, vmax=40, aspect="auto")
axes[0].set_xticks(range(len(metrics))); axes[0].set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
axes[0].set_yticks(range(len(sids))); axes[0].set_yticklabels(sids, fontsize=9)
for i in range(len(sids)):
    for j in range(len(metrics)):
        v = zmat[i, j]
        if np.isfinite(v):
            axes[0].text(j, i, f"{v:+.1f}", ha="center", va="center", fontsize=7,
                         color="white" if abs(v) > 20 else "black")
axes[0].set_title("Z-score of real vs spatially shuffled\n(higher = more sensitive to arrangement)")
axes[0].axvline(0.5, color="white", lw=2)
plt.colorbar(im, ax=axes[0], fraction=0.03)

# 右：B_cell 空间场与 CAF / T_NK 的相关
caf_rho = [d[s]["rho_Bfield_vs_CAF"] for s in sids]
tnk_rho = [d[s]["rho_Bfield_vs_TNK"] for s in sids]
x = np.arange(len(sids)); w = 0.38
axes[1].bar(x - w/2, caf_rho, w, label="vs CAF score", color="#E67E22", alpha=0.85)
axes[1].bar(x + w/2, tnk_rho, w, label="vs T/NK score", color="#2980B9", alpha=0.85)
axes[1].axhline(0, color="0.5", lw=1)
axes[1].set_xticks(x); axes[1].set_xticklabels(sids, fontsize=9)
axes[1].set_ylabel("Spearman ρ (B_cell field)")
axes[1].set_title("B_cell field is NOT a T-cell density rewrite\n(ρ≈0 vs T/NK)")
axes[1].legend(fontsize=8)
for i, (a, b) in enumerate(zip(caf_rho, tnk_rho)):
    axes[1].text(i - w/2, a + (0.02 if a >= 0 else -0.06), f"{a:.2f}", ha="center", fontsize=7)
    axes[1].text(i + w/2, b + (0.02 if b >= 0 else -0.06), f"{b:.2f}", ha="center", fontsize=7)

fig.tight_layout()
p = P.figure("benchmark_ext_dual_cohort.png")
fig.savefig(p, dpi=150); plt.close(fig)
print("saved", p)
