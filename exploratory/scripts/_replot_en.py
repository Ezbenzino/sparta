# 英文版关键图重出：S2 跨队列 + B_mAb 敏感性热图
import sys, json
sys.path.insert(0, r"D:\sparta")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sparta.io_ import Paths

P = Paths(__import__("sparta.io_", fromlist=["load_config"]).load_config(None))

# ===== S2 跨队列 =====
slides = ["MEL01","MEL02","MEL03","MEL04","CSCC01","CSCC02","CSCC03","CSCC04"]
fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
colors = {"MEL": "#C0392B", "CSCC": "#2E8B57"}
mel_r, cscc_r = [], []
for sid in slides:
    d = json.load(open(rf"D:\sparta\results\counterfactual\{sid}.json", encoding="utf-8"))
    s2 = d.get("s2") or {}
    for k, v in s2.get("per_k", {}).items():
        r = v.get("ratio_vs_in_cut"); p = v.get("p_vs_in_cut")
        if r is None: continue
        q = "MEL" if sid.startswith("MEL") else "CSCC"
        ax[0].scatter(k, r, color=colors[q], s=70, alpha=0.85,
                      marker="o" if p < 0.05 else "^")
        (mel_r if q=="MEL" else cscc_r).append(r)
ax[0].axhline(1.0, color="gray", ls="--", lw=1)
ax[0].set_xlabel("number of removed spots k")
ax[0].set_ylabel("residual barrier ratio (dispersed / contiguous)")
ax[0].set_title("S2 ring-breaking: contiguous gap vs equal dispersed removal\n(filled = p<0.05)")
ax[0].legend([f"Melanoma meta. (n={len(mel_r)})", f"cSCC primary (n={len(cscc_r)})"],
             [plt.Line2D([0],[0],marker='o',color='w',markerfacecolor=colors['MEL']),
              plt.Line2D([0],[0],marker='o',color='w',markerfacecolor=colors['CSCC'])])
ax[1].bar(["Melanoma\n(metastasis)", "cSCC\n(primary)"],
          [np.median(mel_r), np.median(cscc_r)], color=[colors["MEL"], colors["CSCC"]])
ax[1].set_ylabel("median ratio vs in-cut (dispersed/contiguous)")
ax[1].set_title("Effect size of barrier continuity by cohort")
for i, v in enumerate([np.median(mel_r), np.median(cscc_r)]):
    ax[1].text(i, v+0.01, f"{v:.2f}x", ha="center")
fig.tight_layout()
fig.savefig(P.figure("s2_cross_cohort.png"), dpi=150); plt.close(fig)
print("saved s2_cross_cohort.png")

# ===== B_mAb 敏感性热图（英文标签）=====
sens = json.load(open(r"D:\sparta\results\validation\bmab_sensitivity.json", encoding="utf-8"))
for gname, gd in sens["grids"].items():
    g = gd["grid"]; M = np.array(gd["crosslink_pct"]); Pot = np.array(gd["potential"])
    xs, ys = g["xs"], g["ys"]
    xlab = {"lam_x_beta": "lambda (ECM decay)",
            "xi0_x_beta": "xi0_nm (baseline mesh size)",
            "kappa_x_kd": "kappa_w (binding-site weight)"}[gname]
    ylab = {"lam_x_beta": "beta (crosslink contraction)",
            "xi0_x_beta": "beta (crosslink contraction)",
            "kappa_x_kd": "kd_eff (dissociation constant)"}[gname]
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 4.8))
    im0 = ax[0].imshow(M, origin="lower", aspect="auto", cmap="RdBu_r", vmin=0, vmax=100)
    ax[0].set_xticks(range(len(xs))); ax[0].set_xticklabels(xs, fontsize=8)
    ax[0].set_yticks(range(len(ys))); ax[0].set_yticklabels(ys, fontsize=8)
    ax[0].set_xlabel(xlab); ax[0].set_ylabel(ylab)
    ax[0].set_title(f"crosslink contribution %  ({gname})")
    fig.colorbar(im0, ax=ax[0])
    for i in range(len(ys)):
        for j in range(len(xs)):
            ax[0].text(j, i, f"{M[i,j]:.0f}", ha="center", va="center", fontsize=6)
    im1 = ax[1].imshow(np.log10(np.maximum(Pot, 0.01)), origin="lower", aspect="auto", cmap="viridis")
    ax[1].set_xticks(range(len(xs))); ax[1].set_xticklabels(xs, fontsize=8)
    ax[1].set_yticks(range(len(ys))); ax[1].set_yticklabels(ys, fontsize=8)
    ax[1].set_xlabel(xlab); ax[1].set_ylabel(ylab)
    ax[1].set_title(f"log10(dissociation potential)  ({gname})")
    fig.colorbar(im1, ax=ax[1])
    fig.tight_layout()
    fig.savefig(P.figure(f"bmab_sensitivity_{gname}.png"), dpi=150); plt.close(fig)
    print(f"saved bmab_sensitivity_{gname}.png")
