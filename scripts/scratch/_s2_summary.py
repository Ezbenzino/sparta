# S2 跨队列汇总 + 屏障带厚度分析（论文补充材料）
import sys, json, glob
sys.path.insert(0, r"D:\sparta")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

from sparta.io_ import Paths, load_config, load_graph
from sparta.barrier import compute_b_cell, scores_from_adata

cfg = load_config(None)
P = Paths(cfg)
slides = ["MEL01","MEL02","MEL03","MEL04","CSCC01","CSCC02","CSCC03","CSCC04"]
cfg_cell = cfg["barrier"]["b_cell"]

# 1) 屏障带厚度：cut_nodes 的最近邻间距与"带内密度"作为厚度代理
thick = {}
for sid in slides:
    adata = sc.read_h5ad(P.scored(sid))
    A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
    S = scores_from_adata(adata)
    coords = np.asarray(adata.obsm["spatial_um"])
    base = compute_b_cell(A, S["ecm"], S["caf"], source, sink, **cfg_cell)
    cut = np.asarray(base["cut_nodes"], int)
    from scipy.spatial import cKDTree
    dnn = cKDTree(coords[cut]).query(coords[cut], k=2)[0][:,1]
    # 厚度代理 = 带内中位最近邻距离 × 带内邻接度（越高带越"弥散"）
    # 简单稳健版本：带内节点的中位 NN 距离（密=薄带，疏=厚带）
    thick[sid] = dict(n_cut=len(cut), nn_med=float(np.median(dnn)),
                      nn_p90=float(np.quantile(dnn, 0.9)))

print(f"{'切片':<7}{'割节点':>8}{'带内中位NN(μm)':>14}{'NN P90':>9}")
for sid in slides:
    t = thick[sid]
    print(f"{sid:<7}{t['n_cut']:>8}{t['nn_med']:>14.1f}{t['nn_p90']:>9.1f}")

# 2) 汇总图：每个切片在每个 k 的 ratio，按队列着色
fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
colors = {"MEL": "#C0392B", "CSCC": "#2E8B57"}
mel_x, cscc_x = [], []
for sid in slides:
    d = json.load(open(rf"D:\sparta\results\counterfactual\{sid}.json", encoding="utf-8"))
    s2 = d.get("s2") or d.get("S2") or {}
    for k, v in s2.get("per_k", {}).items():
        r = v.get("ratio_vs_in_cut")
        p = v.get("p_vs_in_cut")
        if r is None: continue
        q = "MEL" if sid.startswith("MEL") else "CSCC"
        c = colors[q]
        ax[0].scatter(k, r, color=c, s=60,
                      marker="o" if p < 0.05 else "^", alpha=0.8)
        (mel_x if q=="MEL" else cscc_x).append(r)
ax[0].axhline(1.0, color="gray", ls="--", lw=1)
ax[0].set_xlabel("移除 spot 数 k")
ax[0].set_ylabel("剩余屏障比 (分散/连续)")
ax[0].set_title("S2: 连续缺口 vs 等量分散移除 (实心=p<0.05)")
ax[0].legend([f"MEL (n={len(mel_x)})", f"CSCC (n={len(cscc_x)})"],
             [plt.Line2D([0],[0],marker='o',color='w',markerfacecolor=colors['MEL']),
              plt.Line2D([0],[0],marker='o',color='w',markerfacecolor=colors['CSCC'])])

ax[1].bar(["MEL (转移灶)", "CSCC (原发)"],
          [np.median(mel_x), np.median(cscc_x)],
          color=[colors["MEL"], colors["CSCC"]])
ax[1].set_ylabel("中位 ratio_vs_in_cut")
ax[1].set_title("连续性的效应量（队列对比）")
for i, v in enumerate([np.median(mel_x), np.median(cscc_x)]):
    ax[1].text(i, v+0.01, f"{v:.2f}x", ha="center")
fig.tight_layout()
p = P.figure("s2_cross_cohort.png")
fig.savefig(p, dpi=150); plt.close(fig)
print(f"已保存 {p}")

json.dump(thick, open(r"D:\sparta\results\validation\s2_thickness.json","w"), indent=2)
print("已写出 s2_thickness.json")
