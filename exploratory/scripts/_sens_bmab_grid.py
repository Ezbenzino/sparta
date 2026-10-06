# B_mAb 二维参数网格：交联主导是否稳健？产出 JSON + 热图
import sys, json, itertools
sys.path.insert(0, r"D:\sparta")
import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sparta.io_ import Paths, load_config, load_graph, save_json
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


def run(mab_cfg):
    dd = dissociation_drivers(A, S["ecm"], S["caf"], S["crosslink"], S["ag_target"],
                              source, vessel, cfg_cell=cfg_cell, cfg_mab=mab_cfg, r_nm=r_nm)
    return dict(potential=dd["dissociation_potential"],
                ecm=dd["frac_shared_ecm"], antigen=dd["frac_antigen"],
                crosslink=dd["frac_crosslink"])


grids = {
    "lam_x_beta": dict(
        x_name="lam (ECM 衰减)", y_name="beta (交联收缩)",
        xs=[0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0],
        ys=[0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.5, 6.0, 9.0, 12.0],
        key=("lam", "beta")),
    "xi0_x_beta": dict(
        x_name="xi0_nm (基线网孔)", y_name="beta (交联收缩)",
        xs=[8, 12, 16, 20, 30, 40, 60, 80],
        ys=[0.5, 1.0, 1.5, 2.0, 3.0, 4.5, 6.0, 9.0, 12.0],
        key=("xi0_nm", "beta")),
    "kappa_x_kd": dict(
        x_name="kappa_w (BSB权重)", y_name="kd_eff (解离常数)",
        xs=[0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 6.0, 10.0],
        ys=[0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0],
        key=("kappa_w", "kd_eff")),
}

results = {}
for gname, g in grids.items():
    kx, ky = g["key"]
    M = np.zeros((len(g["ys"]), len(g["xs"])))   # 交联%
    Pot = np.zeros_like(M)
    for i, yv in enumerate(g["ys"]):
        for j, xv in enumerate(g["xs"]):
            m = dict(base_mab)
            m[kx] = xv; m[ky] = yv
            r = run(m)
            M[i, j] = r["crosslink"] * 100
            Pot[i, j] = r["potential"]
    results[gname] = dict(grid=g, crosslink_pct=M.tolist(), potential=Pot.tolist())

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    im0 = ax[0].imshow(M, origin="lower", aspect="auto", cmap="RdBu_r",
                       vmin=0, vmax=100)
    ax[0].set_xticks(range(len(g["xs"]))); ax[0].set_xticklabels(g["xs"], fontsize=8)
    ax[0].set_yticks(range(len(g["ys"]))); ax[0].set_yticklabels(g["ys"], fontsize=8)
    ax[0].set_xlabel(g["x_name"]); ax[0].set_ylabel(g["y_name"])
    ax[0].set_title(f"交联贡献 %  ({gname})")
    fig.colorbar(im0, ax=ax[0])
    for i in range(len(g["ys"])):
        for j in range(len(g["xs"])):
            ax[0].text(j, i, f"{M[i,j]:.0f}", ha="center", va="center", fontsize=6)

    im1 = ax[1].imshow(np.log10(np.maximum(Pot, 0.01)), origin="lower",
                       aspect="auto", cmap="viridis")
    ax[1].set_xticks(range(len(g["xs"]))); ax[1].set_xticklabels(g["xs"], fontsize=8)
    ax[1].set_yticks(range(len(g["ys"]))); ax[1].set_yticklabels(g["ys"], fontsize=8)
    ax[1].set_xlabel(g["x_name"]); ax[1].set_ylabel(g["y_name"])
    ax[1].set_title(f"log10(解离潜力)  ({gname})")
    fig.colorbar(im1, ax=ax[1])
    fig.tight_layout()
    fig.savefig(P.figure(f"bmab_sensitivity_{gname}.png"), dpi=150)
    plt.close(fig)
    print(f"已保存 bmab_sensitivity_{gname}.png")

# 摘要：默认参数邻域内交联%的最小值（稳健性区间）
# lam in [1,8], beta in [1.5,6], xi0 in [10,40]
m = dict(base_mab)
summary = {}
for gname in ["lam_x_beta", "xi0_x_beta"]:
    g = results[gname]["grid"]
    M = np.array(results[gname]["crosslink_pct"])
    Pot = np.array(results[gname]["potential"])
    # 提取"中段"参数的交联%范围
    xs = g["xs"]; ys = g["ys"]
    jj = [j for j in range(len(xs)) if 1.0 <= xs[j] <= 8.0] if gname == "lam_x_beta" else \
         [j for j in range(len(xs)) if 10 <= xs[j] <= 40]
    ii = [i for i in range(len(ys)) if 1.5 <= ys[i] <= 6.0]
    sub = M[np.ix_(ii, jj)]
    potsub = Pot[np.ix_(ii, jj)]
    summary[gname] = dict(
        crosslink_pct_min=float(sub.min()), crosslink_pct_max=float(sub.max()),
        potential_min=float(potsub.min()), potential_max=float(potsub.max()),
        n_cells=int(sub.size),
        n_above_GO= int((potsub > 1.0).sum()))

out = dict(slide=SID, base_params=base_mab, grids=results, summary=summary)
p = P.validation("bmab_sensitivity.json")
save_json(p, out)
print(f"\n已写出 {p}")
print(f"\n稳健性摘要（合理参数邻域）:")
for k, v in summary.items():
    print(f"  {k}: 交联%范围 {v['crosslink_pct_min']:.1f}–{v['crosslink_pct_max']:.1f}  "
          f"潜力范围 {v['potential_min']:.2f}–{v['potential_max']:.2f}  "
          f"GO占比 {v['n_above_GO']}/{v['n_cells']}")
