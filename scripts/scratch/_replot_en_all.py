# 英文重出全部 barrier + counterfactual 图（8 张 × 2）
import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.stats import spearmanr
sys.path.insert(0, r"D:\sparta")

from sparta.io_ import Paths, load_config, load_graph
from sparta.barrier import scores_from_adata
from sparta.counterfactual import s1_spatial_permutation

cfg = load_config(None)
P = Paths(cfg)
slides = ["MEL01", "MEL02", "MEL03", "MEL04", "CSCC01", "CSCC02", "CSCC03", "CSCC04"]

import scanpy as sc

# ---------------- 英文绘图函数（内联） ----------------
def barrier_line(coords, cut_edges, node_color, source, sink, ax, title):
    ax.scatter(coords[:, 0], coords[:, 1], c=node_color, s=6, cmap="Greys", alpha=0.55, linewidths=0)
    if len(source):
        ax.scatter(coords[source, 0], coords[source, 1], s=14, c="#2E86C1", label="Source: immune entry", linewidths=0)
    if len(sink):
        ax.scatter(coords[sink, 0], coords[sink, 1], s=14, c="#C0392B", label="Sink: tumor core", linewidths=0)
    if len(cut_edges):
        segs = [[coords[u], coords[v]] for u, v in cut_edges]
        ax.add_collection(LineCollection(segs, colors="#F1C40F", linewidths=2.0, alpha=0.95, zorder=5))
        ax.plot([], [], color="#F1C40F", lw=2.0, label="Min-cut barrier line")
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=10); ax.legend(loc="upper right", fontsize=7, framealpha=0.85)

def barrier_field(coords, values, ax, title, cbar_label="B_mAb"):
    finite = np.isfinite(values)
    sc_ = ax.scatter(coords[finite, 0], coords[finite, 1], c=np.asarray(values)[finite],
                     cmap="magma_r", s=10, linewidths=0)
    if (~finite).any():
        ax.scatter(coords[~finite, 0], coords[~finite, 1], s=10, c="0.85", linewidths=0, label="unreachable")
        ax.legend(loc="upper right", fontsize=7)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=10)
    cb = plt.colorbar(sc_, ax=ax, fraction=0.045, pad=0.02); cb.set_label(cbar_label)

def decoupling(ecm_caf, b_mab, ax):
    x, y = np.asarray(ecm_caf, float), np.asarray(b_mab, float)
    ok = np.isfinite(x) & np.isfinite(y)
    ax.scatter(x[ok], y[ok], s=8, c="#34495E", alpha=0.5, linewidths=0)
    rho, p = spearmanr(x[ok], y[ok]) if ok.sum() > 3 else (np.nan, np.nan)
    mx, my = np.nanmedian(x[ok]), np.nanmedian(y[ok])
    ax.axvline(mx, color="0.8", lw=1); ax.axhline(my, color="0.8", lw=1)
    ax.text(0.98, 0.03, f"Spearman ρ = {rho:.3f}\np = {p:.2e}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8)
    ax.text(0.02, 0.97, "mAb blocked\nT cells pass", transform=ax.transAxes,
            ha="left", va="top", fontsize=7.5, color="#C0392B")
    ax.set_xlabel("B_cell (T-cell migration barrier)"); ax.set_ylabel("B_mAb (antibody transport barrier)")
    ax.set_title("Dual-barrier decoupling", fontsize=10)

def perm_null(res, ax):
    null = np.asarray(res["null"]); real = res["b_real"]
    ax.hist(null, bins=40, color="0.75", edgecolor="white", label="Spatial permutation null")
    ax.axvline(real, color="#C0392B", lw=2.2, label=f"Observed = {real:.3f}")
    ax.set_xlabel("B_cell"); ax.set_ylabel("Frequency")
    ax.set_title(f"S1 spatial shuffle (mode={res['mode']}, z={res['z']:.1f}, p={res['p_emp']:.4f})", fontsize=9.5)
    ax.legend(fontsize=8)

def ring_break(res, ax):
    ks = sorted(res["per_k"]); b0 = res["b0"]
    tgt = [res["per_k"][k]["targeted"] / b0 for k in ks]
    rg = [res["per_k"][k]["rand_global_mean"] / b0 for k in ks]
    rc = [res["per_k"][k]["rand_in_cut_mean"] / b0 for k in ks]
    ax.plot(ks, tgt, "o-", color="#C0392B", label="Continuous gap (targeted)")
    ax.plot(ks, rc, "s--", color="#E67E22", label="Within-barrier diffuse (key ctrl)")
    ax.plot(ks, rg, "^:", color="0.55", label="Global random")
    ax.axhline(1.0, color="0.8", lw=1, ls="-")
    ax.set_xlabel("Removed spots k"); ax.set_ylabel("Remaining / baseline barrier")
    ax.set_title("S2 ring breaking", fontsize=10); ax.legend(fontsize=8)

def size_scan(res, ax):
    r = res["radii_nm"]
    ax.plot(r, res["mean_all"], "o-", color="0.45", label="Whole-section mean")
    if np.isfinite(res["mean_core"]).any():
        ax.plot(r, res["mean_core"], "s-", color="#C0392B", label="Tumor-core mean")
    ax.axvline(5.5, color="#2E86C1", ls="--", lw=1.2)
    ax.text(5.6, ax.get_ylim()[1] * 0.92, "IgG / ADC (5.5 nm)", fontsize=8, color="#2E86C1")
    ax.axvline(0.5, color="#27AE60", ls="--", lw=1.2)
    ax.text(0.6, ax.get_ylim()[1] * 0.78, "small mol.", fontsize=8, color="#27AE60")
    ax.set_xlabel("Hydrodynamic radius (nm)"); ax.set_ylabel("B_mAb")
    ax.set_title("S3 size scan", fontsize=10); ax.legend(fontsize=8)

for sid in slides:
    # ---- barrier ----
    adata = sc.read_h5ad(P.scored(sid))
    coords = np.asarray(adata.obsm["spatial_um"], float)
    S = scores_from_adata(adata)
    bn = np.load(P.barrier(sid))
    mincut = json.load(open(P.mincut(sid), encoding="utf-8"))
    A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
    ecm_caf = 0.5 * (S["ecm"] + S["caf"])

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.2))
    barrier_line(coords, mincut["cut_edges"], S["caf"], source, sink, ax[0],
                 f"{sid} min-cut barrier (B_cell={bn['b_cell'][0]:.3f})")
    barrier_field(coords, bn["b_mab"], ax[1], f"B_mAb (r=5.5 nm)")
    decoupling(ecm_caf, bn["b_mab"], ax[2])
    fig.tight_layout()
    fig.savefig(P.figure(f"barrier_{sid}.png"), dpi=140); plt.close(fig)

    # ---- counterfactual ----
    cf = json.load(open(P.counterfactual(sid), encoding="utf-8"))
    cfg_cell = cfg["barrier"]["b_cell"]
    # S1: 重算 null（seed 固定，n=500 足够直方图）
    mode0 = cf["s1"]["fixed"]["mode"] if "fixed" in cf["s1"] else list(cf["s1"].keys())[0]
    s1mode = "fixed" if "fixed" in cf["s1"] else mode0
    kw = {}
    if s1mode == "follow":
        kw = dict(endothelial=S.get("ecm"), t_nk=S.get("caf"), malignant=S.get("ag_target"),
                  cfg_source_sink={k: cfg["source_sink"][k] for k in
                                   ("q_vessel", "q_immune_nbr", "q_malig", "q_core")})
    s1r = s1_spatial_permutation(A, S["ecm"], S["caf"], source, sink, cfg_cell,
                                 n_perm=500, seed=cfg["seed"], mode=s1mode, **kw)
    s1plot = dict(s1r)
    s1plot["mode"] = s1mode
    s1plot["p_emp"] = cf["s1"][s1mode]["p_emp"]
    s1plot["z"] = cf["s1"][s1mode]["z"]
    s1plot["b_real"] = cf["s1"][s1mode]["b_real"]

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
    perm_null(s1plot, ax[0])
    ring_break(cf["s2"], ax[1])
    size_scan(cf["s3"], ax[2])
    fig.tight_layout()
    fig.savefig(P.figure(f"counterfactual_{sid}.png"), dpi=140); plt.close(fig)
    print(f"replotted {sid}")

print("done: barrier + counterfactual EN for", len(slides), "slides")
