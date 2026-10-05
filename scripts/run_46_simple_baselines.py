#!/usr/bin/env python
"""
run_46_simple_baselines.py —— 「为什么要搞这么复杂？」的头对头对比
==================================================================
审稿人必问：*"你这个图论 / 最小割 / 屏蔽泊松的方法，比直接算欧氏距离或
基质密度好在哪里？"* 本脚本给出可量化的答案。

对比的**不是**拟合优度，而是**决策质量**：
    给你一个干预预算 k（= 切片 spot 数的 1%–10%），
    让你挑 k 个 spot 去做基质正常化，谁挑的位置真能把屏障打下来？

真值（held-out，昂贵）
----------------------
``results/intervention/{sid}.intervention.npz`` 里的 ``delta_b_cell`` /
``delta_b_mab``——run_44 对**每一个 spot** 都完整重算了一遍屏障得到的
反事实效应。它不参与任何候选者的构造，因此可以作为公平竞争的真值：
    · 真值的获得需要 n 次最小割 / n 次稀疏线性求解（贵）
    · 所有候选者（含 SPARTA 自身的廉价场）都是 O(n log n) 或单次求解（便宜）

候选者（全部是"简单方法"，除非标注 SPARTA）
--------------------------------------------
  d_boundary      有符号欧氏距离到肿瘤边界（负 = 瘤内深度，正 = 瘤外距离），
                  取 -|d| —— 越靠近边界分越高。纯坐标 + 一个恶性签名，无图。
  d_vessel        欧氏距离到最近血管 spot。经典的"离血管多远"启发式，无图。
  stromal_density 半径球内 ECM+CAF（B_cell 臂）/ ECM+交联（B_mAb 臂）的局部均值。
                  纯密度，无图、无传输。
  naive_score     0.5*(ECM+CAF) 或 0.5*(ECM+交联)，单点无邻域平滑（run_24 的基线）。
  cv_ridge        上述全部简单特征 + 恶性/缺氧等签名的 **5 折 CV 岭回归**。
                  代表"简单特征族能达到的最好线性模型"，符号由数据学（交叉验证），
                  不受人为定向的影响。
  random          随机挑 k 个（10 次平均）—— 下界刻度。
  sparta_field    SPARTA 的逐点场：B_cell 用 Dijkstra 迁移代价场，B_mAb 用屏蔽
                  泊松场。一次模型求值，不做逐 spot 重算。
  sparta_cut      SPARTA 最小割封锁带成员（0/1）。一次最小割。
  sparta_delta    SPARTA 逐 spot 反事实排序本身 —— **上界参考**，不是参赛者
                  （它用了真值来排序，画出来是给其他曲线当天花板的）。

指标
----
  rho             候选者分数与真值 Δ 的 Spearman（自然定向下的有符号值；
                  同时报 |rho|，便于判断人为定向是否吃亏）
  recall@k        候选者 top-k 命中真值 top-k 的比例
  achieved@k      按候选者的 top-k **真的做一次联合干预**，屏障下降量 / 锚点
                  （B_cell 锚点 = full_cut_breach，B_mAb 锚点 = all_removed_mab）
                  —— 这是主指标，因为它就是"决策质量"本身。

无 Scanpy 路线：全部输入来自 data/interim/{sid}.nodes.npz + {sid}.graph.npz。
BRCA01/02/LN01 为外部验证臂，与主队列分开汇总（其 sink 依赖的 B_cell 读数为探索性）。

输出
----
  results/validation/simple_baselines.json
  results/figures/baseline_comparison.png / .pdf
  results/figures/baseline_comparison_rho.png / .pdf
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from scipy import stats  # noqa: E402

from sparta.io_ import (Paths, load_config, load_graph, load_json,  # noqa: E402
                        save_json, set_seed, stamp_run)
from sparta.node_tables import coords_from_nodes, load_nodes, scores_from_nodes  # noqa: E402

K_FRACS = (0.01, 0.02, 0.05, 0.10)   # 干预预算，占切片 spot 数的比例
N_RANDOM = 10                          # 随机基线的重复次数
RIDGE_ALPHA = 1.0
RIDGE_FOLDS = 5

# 方法顺序（图例顺序）与配色：色相只承载"哪个方法"，两个面板共用
METHODS = [
    ("random",          "Random",                     "#9e9e9e", "--"),
    ("d_boundary",      "Euclidean dist. to tumour boundary", "#2a78d6", "-"),
    ("d_vessel",        "Euclidean dist. to nearest vessel",  "#7b4fa8", "-"),
    ("stromal_density", "Local stromal density",      "#2f9e6f", "-"),
    ("naive_score",     "Naive score (spot-wise)",    "#c9972f", "-"),
    ("cv_ridge",        "CV ridge (all simple feats)", "#8c6d3f", "-"),
    ("sparta_field",    "SPARTA per-spot field",      "#1fa8a8", "-"),
    ("sparta_cut",      "SPARTA min-cut band",        "#0b4f9e", "-"),
    ("sparta_delta",    "SPARTA counterfactual (reference)", "#eb6834", "--"),
]

MCOLOR = {m[0]: m[2] for m in METHODS}
MLABEL = {m[0]: m[1] for m in METHODS}


# --------------------------------------------------------------------------
# 简单特征：全部只吃坐标 + 签名，不碰图、不做传输
# --------------------------------------------------------------------------
def signed_dist_to_boundary(xy: np.ndarray, malig_n: np.ndarray,
                            q_malig: float = 0.70) -> np.ndarray | None:
    """有符号欧氏距离到肿瘤边界。负 = 瘤内深度，正 = 瘤外距离。

    纯几何：只用到坐标与一个恶性签名，不涉及图、不涉及传输。
    返回值中 |d| 小 = 位于肿瘤-基质交界面上。
    """
    n = len(xy)
    tum = malig_n >= np.quantile(malig_n, q_malig)
    if tum.sum() == 0 or (~tum).sum() == 0:
        return None
    t_tum = cKDTree(xy[tum])
    t_non = cKDTree(xy[~tum])
    d_to_tum = t_tum.query(xy, k=1)[0]
    d_to_non = t_non.query(xy, k=1)[0]
    return np.where(tum, -d_to_non, d_to_tum)


def dist_to_set(xy: np.ndarray, idx: np.ndarray) -> np.ndarray | None:
    """欧氏距离到集合 idx 中最近的点。集合为空时返回 None（不要静默填 0）。"""
    idx = np.asarray(idx, int)
    if len(idx) == 0:
        return None
    return cKDTree(xy[idx]).query(xy, k=1)[0]


def local_mean(xy: np.ndarray, x: np.ndarray, radius: float) -> np.ndarray:
    """半径球内的局部均值（含自身）。纯欧氏邻域，不用图。"""
    n = len(xy)
    tree = cKDTree(xy)
    pairs = tree.query_pairs(r=radius, output_type="ndarray")
    if len(pairs) == 0:
        return x.copy()
    A = sp.coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])),
                      shape=(n, n)).tocsr()
    A = A + A.T + sp.eye(n)
    deg = np.asarray(A.sum(axis=1)).ravel()
    return np.asarray(A @ x).ravel() / deg


def _cv_ridge_predict(X: np.ndarray, y: np.ndarray, alpha: float = RIDGE_ALPHA,
                      k: int = RIDGE_FOLDS, seed: int = 0) -> np.ndarray:
    """k 折交叉验证的岭回归预测（折内标准化、折内中位数填补）。

    用 CV 而不是样本内拟合，否则"简单特征族的最好线性模型"会靠过拟合
    拿到不公平的优势。
    """
    n = len(y)
    pred = np.full(n, np.nan)
    idx = np.arange(n)
    rng = np.random.default_rng(seed)
    rng.shuffle(idx)
    for fold in np.array_split(idx, k):
        tr = np.setdiff1d(idx, fold)
        Xtr = np.nan_to_num(X[tr], nan=np.nan)
        med = np.nanmedian(Xtr, axis=0)
        med = np.where(np.isfinite(med), med, 0.0)
        Xtr = np.where(np.isfinite(Xtr), Xtr, med)
        sd = Xtr.std(axis=0)
        sd = np.where(sd < 1e-12, 1.0, sd)
        mu = Xtr.mean(axis=0)
        Xtr_s = (Xtr - mu) / sd
        Xte_s = (np.where(np.isfinite(X[fold]), X[fold], med) - mu) / sd
        ytr = y[tr]
        ymu = float(ytr.mean())
        p = Xtr_s.shape[1]
        w = np.linalg.solve(Xtr_s.T @ Xtr_s + alpha * np.eye(p),
                            Xtr_s.T @ (ytr - ymu))
        pred[fold] = Xte_s @ w + ymu
    return pred


# --------------------------------------------------------------------------
# 工具
# --------------------------------------------------------------------------
def _finite_score(v: np.ndarray | None, n: int) -> np.ndarray | None:
    """把含 inf/NaN 的分数变成可用于排序的有限分数。None 表示该特征不可用。"""
    if v is None:
        return None
    v = np.asarray(v, float).copy()
    if v.shape[0] != n:
        return None
    bad = ~np.isfinite(v)
    if bad.all():
        return None
    if bad.any():
        # 不可达/无穷处填最大值（保守：排在最后，不会被当成高价值位点）
        v[bad] = np.nanmax(v[~bad]) if (~bad).any() else 0.0
    return v


def _topk(score: np.ndarray, k: int) -> np.ndarray:
    """取分数最高的 k 个索引（降序、稳定）。"""
    k = int(min(max(k, 1), len(score)))
    return np.argsort(-score, kind="stable")[:k]


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(stats.spearmanr(a, b).statistic)


def _cohort(sid: str) -> str:
    if sid.startswith("CSCC"):
        return "CSCC"
    if sid.startswith("MEL"):
        return "MEL"
    return "EXT"


# --------------------------------------------------------------------------
# 单切片
# --------------------------------------------------------------------------
def process_section(sid: str, cfg: dict, P: Paths, outdir: Path) -> tuple[str, dict]:
    arr_path = outdir / f"{sid}.intervention.npz"
    if not arr_path.exists():
        return sid, dict(status="no_intervention_npz")

    ranking = load_json(P.validation("intervention_ranking.json"))
    row44 = ranking.get("per_section", {}).get(sid)
    if row44 is None:
        return sid, dict(status="no_ranking_row")

    z = np.load(arr_path, allow_pickle=True)
    delta = {"cell": np.asarray(z["delta_b_cell"], float),
             "mab": np.asarray(z["delta_b_mab"], float)}
    cut_mask = np.asarray(z["cut_mask"], bool)
    xy = np.asarray(z["coords_um"], float)
    n = len(xy)

    A, _D, source, sink, vessel, _ = load_graph(P.graph(sid))
    nodes = load_nodes(Path(P.interim) / f"{sid}.nodes.npz")
    S = scores_from_nodes(nodes)
    xy_n = coords_from_nodes(nodes)
    if xy_n.shape == xy.shape and np.allclose(xy_n, xy):
        pass                      # 坐标一致，放心用
    else:
        xy = xy_n                 # 以节点表为准（两者行序同源，理论上必然一致）

    cfg_cell = cfg["barrier"]["b_cell"]
    cfg_mab = cfg["barrier"]["b_mab"]
    q_malig = cfg["source_sink"]["q_malig"]
    radius = float(cfg["graph"]["radius_um"])
    low_q = float(row44.get("low_q", 0.05))

    malig_n = np.asarray(nodes.get("obs__Malignant_n"), float)
    if malig_n is None or malig_n.shape[0] != n:
        malig_n = S.get("malignant", np.full(n, 0.5))

    # ---- 简单特征（不含任何 SPARTA 输出）----
    d_bnd_raw = signed_dist_to_boundary(xy, malig_n, q_malig)
    d_boundary = None if d_bnd_raw is None else -np.abs(d_bnd_raw)   # 越近边界分越高
    d_vessel = dist_to_set(xy, vessel)
    ecm, caf, xl = S["ecm"], S["caf"], S["crosslink"]
    ag = S["ag_target"]

    dens_cell = local_mean(xy, 0.5 * (ecm + caf), radius)
    dens_mab = local_mean(xy, 0.5 * (ecm + xl), radius)
    naive_cell = 0.5 * (ecm + caf)
    naive_mab = 0.5 * (ecm + xl)

    # 简单特征矩阵（供 CV 岭回归）
    simple_cols = [
        ("d_boundary", d_boundary),
        ("d_vessel", d_vessel),
        ("abs_d_boundary", None if d_bnd_raw is None else np.abs(d_bnd_raw)),
        ("malignant_n", malig_n),
        ("ecm_n", ecm), ("caf_n", caf), ("xlink_n", xl), ("ag_n", ag),
        ("hypoxia_n", S.get("hypoxia", np.full(n, np.nan))),
    ]
    simple_cols.append(("dens_cell", dens_cell))
    simple_cols.append(("dens_mab", dens_mab))
    X = np.column_stack([c for _, c in simple_cols if c is not None
                         and np.asarray(c).shape[0] == n])
    keep = [i for i, c in enumerate(X.T)
            if np.isfinite(c).any() and np.nanstd(c) > 1e-12]
    X = X[:, keep] if len(keep) else np.zeros((n, 0))

    # ---- SPARTA 的廉价读数 ----
    from sparta.barrier import (compute_b_cell, compute_b_cell_field,
                                compute_b_mab)
    base_cell = compute_b_cell(A, ecm, caf, source, sink, **cfg_cell)
    b_cell0 = float(base_cell["b_cell"])
    cut_nodes_rt = np.asarray(base_cell.get("cut_nodes", []), int)

    base_mab = compute_b_mab(A, ecm, xl, ag, vessel, **cfg_mab)
    reach = np.asarray(base_mab["reachable"], bool)
    b_mab0 = float(np.nanmean(base_mab["b_mab"][reach])) if reach.any() else np.nan

    cell_field = _finite_score(
        compute_b_cell_field(A, ecm, caf, source, **cfg_cell)["b_cell_field"], n)
    mab_field = _finite_score(base_mab["b_mab"], n)
    cut_score = cut_mask.astype(float)          # run_44 落盘的封锁带
    if len(cut_nodes_rt):                        # 兜底：与当前参数重算一致
        cm = np.zeros(n, bool)
        cm[cut_nodes_rt] = True
        if not cm.any() or cm.sum() == cut_mask.sum():
            cut_score = cm.astype(float)

    q_ecm = float(np.quantile(ecm, low_q))
    q_caf = float(np.quantile(caf, low_q))
    q_xl = float(np.quantile(xl, low_q))

    anchors = row44.get("anchors", {})
    anchor = {"cell": float(anchors.get("full_cut_breach", np.nan)),
              "mab": float(anchors.get("all_removed_mab", np.nan))}

    def _achieved(idx_set: np.ndarray, op: str) -> float:
        """按 idx_set 真的做一次联合干预，返回 下降量 / 锚点。"""
        if not np.isfinite(anchor[op]) or anchor[op] <= 0:
            return float("nan")
        if op == "cell":
            e, c = ecm.copy(), caf.copy()
            e[idx_set] = q_ecm
            c[idx_set] = q_caf
            b = compute_b_cell(A, e, c, source, sink, **cfg_cell)["b_cell"]
            if not (np.isfinite(b) and np.isfinite(b_cell0)):
                return float("nan")
            return float((b_cell0 - b) / anchor[op])
        e, x = ecm.copy(), xl.copy()
        e[idx_set] = q_ecm
        x[idx_set] = q_xl
        bm = compute_b_mab(A, e, x, ag, vessel, **cfg_mab)
        m = float(np.nanmean(bm["b_mab"][reach]))
        if not (np.isfinite(m) and np.isfinite(b_mab0)):
            return float("nan")
        return float((b_mab0 - m) / anchor[op])

    rng = np.random.default_rng(cfg["seed"])
    out_ops = {}
    for op, gt in (("cell", delta["cell"]), ("mab", delta["mab"])):
        gt = _finite_score(gt, n)
        if gt is None:
            out_ops[op] = dict(status="no_ground_truth")
            continue

        if op == "cell":
            dens, naive, field = dens_cell, naive_cell, cell_field
        else:
            dens, naive, field = dens_mab, naive_mab, mab_field

        ridge_pred = _cv_ridge_predict(X, gt, seed=cfg["seed"]) if X.shape[1] else None

        cands = {
            "d_boundary": d_boundary,
            "d_vessel": d_vessel,
            "stromal_density": dens,
            "naive_score": naive,
            "cv_ridge": ridge_pred,
            "sparta_field": field,
            "sparta_cut": cut_score,
            "sparta_delta": gt,
        }
        # 方向约定：分数越高 = 越该干预。d_vessel / 密度 / naive 取原符号
        # （远处更缺血、基质越密越挡、签名越高越挡），这是事前的生物学定向，
        # 不是看着真值调的。同时报 |rho| 以便读者判断定向有没有吃亏。

        per_method = {}
        gt_order = _topk(gt, max(1, int(round(0.05 * n))))
        gt_top5 = set(gt_order.tolist())

        for name, sc in cands.items():
            sc = _finite_score(sc, n)
            if sc is None:
                per_method[name] = dict(available=False)
                continue
            rho = _spearman(sc, gt)
            rec = {}
            ach = {}
            for f in K_FRACS:
                k = max(1, int(round(f * n)))
                top = _topk(sc, k)
                rec[f"recall_{int(f*100)}pct"] = float(
                    len(set(top.tolist()) & gt_top5) / max(k, 1))
                ach[f"achieved_{int(f*100)}pct"] = _achieved(top, op)
            per_method[name] = dict(
                available=True,
                rho=rho,
                rho_abs=abs(rho) if np.isfinite(rho) else float("nan"),
                **rec, **ach,
            )

        # 随机基线：多次抽样取均值（下界刻度）
        rec_r, ach_r = {}, {}
        rand_curves = []
        for f in K_FRACS:
            k = max(1, int(round(f * n)))
            hits, drops, curves = [], [], []
            for _ in range(N_RANDOM):
                top = rng.choice(n, size=k, replace=False)
                hits.append(len(set(top.tolist()) & gt_top5) / max(k, 1))
                d = _achieved(np.asarray(top, int), op)
                drops.append(d)
                curves.append(d)
            rec_r[f"recall_{int(f*100)}pct"] = float(np.mean(hits))
            ach_r[f"achieved_{int(f*100)}pct"] = float(np.nanmean(drops))
            rand_curves.append(float(np.nanstd(drops)))
        per_method["random"] = dict(available=True, rho=0.0, rho_abs=0.0,
                                    **rec_r, **ach_r)

        out_ops[op] = dict(
            status="ok",
            baseline=float(b_cell0 if op == "cell" else b_mab0),
            anchor=anchor[op],
            n=n,
            n_cut=int(cut_mask.sum()),
            methods=per_method,
        )

    return sid, dict(status="ok", cohort=_cohort(sid), n=n,
                     low_q=low_q, operators=out_ops)


# --------------------------------------------------------------------------
# 汇总与出图
# --------------------------------------------------------------------------
def _collect(rows: dict, op: str, key: str) -> dict[str, np.ndarray]:
    out: dict[str, list] = {}
    for r in rows.values():
        o = r.get("operators", {}).get(op, {})
        if o.get("status") != "ok":
            continue
        for mname, m in o.get("methods", {}).items():
            if not m.get("available"):
                continue
            v = m.get(key)
            if v is not None and np.isfinite(v):
                out.setdefault(mname, []).append(float(v))
    return {k: np.asarray(v) for k, v in out.items()}


def summarize(rows: dict, op: str) -> dict:
    agg = {}
    for key in ("rho", "rho_abs",
                *(f"recall_{int(f*100)}pct" for f in K_FRACS),
                *(f"achieved_{int(f*100)}pct" for f in K_FRACS)):
        d = _collect(rows, op, key)
        agg[key] = {m: dict(median=float(np.median(v)), mean=float(np.mean(v)),
                            n=int(len(v)),
                            q25=float(np.percentile(v, 25)),
                            q75=float(np.percentile(v, 75)))
                    for m, v in d.items()}
    return agg


def _make_figures(agg_main: dict, n_slides: int, P: Paths, cfg: dict) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    written = []

    # ---- 图 1：实际取得的干预效果 vs 预算 ----
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.9))
    for ax, op, ttl in ((axes[0], "cell",
                         "B_cell  —  fraction of a full cut-breach achieved"),
                        (axes[1], "mab",
                         "B_mAb  —  fraction of the matrix-free limit achieved")):
        for mname, label, col, ls in METHODS:
            xs, ys, es = [], [], []
            for f in K_FRACS:
                key = f"achieved_{int(f*100)}pct"
                m = agg_main[op].get(key, {}).get(mname)
                if m is None:
                    continue
                xs.append(f * 100)
                ys.append(m["mean"])
                es.append((m["q75"] - m["q25"]) / 2.0)
            if not xs:
                continue
            ax.errorbar(xs, ys, yerr=es, color=col, ls=ls,
                        lw=2.4 if mname.startswith("sparta") else 1.6,
                        marker="o", ms=4.5, capsize=2.5,
                        alpha=0.95, label=label)
        ax.set_xlabel("intervention budget  (k as % of spots)")
        ax.set_ylabel("achieved barrier reduction\n(fraction of anchor)")
        ax.set_title(ttl, fontsize=10.5)
        ax.set_xticks([f * 100 for f in K_FRACS])
        ax.axhline(0.0, color="0.8", lw=0.8, zorder=0)
        ax.grid(alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
    axes[0].legend(fontsize=7.6, loc="upper left", framealpha=0.94, ncol=1)
    fig.suptitle(f"Does the graph earn its keep?  "
                 f"main cohort, n = {n_slides} sections; "
                 f"error bars = IQR across sections", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for ext in ("png", "pdf"):
        p = P.figure(f"baseline_comparison.{ext}")
        fig.savefig(p, dpi=200 if ext == "png" else None)
        written.append(str(p))
    plt.close(fig)

    # ---- 图 2：排序质量（Spearman）----
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.4))
    order = [m[0] for m in METHODS if m[0] != "random"]
    for ax, op, ttl in ((axes[0], "cell", "B_cell:  Spearman with per-spot ground-truth Δ"),
                        (axes[1], "mab", "B_mAb:  Spearman with per-spot ground-truth Δ")):
        med, lo, hi, cols, labels = [], [], [], [], []
        for mname in order:
            m = agg_main[op].get("rho_abs", {}).get(mname)
            if m is None:
                continue
            med.append(m["median"]); lo.append(m["q25"]); hi.append(m["q75"])
            cols.append(MCOLOR[mname])
            labels.append(MLABEL[mname])
        ypos = np.arange(len(med))
        ax.barh(ypos, med, color=cols, alpha=0.9, height=0.68)
        ax.errorbar(med, ypos, xerr=[np.subtract(med, lo), np.subtract(hi, med)],
                    fmt="none", ecolor="0.25", elinewidth=1.0, capsize=2.5)
        ax.set_yticks(ypos)
        ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("|Spearman rho| with ground-truth Δ  (median, IQR)")
        ax.set_title(ttl, fontsize=10.5)
        ax.grid(axis="x", alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
    fig.suptitle("Ranking quality: how well each score orders the true "
                 "counterfactual effect (|rho| — sign-agnostic)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    for ext in ("png", "pdf"):
        p = P.figure(f"baseline_comparison_rho.{ext}")
        fig.savefig(p, dpi=200 if ext == "png" else None)
        written.append(str(p))
    plt.close(fig)

    return written


# --------------------------------------------------------------------------
def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--slides", nargs="*", default=None,
                    help="只跑这些切片（默认：intervention 目录下有 npz 的全部）")
    ap.add_argument("--figures-only", action="store_true",
                    help="不重算，用已有 JSON 重绘图")
    args = ap.parse_args()

    cfg = load_config(None)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    outdir = Path(P.results) / "intervention"
    json_path = Path(P.validation("simple_baselines.json"))

    if args.figures_only:
        old = load_json(str(json_path))
        _make_figures(old["main_cohort"]["aggregate"],
                      old["main_cohort"]["n_sections"], P, cfg)
        print("figures redrawn")
        return

    if args.slides:
        sids = list(args.slides)
    else:
        sids = sorted(p.stem.split(".")[0] for p in outdir.glob("*.intervention.npz"))
    print(f"run_46 simple baselines: {len(sids)} sections\n")

    rows = {}
    t_all = time.perf_counter()
    for sid in sids:
        t0 = time.perf_counter()
        sid, r = process_section(sid, cfg, P, outdir)
        rows[sid] = r
        if r.get("status") != "ok":
            print(f"  [skip] {sid}: {r.get('status')}")
            continue
        oc = r["operators"].get("cell", {})
        om = r["operators"].get("mab", {})
        a5c = oc.get("methods", {}).get("sparta_delta", {}).get("achieved_5pct")
        b5c = oc.get("methods", {}).get("d_boundary", {}).get("achieved_5pct")
        a5m = om.get("methods", {}).get("sparta_delta", {}).get("achieved_5pct")
        b5m = om.get("methods", {}).get("stromal_density", {}).get("achieved_5pct")
        print(f"  [done] {sid:<8} n={r['n']:<5} "
              f"B_cell@5%: SPARTA {a5c if a5c is not None else float('nan'):.3f} "
              f"vs dist {b5c if b5c is not None else float('nan'):.3f} | "
              f"B_mAb@5%: SPARTA {a5m if a5m is not None else float('nan'):.3f} "
              f"vs density {b5m if b5m is not None else float('nan'):.3f}   "
              f"({time.perf_counter()-t0:.1f}s)")

    main_rows = {k: v for k, v in rows.items()
                 if v.get("status") == "ok" and v.get("cohort") in ("CSCC", "MEL")}
    ext_rows = {k: v for k, v in rows.items()
                if v.get("status") == "ok" and v.get("cohort") == "EXT"}

    result = dict(
        main_cohort=dict(n_sections=len(main_rows),
                         aggregate={op: summarize(main_rows, op)
                                    for op in ("cell", "mab")}),
        external_arm=dict(n_sections=len(ext_rows),
                          sections=sorted(ext_rows),
                          aggregate={op: summarize(ext_rows, op)
                                     for op in ("cell", "mab")}),
        per_section=rows,
        settings=dict(k_fracs=list(K_FRACS), n_random=N_RANDOM,
                      ridge_alpha=RIDGE_ALPHA, ridge_folds=RIDGE_FOLDS),
        semantics=(
            "Ground truth = per-spot counterfactual delta from run_44 (full "
            "recompute for every spot). It is NOT used to build any candidate "
            "ranking except 'sparta_delta', which is drawn as a reference "
            "ceiling, not as a competitor. 'achieved_x pct' = the barrier "
            "reduction obtained by REALLY ablating that method's top-k spots as "
            "one set, divided by the scale anchor (full_cut_breach for B_cell, "
            "all_removed_mab for B_mAb). Higher = better target selection. "
            "Sign convention for simple baselines is fixed a priori on "
            "biological grounds (near boundary / far from vessel / dense matrix "
            "= more barrier), never chosen by looking at the ground truth."),
        meta=stamp_run(cfg, {"module": "M46-simple-baselines"}),
    )
    save_json(str(json_path), result)
    print(f"\nWrote {json_path}")

    figs = _make_figures(result["main_cohort"]["aggregate"],
                         len(main_rows), P, cfg)
    for f in figs:
        print(f"Wrote {f}")

    # 控制台摘要
    print("\n=== main cohort: achieved barrier reduction at 5% budget "
          "(median, fraction of anchor) ===")
    for op, nm in (("cell", "B_cell"), ("mab", "B_mAb")):
        agg = result["main_cohort"]["aggregate"][op]["achieved_5pct"]
        print(f"  {nm}:")
        for mname, label, _c, _l in METHODS:
            if mname in agg:
                print(f"    {label:<42} {agg[mname]['median']:+.3f}")

    print(f"\ntotal {time.perf_counter()-t_all:.0f}s")


if __name__ == "__main__":
    main()
