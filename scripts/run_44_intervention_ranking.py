#!/usr/bin/env python
"""
run_44_intervention_ranking.py —— 逐 spot 反事实干预排序（定位屏障）
====================================================================
对每张切片的**每一个 spot**做"完全通透"反事实（把该处 ecm/caf/crosslink
降到全片 5% 分位，与 S2 同口径），重算 B_cell 与 B_mAb，得到两张按
干预价值排序的表。回答："在这个模型里，打通哪些位置对屏障影响最大？"

读数（每个 spot 两列）：
  delta_b_cell : 干预后切片级最小割屏障的下降（绝对值）
  delta_b_mab  : 干预后 b_mab 场（reachable 均值）的下降（绝对值）

刻度锚点（让相对效应可跨切片比较）：
  full_cut_breach      : 整条最小割封锁带物质全去掉后的 B_cell 下降
  all_removed_mab      : 全片阻力物质去掉后的 B_mAb 下降

统计汇总：
  · top-5% 干预位点的空间富集（落在最小割上的比例 vs 本底）
  · 两个算子对"哪里最该打通"的秩相关（Spearman）与 top-5% 重叠（Jaccard）
  —— 后者直接把干预排序与耦合分解叙事接上：如果两个屏障共享封锁位置，
    联合干预的靶点就是同一个。

无 Scanpy 路线：全部输入来自 data/interim/{sid}.nodes.npz + {sid}.graph.npz。
BRCA01/02/LN01 为外部验证臂，汇总与主队列（15 cSCC + 4 MEL）分开报告，
其 sink 依赖的 B_cell 读数在稿件口径中为探索性（见 configs/default.yaml
admission_overrides.brca_10x_vis 的审查注记）。

输出：
  results/intervention/{sid}.intervention.npz   逐 spot 数组（含坐标）
  results/validation/intervention_ranking.json  汇总
  results/figures/intervention_ranking_{sid}.png  逐切片效应空间图
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

from sparta.io_ import (Paths, load_config, load_graph, load_json,  # noqa: E402
                        patient_map, save_json, stamp_run)
from sparta.node_tables import coords_from_nodes, load_nodes, scores_from_nodes  # noqa: E402
from sparta.intervention import (intervention_anchors, joint_intervention_curve,  # noqa: E402
                                spot_intervention_ranking)

TOP_LIST_K = 20          # JSON 里逐切片保留的 top 名单长度
FRACS = (0.01, 0.05)     # top-k% 汇总口径


def _cohort(sid: str) -> str:
    if sid.startswith("MEL_THR"):
        return "REP"  # melanoma replication cohort (Thrane et al. 2018)
    if sid.startswith("CSCC"):
        return "CSCC"
    if sid.startswith("MEL"):
        return "MEL"
    return "EXT"  # BRCA01 / BRCA02 / LN01（外部验证臂）


def _top_stats(delta: np.ndarray, anchor: float, cut_mask: np.ndarray,
               n_cut: int, n: int) -> dict:
    """单个算子列的排序统计。anchor 是该算子的刻度锚点下降量。"""
    d = np.asarray(delta, float)
    order = np.argsort(-d)
    out = dict(
        max=float(d[order[0]]),
        n_zero_effect=int(np.sum(d <= 0)),
        frac_zero_effect=float(np.mean(d <= 0)),
    )
    if np.isfinite(anchor) and anchor > 0:
        out["max_over_anchor"] = float(d[order[0]] / anchor)
    for f in FRACS:
        k = max(1, int(round(f * n)))
        topk = order[:k]
        out[f"top{int(f*100)}pct"] = dict(
            k=int(k),
            mean=float(np.mean(d[topk])),
            sum_indep=float(np.sum(d[topk])),
            sum_indep_over_anchor=(float(np.sum(d[topk]) / anchor)
                                   if np.isfinite(anchor) and anchor > 0 else float("nan")),
            in_cut_frac=float(np.mean(cut_mask[topk])) if n_cut else float("nan"),
            base_cut_frac=float(n_cut / n) if n else float("nan"),
        )
    return out, order


def _plot_section(sid: str, xy: np.ndarray, dcell: np.ndarray, dmab: np.ndarray,
                  cut_nodes: np.ndarray, n: int, path: Path,
                  joint: dict | None = None) -> None:
    """逐切片干预效应空间图：左 B_cell，右 B_mAb。

    配色标尺用 97 分位而非最大值：大多数 spot 的 Δ 接近 0（不在封锁线上），
    若按最大值归一化，热点会被压成一片黑，图就白画了。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    in_cut = np.zeros(n, bool)
    if len(cut_nodes):
        in_cut[cut_nodes] = True

    ann = ""
    if joint and "b_cell" in joint:
        pts = {p["k"]: p["joint_frac_of_anchor"] for p in joint["b_cell"]["points"]}
        if 20 in pts:
            ann = f"  |  joint top-20 spots = {pts[20]*100:.0f}% of a full breach"

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    for ax, d, ttl in ((axes[0], dcell, "ΔB_cell  (topological: opens the blockade)"),
                       (axes[1], dmab, "ΔB_mAb  (diffusive: distributed)")):
        pos = d[d > 0]
        vmax = float(np.percentile(pos, 97)) if len(pos) else 1.0
        vmax = max(vmax, 1e-12)
        s = ax.scatter(xy[:, 0], xy[:, 1], c=np.clip(d, 0, vmax), s=6,
                       cmap="magma", norm=Normalize(0, vmax), rasterized=True,
                       linewidths=0)
        if len(cut_nodes):
            ax.scatter(xy[in_cut, 0], xy[in_cut, 1], s=8, facecolors="none",
                       edgecolors="#4dd0e1", linewidths=0.5, alpha=0.85,
                       zorder=3, label=f"min-cut band (n={len(cut_nodes)})")
            ax.legend(loc="upper right", fontsize=7, frameon=False)
        ax.set_title(f"{ttl}\ncolour clipped at 97th pct", fontsize=9)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(s, ax=ax, shrink=0.8)
    fig.suptitle(f"{sid} — counterfactual intervention ranking (cyan = min-cut band){ann}",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _section_sids(P: Paths, include_replication: bool = False) -> list[str]:
    """Sections with a node table whose ledger status is "ingested" (primary + external arm).

    2026-10-05: the melanoma replication cohort (ledger status "replication") also has node tables;
    globbing every *.nodes.npz would silently pull it into the main-cohort summary. It is added
    only with ``include_replication`` and is then summarised as its own cohort ("REP").
    """
    from sparta.io_ import load_ledger
    wanted = {"ingested", "replication"} if include_replication else {"ingested"}
    ok = {r["slide_id"] for r in load_ledger(P) if str(r.get("status", "")).strip() in wanted}
    return sorted(s for s in (p.stem.split(".")[0] for p in Path(P.interim).glob("*.nodes.npz")) if s in ok)


def process_section(args) -> tuple[str, dict]:
    """处理单张切片：排序、落盘数组、出图。供串行/多进程复用。

    各切片彼此独立且完全确定，因此并行不改变任何数值结果。
    """
    sid, outdir_s, force = args[:3]
    refresh_cut = bool(args[3]) if len(args) > 3 else False
    outdir = Path(outdir_s)
    cfg = load_config(None)
    P = Paths(cfg)
    pmap = patient_map(P)
    cfg_cell = cfg["barrier"]["b_cell"]
    cfg_mab = cfg["barrier"]["b_mab"]

    arr_path = outdir / f"{sid}.intervention.npz"
    A_, S_, SRC_ = None, None, None
    if arr_path.exists() and not force and refresh_cut:
        # 2026-10-05: per-spot deltas depend only on max-flow values, which the exact-cut fix did not
        # change (< 1e-12); the cut geometry, the anchors and the joint curves are recomputed.
        z = np.load(arr_path, allow_pickle=True)
        dcell, dmab = z["delta_b_cell"], z["delta_b_mab"]
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        nodes = load_nodes(Path(P.interim) / f"{sid}.nodes.npz")
        S = scores_from_nodes(nodes)
        xy = coords_from_nodes(nodes)
        n = A.shape[0]
        anc = intervention_anchors(A, S["ecm"], S["caf"], S["crosslink"], S["ag_target"], source, sink,
                                   vessel, cfg_cell, cfg_mab, low_q=0.05)
        cut_nodes = anc["cut_nodes"]
        cut_mask_new = np.zeros(n, bool)
        cut_mask_new[cut_nodes] = True
        np.savez_compressed(arr_path, delta_b_cell=dcell, delta_b_mab=dmab, cut_mask=cut_mask_new,
                            coords_um=z["coords_um"], obs_names=z["obs_names"])
        res = dict(b_cell_base=anc["b_cell0"], b_mab_base=anc["b_mab0"],
                   full_cut_breach=anc["full_cut_breach"], all_removed_cell=anc["all_removed_cell"],
                   all_removed_mab=anc["all_removed_mab"], low_q=0.05)
        cached = False
        A_, S_, SRC_ = A, S, (source, sink, vessel)
    elif arr_path.exists() and not force:
        # 断点续跑：数组已存在则跳过重算，但仍返回该切片的汇总
        z = np.load(arr_path, allow_pickle=True)
        dcell, dmab = z["delta_b_cell"], z["delta_b_mab"]
        cut_nodes = np.flatnonzero(z["cut_mask"])
        n = len(dcell)
        xy = z["coords_um"]
        # 锚点无法从数组恢复 -> 若需完整汇总则用 --force 重算
        res = dict(b_cell_base=float("nan"), b_mab_base=float("nan"),
                   full_cut_breach=float("nan"), all_removed_cell=float("nan"),
                   all_removed_mab=float("nan"), low_q=float("nan"))
        cached = True
    else:
        t0 = time.perf_counter()
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        nodes = load_nodes(Path(P.interim) / f"{sid}.nodes.npz")
        S = scores_from_nodes(nodes)
        xy = coords_from_nodes(nodes)
        n = A.shape[0]
        res = spot_intervention_ranking(
            A, S["ecm"], S["caf"], S["crosslink"], S["ag_target"],
            source, sink, vessel, cfg_cell, cfg_mab)
        dcell = np.asarray(res["delta_b_cell"], float)
        dmab = np.asarray(res["delta_b_mab"], float)
        cut_nodes = np.asarray(res["cut_nodes"], int)
        cut_mask = np.zeros(n, bool)
        if len(cut_nodes):
            cut_mask[cut_nodes] = True
        np.savez_compressed(
            arr_path,
            delta_b_cell=dcell, delta_b_mab=dmab, cut_mask=cut_mask, coords_um=xy,
            obs_names=np.asarray([str(x) for x in nodes["obs_names"]]),
        )
        cached = False
        A_, S_, SRC_ = A, S, (source, sink, vessel)

    cut_nodes = np.asarray(cut_nodes, int)
    cut_mask = np.zeros(n, bool)
    if len(cut_nodes):
        cut_mask[cut_nodes] = True

    stat_cell, order_cell = _top_stats(dcell, res["full_cut_breach"], cut_mask,
                                       len(cut_nodes), n)
    stat_mab, order_mab = _top_stats(dmab, res["all_removed_mab"], cut_mask,
                                     len(cut_nodes), n)

    # 联合干预曲线：按排序同时打通 top-k（回答"打到第 k 个共撕开多大口子"）
    joint = {}
    if A_ is not None and SRC_ is not None:
        source_, sink_, vessel_ = SRC_
        ks = tuple(sorted({1, 5, 20, 50} | {max(1, int(round(f * n))) for f in (0.01, 0.02, 0.05, 0.10)}))
        joint["b_cell"] = joint_intervention_curve(
            A_, S_["ecm"], S_["caf"], S_["crosslink"], S_["ag_target"],
            source_, sink_, vessel_, cfg_cell, cfg_mab, order_cell,
            ks=ks, low_q=res["low_q"], operator="cell")
        joint["b_mab"] = joint_intervention_curve(
            A_, S_["ecm"], S_["caf"], S_["crosslink"], S_["ag_target"],
            source_, sink_, vessel_, cfg_cell, cfg_mab, order_mab,
            ks=ks, low_q=res["low_q"], operator="mab")

    if np.std(dcell) > 0 and np.std(dmab) > 0:
        rho = float(stats.spearmanr(dcell, dmab).statistic)
    else:
        rho = float("nan")
    k5 = max(1, int(round(0.05 * n)))
    t5c, t5m = set(order_cell[:k5].tolist()), set(order_mab[:k5].tolist())
    jacc = len(t5c & t5m) / max(len(t5c | t5m), 1)

    if not cached:
        try:
            _plot_section(sid, xy, dcell, dmab, cut_nodes, n,
                          Path(P.figure(f"intervention_ranking_{sid}.png")), joint)
        except Exception as e:  # noqa: BLE001 —— 图失败不应中断主统计
            print(f"  [warn] figure failed for {sid}: {e}")

    nodes_names = np.load(arr_path, allow_pickle=True)["obs_names"]

    def _top_list(order, d):
        return [dict(spot=str(nodes_names[i]), idx=int(i),
                     x_um=float(xy[i, 0]), y_um=float(xy[i, 1]),
                     delta=float(d[i]), in_cut=bool(cut_mask[i]))
                for i in order[:TOP_LIST_K]]

    row = dict(
        patient=pmap.get(sid, "?"),
        cohort=_cohort(sid),
        n=n,
        n_cut_nodes=int(len(cut_nodes)),
        b_cell_base=res["b_cell_base"],
        b_mab_base=res["b_mab_base"],
        anchors=dict(
            full_cut_breach=res["full_cut_breach"],
            all_removed_cell=res["all_removed_cell"],
            all_removed_mab=res["all_removed_mab"],
        ),
        b_cell=stat_cell,
        b_mab=stat_mab,
        operator_agreement=dict(
            spearman_delta_cell_vs_delta_mab=rho,
            top5pct_jaccard=float(jacc),
        ),
        top_b_cell=_top_list(order_cell, dcell),
        top_b_mab=_top_list(order_mab, dmab),
        joint_intervention=joint,
        low_q=res["low_q"],
        from_cache=bool(cached),
    )
    return sid, row


def main() -> None:
    import argparse
    import multiprocessing as mp

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--jobs", type=int, default=min(8, mp.cpu_count()),
                    help="并行进程数（切片彼此独立，并行不改变结果）")
    ap.add_argument("--force", action="store_true",
                    help="忽略已存在的逐切片 npz，全部重算")
    ap.add_argument("--refresh-cut", action="store_true",
                    help="keep cached per-spot deltas; recompute the exact cut, anchors and joint curves")
    ap.add_argument("--include-replication", action="store_true",
                    help="also rank the replication-cohort sections (summarised separately as REP)")
    ap.add_argument("--figures-only", action="store_true",
                    help="不重算，只用已保存的 npz + JSON 重绘逐切片图")
    args = ap.parse_args()

    cfg = load_config(None)
    P = Paths(cfg)
    outdir = Path(P.results) / "intervention"
    outdir.mkdir(parents=True, exist_ok=True)

    if args.figures_only:
        ranking = load_json(P.validation("intervention_ranking.json"))
        n_done = 0
        for sid, row in sorted(ranking["per_section"].items()):
            arr = outdir / f"{sid}.intervention.npz"
            if not arr.exists():
                continue
            z = np.load(arr, allow_pickle=True)
            try:
                _plot_section(sid, z["coords_um"], z["delta_b_cell"], z["delta_b_mab"],
                              np.flatnonzero(z["cut_mask"]), len(z["delta_b_cell"]),
                              Path(P.figure(f"intervention_ranking_{sid}.png")),
                              row.get("joint_intervention"))
                n_done += 1
            except Exception as e:  # noqa: BLE001
                print(f"  [warn] figure failed for {sid}: {e}")
        print(f"Redrew {n_done} figures from cached arrays")
        return

    sids = _section_sids(P, args.include_replication)
    todo = [(sid, str(outdir), args.force, args.refresh_cut) for sid in sids]
    print(f"Intervention ranking on {len(sids)} sections, jobs={args.jobs}: {sids}\n")

    rows = {}
    t_start = time.perf_counter()
    if args.jobs > 1:
        with mp.Pool(processes=args.jobs) as pool:
            for sid, row in pool.imap_unordered(process_section, todo):
                rows[sid] = row
                print(f"[done] {sid:<8} n={row['n']:<5} cut={row['n_cut_nodes']:<4} "
                      f"maxΔB_cell/full={row['b_cell'].get('max_over_anchor', float('nan')):.3f} "
                      f"ρ(Δ,Δ)={row['operator_agreement']['spearman_delta_cell_vs_delta_mab']:+.3f}",
                      flush=True)
    else:
        for job in todo:
            sid, row = process_section(job)
            rows[sid] = row
            print(f"[done] {sid:<8} n={row['n']:<5} cut={row['n_cut_nodes']:<4} "
                  f"maxΔB_cell/full={row['b_cell'].get('max_over_anchor', float('nan')):.3f} "
                  f"ρ(Δ,Δ)={row['operator_agreement']['spearman_delta_cell_vs_delta_mab']:+.3f}",
                  flush=True)
    print(f"\nAll sections done in {time.perf_counter()-t_start:.0f}s")

    # ---- 汇总：主队列 vs 外部臂 ----
    def _med(key, sidsel):
        vals = [rows[s][key] for s in sidsel
                if np.isfinite(rows[s].get(key, np.nan))]
        return float(np.median(vals)) if vals else float("nan")

    def _med_nested(path_keys, sidsel):
        vals, obj = [], rows
        for s in sidsel:
            try:
                v = rows[s]
                for k in path_keys:
                    v = v[k]
                vals.append(float(v))
            except (KeyError, TypeError, ValueError):
                pass
        return float(np.median(vals)) if vals else float("nan")

    main_sids = [s for s in rows if rows[s]["cohort"] in ("CSCC", "MEL")]
    ext_sids = [s for s in rows if rows[s]["cohort"] == "EXT"]

    def _med_joint(op, k_sel, sidsel):
        """联合曲线上某个 k 的 median(joint_frac_of_anchor)。k_sel 为 '5pct' 或整数。"""
        vals = []
        for s in sidsel:
            j = rows[s].get("joint_intervention", {}).get(op)
            if not j:
                continue
            if isinstance(k_sel, str) and k_sel.endswith("pct"):
                k = max(1, int(round(float(k_sel[:-3]) / 100 * rows[s]["n"])))
            else:
                k = int(k_sel)
            cand = [p for p in j.get("points", []) if p["k"] == k]
            if cand and np.isfinite(cand[0]["joint_frac_of_anchor"]):
                vals.append(float(cand[0]["joint_frac_of_anchor"]))
        return float(np.median(vals)) if vals else float("nan")

    def _max_joint(op, k_sel, sidsel):
        vals = []
        for s in sidsel:
            j = rows[s].get("joint_intervention", {}).get(op)
            if not j:
                continue
            k = max(1, int(round(float(k_sel[:-3]) / 100 * rows[s]["n"])))
            cand = [p for p in j.get("points", []) if p["k"] == k]
            if cand and np.isfinite(cand[0]["joint_frac_of_anchor"]):
                vals.append(float(cand[0]["joint_frac_of_anchor"]))
        return float(np.max(vals)) if vals else float("nan")

    def _cohort_summary(sidsel):
        return dict(
            n_sections=len(sidsel),
            n_patients=len({rows[s]["patient"] for s in sidsel}),
            median_n_cut_nodes=_med("n_cut_nodes", sidsel),
            median_max_delta_cell_over_full_breach=_med_nested(
                ["b_cell", "max_over_anchor"], sidsel),
            median_top5pct_sum_cell_over_full_breach=_med_nested(
                ["b_cell", "top5pct", "sum_indep_over_anchor"], sidsel),
            median_top5pct_in_cut_frac=_med_nested(
                ["b_cell", "top5pct", "in_cut_frac"], sidsel),
            median_base_cut_frac=_med_nested(
                ["b_cell", "top5pct", "base_cut_frac"], sidsel),
            median_joint_frac_cell_at_k20=_med_joint("b_cell", 20, sidsel),
            median_joint_frac_cell_at_top1pct=_med_joint("b_cell", "1pct", sidsel),
            max_joint_frac_cell_at_top1pct=_max_joint("b_cell", "1pct", sidsel),
            median_joint_frac_cell_at_top5pct=_med_joint("b_cell", "5pct", sidsel),
            median_max_delta_mab_over_all_removed=_med_nested(
                ["b_mab", "max_over_anchor"], sidsel),
            median_top5pct_sum_mab_over_all_removed=_med_nested(
                ["b_mab", "top5pct", "sum_indep_over_anchor"], sidsel),
            median_top5pct_in_cut_frac_mab=_med_nested(
                ["b_mab", "top5pct", "in_cut_frac"], sidsel),
            median_joint_frac_mab_at_k20=_med_joint("b_mab", 20, sidsel),
            median_joint_frac_mab_at_top1pct=_med_joint("b_mab", "1pct", sidsel),
            max_joint_frac_mab_at_top1pct=_max_joint("b_mab", "1pct", sidsel),
            median_joint_frac_mab_at_top5pct=_med_joint("b_mab", "5pct", sidsel),
            median_joint_frac_mab_at_top10pct=_med_joint("b_mab", "10pct", sidsel),
            median_frac_zero_effect_cell=_med_nested(
                ["b_cell", "frac_zero_effect"], sidsel),
            median_spearman_operator_agreement=_med_nested(
                ["operator_agreement", "spearman_delta_cell_vs_delta_mab"], sidsel),
            median_top5pct_jaccard=_med_nested(
                ["operator_agreement", "top5pct_jaccard"], sidsel),
        )

    rep_sids = [s for s in rows if rows[s]["cohort"] == "REP"]
    summary = dict(
        main_cohort=_cohort_summary(main_sids),
        external_arm=_cohort_summary(ext_sids),
        per_cohort={c: _cohort_summary([s for s in rows if rows[s]["cohort"] == c])
                    for c in ("CSCC", "MEL", "EXT", "REP") if any(rows[s]["cohort"] == c for s in rows)},
    )
    if rep_sids:
        summary["replication_cohort"] = _cohort_summary(rep_sids)

    out = dict(
        per_section=rows,
        summary=summary,
        semantics=(
            "Single-spot counterfactual within the SPARTA transport model: a spot's "
            "ecm/caf/crosslink are set to the section's 5th percentile "
            "(matrix-normalisation-like ablation, same convention as S2); "
            "ag_target (absorption) is untouched. delta_b_cell is the drop of the "
            "section-level min-cut barrier; delta_b_mab is the drop of the "
            "reachable-mean of the b_mab field. Rankings are model-defined "
            "counterfactuals, not measured permeability. "
            "READING THE NUMBERS: (1) every spot is tested, so the per-spot delta "
            "distribution IS the empirical null ('what a random spot achieves'), and no "
            "separate random control is needed. (2) 'sum_indep' inside top*k*pct sums "
            "INDEPENDENT single-spot drops and is therefore an upper bound (barriers are "
            "sub-additive); use joint_intervention (top-k ablated as one set) for the "
            "interpretable quantity. (3) full_cut_breach (B_cell) and all_removed_mab "
            "(B_mAb) are scale anchors, not experiments: they make the effect size "
            "comparable across sections and parameters. (4) B_mAb deltas can be very "
            "slightly negative for individual spots (absorption redistribution); see "
            "sparta/intervention.py."
        ),
        meta=stamp_run(cfg, {"module": "M44-intervention-ranking"}),
    )
    p = P.validation("intervention_ranking.json")
    save_json(p, out)
    print(f"\nWrote {p}")
    print(f"Arrays in {outdir}/; figures as results/figures/intervention_ranking_*.png")

    s = summary["main_cohort"]
    print("\n=== Main cohort (median across sections) ===")
    print(f"  top single spot        achieves {s['median_max_delta_cell_over_full_breach']*100:.1f}% "
          f"of a full min-cut breach (B_cell)")
    print(f"  joint top-20 spots     achieves {s['median_joint_frac_cell_at_k20']*100:.1f}% "
          f"of a full min-cut breach (B_cell)")
    print(f"  joint top-5% of spots  achieves {s['median_joint_frac_cell_at_top5pct']*100:.1f}% "
          f"of a full min-cut breach (B_cell)")
    print(f"  [upper bound] sum of independent top-5% deltas = "
          f"{s['median_top5pct_sum_cell_over_full_breach']*100:.0f}% (sub-additive; see semantics)")
    print(f"  top-5% B_cell spots inside min-cut band: "
          f"{s['median_top5pct_in_cut_frac']*100:.0f}% vs base "
          f"{s['median_base_cut_frac']*100:.0f}%")
    print(f"  B_mAb: top single spot   achieves {s['median_max_delta_mab_over_all_removed']*100:.2f}% "
          f"of the no-barrier limit; joint top-20 = "
          f"{s['median_joint_frac_mab_at_k20']*100:.1f}%")
    print(f"  top-5% B_mAb spots inside min-cut band: "
          f"{s['median_top5pct_in_cut_frac_mab']*100:.0f}% "
          f"(diffusive barrier is distributed, not a blockade line)")
    print(f"  zero-effect spots (B_cell): {s['median_frac_zero_effect_cell']*100:.1f}%")
    print(f"  operator agreement: rho(delta)={s['median_spearman_operator_agreement']:+.3f}, "
          f"top-5% Jaccard={s['median_top5pct_jaccard']:.2f}")


if __name__ == "__main__":
    main()
