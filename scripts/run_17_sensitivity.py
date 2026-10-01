#!/usr/bin/env python
"""
run_17_sensitivity.py —— B_mAb 参数敏感性网格（多切片版，收编自 scratch/_sens_bmab_grid.py）
================================================================================================

输入：{sid}.scored.h5ad + {sid}.graph.npz
输出：results/validation/bmab_sensitivity.json（多切片）
      results/figures/bmab_sensitivity_{grid}_{sid}.png（逐片热图）
      results/figures/bmab_sensitivity_beta_curves.png（跨切片汇总曲线）
上游模块：run_03_graph.py
下游模块：run_16_figures.py（Fig 2b 可改用跨切片曲线）

相对 scratch 版改了四处
-----------------------
1. **路径走 io_.Paths**，不再硬编码 `D:\\sparta`。
2. **支持多张切片**。旧版只跑 MEL01 一张，骨架却把它表述成全局稳健性——
   这是审查报告里点名的短板。现在逐片跑并给出跨片范围。
3. **"合理参数邻域"变成显式参数并打印出来**。旧版把这个窗口
   （lam∈[1,8]、beta∈[1.5,6]、xi0∈[10,40]）埋在注释里，
   而 "88–92% 的网格点判 GO" 这个数字完全取决于窗口怎么划——
   窗口是事后划的还是预先定的，审稿人一定会问。现在它必须被显式写出来。
4. **色标改对**。旧版交联占比（0–100% 的量级）用了发散色标 RdBu_r，
   解离潜力（真正有阈值 1.0 的极性量）反而用了 viridis。两个都反了：
     · 量级 -> 单色相 light→dark；
     · 极性 -> 双色相 + 中性灰中点（这里中点是 log10(1)=0，即 GO 阈值）。

用法
----
    python scripts/run_17_sensitivity.py --slides MEL01 MEL03 CSCC01 CSCC03
    python scripts/run_17_sensitivity.py --slides MEL01 --grids lam_x_beta
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import Paths, load_config, load_graph, save_json, set_seed, stamp_run  # noqa: E402

MAB, INK, INK2, SURFACE = "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"
CELL = "#2a78d6"

GRIDS = {
    "lam_x_beta": dict(
        x_name=r"$\lambda$  (ECM decay)", y_name=r"$\beta$  (mesh contraction)",
        xs=[0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0],
        ys=[0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.5, 6.0, 9.0, 12.0],
        key=("lam", "beta")),
    "xi0_x_beta": dict(
        x_name=r"$\xi_0$  (baseline mesh, nm)", y_name=r"$\beta$  (mesh contraction)",
        xs=[8, 12, 16, 20, 30, 40, 60, 80],
        ys=[0.5, 1.0, 1.5, 2.0, 3.0, 4.5, 6.0, 9.0, 12.0],
        key=("xi0_nm", "beta")),
    "kappa_x_kd": dict(
        x_name=r"$\kappa_w$  (binding-site weight)", y_name=r"$K_{d,eff}$",
        xs=[0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 6.0, 10.0],
        ys=[0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0],
        key=("kappa_w", "kd_eff")),
}


def _seq(hex_hi):
    """单色相 light->dark（量级用）。"""
    from matplotlib.colors import LinearSegmentedColormap, to_rgb
    r, g, b = to_rgb(hex_hi)
    return LinearSegmentedColormap.from_list(
        "seq", [(1 - .96 * (1 - r), 1 - .96 * (1 - g), 1 - .96 * (1 - b)),
                (r, g, b), (r * .42, g * .42, b * .42)])


def _div():
    """双色相 + 中性灰中点（极性用；中点 = GO 阈值 log10(1)=0）。"""
    from matplotlib.colors import LinearSegmentedColormap, to_rgb
    lo, hi = to_rgb(CELL), to_rgb(MAB)
    return LinearSegmentedColormap.from_list(
        "div", [(lo[0] * .45, lo[1] * .45, lo[2] * .45), lo,
                (0.90, 0.895, 0.88), hi, (hi[0] * .45, hi[1] * .45, hi[2] * .45)])


def main():
    ap = argparse.ArgumentParser(
        description="B_mAb 参数敏感性网格（多切片）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slides", nargs="+", required=True,
                    help="**至少给 3 张、且跨两个队列**，否则不能称为全局稳健性")
    ap.add_argument("--grids", nargs="+", default=list(GRIDS), choices=list(GRIDS))
    ap.add_argument("--lam-window", nargs=2, type=float, default=[1.0, 8.0],
                    help="'合理邻域'的 lam 范围。**预先指定，不要看完结果再划**")
    ap.add_argument("--beta-window", nargs=2, type=float, default=[1.5, 6.0])
    ap.add_argument("--xi0-window", nargs=2, type=float, default=[10.0, 40.0])
    ap.add_argument("--config", default=None)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    try:
        import scanpy as sc
    except ImportError:
        sys.exit("需要 scanpy 读取签名分数。")

    from sparta.barrier import scores_from_adata
    from sparta.validate import dissociation_drivers

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    cfg_cell = cfg["barrier"]["b_cell"]
    base_mab = {k: v for k, v in cfg["barrier"]["b_mab"].items() if k != "r_nm"}
    r_nm = cfg["barrier"]["b_mab"]["r_nm"]

    print("=" * 78)
    print(f"B_mAb 敏感性网格　{len(args.slides)} 张切片　网格 {args.grids}")
    print(f"'合理邻域'（预先指定）：lam {args.lam_window}　beta {args.beta_window}"
          f"　xi0 {args.xi0_window}")
    print("  ⚠ 这个窗口决定了'X% 的网格点判 GO'这个数字。窗口必须写进 Methods，")
    print("    并说明它是预先定的——事后按结果划窗口是隐性的挑选。")
    print("=" * 78)

    if len(args.slides) < 3:
        print(f"\n[warn] 只给了 {len(args.slides)} 张切片。少于 3 张时，"
              f"正文只能写'以 {args.slides[0]} 为例'，不能称为全局稳健性。\n")

    out = {"slides": {}, "base_params": base_mab,
           "windows": dict(lam=args.lam_window, beta=args.beta_window,
                           xi0_nm=args.xi0_window),
           "meta": stamp_run(cfg, {"module": "M17-sensitivity"})}

    for sid in args.slides:
        try:
            adata = sc.read_h5ad(P.scored(sid))
            A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        except FileNotFoundError as e:
            print(f"{sid}: 跳过（{e}）")
            continue
        S = scores_from_adata(adata)

        def _run(mab_cfg):
            dd = dissociation_drivers(A, S["ecm"], S["caf"], S["crosslink"],
                                      S["ag_target"], source, vessel,
                                      cfg_cell=cfg_cell, cfg_mab=mab_cfg, r_nm=r_nm)
            return dd["frac_crosslink"] * 100, dd["dissociation_potential"]

        res = {}
        for gname in args.grids:
            g = GRIDS[gname]
            kx, ky = g["key"]
            M = np.zeros((len(g["ys"]), len(g["xs"])))
            Pot = np.zeros_like(M)
            t0 = time.time()
            for i, yv in enumerate(g["ys"]):
                for j, xv in enumerate(g["xs"]):
                    m = dict(base_mab); m[kx] = xv; m[ky] = yv
                    M[i, j], Pot[i, j] = _run(m)
                print(f"  {sid} {gname}: {i+1}/{len(g['ys'])} 行"
                      f"（{time.time()-t0:.0f}s）", end="\r", flush=True)
            print(f"  {sid} {gname}: {M.size} 个网格点，{time.time()-t0:.0f}s" + " " * 20)

            # 预先指定的窗口
            win = {"lam": args.lam_window, "beta": args.beta_window,
                   "xi0_nm": args.xi0_window}
            jj = [j for j, v in enumerate(g["xs"])
                  if kx not in win or win[kx][0] <= v <= win[kx][1]]
            ii = [i for i, v in enumerate(g["ys"])
                  if ky not in win or win[ky][0] <= v <= win[ky][1]]
            sub, psub = M[np.ix_(ii, jj)], Pot[np.ix_(ii, jj)]
            res[gname] = dict(
                grid=g, crosslink_pct=M.tolist(), potential=Pot.tolist(),
                window=dict(crosslink_pct_min=float(sub.min()),
                            crosslink_pct_max=float(sub.max()),
                            potential_min=float(psub.min()),
                            potential_max=float(psub.max()),
                            n_cells=int(sub.size),
                            n_above_GO=int((psub > 1.0).sum())),
                full=dict(crosslink_pct_min=float(M.min()),
                          crosslink_pct_max=float(M.max()),
                          n_cells=int(M.size),
                          n_above_GO=int((Pot > 1.0).sum())))
            w, f = res[gname]["window"], res[gname]["full"]
            print(f"      预设窗口内：交联% {w['crosslink_pct_min']:.1f}–"
                  f"{w['crosslink_pct_max']:.1f}　GO {w['n_above_GO']}/{w['n_cells']}")
            print(f"      整个网格上：交联% {f['crosslink_pct_min']:.1f}–"
                  f"{f['crosslink_pct_max']:.1f}　GO {f['n_above_GO']}/{f['n_cells']}"
                  f"　← 正文必须同时报这一行")

            if not args.no_plot:
                _heatmaps(P, sid, gname, g, M, Pot)
        out["slides"][sid] = res

    if not out["slides"]:
        sys.exit("没有任何切片可跑。")

    # ---- 跨切片汇总 ----
    print("\n" + "=" * 78)
    print("跨切片汇总（预设窗口内）")
    print("=" * 78)
    for gname in args.grids:
        vals = [out["slides"][s][gname]["window"] for s in out["slides"]
                if gname in out["slides"][s]]
        if not vals:
            continue
        go = sum(v["n_above_GO"] for v in vals); tot = sum(v["n_cells"] for v in vals)
        print(f"  {gname:<14} 交联% {min(v['crosslink_pct_min'] for v in vals):.1f}–"
              f"{max(v['crosslink_pct_max'] for v in vals):.1f}"
              f"　GO {go}/{tot}（{go/tot*100:.0f}%）")

    if not args.no_plot and "lam_x_beta" in args.grids and len(out["slides"]) > 1:
        _beta_curves(P, out, args.grids)

    p = P.validation("bmab_sensitivity.json")
    save_json(p, out)
    print(f"\n已写出 {p}")
    print("提醒：run_16_figures.py 的 Fig 2b 目前读的是单切片曲线，"
          "跑完多切片后可以改成读 bmab_sensitivity_beta_curves 那一版。")


def _heatmaps(P, sid, gname, g, M, Pot):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), facecolor=SURFACE)
    # 左：量级 -> 单色相
    im0 = ax[0].imshow(M, origin="lower", aspect="auto", cmap=_seq(MAB),
                       vmin=0, vmax=100)
    ax[0].set_title(f"crosslinking share of B$_{{mAb}}$ variance (%)  —  {sid}",
                    fontsize=10, color=INK, loc="left")
    fig.colorbar(im0, ax=ax[0])
    # 右：极性（阈值 1.0）-> 双色相 + 中性中点
    L = np.log10(np.maximum(Pot, 1e-2))
    lim = float(max(abs(L.min()), abs(L.max()), 0.5))
    im1 = ax[1].imshow(L, origin="lower", aspect="auto", cmap=_div(),
                       norm=TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim))
    ax[1].set_title(f"log$_{{10}}$ dissociation potential  (0 = GO threshold)  —  {sid}",
                    fontsize=10, color=INK, loc="left")
    fig.colorbar(im1, ax=ax[1])
    for a, Z, fmt in ((ax[0], M, "{:.0f}"), (ax[1], Pot, "{:.1f}")):
        a.set_xticks(range(len(g["xs"]))); a.set_xticklabels(g["xs"], fontsize=7.4)
        a.set_yticks(range(len(g["ys"]))); a.set_yticklabels(g["ys"], fontsize=7.4)
        a.set_xlabel(g["x_name"], fontsize=8.6, color=INK2)
        a.set_ylabel(g["y_name"], fontsize=8.6, color=INK2)
        for i in range(Z.shape[0]):
            for j in range(Z.shape[1]):
                a.text(j, i, fmt.format(Z[i, j]), ha="center", va="center",
                       fontsize=5.8, color=INK)
    fig.tight_layout()
    p = P.figure(f"bmab_sensitivity_{gname}_{sid}.png")
    fig.savefig(p, dpi=150, facecolor=SURFACE); plt.close(fig)
    print(f"      已保存 {p.name}")


def _beta_curves(P, out, grids):
    """跨切片的 交联% vs beta 曲线——Fig 2b 的多切片版。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, a = plt.subplots(figsize=(6.4, 4.4), facecolor=SURFACE)
    a.set_facecolor(SURFACE)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
    a.grid(axis="y", color="#eeece6", lw=0.8); a.set_axisbelow(True)
    for sid, res in out["slides"].items():
        g = res["lam_x_beta"]["grid"]
        M = np.asarray(res["lam_x_beta"]["crosslink_pct"], float)
        j = g["xs"].index(3.0)
        a.plot(g["ys"], M[:, j], lw=1.8, alpha=0.9, color=MAB,
               marker="o" if sid.startswith("CSCC") else "s", ms=4.5,
               markerfacecolor=MAB if sid.startswith("CSCC") else SURFACE,
               markeredgecolor=MAB, markeredgewidth=1.2)
        a.annotate(sid, xy=(g["ys"][-1], M[-1, j]), xytext=(5, 0),
                   textcoords="offset points", fontsize=7, color=INK2, va="center")
    a.axvline(3.0, color="#8f8e88", lw=0.9, ls=(0, (3, 3)))
    a.set_xscale("log"); a.set_xticks(g["ys"])
    a.set_xticklabels([f"{y:g}" for y in g["ys"]], fontsize=7.4)
    a.set_xlabel(r"$\beta$  (mesh-contraction steepness)", fontsize=8.6, color=INK2)
    a.set_ylabel("crosslinking share (%)", fontsize=8.6, color=INK2)
    a.tick_params(labelsize=7.6, colors=INK2, length=3)
    a.set_title("The crosslinking share is a function of $\\beta$ in every section",
                fontsize=10.5, color=INK, loc="left", pad=10)
    fig.tight_layout()
    p = P.figure("bmab_sensitivity_beta_curves.png")
    fig.savefig(p, dpi=300, facecolor=SURFACE); plt.close(fig)
    print(f"  已保存 {p.name}")


if __name__ == "__main__":
    main()
