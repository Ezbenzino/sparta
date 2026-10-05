"""Online Resource 1 figures (S1 maps, S3 domain scans, S4 radius, S5 reproduction, S6 runtime)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from scipy.stats import rankdata  # noqa: E402

import figstyle as S  # noqa: E402
from isdata import EXTERNAL_ORDER, PRIMARY_ORDER, VAL, j  # noqa: E402
from sparta.barrier import compute_b_cell_field, compute_b_mab  # noqa: E402
from sparta.io_ import Paths, load_config, load_graph  # noqa: E402
from sparta.node_tables import coords_from_nodes, load_nodes, scores_from_nodes  # noqa: E402

OUT = Path(__file__).resolve().parents[2] / "results" / "figures" / "is" / "supplement"
CMAP_CELL = LinearSegmentedColormap.from_list("cell", ["#eef5fd", "#b7d3f6", "#5598e7", "#1c5cab", "#0d366b"])
CMAP_MAB = LinearSegmentedColormap.from_list("mab", ["#fdf1ea", "#f6b896", "#eb6834", "#b8461b", "#6e2a10"])
CMAP_GREY = LinearSegmentedColormap.from_list("grey", ["#f4f3f0", "#c3c2b7", "#898781", "#52514e", "#1f1f1e"])


def fence(ax, xy, edges, half=50.0, lw=0.5):
    """Short bar across each cut edge; half-length in micrometres (fixed, so long diagonal
    edges of the first-generation ST lattice do not produce crossing bars)."""
    segs = []
    for (u, v) in edges:
        mid = xy[[u, v]].mean(axis=0)
        d = xy[v] - xy[u]
        nrm = np.array([-d[1], d[0]]) / (np.linalg.norm(d) + 1e-12) * half
        segs.append([mid - nrm, mid + nrm])
    ax.add_collection(LineCollection(segs, colors=S.INK, linewidths=lw, zorder=4))


def maps(sids, name, rows_per_page):
    cfg = load_config(None)
    P = Paths(cfg)
    sn = j("spatial_null_check.json")["per_slide"]
    ext = j("ext_validation.json")["per_slide"]
    ncol = 6
    nrow = int(np.ceil(len(sids) / 2))
    fig = plt.figure(figsize=(S.FULL_W, min(S.MAX_H, nrow * 37 * S.MM + 10 * S.MM)))
    for i, sid in enumerate(sids):
        nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
        A, D, src, snk, ves, _ = load_graph(P.graph(sid))
        Sx = scores_from_nodes(nodes)
        xy = coords_from_nodes(nodes)
        xy = xy - xy.min(axis=0)
        xy[:, 1] = xy[:, 1].max() - xy[:, 1]
        bcf = compute_b_cell_field(A, Sx["ecm"], Sx["caf"], src, **cfg["barrier"]["b_cell"])["b_cell_field"]
        bmf = compute_b_mab(A, Sx["ecm"], Sx["crosslink"], Sx["ag_target"], ves, **cfg["barrier"]["b_mab"])["b_mab"]
        cut = json.load(open(P.mincut(sid), encoding="utf-8"))["cut_edges"]
        rho = sn[sid]["real_rho_partial"] if sid in sn else ext[sid]["spatial_null"]["real_rho_partial"]
        pitch = 200.0 if (sid.startswith("CSCC") and sid not in ("CSCC01", "CSCC02", "CSCC03", "CSCC04")) else 100.0
        r_, c_ = divmod(i, 2)
        for k, (val, cm) in enumerate([(Sx["ecm"], CMAP_GREY), (rankdata(bcf) / len(bcf), CMAP_CELL),
                                       (rankdata(bmf) / len(bmf), CMAP_MAB)]):
            ax = fig.add_axes([0.005 + (c_ * 3 + k) * (0.99 / ncol), 1 - (r_ + 1) / nrow * 0.975,
                               0.99 / ncol - 0.006, 0.975 / nrow - 0.03])
            ax.set_aspect("equal")
            ax.axis("off")
            span = max(np.ptp(xy[:, 0]), np.ptp(xy[:, 1])) + 2 * pitch
            cx, cy = 0.5 * (xy[:, 0].min() + xy[:, 0].max()), 0.5 * (xy[:, 1].min() + xy[:, 1].max())
            ax.set_xlim(cx - span / 2, cx + span / 2)
            ax.set_ylim(cy - span / 2, cy + span / 2)
            fig.canvas.draw()
            bb = ax.get_window_extent()
            d_pt = pitch / span * bb.width / fig.dpi * 72 * 0.95
            ax.scatter(xy[:, 0], xy[:, 1], c=val, s=d_pt ** 2, cmap=cm, vmin=0, vmax=1, edgecolors="none",
                       rasterized=True)
            if k == 1:
                fence(ax, xy, cut, half=0.5 * pitch, lw=0.45)
            if k == 0:
                rtxt = (f"{rho:.2f}" if abs(rho) >= 0.005 else f"{rho:.4f}").replace("-", "\u2212")
                ax.text(0.0, 1.02, f"{sid}  ρ = {rtxt}", transform=ax.transAxes, fontsize=6.4, va="bottom", ha="left")
            elif r_ == 0:
                ax.text(0.5, 1.02, ["", "B$_{\\rm cell}$ field", "B$_{\\rm mAb}$ field"][k], transform=ax.transAxes,
                        fontsize=6.4, va="bottom", ha="center", color=S.INK2)
    S.save(fig, OUT, name, eps=False)


def domain_scans():
    nd = j("ndomains_sensitivity.json")["per_slide"]
    lam = j("benchmark_lambda_sensitivity.json")["per_slide"]
    fig, axs = plt.subplots(1, 4, figsize=(S.FULL_W, 52 * S.MM))
    grid = [4, 6, 8, 10, 12]
    for metric, axe, axp in [("n", axs[0], axs[1])]:
        pass
    for k, (dat, xs, xlabel, key_fmt) in enumerate([(nd, grid, "Number of domains", str),
                                                    (lam, [0.1, 0.3, 0.5, 0.7], "Neighbourhood weight λ", str)]):
        for sid in PRIMARY_ORDER:
            runs = dat[sid]["runs"]
            e = [runs[key_fmt(x)]["enrichment"] for x in xs]
            p = [100 * runs[key_fmt(x)]["precision_of_boundary_for_cut"] for x in xs]
            mk = S.SHAPE[S.platform_key(sid)]
            axs[2 * k].plot(xs, e, color=S.LIGHT, lw=0.6, marker=mk, ms=2, mfc=S.MUTED, mec="none")
            axs[2 * k + 1].plot(xs, p, color=S.LIGHT, lw=0.6, marker=mk, ms=2, mfc=S.MUTED, mec="none")
        e_med = [np.median([dat[s]["runs"][key_fmt(x)]["enrichment"] for s in PRIMARY_ORDER]) for x in xs]
        p_med = [np.median([100 * dat[s]["runs"][key_fmt(x)]["precision_of_boundary_for_cut"] for s in PRIMARY_ORDER]) for x in xs]
        axs[2 * k].plot(xs, e_med, color=S.INK, lw=1.4)
        axs[2 * k + 1].plot(xs, p_med, color=S.INK, lw=1.4)
        axs[2 * k].axhline(1.0, color=S.INK2, lw=0.6)
        axs[2 * k].set_xlabel(xlabel)
        axs[2 * k + 1].set_xlabel(xlabel)
        axs[2 * k].set_ylabel("Cut-edge enrichment\namong domain boundaries")
        axs[2 * k + 1].set_ylabel("Boundary edges that\nare cut edges (%)")
        axs[2 * k + 1].set_ylim(0, 40)
    for ax, l in zip(axs, "abcd"):
        S.panel(ax, l, x=-0.42, y=1.02)
    fig.tight_layout(w_pad=1.2)
    S.save(fig, OUT, "FigS3_domain_scans", eps=False)


def radius():
    rs = j("radius_sensitivity.json")["per_slide"]
    fig, ax = plt.subplots(figsize=(S.COL_W, 62 * S.MM))
    for sid in PRIMARY_ORDER:
        d = rs[sid]
        xs = sorted(float(k) for k in d)
        ys = [d[f"{x:.1f}"]["rho_partial"] if f"{x:.1f}" in d else d[str(x)]["rho_partial"] for x in xs]
        ax.plot(xs, ys, color=S.LIGHT, lw=0.7, marker=S.SHAPE[S.platform_key(sid)], ms=2.4, mfc=S.MUTED, mec="none")
    ax.axhline(0, color=S.INK2, lw=0.6)
    ax.set_xlabel("Graph radius (μm)")
    ax.set_ylabel("Partial Spearman ρ")
    S.save(fig, OUT, "FigS4_radius", eps=False)


def reproduction():
    rc = j("reproduction_check.json")["per_slide"]
    ns = j("spatial_null_nscore.json")["per_slide"] if (VAL / "spatial_null_nscore.json").exists() else None
    sn = j("spatial_null_check.json")["per_slide"]
    fig, axs = plt.subplots(1, 3, figsize=(S.FULL_W, 56 * S.MM))
    a = np.array([[rc[s]["locked_real_rho"], rc[s]["real_rho"]] for s in PRIMARY_ORDER])
    axs[0].plot([-0.05, 0.45], [-0.05, 0.45], color=S.LIGHT, lw=0.8)
    axs[0].plot(a[:, 0], a[:, 1], "o", ms=3, color=S.INK)
    axs[0].set_xlabel("Stored ρ (original analysis)")
    axs[0].set_ylabel("Recomputed ρ (node tables)")
    b = np.array([[rc[s]["locked_null_std"], rc[s]["null_std"]] for s in PRIMARY_ORDER])
    lim = [0.03, 0.27]
    axs[1].plot(lim, lim, color=S.LIGHT, lw=0.8)
    axs[1].plot(b[:, 0], b[:, 1], "o", ms=3, color=S.INK)
    axs[1].set_xlabel("Stored surrogate-null SD")
    axs[1].set_ylabel("Recomputed surrogate-null SD")
    if ns:
        c = np.array([[sn[s]["empirical_p_one_sided"], ns[s]["empirical_p_one_sided"]] for s in PRIMARY_ORDER])
        axs[2].plot([1e-3, 1], [1e-3, 1], color=S.LIGHT, lw=0.8)
        axs[2].plot(c[:, 0], c[:, 1], "o", ms=3, color=S.INK)
        axs[2].set_xscale("log")
        axs[2].set_yscale("log")
        axs[2].axvline(0.05, color=S.LIGHT, lw=0.6)
        axs[2].axhline(0.05, color=S.LIGHT, lw=0.6)
        axs[2].set_xlabel("p, raw-spectrum surrogates")
        axs[2].set_ylabel("p, normal-score surrogates")
    for ax, l in zip(axs, "abc"):
        S.panel(ax, l, x=-0.35, y=1.02)
    fig.tight_layout(w_pad=1.5)
    S.save(fig, OUT, "FigS5_reproduction", eps=False)


def runtime():
    rt = j("runtime_benchmark.json")["per_slide"]
    fig, ax = plt.subplots(figsize=(S.COL_W, 62 * S.MM))
    n = np.array([rt[s]["n_nodes"] for s in PRIMARY_ORDER])
    for key, lab, mk, col in [("t_sparta_core_s", "SPARTA core", "o", S.INK),
                              ("t_banksy_style_s", "BANKSY-style domains", "s", S.MUTED),
                              ("t_squidpy_s", "Squidpy enrichment", "^", S.LIGHT)]:
        t = np.array([rt[s][key] for s in PRIMARY_ORDER])
        ax.plot(n, t, mk, ms=3.4, color=col, label=lab, ls="none")
    ax.set_yscale("log")
    ax.set_xlabel("Spots per section")
    ax.set_ylabel("Wall time (s)")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.0), fontsize=6.6)
    S.save(fig, OUT, "FigS6_runtime", eps=False)


if __name__ == "__main__":
    S.apply()
    which = sys.argv[1:] or ["maps", "domains", "radius", "repro", "runtime"]
    if "maps" in which:
        maps(PRIMARY_ORDER[:12], "FigS1a_maps", 6)
        maps(PRIMARY_ORDER[12:] + EXTERNAL_ORDER, "FigS1b_maps", 5)
    if "domains" in which:
        domain_scans()
    if "radius" in which:
        radius()
    if "repro" in which:
        reproduction()
    if "runtime" in which:
        runtime()
