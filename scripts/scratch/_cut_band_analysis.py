# 割带特征正式产物：S1 不显著切片的机制解释
# 2026-08-28 扩到 19 张（8 Visium + 11 legacy）：
#   - S1 z 不再硬编码，从 results/counterfactual/{sid}.json 的 s1.fixed.z 读；
#   - Ripley 半径按平台缩放（Visium 250 μm / legacy 500 μm，同为 2.5 个点间距）；
#   - 队列改用形状 + 分组 + 直接标注（项目出图约定：颜色留给屏障模态）。
import json
import sys

sys.path.insert(0, r"D:\sparta")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.spatial import cKDTree

from sparta.io_ import Paths, load_config, load_graph
from sparta.barrier import scores_from_adata

cfg = load_config(None)
P = Paths(cfg)
VISIUM = ["MEL01", "MEL02", "MEL03", "MEL04",
          "CSCC01", "CSCC02", "CSCC03", "CSCC04"]
LEGACY = ["CSCC05", "CSCC06", "CSCC07", "CSCC08", "CSCC09", "CSCC10",
          "CSCC11", "CSCC12", "CSCC14", "CSCC15", "CSCC16"]
slides = VISIUM + LEGACY
PITCH = {**{s: 100.0 for s in VISIUM}, **{s: 200.0 for s in LEGACY}}

def ripley(coords, idx, r):
    pts = coords[idx]; n = len(pts)
    if n < 2: return 0.0
    xmin, ymin = coords.min(0); xmax, ymax = coords.max(0)
    area = max((xmax-xmin)*(ymax-ymin), 1.0)
    tree = cKDTree(pts)
    K = area * tree.count_neighbors(tree, r=r) / (n*(n-1))
    return float(np.sqrt(max(K/np.pi, 0)) - r)

out = {}
for sid in slides:
    adata = sc.read_h5ad(P.scored(sid))
    A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
    S = scores_from_adata(adata)
    coords = np.asarray(adata.obsm["spatial_um"])
    npz = np.load(P.interim / f"{sid}.barrier.npz", allow_pickle=True)
    cut = np.asarray(npz["cut_nodes"], int)
    n = A.shape[0]
    caf = S["caf"]
    is_high = caf >= np.quantile(caf, 0.85)
    z = float(json.load(open(P.counterfactual(sid), encoding="utf-8"))["s1"]["fixed"]["z"])
    out[sid] = dict(
        s1_z=z, n_cut=int(len(cut)), frac_cut=float(len(cut)/n),
        caf_on_cut=float(np.mean(caf[cut])),
        high_caf_on_cut_pct=float(is_high[cut].mean()*100),
        enrichment=float(is_high[cut].mean()/max(is_high.mean(), 1e-9)),
        compact_um=float(np.std(coords[cut], axis=0).mean()),
        ripley_r_um=PITCH[sid]*2.5,
        high_caf_ripleyL=ripley(coords, np.where(is_high)[0], PITCH[sid]*2.5),
        platform=("visium" if sid in VISIUM else "legacy_st"),
        mechanism=("signal-sparse narrow band" if len(cut)/n < 0.15
                   else "non-CAF-dominated diffuse band" if is_high[cut].mean() < 0.20
                   else "continuous CAF band"))

json.dump(out, open(P.validation("cut_band_analysis.json"), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)

fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
MARK = {"visium": "o", "legacy_st": "^"}
for grp, lab in [(VISIUM, "Visium"), (LEGACY, "first-gen ST")]:
    for sid in grp:
        r = out[sid]
        ax[0].scatter(r["frac_cut"]*100, r["s1_z"], marker=MARK[r["platform"]],
                      color="#3E6B8A", edgecolor="#0b0b0b", linewidth=.4,
                      s=70, alpha=.85)
        ax[0].annotate(sid, (r["frac_cut"]*100, r["s1_z"]), fontsize=7,
                       textcoords="offset points", xytext=(4, 3))
    ax[0].scatter([], [], marker=MARK[out[grp[0]]["platform"]], color="#3E6B8A",
                  edgecolor="#0b0b0b", linewidth=.4, label=f"{lab} (n={len(grp)})")
ax[0].axhline(2.0, color="gray", ls="--", lw=1, label="z = 2")
ax[0].set_xlabel("cut-band size (% of spots)")
ax[0].set_ylabel("S1-fixed z-score")
ax[0].set_title("S1 significance vs cut-band width (19 sections)")
ax[0].legend(fontsize=8)

mechs = {}
for sid in slides:
    mechs.setdefault(out[sid]["mechanism"], []).append(sid)
ax[1].axis("off")
y = 0.92
ax[1].set_title("Mechanism of S1-insensitive barriers", fontsize=10)
for m, sl in mechs.items():
    zs = [out[s]["s1_z"] for s in sl]
    ax[1].text(0.05, y, f"[{', '.join(sl)}]\n{m}  (z={min(zs):+.1f}..{max(zs):+.1f})",
               fontsize=7.5, va="top", transform=ax[1].transAxes)
    y -= 0.12 + 0.05 * (len(sl) > 8)
ax[1].text(0.05, y-0.10, "S1-fixed insensitive does NOT mean no barrier;\nit means the barrier is NOT a strict\ncontinuous CAF ring in these slides.",
           fontsize=8, style="italic", color="#555", va="top", transform=ax[1].transAxes)
fig.tight_layout()
fig.savefig(P.figure("cut_band_analysis.png"), dpi=150); plt.close(fig)
print("saved cut_band_analysis.png + cut_band_analysis.json")
for sid in slides:
    r = out[sid]
    print(f"{sid}: frac_cut={r['frac_cut']*100:.1f}%  highCAF_on_cut={r['high_caf_on_cut_pct']:.1f}%  "
          f"enrich={r['enrichment']:.2f}  z={r['s1_z']:+.1f}  -> {r['mechanism']}")
