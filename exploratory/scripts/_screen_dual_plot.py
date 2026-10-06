# 英文双队列决策图（文章 Fig 1）
import sys
sys.path.insert(0, r"D:\sparta")
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sparta.io_ import Paths

P = Paths(__import__("sparta.io_", fromlist=["load_config"]).load_config(None))
d = json.load(open(r"D:\sparta\results\validation\screen_decision.json", encoding="utf-8"))
rows = d["per_slide"]
sids = list(rows)
cohort = {s: ("Melanoma\n(metastasis)" if s.startswith("MEL") else "cSCC\n(primary)")
          for s in sids}
cmap = {"Melanoma\n(metastasis)": "#C0392B", "cSCC\n(primary)": "#2E8B57"}

fig, ax = plt.subplots(1, 2, figsize=(13.5, 4.8), gridspec_kw={"width_ratios": [1.6, 1]})
# 左：驱动分解堆叠
bottom = np.zeros(len(sids))
for key, lab, c in (("frac_shared_ecm", "shared ECM (both barriers)", "#8FA3B8"),
                    ("frac_antigen", "antigen BSB (mAb only)", "#C0392B"),
                    ("frac_crosslink", "crosslink/size (mAb only)", "#E67E22")):
    v = np.array([rows[s][key] for s in sids])
    ax[0].bar(sids, v, bottom=bottom, label=lab, color=c)
    bottom += v
ax[0].set_ylabel("fraction of B_mAb variance")
ax[0].set_title("Driver decomposition of B_mAb variance")
ax[0].legend(fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=3)
for i, s in enumerate(sids):
    pass  # 潜力值已由右图展示，避免柱内文字重叠
# 队列分隔线
for x in [3.5]:
    ax[0].axvline(x + 0.5, color="gray", ls="--", lw=1)
ax[0].text(1.5, 1.02, "Melanoma meta.", ha="center", fontsize=9, color="#C0392B", weight="bold")
ax[0].text(5.5, 1.02, "cSCC primary", ha="center", fontsize=9, color="#2E8B57", weight="bold")
ax[0].set_ylim(0, 1.15)
ax[0].tick_params(axis="x", rotation=45)

# 右：解离潜力
pots = np.array([rows[s]["dissociation_potential"] for s in sids])
cols = [cmap[cohort[s]] for s in sids]
ax[1].bar(sids, pots, color=cols)
ax[1].axhline(1.0, color="#2E8B57", ls="--", lw=1.2, label="GO threshold = 1.0")
ax[1].axhline(0.3, color="#C0392B", ls=":", lw=1.2, label="PIVOT threshold = 0.3")
ax[1].set_ylabel("dissociation potential ratio")
ax[1].set_title(f"Decision: GO  (median {d['median_potential']:.1f}, n={len(sids)})")
ax[1].legend(fontsize=8)
ax[1].tick_params(axis="x", rotation=45)
fig.tight_layout()
p = P.figure("screen_decision_dual_cohort.png")
fig.savefig(p, dpi=150); plt.close(fig)
print("saved", p)
