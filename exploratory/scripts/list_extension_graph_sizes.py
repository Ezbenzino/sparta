import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
rows = list(csv.DictReader((ROOT / "data/ledger.csv").open(encoding="utf-8")))
out = []
for row in rows:
    if row["status"] != "extension":
        continue
    path = ROOT / "data/interim" / f"{row['slide_id']}.graph.npz"
    with np.load(path, allow_pickle=True) as z:
        meta = json.loads(str(z["meta"][0]))
    out.append((
        meta["n_nodes"],
        row["slide_id"],
        row["platform"],
        meta["largest_component_frac"],
    ))
for n, slide, platform, lcc in sorted(out, reverse=True):
    print(f"{slide:<20} {platform:<10} n={n:6d} LCC={lcc:.3f}")
