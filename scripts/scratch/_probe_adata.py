# -*- coding: utf-8 -*-
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import scanpy as sc
from sparta.io_ import Paths, load_config
cfg = load_config(None); P = Paths(cfg)
adata = sc.read_h5ad(P.scored("MEL01"))
print("obs cols:", [c for c in adata.obs.columns if not c.startswith("_")][:40])
print("---")
for g in ["CD3D","CD8A","GZMB","CD4","MKI67","CD3E","NKG7","PRF1","PECAM1","VWF"]:
    print(f"  {g}: {'YES' if g in adata.var_names else 'no'}")
