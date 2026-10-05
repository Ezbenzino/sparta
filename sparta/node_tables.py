"""
Node tables: scanpy-free access to the per-spot model inputs
============================================================
``scripts/run_35_export_node_tables.py`` writes ``data/interim/{sid}.nodes.npz``
from the M2 output.  This module reads them back into exactly the dictionaries
that ``sparta.barrier`` consumes, so every downstream statistic can be recomputed
with numpy / scipy / networkx only.

Row order equals the node order of ``{sid}.graph.npz`` (both follow adata.obs).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .barrier import _SCORE_KEYS

__all__ = ["load_nodes", "scores_from_nodes", "coords_from_nodes", "gene_from_nodes"]


def load_nodes(path) -> dict:
    """Load a node table into a plain dict (arrays + parsed meta)."""
    z = np.load(Path(path), allow_pickle=True)
    out = {k: z[k] for k in z.files if k != "meta"}
    out["meta"] = json.loads(str(z["meta"][0])) if "meta" in z.files else {}
    return out


def scores_from_nodes(nodes: dict, keys: dict | None = None,
                      missing_out: list | None = None) -> dict:
    """Same contract as ``barrier.scores_from_adata`` (neutral 0.5 fill + record)."""
    keys = keys or _SCORE_KEYS
    n = len(nodes["obs_names"])
    out = {}
    for k, col in keys.items():
        arr = nodes.get(f"obs__{col}")
        if arr is not None:
            out[k] = np.asarray(arr, float)
        else:
            out[k] = np.full(n, 0.5)
            if missing_out is not None:
                missing_out.append(k)
    return out


def coords_from_nodes(nodes: dict) -> np.ndarray:
    """Micrometre coordinates (obsm['spatial_um'])."""
    return np.asarray(nodes["obsm__spatial_um"], float)


def gene_from_nodes(nodes: dict, gene: str, layer: str = "logcp10k"):
    """One marker gene from the exported panel (None if not measured)."""
    genes = [str(g) for g in nodes.get("genes", [])]
    if gene not in genes:
        return None
    arr = nodes["gene_logcp10k"] if layer == "logcp10k" else nodes["gene_counts"]
    return np.asarray(arr[:, genes.index(gene)], float)
