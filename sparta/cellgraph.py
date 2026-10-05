"""
Cell-resolution graphs for SPARTA (single-cell spatial proteomics / imaging data)
================================================================================
Input  : per-cell centroids (µm), cell-type labels and per-cell marker intensities
Output : the same arrays that ``sparta.barrier`` consumes (A, D, source, sink, inputs)
Upstream  : any segmented multiplexed image (here: CODEX, Schürch et al. 2020)
Downstream: ``barrier.compute_b_cell`` / ``compute_b_cell_field`` (unchanged operators),
            ``scripts/run_47_codex_validation.py``

Why a separate module
---------------------
The section pipeline builds a radius graph over spots whose inputs are expression scores.
At single-cell resolution two things change and nothing else:

1. the graph is a Delaunay triangulation of cell centroids (cells are not on a lattice);
2. some cell types are *withheld* from the graph, so that their positions can serve as a
   held-out outcome. Immune cells are withheld in the CODEX validation: the barrier is
   computed from the stromal/tumour scaffold only, and the measured CD8+ T-cell positions are
   compared with it afterwards. Nothing about immune cells can leak into the barrier.

The operators themselves are not modified. Only NumPy and SciPy are imported here; the cut
partition reuses ``barrier.exact_min_cut``.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import Delaunay, cKDTree
from scipy.stats import rankdata

__all__ = [
    "delaunay_graph", "largest_component", "rank01", "core_sinks",
    "cut_partition", "uniform_max_flow", "tissue_raster", "assign_nearest",
    "region_density_ratio", "contact_enrichment_z", "ripley_l",
]


def delaunay_graph(coords_um: np.ndarray, max_len_um: float = 50.0):
    """Delaunay triangulation of cell centroids, edges longer than ``max_len_um`` removed.

    Returns
    -------
    A : (n, n) CSR binary symmetric adjacency
    D : (n, n) CSR symmetric edge lengths (µm)
    """
    xy = np.asarray(coords_um, float)
    n = len(xy)
    if n < 3:
        raise ValueError("Delaunay needs at least 3 cells")
    tri = Delaunay(xy)
    s = tri.simplices
    e = np.vstack([s[:, [0, 1]], s[:, [1, 2]], s[:, [0, 2]]])
    e.sort(axis=1)
    e = np.unique(e, axis=0)
    L = np.linalg.norm(xy[e[:, 0]] - xy[e[:, 1]], axis=1)
    keep = L <= max_len_um
    e, L = e[keep], L[keep]
    rows = np.concatenate([e[:, 0], e[:, 1]])
    cols = np.concatenate([e[:, 1], e[:, 0]])
    A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    D = sp.csr_matrix((np.concatenate([L, L]), (rows, cols)), shape=(n, n))
    return A, D


def largest_component(A: sp.spmatrix) -> np.ndarray:
    """Boolean mask of the largest connected component."""
    _, lab = connected_components(A, directed=False)
    big = np.bincount(lab).argmax()
    return lab == big


def rank01(x: np.ndarray) -> np.ndarray:
    """(rank - 1) / (n - 1) with average ties, as for the section scores."""
    x = np.asarray(x, float)
    if len(x) <= 1:
        return np.full(len(x), 0.5)
    return (rankdata(x) - 1.0) / (len(x) - 1.0)


def core_sinks(A: sp.spmatrix, tumour: np.ndarray, q_core: float = 0.5) -> np.ndarray:
    """Tumour nodes at or beyond the ``q_core`` quantile of hop distance to non-tumour nodes.

    Same rule as ``graph.define_source_sink`` (unweighted Dijkstra from all non-tumour nodes).
    """
    tumour = np.asarray(tumour, bool)
    t_idx = np.flatnonzero(tumour)
    nt_idx = np.flatnonzero(~tumour)
    if len(t_idx) == 0 or len(nt_idx) == 0:
        return t_idx
    d = dijkstra(sp.csr_matrix(A), directed=False, indices=nt_idx, min_only=True)
    d = np.where(np.isfinite(d), d, 0.0)
    thr = np.quantile(d[t_idx], q_core)
    return t_idx[d[t_idx] >= thr]


def cut_partition(pairs: np.ndarray, cap: np.ndarray, nodes: Sequence[int],
                  source: Sequence[int], sink: Sequence[int]):
    """Max flow and the source-side node set of an exact minimum cut.

    Uses ``barrier.exact_min_cut`` (integer-scaled capacities), the same routine as
    ``barrier.compute_b_cell``. Returns ``(max_flow, source_side)`` where ``source_side`` is a
    boolean mask of length max(nodes) + 1.
    """
    from .barrier import exact_min_cut
    f, reach, _, _ = exact_min_cut(nodes, pairs, cap, source, sink)
    n = int(max(nodes)) + 1
    side = np.zeros(n, bool)
    side[[int(u) for u in reach if u != "s"]] = True
    return float(f), side


def uniform_max_flow(pairs: np.ndarray, nodes: Sequence[int], source: Sequence[int],
                     sink: Sequence[int], cap_value: float):
    """Max flow (and source side) with every edge at the same capacity: geometry only."""
    cap = np.full(len(pairs), float(cap_value))
    return cut_partition(pairs, cap, nodes, source, sink)


def tissue_raster(coords_um: np.ndarray, step_um: float = 5.0, radius_um: float = 10.0):
    """Raster points (µm) within ``radius_um`` of any cell centroid: an area estimate."""
    xy = np.asarray(coords_um, float)
    lo, hi = xy.min(axis=0) - radius_um, xy.max(axis=0) + radius_um
    gx = np.arange(lo[0], hi[0] + step_um, step_um)
    gy = np.arange(lo[1], hi[1] + step_um, step_um)
    P = np.column_stack([a.ravel() for a in np.meshgrid(gx, gy)])
    d, _ = cKDTree(xy).query(P, k=1)
    return P[d <= radius_um]


def assign_nearest(points: np.ndarray, ref_coords: np.ndarray) -> np.ndarray:
    """Index of the nearest reference point for every point."""
    if len(points) == 0:
        return np.zeros(0, int)
    _, j = cKDTree(np.asarray(ref_coords, float)).query(np.asarray(points, float), k=1)
    return j


def region_density_ratio(cell_in: np.ndarray, area_in: float, cell_out: np.ndarray,
                         area_out: float, pseudo: float = 0.5) -> float:
    """log2 of (cells per area inside) / (cells per area outside), with a pseudocount."""
    n_in = float(np.sum(cell_in))
    n_out = float(np.sum(cell_out))
    if area_in <= 0 or area_out <= 0:
        return float("nan")
    return float(np.log2(((n_in + pseudo) / area_in) / ((n_out + pseudo) / area_out)))


def contact_enrichment_z(A: sp.spmatrix, lab_a: np.ndarray, lab_b: np.ndarray,
                         n_perm: int = 1000, seed: int = 0) -> float:
    """z-score of the number of a–b contacts against label permutation (Squidpy-style)."""
    A = sp.triu(sp.csr_matrix(A), k=1).tocoo()
    u, v = A.row, A.col
    la = np.asarray(lab_a, bool)
    lb = np.asarray(lab_b, bool)
    obs = np.sum((la[u] & lb[v]) | (lb[u] & la[v]))
    rng = np.random.default_rng(seed)
    n = len(la)
    # permute the joint label vector (a, b, other) so that the counts are preserved
    lab = np.where(la, 1, np.where(lb, 2, 0))
    null = np.empty(n_perm)
    for k in range(n_perm):
        p = lab[rng.permutation(n)]
        null[k] = np.sum(((p[u] == 1) & (p[v] == 2)) | ((p[u] == 2) & (p[v] == 1)))
    sd = null.std()
    return float((obs - null.mean()) / sd) if sd > 0 else float("nan")


def ripley_l(points: np.ndarray, r_um: float, area_um2: float) -> float:
    """Ripley's L(r) − r without edge correction (descriptive comparator only)."""
    P = np.asarray(points, float)
    n = len(P)
    if n < 2 or area_um2 <= 0:
        return float("nan")
    t = cKDTree(P)
    pairs = t.count_neighbors(t, r_um) - n          # ordered pairs within r, self excluded
    K = area_um2 * pairs / (n * (n - 1))
    return float(np.sqrt(K / np.pi) - r_um)
