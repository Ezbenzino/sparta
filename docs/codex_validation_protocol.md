# Pre-specified protocol — SPARTA cellular barrier against measured CD8⁺ T-cell positions (CODEX, colorectal cancer)

Written 2026-10-05 00:20 (Asia/Shanghai; 2026-10-04 16:20 UTC) **before any barrier or outcome
quantity was computed on these data**. The only inspection done beforehand was a tabulation of
cell-type counts, patients, groups and core sizes (to define the node sets below). The SHA-256 of
this file is recorded in `results/validation/codex_validation.json`; any later deviation is listed
at the end of this file under *Deviations*, with its reason.

## Question

Does the SPARTA minimum-cut barrier, computed only from the stromal scaffold of a tissue section,
predict where CD8⁺ T cells are actually found? The spatial-transcriptomics analyses in the
manuscript have no measured transport outcome. Multiplexed protein imaging at single-cell
resolution measures every T cell's position, so the T cells can be withheld from the barrier
computation and used as a held-out outcome.

## Data

Schürch et al. (2020) *Cell* 182:1341–1359, CODEX imaging of colorectal-cancer tissue microarrays:
258,385 segmented cells, 140 cores (4 per patient), 35 patients (17 Crohn's-like reaction, CLR;
18 diffuse inflammatory infiltration, DII). Published cell-type labels (`ClusterName`) and
per-cell marker intensities are used as distributed; nothing is re-clustered. Pixel size
0.377 µm (20× objective, as reported by the authors).

## Node sets

- Excluded: `dirt`.
- **Scaffold** (the only cells that define the graph and the barrier): `tumor cells`, `stroma`,
  `smooth muscle`, `vasculature`, `lymphatics`, `adipocytes`, `nerves`, `undefined`.
- **Withheld** (never enter the barrier computation): every immune or mixed immune cluster
  (all T-cell, B-cell, plasma-cell, NK, macrophage, monocyte, dendritic-cell and granulocyte
  clusters, `immune cells`, `immune cells / vasculature`, `tumor cells / immune cells`).
- **Outcome cells**: `CD8+ T cells` (primary); all T-cell clusters (secondary).

## Graph, inputs and barrier (operator unchanged)

- Delaunay triangulation of scaffold-cell centroids in µm; edges longer than 50 µm removed;
  analysis on the largest connected component (LCC), as for the section graphs.
- Inputs, rank-normalised within the core over LCC scaffold nodes as $(\mathrm{rank}-1)/(n-1)$:
  matrix $E$ = Collagen IV (the only core-matrisome protein in the panel);
  fibroblast $F$ = rank of the mean of the within-core ranks of αSMA and vimentin.
- Edge capacity Eq. (1) with the manuscript's parameters $a=3$, $b=8$, $c=4$ — **no tuning**.
- Sources: `vasculature` cells in the LCC (T-cell neighbourhoods are *not* used, unlike the
  section analyses, because T cells are the outcome).
- Sinks: `tumor cells` in the LCC at or beyond the median graph (hop) distance to the nearest
  non-tumour scaffold node (same rule as the sections, `graph.define_source_sink`).
- $B_{\mathrm{cell}} = 1/(\text{max flow}+\epsilon)$ (unchanged operator).
- Geometry-normalised barrier $B_{\mathrm{rel}} = F_{\mathrm{open}}/F$, where $F$ is the max flow
  and $F_{\mathrm{open}}$ the max flow on the same graph with every capacity set to $\sigma(a)$
  (matrix-free tissue of identical geometry). $\log B_{\mathrm{rel}} = \log B_{\mathrm{cell}} - \log B_{\mathrm{open}}$
  removes the number of vessels and the size of the tumour interface and keeps only the
  contribution of matrix arrangement.

## Inclusion (decided from scaffold composition only)

A core is analysed if its LCC has ≥ 100 scaffold nodes, ≥ 5 sources, ≥ 25 tumour cells and
≥ 10 sinks, and ≥ 10 CD8⁺ T cells lie in its tissue area.

## Outcome (held out)

Tissue area: points of a 5-µm raster within 10 µm of any non-artefact cell centroid. Every CD8⁺ T
cell and every tissue raster point is assigned to its nearest LCC scaffold node. The tumour-core
region is the set assigned to sink nodes; the rest is the remainder.

**Primary outcome**: core infiltration ratio
$\mathrm{IR} = \log_2\!\left[\dfrac{(n_{\mathrm{core}}+0.5)/A_{\mathrm{core}}}{(n_{\mathrm{rest}}+0.5)/A_{\mathrm{rest}}}\right]$,
the relative CD8⁺ density in the tumour core (negative = depleted, i.e. excluded).

## Primary analysis

Hypothesis: higher $B_{\mathrm{rel}}$ goes with lower IR. Unit of inference = patient.
Primary test: one-sided Spearman correlation between patient means of $\log B_{\mathrm{rel}}$ and of
IR over included cores (35 patients if all have ≥ 1 included core). Supporting: core-level Spearman
with a patient-cluster bootstrap 95% CI (2,000 resamples of patients).

## Comparators (same cores, same outcome; oriented so that larger = more barrier)

global mean of $(8E+4F)/12$ over scaffold nodes; peritumoural mean of the same quantity over
non-tumour scaffold nodes within 30 µm of a tumour cell; stromal fraction (stroma + smooth muscle
among scaffold nodes); tumour–stroma contact enrichment z-score (Delaunay contacts against
1,000 label permutations); Ripley's $L$ of matrix-rich scaffold cells ($E \ge 0.7$) at 50 µm;
SPARTA $B_{\mathrm{cell}}$ field at the core (median shortest-path cost of the sinks);
absolute $B_{\mathrm{cell}}$. For each: core-level and patient-level Spearman with IR.
SPARTA versus the best comparator: patient-cluster bootstrap of the difference in |ρ|.
Added value: partial Spearman of $\log B_{\mathrm{rel}}$ with IR given tumour fraction, stromal
fraction, global matrix and peritumoural matrix (rank-regression residuals), patient-cluster
bootstrap CI.

## Secondary analyses

S1 Absolute access: $B_{\mathrm{cell}}$ against CD8⁺ density in the core (cells/mm²).
S2 Blockade line: per core, CD8⁺ density on the sink side of the SPARTA cut divided by that on the
source side, compared with the same ratio for the geometry-only cut (uniform capacities); paired
Wilcoxon signed-rank test on patient means of the log-ratio difference.
S3 CLR versus DII: patient-mean $\log B_{\mathrm{rel}}$ (two-sided Mann–Whitney; exploratory).
Sensitivity: edge cut-off 30 and 100 µm; $F$ = αSMA only; sinks = all tumour cells; outcome =
all T cells.

## Reporting

All results are reported whatever their direction. Cores, patients and exclusions are counted.
Outputs: `results/validation/codex_validation.json`, `results/validation/codex_cores.csv`.

## Deviations

1. Edited before any computation (2026-10-05 00:35 Asia/Shanghai): the sink rule was first written
   with a µm-weighted distance; it now uses the hop distance, which is what the section pipeline
   (`graph.define_source_sink`) uses. No result existed at the time of the edit.
2. Found when the results were first summarised (2026-10-05 01:00 Asia/Shanghai): the comparator
   "global mean of $(8E+4F)/12$" is constant (0.5) in every core, because $E$ and $F$ are within-core
   ranks. It is reported as degenerate rather than as a result, and the pre-specified partial
   correlation therefore adjusts for tumour fraction, stromal fraction and peritumoural matrix only.
3. Before the final run (2026-10-05): every minimum cut is now computed with integer-scaled capacities
   (`barrier.exact_min_cut`), because the floating-point partition returned by networkx was not always a
   minimum cut. Max-flow values, and therefore every barrier value, are unchanged; the source and sink
   sides used by secondary analysis S2 can change.

## Post hoc analyses (not pre-specified; reported as such)

- Tumour–stroma contact enrichment (intermixing) was the strongest single correlate of the outcome,
  in the direction opposite to a barrier reading. The partial association of $\log B_{\mathrm{rel}}$
  with IR is therefore also reported given intermixing, and given intermixing plus composition;
  intermixing and peritumoural matrix are also reported given all other summaries.
- Sensitivity of the primary analysis to the Eq. (1) weights: halved ($b=4$, $c=2$), increased by half
  ($b=12$, $c=6$), matrix only ($b=12$, $c=0$) and fibroblast only ($b=0$, $c=12$).
- Geometry-only comparators (added 2026-10-05 10:00 Asia/Shanghai, after the results above had been seen,
  to answer the question of what the operator adds over simple summaries): the median distance from
  vascular cells to the nearest tumour cell (`vessel_tumour_distance`) and the median distance from sinks to
  the nearest non-tumour cell (`core_depth`). Their correlations with IR are reported, and the added value
  of $\log B_{\mathrm{rel}}$ is also reported given the two distances and given composition plus distances.
  The core-depth comparator is partly tied to the outcome by construction, because IR is measured in the
  region assigned to the sinks whose depth it describes.
  The partial association of $\log B_{\mathrm{rel}}$ with IR is also reported given each of peritumoural matrix,
  vessel-to-tumour distance, core depth and intermixing alone, and given all four together.

