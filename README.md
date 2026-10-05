# SPARTA

**SPA**tial **R**esistance to **T**herapeutic **A**gents — a graph-transport
framework for the two physically distinct barriers that stand between a drug and
a tumour nest in spatial transcriptomics data.

> 中文版说明见 [`README.zh.md`](README.zh.md)（更详细，含开发笔记与踩坑记录）。

## SPARTA in three points

- It builds one spatial graph from each transcriptomics section and uses it to
  represent cellular migration and antibody diffusion as separate transport
  problems.
- The cellular operator is a source–sink minimum cut on the largest connected
  component; the antibody operator is a screened diffusion–absorption field
  with an explicit molecular-size exclusion term.
- It is a deterministic, CPU-only analysis with no learned parameters. The
  manuscript analysis covers 19 primary sections from 7 patients, an independent
  melanoma replication cohort (8 sections, 4 patients) and, for validation against
  measured CD8⁺ T-cell positions, 140 CODEX cores from 35 colorectal-cancer patients.

## Results reported in the manuscript

The following values are the manuscript’s locked results and should not be
changed without the author’s review: partial Spearman correlation is positive
in 18/19 sections (median ρ = +0.217); at β = 3, size exclusion explains a
median 97.5% of `B_mAb` variance; removing the shared matrix input leaves
significant coupling in 12/15 cSCC sections (median retention 81%), while the
four melanoma sections from one patient remain unresolved; and 60–65% of
edges exclude IgG. The accepted cohort contains 15 cSCC sections from 6
patients and 4 metastatic melanoma sections from 1 patient. CSCC13 remains
excluded (median UMI 289.5, below the 300 threshold).

External-validation arm (R7, added 2026-10-03): three 10x public Visium
sections (two breast cancer, one lymph node) reproduce the coupling in 3/3
(median partial ρ = +0.281), all three survive the spectral phase-randomisation
null, and the S2 matched-selection and tool-benchmark results transfer. These
sections are analysed separately and never enter the main-cohort numbers above.
The BRCA02 spatial-null result is at the Monte Carlo boundary (raw and BH-adjusted
p = 0.0499 with 500 surrogates); it should be described as borderline, not as a
robust third replication.


## v2.1.0: Interdisciplinary Sciences submission

Release v2.1.0 accompanies the manuscript submitted to *Interdisciplinary Sciences:
Computational Life Sciences*. The v2.0.0 values above remain correct as raw summaries, but the
analyses below change how they should be read. Every number quoted in the manuscript is written
by `scripts/is_figures/facts.py` to `results/validation/is_manuscript_facts.json`.

- **Simulation benchmark with ground truth** (`run_38`, `run_38b`): 540 simulated
  tissues in 9 planted-barrier geometries, with access measured by agent-based random
  walks. The minimum cut tracked lost access (Spearman ρ = 0.87) and was the only
  summary that separated closed capsules from capsules with a 5% gap (AUC 0.85;
  best alternative 0.72). The selection-matched gap experiment (S2) gave ratios
  near 1 on planted closed capsules, so it is reported as descriptive only.
- **Calibration on every real graph** (`run_37`): the graph-spectral surrogate test held the
  false-positive rate near 5% on all 22 section graphs (pooled 0.047), whereas a
  spot-level test reached 0.246. For strongly skewed fields use the normal-score variant
  (`sparta.spatial_stats.normal_scores`; pooled 0.064 against 0.150 for the raw spectrum).
- **Patient-level inference** (`run_36`): nested random-effects model (sections within patients),
  exact patient sign-flip test and leave-one-patient-out. Pooled ρ = 0.20 (95% CI 0.14–0.27);
  positive in 7/7 patients.
- **Structural nulls** (`run_39`): with the real ECM score kept and the crosslinking and ligand
  inputs replaced by graph-spectral surrogates, the model reproduces a median 60% of
  the observed coupling. The tissue-specific excess is small (pooled 0.07, 95% CI
  0.01–0.14).
- **Robustness** (`run_41`, `run_42`, `run_43`): alternative vessel-distance adjustments (median ρ
  0.170–0.246), a size-exclusion scan over β and a normal-score re-test of the primary
  association (15/19 sections BH q < 0.05).
- **Reproduction without Scanpy** (`run_35`, `run_35b`, `sparta/node_tables.py`): the per-spot node
  tables in `data/interim/*.nodes.npz` reproduce the stored field associations (max |Δρ| =
  0.0009) without Scanpy, AnnData or h5py.

New statistics live in `sparta/spatial_stats.py` and are covered by `tests/test_is_revision.py`.
The two fields co-vary in 18/19 sections and in every patient, but most of that coupling is
produced by shared inputs and graph geometry. Outputs are model-defined; they are not measurements
of antibody exposure, T-cell passage or treatment response.

| Script | Purpose |
|---|---|
| `run_35_export_node_tables.py` | export Scanpy-free per-spot node tables (run once on the machine that holds the `.h5ad` files) |
| `run_35b_verify_reproduction.py` | recompute the field associations from the node tables and compare with the stored results |
| `run_36_patient_level.py` | nested random-effects model, patient sign-flip test, leave-one-patient-out |
| `run_37_null_calibration.py` | false-positive rate of the surrogate tests on every real section graph |
| `run_38_synthetic_benchmark.py`, `run_38b_synthetic_s2.py` | planted-barrier simulation benchmark with agent-based access |
| `run_39_coupling_decomposition.py` | geometry and construction null models for the field coupling |
| `run_40_graph_coverage_table.py` | per-section graph coverage, components and stranded spots |
| `run_41_adjustment_robustness.py` | alternative adjustments for vessel distance |
| `run_42_size_exclusion_scan.py` | fraction of IgG-excluding edges as a function of β |
| `run_43_nscore_spatial_null.py` | primary association re-tested with normal-score surrogates |
| `run_44_intervention_ranking.py` | per-spot counterfactual intervention ranking for both barriers (`--jobs`, `--force`, `--figures-only`) |
| `run_45_tcga_survival.py` | optional exploratory bulk concept check: TCGA-SKCM signatures vs overall survival (requires network on first run) |
| `run_46_simple_baselines.py` | exploratory internal-consistency check of target rankings against SPARTA's own counterfactual (not a validation) |
| `run_47_codex_validation.py` | v2.2: cellular barrier against measured CD8⁺ T-cell positions in CODEX images (pre-specified protocol) |
| `run_48_replication_cohort.py` | v2.2: melanoma replication cohort, Scanpy-free processing (`--validate` checks equivalence) and statistics |
| `run_49_intervention_targeting.py` | v2.2: joint ablation curves, half-breach budgets and targeting comparisons for the intervention maps |
| `run_50_exact_mincut_refresh.py` | v2.2: recompute every stored minimum cut exactly and document old versus new |
| `run_51_domain_comparison_exact_cut.py` | v2.2: BANKSY-style domain comparison with the exact cuts (Scanpy optional; reproduces the archived Scanpy numbers exactly) |
| `run_52_simple_baselines_spearman.py` | v2.2: both operators against stromal density, distance to the tumour and neighbourhood enrichment (30 sections) |
| `run_53_intervention_domain_enrichment.py` | v2.2: tumour / stromal / immune composition of the top 1% and 5% intervention spots against random placement |
| `scripts/is_figures/*.py` | all IS figures (`fig1`–`fig8`, `figS_supplement`, `figS7_intervention`) and `facts.py` |

---

## v2.2.0 (2026-10-05): measured-outcome validation, replication cohort, exact cuts, intervention maps

None of these additions changes a locked primary-cohort association number above. Every number quoted in
the manuscript is still written by `scripts/is_figures/facts.py` to `results/validation/is_manuscript_facts.json`.

### 1. The cellular operator against measured CD8⁺ T-cell positions (CODEX, colorectal cancer)

`sparta/cellgraph.py`, `scripts/run_47_codex_validation.py`, protocol `docs/codex_validation_protocol.md`
(written before any barrier or outcome was computed; its SHA-256 is stored with the results).
In the CODEX images of Schürch et al. (2020; 140 cores, 35 patients), immune cells are withheld from the
graph: the barrier is computed from the tumour/stromal scaffold only (Delaunay graph, collagen IV and
αSMA/vimentin protein in place of the ECM and CAF scores, Eq. 1 unchanged), and the measured CD8⁺ T-cell
positions are the held-out outcome. The geometry-normalised barrier `B_rel = F_open / F` predicted CD8⁺
depletion from tumour cores (114 cores, 35 patients; patient-level ρ = −0.53, one-sided p = 0.0006; core level
−0.43, 95% CI −0.57 to −0.26), in all five sensitivity variants. It was **not uniquely informative**: peritumoural
matrix density did as well (−0.54) and tumour–stroma intermixing better (+0.58, opposite direction); the partial
association given composition was −0.19 (−0.33 to −0.03) and fell to −0.14 (−0.33 to 0.03) once intermixing
was added (post hoc). The cut was not the line beyond which T cells were depleted (pre-specified secondary
analysis negative). Outputs: `results/validation/codex_validation.json`, `codex_cores.csv`.

### 2. Independent melanoma replication cohort (Thrane et al. 2018)

`scripts/run_48_replication_cohort.py`, `sparta/lite.py`, plan `docs/replication_melanoma_protocol.md`.
Eight first-generation ST sections of lymph-node metastases from four stage III melanoma patients, processed
with the locked pipeline and parameters (admission entry `mel_thrane2018`: 200-spot minimum for the
1,007-spot array; C3/C7 waived and disclosed). The earlier ingestion had kept "SYMBOL ENSG…" gene names, so
no endothelial marker matched; names are now cleaned. Scanpy is not needed: `sparta/lite.py` re-implements
QC, normalisation and `score_genes`, and `--validate` shows that it reproduces the archived scores of three
primary first-generation sections to 2.7 × 10⁻⁷ (ranks identical). Result: 8/8 sections positive (median ρ
0.31), 8/8 BH q < 0.05, 4/4 patients positive; pooled ρ = 0.26 (95% CI 0.09–0.42), meeting the pre-specified
criterion; all 11 patients pooled 0.21 (0.16–0.27). The construction null reproduced a median 34% of the
coupling (60% in the primary cohort). The ledger status of these sections is `replication`, so
`io_.admitted_slides()` and every primary-cohort script ignore them.

### 3. Exact minimum cuts (bug fix)

`barrier.exact_min_cut`, `scripts/run_50_exact_mincut_refresh.py`, `tests/test_v22_additions.py`.
networkx's preflow-push on floating-point capacities returned the correct **flow value** (|Δ| < 5 × 10⁻¹² relative)
but its residual-graph partition was **not always a minimum cut**: in 14 of 30 section graphs the returned cut
edges added up to 0.03–9.7% more or less than the maximum flow. Cuts are now computed with capacities in
integer units of 2⁻⁴⁰ (cut capacity = max flow to 4 × 10⁻¹²). `B_cell` values are unchanged; everything that uses
cut geometry was recomputed (`mincut.json`, `barrier.npz` cut fields, S2 gap experiments, intervention maps,
CODEX blockade-line analysis, figures). The BANKSY-style domain comparison (`run_13b`, `run_34`) was recomputed
with `run_51_domain_comparison_exact_cut.py`, which re-implements the Scanpy preprocessing (Seurat-flavour HVG,
scale, ARPACK PCA) in `sparta/lite.py` and reads `.h5ad` through libhdf5 when h5py is absent; given the old
cuts it reproduced the archived numbers exactly in 19/19 sections, and with the exact cuts the median
enrichment moved from 1.15 to 1.14 (precision unchanged). Six stored `barrier.npz` files also carried max-flow values computed before the
largest-component rule; they were refreshed (no reported statistic had read them).

### 4. Counterfactual intervention maps (modality-specific; Online Resource 1, Fig. S7)

`sparta/intervention.py`, `scripts/run_44_intervention_ranking.py` (`--refresh-cut`, `--include-replication`),
`scripts/run_49_intervention_targeting.py`, `scripts/run_53_intervention_domain_enrichment.py`.
Every spot of all 30 sections (33,392 spots) is ablated once (ECM/CAF/crosslinking set to the section's 5th
percentile, the S2 convention). Two properties are **mathematical, not findings**: by max-flow/min-cut duality,
single-spot effects on `B_cell` are zero away from the (now exact) cut, so the top-ranked spots lie on the cut.
The informative quantities are joint. Cellular barrier, focal: the best single spot removes a median 19% of a
full cut breach (up to 46%, CSCC07) and the top 1% of spots together a median 62% (33–86%); 20 map-chosen spots
recover 77% (random cut spots 36%, densest matrix 14%). Antibody barrier, diffuse: the top 1% removes 12%, the
top 5% 31%, the top 10% 44% of the matrix-free limit. The two maps barely overlap (top-5% Jaccard median 5.7%,
range 2–20%). High-impact spots are under-represented in immune-dominated spots for both operators (top 5%:
29% and 28% vs 36%; patient-level sign-flip p = 0.027 and 0.012), modestly stromal for B_cell and tumour for
B_mAb. All are in-model counterfactuals.

### 5. What the operators add over simple summaries (main-text Fig. 8)

`scripts/run_52_simple_baselines_spearman.py`; synthetic benchmark and CODEX analysis extended with
geometry-only summaries (`run_38`: vessel–nest distance; `run_47`: vessel–tumour distance and tumour-core
depth, post hoc). Within sections, the B_cell field shares about a quarter of its rank variance with local
stromal density (median ρ = 0.50) and little with distance to the tumour (0.07) or neighbourhood enrichment
(0.34); the B_mAb field is weakly related to all three (−0.14 to −0.16). Across 29 tumour sections the
matrix-dependent cellular barrier `log B_rel` is unrelated to all three section-level summaries. In simulations
no simple summary separates closed from gapped capsules (AUC ≤ 0.72 vs 0.85). In CODEX, peritumoural matrix
density matches the cut and tumour-core depth (post hoc) is the strongest single correlate (−0.60), but the cut
keeps a partial association given composition and both distances (−0.20, −0.35 to −0.04); not given all four
simple summaries together (−0.11, −0.32 to 0.06).

### Exploratory, not in the manuscript

`run_45_tcga_survival.py` (bulk TCGA-SKCM survival; barrier summary null, immune composite positive control)
and `run_46_simple_baselines.py` are kept as exploratory scripts. `run_46` scores candidate rankings against
SPARTA's own in-model counterfactual, so it is an internal-consistency check, not a validation.

## What it does

Most spatial analyses of immunotherapy resistance ask one question: *can a T cell
reach the tumour?* An anti-PD-1 antibody is not a cell. It is a ~150 kDa
macromolecule with a hydrodynamic radius near 5.5 nm that must diffuse through
interstitial matrix. SPARTA models an IgG-sized molecule, using a PD-L1/PD-L2
ligand-expression score as an uncalibrated absorption proxy. This is not a
pharmacologically specific anti-PD-1 model, and the effective mesh is not
measured in the transcriptomic data. That is a distinct, model-defined transport
problem with different physics.

SPARTA poses both on **the same spatial graph**, changing only the edge-weight
semantics and the operator:

| Component | Obstructs | Operator | What you get back |
|---|---|---|---|
| `B_cell` | CD8⁺ T cells (~10 µm, actively migrating) | **source–sink minimum cut** | a barrier magnitude *and* an explicit blockade geometry |
| `B_mAb` | IgG / ADC (~5.5 nm, passively diffusing) | **screened Poisson** `(L + diag(κ))φ = 0`, `B = −log φ` | a per-spot exposure deficit field |
| `B_meta` | small molecules (delivered but ineffective) | graph distance from vasculature × metabolic state | a per-spot refuge score |

When absorption vanishes the second operator reduces to the harmonic /
effective-resistance limit, so the graph-transport reading is exact rather than
metaphorical.

**CPU only. No GPU, no deep learning.** The core operators (`barrier`, `graph`,
`counterfactual`, `synthetic`) depend on NumPy, SciPy and NetworkX alone.

---

## Install

```powershell
py -3.13 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m pip install -e ".[spatial]"
```

`requirements.txt` records the exact package versions from the analysis
environment. `environment.yml` is a Conda alternative. Minimal install — enough
to run the core operators and test scripts without Scanpy:

```bash
pip install numpy scipy networkx matplotlib pyyaml
```

## Five-minute check

```bash
python scripts/run_00_demo.py --fast   # end-to-end on a synthetic section
python tests/test_barrier.py           # 9 tests
python tests/test_counterfactual.py    # 4 tests
python tests/test_loaders.py           # 5 tests
python tests/test_statistics.py        # 4 tests
python tests/test_is_revision.py        # 8 tests (IS revision statistics; no Scanpy)
python tests/test_intervention.py      # 4 tests (v2.2 intervention ranking; no Scanpy)
python tests/test_v22_additions.py     # 5 tests (exact min cut, cell graphs, Scanpy-free scoring)
python tests/test_new_experiments.py   # 3 result smoke checks (skip before real-data outputs exist)
```

All eight scripts run without pytest. The seven synthetic/core scripts contain
39 checks; the three real-result smoke checks run when their JSON outputs are
present. `run_00_demo.py` is also the fastest way to understand the framework
— read it once before touching real data.

---

## Pipeline

```
data/raw/{slide}/
      │  run_00b_ingest.py     ingest + data ledger (patient / replicate)
      │  run_00c_admission_audit.py   auditable C1–C7 admission record
      ▼  run_01_qc.py          QC + admission checks
{slide}.qc.h5ad
      ▼  run_02_score.py       signature scoring + rank normalisation
{slide}.scored.h5ad
      ▼  run_03_graph.py       radius graph + source / sink / vessel sets
{slide}.graph.npz
      ▼  run_04_barrier.py     three barrier components + blockade overlay
{slide}.barrier.npz
      ├─▶ run_05_counterfactual.py    S1 rearrangement / S2 continuity / S3 size
      └─▶ run_06_validate.py          consistency, decoupling, morphology
```

Analysis and reporting on top of that:

| Script | Purpose |
|---|---|
| `run_12_paper_stats.py` | section-level manuscript summaries and multiple-testing corrections |
| `run_13_benchmark_tools.py` | comparison against a BANKSY-style spatial-domain segmentation |
| `run_14_shared_ecm_check.py` | robustness check: recompute with the shared matrix term removed from both operators |
| `run_15_figure1.py`, `run_16_figures.py` | all manuscript figures |
| `run_17_sensitivity.py` | `B_mAb` parameter grids across sections |
| `run_20_pending_figures.py` | Figures 3–5 and the graphical abstract from saved outputs |
| `run_21_mesh_stats.py` | per-edge effective mesh size and IgG size-exclusion statistics |
| `run_22_null_crosslink.py` | permutation check for rank-normalised crosslinking (all 19 primary sections) |
| `run_23_stromal_intervention.py` | in silico score-scaling scenario for both barriers |
| `run_24_naive_baseline.py` | graph-field versus node-score comparison |
| `run_25_biological_validation.py` | exploratory alignment with measured signatures; see the source-definition caveat in `docs/review_for_journal.md` |
| `run_26_radius_sensitivity.py` | graph-radius sensitivity across both platforms |
| `run_27_spatial_null.py`, `run_28_s2_matched.py` | graph-spectral association null and selection-matched S2 counterfactual |
| `run_32_prereadiness_audit.py` | metric-specific missing-input audit and disconnected source/sink counts |
| `run_33_input_connectivity_sensitivity.py` | input-complete-case and largest-component association sensitivities |
| `run_34_benchmark_lambda_sensitivity.py` | neighborhood-weight sensitivity for the BANKSY-style descriptive benchmark |

Core model parameters live in `configs/default.yaml`; analysis-specific grids, null counts and sensitivity subsets are recorded by their scripts and output metadata.

---

## Three things worth knowing before you use it

**1. Sections from one patient are not independent.** The data ledger carries
`patient` and `replicate` columns and the statistics aggregate at both levels.
Report both; never quote a section count as a sample size.

**2. `ξ₀` and `β` are qualitative scale parameters, not measurements.** The
size-exclusion law is written in nanometres but takes a *within-section
rank-normalised* crosslinking score as input. That fixes the nominal mesh range
and the exclusion threshold identically in every section, so those two
parameters govern how sharply the model separates permeable from impermeable
tissue — they do not measure mesh size. Any conclusion that depends on them is
reported with that dependence, and `run_17_sensitivity.py` produces the grid.

**3. The two barriers share an input by construction.** The core matrisome score
enters both the min-cut capacity and the diffusion conductance, so some coupling
between them is structurally guaranteed. `run_14_shared_ecm_check.py` exists to
separate that from any coupling actually present in the tissue, and it should be
run before any claim about the relationship between the two barriers.

---

## Reproducing the analysis

The full analysis can be rerun with one command after the public datasets have
been ingested into `data/interim/` and recorded in `data/ledger.csv`:

```powershell
.\.venv\Scripts\python.exe scripts\run_all.py
```

The runner rebuilds the analysis in dependency order, then regenerates the
manuscript figures and Word file. To rerun the optional manuscript checks in
scripts 22–26 before regenerating the figures and Word file, add
`--include-review-analyses`. It is deterministic and safe to rerun; the
500-permutation counterfactual analysis is the slowest stage and makes the full
run take hours. The core operators themselves take under 0.2 s per section on
CPU. The repository does not include the public raw datasets. For individual
steps and their expected outputs, see
[`docs/runbook_after_review.md`](docs/runbook_after_review.md); for known
method limitations, see [`docs/review_for_journal.md`](docs/review_for_journal.md).

One-time scripts from the exploratory phase are archived under
`scripts/scratch/` with a README explaining which have been folded into the
reproducible pipeline and which have not.

### External-validation arm (R7)

The three external sections (10x Genomics public Visium: human breast cancer
Block A Sections 1/2 and one human lymph node) are **not** part of the
`run_all.py` main cohort. Reproduce them independently:

```powershell
# 1. download (matrix + spatial tarballs, ~350 MB total)
#    https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/
#    ... same pattern for V1_Breast_Cancer_Block_A_Section_2 and V1_Human_Lymph_Node
#    unpack each tar.gz into data/external/brca_vis/<slide>/

# 2. ingest each section (repeat with the corresponding accession and replicate)
.\.venv\Scripts\python.exe scripts\run_00b_ingest.py --input data\external\brca_vis\V1_Breast_Cancer_Block_A_Section_1 --slide BRCA01 --cancer other --platform visium --source "10x Genomics public Visium (breast cancer Block A Section 1)" --accession V1_Breast_Cancer_Block_A_Section_1 --patient BRCA_P1 --replicate s1 --treatment unknown --site primary
# Repeat for BRCA02 (patient BRCA_P1, replicate s2) and LN01
# (patient LN_P1, replicate s1, --site unknown).

# 3. QC with the documented cohort waiver (C7 is unannotated)
.\.venv\Scripts\python.exe scripts\run_01_qc.py --slide BRCA01 --input data\interim\BRCA01.raw.h5ad --platform visium --has-he --cohort brca_10x_vis --force --force-reason "cohort waiver, see configs/default.yaml admission_overrides.brca_10x_vis"
# Repeat for BRCA02; LN01 uses --cohort ln_10x_vis and its recorded waiver.

# 4. score with the external epithelial-marker profile, then build graphs/barriers
# cancer_type=other is not a supported marker set and must not silently fall back.
$slides = 'BRCA01','BRCA02','LN01'
foreach ($slide in $slides) {
  .\.venv\Scripts\python.exe scripts\run_02_score.py --slide $slide --tumor-type brca --plot
  .\.venv\Scripts\python.exe scripts\run_03_graph.py --slide $slide
  .\.venv\Scripts\python.exe scripts\run_04_barrier.py --slide $slide
  .\.venv\Scripts\python.exe scripts\run_05_counterfactual.py --slide $slide --s2-match-groups 100
}

# 5. independent external summary files; main-cohort outputs stay untouched
.\.venv\Scripts\python.exe scripts\run_29_ext_validation.py
.\.venv\Scripts\python.exe scripts\run_30_s2_matched_ext.py
.\.venv\Scripts\python.exe scripts\run_31_benchmark_ext_slides.py
```

The external arm is a cross-dataset reproduction, not independent patient-level
validation: BRCA01/BRCA02 are adjacent sections from one patient and LN01 is a
non-tumour control. The epithelial-marker sink is model-defined and has not been
confirmed against a pathology annotation.

The waiver and the `run_00c_admission_audit.py` accession mapping were
registered on 2026-10-03; see `configs/default.yaml`
(`admission_overrides.brca_10x_vis`, `admission_overrides.ln_10x_vis`).

---

## Data availability

All data are public. Spatial transcriptomics:
[GSE250636](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250636)
(metastatic cutaneous melanoma, Visium) and
[GSE144239](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144239)
(primary cutaneous squamous cell carcinoma; Ji et al., *Cell* 2020).
External validation (R7): [10x Genomics spatial-expression sample repository]
(https://www.10xgenomics.com/resources/datasets), Space Ranger v1.1.0,
`V1_Breast_Cancer_Block_A_Section_1/2` and `V1_Human_Lymph_Node`.
Bulk ICB cohorts used for a concept check only: GSE78220, GSE91061.
Gene sets: MSigDB HALLMARK_HYPOXIA (v7.1); per-signature provenance in
[`docs/gene_set_references.md`](docs/gene_set_references.md).

**This repository contains no raw data.** From v2.1.0 the derived per-spot node tables, spot
graphs, model fields and minimum-cut records of the 22 analysed sections (`data/interim/*.nodes.npz`,
`*.graph.npz`, `*.barrier.npz`, `*.mincut.json`), every result file (`results/validation/`,
`results/counterfactual/`) and the IS figures are versioned, so every manuscript number can be
regenerated without the raw expression matrices.

## Citation

Release v2.1.0 (code, node tables, result files and figure scripts) accompanies the manuscript
submitted to *Interdisciplinary Sciences: Computational Life Sciences*. Its Zenodo DOI will be added
here and in [`CITATION.cff`](CITATION.cff) once the deposit and landing page have been verified.

## License

MIT — see [`LICENSE`](LICENSE).

