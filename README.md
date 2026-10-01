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
  manuscript analysis covers 19 sections from 7 patients.

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

---

## What it does

Most spatial analyses of immunotherapy resistance ask one question: *can a T cell
reach the tumour?* An anti-PD-1 antibody is not a cell. It is a ~150 kDa
macromolecule with a hydrodynamic radius near 5.5 nm that must diffuse through a
matrix whose effective mesh is measured in tens of nanometres, and it can be
consumed en route by binding its target. That is a different transport problem
with different physics.

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
to run the core operators and tests without Scanpy:

```bash
pip install numpy scipy networkx matplotlib pyyaml
```

## Five-minute check

```bash
python scripts/run_00_demo.py --fast   # end-to-end on a synthetic section
python tests/test_barrier.py           # 6 tests
python tests/test_counterfactual.py    # 4 tests
python tests/test_loaders.py           # 4 tests
python tests/test_statistics.py        # 4 tests
```

21/21 should pass. `run_00_demo.py` is also the fastest way to understand the
framework — read it once before touching real data.

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
| `run_12_paper_stats.py` | statistics for the manuscript — BH within per-section families, pre-specified S2 removal fraction, **patient-level aggregation** |
| `run_13_benchmark_tools.py` | comparison against a BANKSY-style spatial-domain segmentation |
| `run_14_shared_ecm_check.py` | robustness check: recompute with the shared matrix term removed from both operators |
| `run_15_figure1.py`, `run_16_figures.py` | all manuscript figures |
| `run_17_sensitivity.py` | `B_mAb` parameter grids across sections |

Every parameter lives in `configs/default.yaml`. Nothing is hard-coded.

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
manuscript figures and Word file. It is deterministic and safe to rerun; the
500-permutation counterfactual analysis is the slowest stage and makes the full
run take hours. The core operators themselves take under 0.2 s per section on
CPU. The repository does not include the public raw datasets. For individual
steps and their expected outputs, see
[`docs/runbook_after_review.md`](docs/runbook_after_review.md); for known
method limitations, see [`docs/review_for_journal.md`](docs/review_for_journal.md).

One-time scripts from the exploratory phase are archived under
`scripts/scratch/` with a README explaining which have been folded into the
reproducible pipeline and which have not.

---

## Data availability

All data are public. Spatial transcriptomics:
[GSE250636](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250636)
(metastatic cutaneous melanoma, Visium) and
[GSE144239](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144239)
(primary cutaneous squamous cell carcinoma; Ji et al., *Cell* 2020).
Bulk ICB cohorts used for a concept check only: GSE78220, GSE91061.
Gene sets: MSigDB HALLMARK_HYPOXIA (v7.1); per-signature provenance in
[`docs/gene_set_references.md`](docs/gene_set_references.md).

**This repository contains no data.** `data/` and `results/` are gitignored.

## Citation

Code archive: [Zenodo DOI 10.5281/zenodo.23086431](https://doi.org/10.5281/zenodo.23086431).
See [`CITATION.cff`](CITATION.cff) for the software citation metadata.

## License

MIT — see [`LICENSE`](LICENSE).
