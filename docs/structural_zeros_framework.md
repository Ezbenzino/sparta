# Construction nulls: separating built-in coupling from tissue signal in multi-operator spatial models

**Status:** methodology note, written 2026-10-04. It reframes the analysis currently
implemented as `scripts/run_39_coupling_decomposition.py` (`results/validation/geometry_null.json`,
`results/validation/shared_input_spatial_null.json`) from *an analysis step inside the SPARTA paper*
into *a reusable framework*, and it states exactly what that framework can and cannot support.

---

## 1. The contribution in one paragraph

Any spatial model that runs **two or more operators over the same graph** and reports that their
outputs co-vary is reporting a quantity it has partly built itself. Two operators that read the same
input score, are anchored to the same source set, and run on the same geometry will correlate even
in tissue with no biology in it. This note gives a small, explicit procedure — three input-level
nulls plus one parameter ablation — that partitions an observed output–output association into
(i) anchoring geometry, (ii) model construction (shared inputs), and (iii) operator-specific
co-arrangement of the inputs in tissue. The residual after (i) and (ii) are removed is the only part
that can be read as tissue signal, and it comes with an explicit null distribution and a calibrated
false-positive rate. The procedure is not specific to sparta: it applies to any pipeline whose
operators accept substitutable per-node inputs on a shared graph.

---

## 2. The failure mode

Given two per-spot model outputs `f_A` and `f_B` computed on the same graph `G` (adjacency `A`,
coordinates, a source/anchoring set), an observed positive association `ρ(f_A, f_B | controls)` can be
produced by three distinct mechanisms:

| # | Mechanism | Why it creates coupling | SPARTA instance |
|---|---|---|---|
| i | **Anchoring geometry** | Both fields are governed by distance from the *same* source/anchor set; smooth fields on the same graph are positively correlated by construction | `B_cell` accumulates from source sets; `B_mAb` is a screened-Poisson field with `φ = 1` on the vessel set |
| ii | **Shared inputs** | Both operators read the same per-node score, so the shared score's spatial pattern appears in both outputs | ECM score enters the min-cut capacity **and** the diffusion conductance; the CAF set overlaps the ECM gene set (median input `ρ(ECM, CAF) = 0.87`) |
| iii | **Co-arrangement in tissue** | Operator-specific inputs really are spatially aligned in the sample | CAF (only in `B_cell`) co-located with crosslinking and ligand absorption (only in `B_mAb`) |

Mechanisms (ii) and (iii) are the ones that matter and are the ones usually not separated.
A correlation driven by (ii) says nothing about the tissue; a correlation that survives (i)+(ii)
is a statement about the sample. Reporting a single observed `ρ` conflates them.

---

## 3. The framework

Four components, all input-level (they perturb the **inputs**, never the labels or the outcomes, so
no outcome-dependent quantity enters the null).

### 3.1 Statistic with a nuisance surface

Use a rank-space partial association: rank-transform both fields, project out a flexible control
surface `[1, z_c, z_c²]` where `z_c` is the rank-transformed nuisance variable, then take the
Spearman correlation of the residuals. Here `z_c` is vessel distance, because both operators are
anchored to the vasculature. This removes the *linear-and-quadratic* component of mechanism (i)
directly, before any null is built. Implemented as `sparta.spatial_stats.partial_spearman` /
`control_projector` (vectorised `partial_spearman_many` for surrogate batches).

### 3.2 Geometry null — isolates (i) beyond the control surface

Recompute the second operator `f_B` from **graph-spectral surrogates of all of its own inputs**, keep
`A` and the *real* first field `f_A`, and re-measure. The surrogates preserve each input's graph
power spectrum (sign randomisation of `|Vᵀx|` on the normalised-Laplacian eigenbasis `V`) and are
rank-normalised back to the same marginal distribution. What remains is coupling produced by the
graph, the anchor set and the operator's functional form alone.

### 3.3 Construction null — isolates (i)+(ii): the part the model builds itself

Recompute `f_B` from the **real shared input** but **surrogate operator-specific inputs**. Everything
that is genuinely shared between the two operators is retained; everything that is specific to `f_B`
is randomised while keeping its spectrum. The mean of this null is the coupling the *construction*
produces. Define:

- **share reproduced by construction** = `E[ρ_construction_null] / ρ_observed`
- **excess over construction** = `ρ_observed − E[ρ_construction_null]` (the tissue-specific remainder)

The excess, not the observed `ρ`, is the quantity that may be interpreted biologically.

### 3.4 Input ablation — rules out "it is just one shared scalar"

Zero out the coupling channel at the *parameter* level rather than the input level: for SPARTA,
`b_ecm = 0` and `c_caf = 12` in the min-cut capacity and `lam = 0` in the conductance, so the shared
ECM term no longer modulates either operator through those channels. Re-test the association with the
same spatial surrogate test used for the primary analysis. **Retained fraction** = `ρ_ablated / ρ_full`.
A high retained fraction means the association does not depend on the specific shared-scalar channel
one removed — it must enter through other routes (or be real).

### 3.5 Surrogate design, calibration and multiplicity

- **Surrogate family.** Graph-spectral sign randomisation is used for the *field-level* spatial null
  (preserves the graph spectrum and hence spatial autocorrelation; `spectral_sign_null`). A
  spot-level test is *not* a substitute: on the 22 real section graphs the graph-spectral surrogate
  held the false-positive rate near nominal (pooled 0.047 at α = 0.05), whereas a spot-level test
  reached 0.246 (`run_37`). For strongly skewed fields use the normal-score variant
  (`sparta.spatial_stats.normal_scores`; pooled 0.064 vs 0.150 for the raw spectrum).
- **Draw counts.** 500 draws for the primary spatial surrogate test; 200 draws per input-level null
  (the input-level null means are stable well before that). Fix the seed.
- **Multiplicity.** Apply Benjamini–Hochberg within the declared family. Report the family
  explicitly — different families are different claims. Report per-section results *and* a pooled
  estimate; never quote a section count as a sample size when sections nest within patients
  (use the patient-level model, `run_36`).

---

## 4. What the framework does and does not license

**Licenses**

- "After removing geometry and construction, X% of the observed coupling remains, with excess
  `ρ = …` (95% CI …), positive in *k/n* sections."
- "The coupling is not an artefact of the single shared scalar channel: retained fraction after
  ablating that channel is … ."
- "Equivalent implementations of the same operator (e.g. with and without Scanpy) reproduce the same
  association to within `Δρ = …`, so the result is not a library artefact."

**Does not license**

- Any claim that a model output measures tissue permeability, drug exposure or treatment response.
  The framework separates *sources of coupling inside a model*; it cannot validate the model.
- Causal or mechanistic statements about the residual excess. A positive excess says the operators'
  inputs are co-arranged in the sample beyond what construction explains — no more.
- A claim that the residual is large. Report it against the observed value, and report the CI.

---

## 5. Worked example (SPARTA, 19-section main cohort)

Shared input: ECM score (in min-cut capacity **and** diffusion conductance). Operators: source–sink
minimum cut (`B_cell` field version) vs screened-Poisson diffusion–absorption (`B_mAb`). Nuisance:
vessel distance. Numbers as stored in `results/validation/*.json` (`run_39`, 2026-10-05 seed):

| Quantity | Value |
|---|---|
| observed partial `ρ` (median across sections) | **0.217** |
| geometry-null mean (median) | 0.046 |
| construction-null mean (median) | 0.116 |
| share reproduced by construction (median) | **0.597 (≈60%)** |
| excess over construction (median) | 0.081 (positive in 16/19; p<0.05 in 4/19; BH q<0.05 in 2/19) |
| pooled tissue-specific excess (patient-level) | 0.07 (95% CI 0.01–0.14) |
| retained fraction after ECM ablation (median) | **0.751** (17/19 positive; 11/19 p<0.05; 7/19 BH q<0.05) |
| median input `ρ(ECM, CAF)` | 0.873 |

Reading: ~60% of the observed field coupling is reproduced by the model's own construction from a
shared input plus graph geometry; the tissue-specific remainder is small but positive and concentrated
in the cSCC sections. The four melanoma sections from one patient remain unresolved (negative/zero
excess), and are reported as such rather than pooled away.

**This is the sentence the framework supports:** *the two model fields co-vary in most sections, most
of that coupling is produced by shared inputs and graph geometry, and a small tissue-specific excess
remains that the construction alone does not reproduce.* It is not the sentence "the two barriers are
independent mechanisms that the tissue couples".

---

## 6. Reuse: applicability checklist and pseudocode

Required:

1. ≥2 operators that run on one shared graph, each accepting explicit per-node inputs.
2. Every input is a per-node vector with a well-defined marginal (so surrogates can be rank-matched).
3. At least one input is *shared* between the operators (otherwise (ii) is trivially zero).
4. A nuisance/anchoring variable measurable per node (distance to source/anchor).
5. Operators must be re-runnable with substituted inputs (no hidden state, no frozen learned encoder
   that cannot accept a substituted input). If an operator is an opaque learned model, the framework
   applies only if the *input layer* can be substituted while keeping the rest fixed; otherwise (ii)
   cannot be isolated and should be declared as a limitation.

```
V, w  = eigendecomposition(normalised_laplacian(A))         # spectrum
x_sur = rank01( V @ (signs * |Vᵀ x|) )                     # spectral surrogate, same marginal

rho_obs = partial_spearman(f_A, f_B, control=anchor)

# geometry null  -> mechanism (i)
rho_geo = mean over draws of partial_spearman(f_A, B(all inputs surrogate), anchor)
# construction null -> mechanisms (i)+(ii)
rho_con = mean over draws of partial_spearman(f_A, B(shared real, other inputs surrogate), anchor)

share_built_in = rho_con / rho_obs
excess         = rho_obs - rho_con
# ablation -> is coupling carried by the one shared scalar channel?
rho_abl = partial_spearman(A(params with shared channel off), B(params with shared channel off), anchor)
retained = rho_abl / rho_obs
# plus: field-level spectral surrogate p, BH within the declared family, per-section + pooled
```

---

## 7. Relation to existing approaches

The framework sits between three established ideas and adds one specific twist.

- **Partial correlation / nuisance conditioning** (3.1) is standard; it decomposes only the part of
  coupling expressible in the chosen control basis.
- **Surrogate-data testing** (Theiler et al. 1992; constrained realisations) supplies the surrogate
  construction used in (3.2)–(3.3), and its requirement that surrogates preserve relevant structure.
- **Variance partitioning / common-factor audits** in ecology and epidemiology share the goal of
  attributing covariation to sources, but partition *observed* variables.
- **The twist:** the nulls are *pushed back through the model*. Rather than permuting outputs or
  labels, the framework substitutes inputs and re-runs the operators, so the null is defined in the
  same functional space as the observed statistic and the "share built in" is a property of the model
  construction, measurable on any dataset — including datasets where the tissue signal is absent.
  This is why it works as an internal validity check that requires no ground truth.

Note the honest boundaries: the decomposition is only as good as the surrogates' structure
preservation, and mechanism (iii) is inferred by *exclusion* — it is reported as "not explained by
(i)+(ii)", not measured directly.

---

## 8. Reporting template

For each section, store and report: `ρ_obs`, `ρ_geo` mean±sd and its p, `ρ_con` mean±sd and its p,
`share_reproduced_by_construction`, `excess_over_construction`, `retained_fraction` and its p, plus
the input cross-correlations. Then a summary block with medians, counts of positive/passed sections,
the pooled/patient-level estimate with CI, and an explicit statement of which sections are
unresolved. Report the null calibration (`run_37`) alongside, so the reader knows the FPR is near
nominal. This is exactly the structure of `results/validation/geometry_null.json` and
`shared_input_spatial_null.json`.

---

## 9. Limitations

1. **Model-internal.** The framework partitions coupling *inside* the model. It cannot show the model
   is right, only that an observed coupling is not entirely self-generated.
2. **Surrogate-dependent.** If surrogates fail to preserve a structurally relevant property, the null
   is biased. The calibration study (`run_37`) bounds this for the field-level test only.
3. **Small residual.** In SPARTA the tissue-specific excess is small (0.07, CI 0.01–0.14) and
   significant in a minority of sections after correction. It must be reported as small.
4. **Unresolved sections.** One patient's melanoma sections do not follow the pattern; the framework
   does not explain why, and they must not be silently pooled.
5. **Not causal.** Excess over the construction null is co-arrangement of inputs, not a mechanism.

---

## 10. Paper framing options

The current manuscript presents this as a robustness analysis ("how much of the coupling is built
in?"). Three ways to raise it to a methodological contribution, in increasing order of restructuring:

- **Option A (lowest cost, framing only).** Move it into the title/abstract and give it a name, e.g.
  *"Construction nulls"*; add one paragraph to Methods formalising (i)–(iii) and a dedicated Results
  subsection. Keeps the current narrative: minimal-risk.
- **Option B (moderate).** Add a short "Applicability" section stating the four requirements and the
  pseudocode, with one sentence on which other multi-operator spatial pipelines could adopt it. Turn
  §6–§7 of this note into that section.
- **Option C (highest cost).** Reframe the paper around the framework with SPARTA as the worked
  example; the two operators then serve as the demonstration that coupling between model outputs is
  largely self-generated unless tested for. Requires rewriting the introduction and re-ordering
  results; the risk is that reviewers ask for a second worked example on a different pipeline.

Recommended for this submission: **A + B**, keeping the SPARTA results as they stand.

---

## 11. Pointers

- Implementation: `scripts/run_39_coupling_decomposition.py`; statistics in `sparta/spatial_stats.py`;
  tests in `tests/test_is_revision.py`.
- Outputs: `results/validation/geometry_null.json` (= construction / geometry nulls),
  `results/validation/shared_input_spatial_null.json` (= ECM ablation + spatial null).
- Null calibration: `scripts/run_37_null_calibration.py` → `results/validation/null_calibration.json`.
- Primary association and its spatial null: `scripts/run_27_spatial_null.py`.
- Patient-level aggregation: `scripts/run_36_patient_level.py`.
