---
title: "Coupled, not separable: the T-cell migration barrier and the antibody mass-transport barrier in cutaneous tumours share a matrix origin"
---

**Working draft, assembled 2026-08-28 from the verified artefacts.**
Author list, affiliations and the alternative tool-forward title are still open.
Every numeric claim carries its source artefact in `monospace brackets`; **strip
those before submission** — they are provenance for the authors, not manuscript text.

# Abstract

## Motivation

Immune checkpoint blockade fails in most patients with advanced cutaneous
malignancy, and spatial transcriptomics has made the microenvironmental
contribution to that failure directly measurable. The field's analytical
vocabulary, however, addresses only one of the two transport problems involved.
Whether a cytotoxic T cell can reach the tumour nest is a cell-migration problem
through a fibre network; whether the anti-PD-1 antibody itself can reach it is a
mass-transport problem for a 5.5 nm macromolecule through a matrix whose
effective mesh is measured in tens of nanometres. The second is almost never
quantified, and whether the two barriers are separable in real tissue has not
been asked.

## Results

We present SPARTA, which poses both questions as transport problems on one
spatial graph, changing only the edge-weight semantics and the operator: T-cell
migration as a source–sink minimum cut, which returns both a barrier magnitude
and an explicit blockade geometry, and antibody delivery as a screened Poisson
diffusion–absorption field that reduces to the harmonic/effective-resistance
limit when absorption vanishes. We applied it to **19 sections from seven
patients** — six primary cutaneous squamous carcinomas profiled across two
platform generations (four Visium, eleven first-generation ST) and one melanoma
patient contributing four extracranial metastatic deposits.

The two barriers prove **not** to be separable. After residualising both on the
weighted graph distance to vasculature the partial correlation is positive in
18/19 sections and significant in 17 (ρ = −0.000 to +0.385, median +0.217), every
one of the seven patients has a positive median, and regions where one barrier is
high while the other is low occupy 1.1–6.9 % of spots against a 6.25 % chance
expectation — fewer discordant spots than independence predicts, in 18/19
sections. The result reproduces across platform generations within the squamous
cohort (median ρ = +0.270 on Visium, +0.207 on first-generation ST), which
separates assay generation from the tumour-type and disease-stage covariates that
were confounded with it in a smaller design.

Removing the matrix term shared by the two operators leaves the association
significant in 15/19 sections — 12/15 squamous sections (a majority in five of
the six squamous patients) with a median 81 % of its original magnitude, and
3/4 melanoma sections with a median 38 % — unchanged whether the geometric control
is a weighted path or a hop count. The single melanoma patient is the one stratum
whose classification depends on that choice, and we report it as unresolved
rather than as a cohort contrast.

Counterfactual experiments support the topological reading: a contiguous gap in
the blockade leaves 1.04–1.39× more residual barrier than removing the same
material scattered within it (16/19 sections, five of seven patients with every
section significant, at a
pre-specified 20 % removal fraction; the effect is larger at 30 % than at 5 % in
15/19), and sweeping molecular radius on fixed tissue raises the tumour-core
barrier monotonically across an order of magnitude. Relative to a spatial-domain
segmentation, just 4–22 % of domain boundaries lie on the cut on either platform
— the transport formulation selects and quantifies the one boundary that limits
flux rather than detecting boundaries.

## Availability and implementation

SPARTA is implemented in Python and runs on CPU; the barrier operators depend
only on NumPy, SciPy and NetworkX. Source code, the analysis pipeline and the
scripts reproducing every figure are available at <GitHub URL> and archived at
<Zenodo DOI>. All data are public (GEO: GSE250636, GSE144239).

---

# Introduction

Immune checkpoint blockade has changed the outlook for cutaneous malignancies,
yet most patients with advanced melanoma or cutaneous squamous cell carcinoma
still fail to obtain durable benefit from anti-PD-1 therapy. Response is only
partly explained by tumour-cell-intrinsic features such as mutational burden or
antigen presentation; a substantial share of variance is attributed to the
organisation of the tumour microenvironment. Spatial transcriptomics has made
that organisation directly measurable, and studies across melanoma, squamous
carcinoma and other solid tumours now routinely describe stromal bands,
immune-excluded phenotypes and fibroblast-rich boundaries around tumour nests.

The dominant framing of this literature is cellular: whether cytotoxic T cells
can reach the tumour nest. Cancer-associated fibroblast bands, dense collagen
and TGF-β-driven programmes have all been implicated, and the analytical
vocabulary is largely that of local spatial statistics — neighbourhood
enrichment, co-localisation scores, Ripley's functions, and spatial-domain
segmentation. These measures answer "what is next to what". They do not answer
"what gets through", and they are, by construction, insensitive to a
distinction that matters physically: a continuous, closed band of fibroblasts
blocks passage, whereas the same number of fibroblasts scattered through the
tissue does not. Connectivity is a topological property, and a local statistic
cannot see it.

A second barrier has received far less attention. An anti-PD-1 antibody is not a
cell: it is a ~150 kDa macromolecule with a hydrodynamic radius near 5.5 nm that
must leave the vasculature and diffuse through interstitial matrix whose
effective mesh size in tumours is measured in tens of nanometres, and that can
be consumed en route by binding to its target. Whether the antibody itself
reaches the tumour nest is a mass-transport problem with different physics from
cell migration — different length scale, different obstruction mechanism,
different dependence on matrix crosslinking — and it is almost never quantified
in spatial-omics analyses. Whether the two barriers are separable in real tissue
has, to our knowledge, not been asked.

We introduce SPARTA, a framework that poses both questions as transport problems
on the same spatial graph, changing only the edge-weight semantics and the
operator. T-cell migration is formulated as a source–sink minimum cut between
immune-entry and tumour-core compartments, which yields both a scalar barrier
strength and an explicit blockade geometry. Antibody delivery is formulated as a
screened Poisson diffusion–absorption problem in which conductance decreases with
matrix density and with size exclusion relative to the local mesh, and in which
target antigen acts as a distributed sink; when absorption vanishes this reduces
to the familiar harmonic/effective-resistance limit, so the graph-transport
reading is preserved. We apply the framework to 19 sections from seven patients
across two independent cutaneous cohorts and two platform generations, decompose
the antibody barrier into its shared and modality-specific components, and test
directly whether the two barriers dissociate.

They do not. After controlling for distance from the vasculature, the two
barriers remain positively associated in every section, and regions where one
barrier is high while the other is low are *rarer* than chance. The dominant
axis of the antibody barrier is a matrix property — size exclusion set by
crosslinking — rather than an antigen property. We therefore report a coupling
rather than a dissociation, and we argue that the coupling is the more
actionable finding: it predicts that matrix-directed intervention should relieve
both modalities at once, whereas modality-specific intervention should not.

---

# Results

## R1. Nineteen sections, seven patients, two tumour types, two platform generations

We assembled 19 spatial transcriptomic sections from two independent, publicly
available cohorts of cutaneous malignancy: four extracranial metastases of
cutaneous melanoma (MEL01–MEL04, GSE250636) and fifteen treatment-naive primary
cutaneous squamous cell carcinomas (CSCC01–CSCC16 minus CSCC13, GSE144239;
Ji et al., 2020).

**These 19 sections come from seven patients, and every count in this paper is
reported at both levels.** Sections from one patient are not independent
observations, so section-level counts are always reported alongside the
patient-level result, and n = 19 is not treated as a sample size anywhere. The
patient assignment was taken from the GEO sample metadata, not inferred, and is
recorded in `data/ledger.csv` (`patient`, `replicate`).

| Cohort | Platform | Sections | Patients | Nodes after QC | Median UMI / spot |
|---|---|---|---|---|---|
| Melanoma, metastatic | Visium | 4 | 1 (Patient B) | 801 – 2 673 | 1 488 – 5 514 |
| cSCC, primary | Visium | 4 | 2 (P4, P6) | 652 – 2 645 | 636 – 16 686 |
| cSCC, primary | First-generation ST | 11 | 4 (P2, P5, P9, P10) | 370 – 672 | 567 – 7 219 |

`[data/ledger.csv; data/interim/*.graph.npz]`

The four melanoma sections are four separate extracranial metastatic deposits
(sternum, cecal nodule, right upper chest wall, ribcage) from a single patient.
GSE250636 contains a second patient, but that patient contributes only
leptomeningeal deposits — a different anatomical compartment, behind the
blood–brain barrier and without dermal stroma — so including them would exchange
n = 1 for an n = 2 in which patient and compartment are perfectly collinear. We
therefore did not include them, and we treat the melanoma cohort as n = 1
patient throughout.
`[GSE250636: GSM7983359/364/365/366 = Patient B extracranial; GSM7983358 =
Patient B leptomeningeal; GSM7983360 = Patient A leptomeningeal]`

**The squamous cohort deliberately spans two platform generations.** CSCC01–04
are early Visium sections from two patients; CSCC05–CSCC16 are first-generation
ST sections (200 µm array pitch, 2016-era chemistry) from four further patients,
three technical replicates each. This is the design feature that lets R3b ask
whether any conclusion is a property of the tissue or of the assay generation,
and it is the only one of the three covariates that moved together in an earlier
eight-section design (tumour type, disease stage, platform generation) that the
extended cohort separates.

**One section was excluded at admission.** CSCC13 (P9, replicate 3) failed
criterion C5 with a median of 289.5 UMI per spot against a cohort threshold of
300. The threshold was fixed, with its date and rationale, before the sections
were ingested; CSCC13 sits alone below a near two-fold gap (the next lowest
section is 567), so the cut does not slice through a dense region. It is recorded
as `status = rejected` in the ledger with the evidence retained in
`data/interim/CSCC13.admission.json`. Exclusion costs one technical replicate and
no patient: P9 is still represented by CSCC11 and CSCC12.
`[data/interim/CSCC13.admission.json; data/ledger.csv]`

### Graph construction differs between platforms only in geometry

Spatial graphs were built by radius adjacency after rescaling coordinates so that
the **median nearest-neighbour distance equals the platform's nominal pitch** —
100 µm for Visium, 200 µm for first-generation ST. Rescaling to the observed
median rather than assuming a coordinate convention matters here, because the
first-generation ST array is **staggered**: array indices satisfy x + y ≡ 0
(mod 2), so nearest neighbours lie on the diagonal and one array index unit
corresponds to 141.4 µm, not 200 µm. Assuming otherwise would inflate every
micrometre-scale quantity by √2, including the distances compared against the
literature-anchored oxygen diffusion limit `d0` = 130 µm.

With `radius_um` = 150 the Visium graphs recover the expected hexagonal
connectivity (mean degree 4.92–5.74 after quality control). With
`radius_um` = 300 the first-generation ST graphs connect both the 200 µm nearest
neighbours and the 282.8 µm second neighbours (mean degree 5.94–7.50 after
quality control; 7.2–7.6 on the full on-tissue set, before QC removes low-count
spots). The choice of 300 over 250
is not a connectivity argument — both values give connected graphs — but a
geometric one: relative to straight-line distance, the hop approximation
`hops × pitch` over-estimates by a median factor of 1.265 at `radius_um` = 250
and under-estimates by 0.894 at 300, and `d0` is anchored to a straight-line
diffusion limit. The analysis reported here does not rely on that approximation
at all (see Methods M4), but the graph geometry is still chosen to make it as
close to unbiased as possible.
`[configs/cscc_legacy_st.yaml; data/interim/*.graph.npz meta]`

Quality control removes spots below 500 UMI, which affects the shallow
first-generation ST sections disproportionately: CSCC11 loses 41 % of its
on-tissue spots (1 145 → 672) and CSCC14 loses 39 % (608 → 370), while the deep
sections lose almost none (CSCC05, 666 → 664).

**Removing spots fragments the graphs, and we report the fragmentation rather
than describing the graphs as connected.** After QC the analysed graphs contain
1–54 connected components, with **90.0–100 % of nodes in the largest** (Visium
1–54 components, single-component in 1/8 sections; first-generation ST 1–16
components, single-component in 6/11). The most fragmented sections are the two
large, shallow Visium ones (CSCC03: 50 components, 94.2 % in the largest;
CSCC04: 54, 90.0 %) and the two shallowest first-generation ST sections
(CSCC11: 15, 94.5 %; CSCC12: 16, 95.1 %).

Fragmentation has a direct consequence for the min cut that we state rather than
hide: a source or sink stranded in a small component contributes no flow to the
main component. In 8 of 19 sections some source or sink nodes lie outside the
largest component — at most 23/187 sources (CSCC04, 12 %) and 17/232 sinks
(MEL02, 7 %). The maximum-flow computation therefore runs effectively on the
largest component plus any self-contained fragments, which biases `B_cell`
(= 1/max-flow) slightly upward in the fragmented sections. The operator already
counts and reports the affected nodes (`n_unreachable` in
`data/interim/*.barrier_meta.json`: 49 and 55 in CSCC03/CSCC04, 0 in the six
singly-connected first-generation ST sections).

Depth is therefore correlated both with final graph size and with connectivity
(Spearman ρ = +0.45 between median UMI and the largest-component fraction), and
all three are checked as possible drivers wherever a result could be sensitive to
them (R3b, R4-S1).
`[data/interim/*.graph.npz; data/interim/*.barrier_meta.json]`

---

## R2. Decomposing the antibody barrier: one axis dominates, and its share is set by a scale parameter

For every section we decomposed the variance of the antibody barrier `B_mAb`
into three ablation-defined components: the matrix channel shared with the
T-cell barrier (`ECM_core`), the size-exclusion channel specific to
macromolecules (`ECM_crosslink`), and the binding-site channel specific to the
target antigen (`CD274`/`PDCD1LG2`).

Across all 19 sections the size-exclusion channel accounts for **95.9–98.9 %**
of `B_mAb` variance (median 97.5 %) and the shared matrix channel for
**0.8–2.8 %** (median 1.8 %). Crosslinking and core matrisome scores were not
collinear in any section (|r| = 0.010–0.512, 0/19 above the 0.80 threshold at
which the estimator conservatively merges the two buckets), so the separation is
not an artefact of the ablation design. The resulting ratio of modality-specific
to shared variance is **35.1–126.5 (median 55.1)**. The decomposition is stable
across both platform generations: no section of either platform falls outside
these ranges.

**The antigen channel could only be measured on 15 of the 19 sections.** In
CSCC10, CSCC14, CSCC15 and CSCC16 — all first-generation ST, all with the lowest
detected gene counts in the study (15 383–17 399) — neither `CD274` nor
`PDCD1LG2` passed the minimum-gene-match threshold, so the `Ag_target` signature
was neutral-filled and its variance contribution is exactly zero. **That zero is
missing data, not a measurement**, and we exclude those four sections from the
antigen range rather than reporting them as 0 %. On the 15 sections where the
signature exists, the antigen channel accounts for **0.47–3.06 % (median 0.96 %)**
of `B_mAb` variance.
`[results/validation/screen_decision.json — var_antigen = 0 exactly in
CSCC10/14/15/16; all 19 sections for the other two channels]`

The mechanism behind that dominance is now measured directly rather than
inferred. At the default parameters the effective mesh size falls below the IgG
hydrodynamic radius (5.5 nm) on the **median edge of every section**: the
median `ξ` is **4.36–4.68 nm** across the 19 sections, and **60.0–65.2 %** of
graph edges (median 62.7 %) are size-excluded for an IgG-sized molecule. The
antibody barrier is therefore not a diffuse gradient but a largely
percolation-limited field, and the small fraction of spots that saturate the
conductance floor (0.0–7.8 %, median 1.0 %) confirms the numerical guard is
rarely binding.
`[results/validation/review_diagnostics.json, run_11_review_diagnostics.py]`

**This number must be reported together with its parameter dependence.**
The size-exclusion term enters through an effective mesh size
`ξ = ξ₀·exp(−β·x)`, where `x` is the *within-section rank-normalised*
crosslinking score. Holding `λ = 3`, the crosslinking share of `B_mAb` variance
moves monotonically with `β`, and it does so almost identically on every section
we scanned — five sections spanning both tumour types and both platform
generations (MEL01, MEL03, CSCC01, CSCC03, CSCC05):

| β | 0.25 | 0.5 | 1.0 | 1.5 | 2.0 | **3.0 (default)** | 6.0 | 12.0 |
|---|---|---|---|---|---|---|---|---|
| crosslinking share of `B_mAb` variance | 0.2–0.3 % | 0.9–1.7 % | 9.1–17.3 % | 50.3–75.0 % | 81.2–85.4 % | **95.9–97.7 %** | 94.5–97.8 % | 81.0–96.0 % |

`[results/validation/bmab_sensitivity.json, grid lam_x_beta, column λ = 3]`

That the five sections agree to within a few percentage points at every value of
β is the point: the share is a property of the scale parameter, not of the
tissue. We therefore do **not** claim that "crosslinking explains 96 % of
antibody barrier variance in skin cancer" as an empirical finding about tumour
biology. The defensible statement is narrower and is the one we make: *given a
size-exclusion model steep enough that mid-range crosslinking already contracts
the mesh below the IgG hydrodynamic radius — which, at the default β, it does on
roughly 63 % of edges — the resulting barrier field is dominated by that axis
rather than by bulk matrix density or by target-antigen availability.*
Methods M2 states exactly what `β` and `ξ₀` are and are not.

The antigen channel is weak under every parameter setting we examined
(0.3–17.9 % across the full grid). We caution against reading this as evidence
that antigen sinks are unimportant in vivo — and the four sections where the
signature could not be scored at all make the point better than the weak values
do. The `Ag_target` signature contains only two genes (`CD274`, `PDCD1LG2`),
which is at the floor of our minimum-gene-match threshold; a two-gene score is
close to the noise floor, and on the shallower first-generation ST sections the
two genes are sometimes not detected at all. The honest reading is that **this
framework cannot resolve the binding-site barrier with the ICB target set**, not
that the binding-site barrier is small. An ADC-style target with higher and more structured
expression (e.g. `ERBB2`, `TACSTD2`) would be the natural test, and the operator
requires no change to run it (Discussion).

---

## R3. The two barriers do not dissociate

The central question of this work is whether the cellular and macromolecular
barriers are separable in real tissue. They are not.

Per-spot barrier fields were computed as the minimum accumulated migration cost
from immune-entry nodes (`B_cell` field) and as the steady-state screened-Poisson
concentration deficit from vessel nodes (`B_mAb`). Because both fields accumulate
away from the vasculature, they share a strong geometric component; we therefore
residualised both on distance to the nearest vessel (quadratic, in rank space)
before correlating them. That distance is a weighted shortest path along the
graph, in micrometres — not a hop count (Methods M4).

**After controlling for vessel distance, the two barriers remain positively
correlated in 18 of 19 sections and significantly so in 17** (Spearman ρ_partial
= −0.000 to +0.385, median **+0.217**). The two exceptions are CSCC15
(ρ_partial = −0.0004, p = 0.99) and MEL02 (+0.055, p = 0.12); neither is
negative. Uncontrolled correlations are positive in all 19 (ρ = +0.052 to +0.429,
median +0.246).
`[results/validation/decoupling.json and results/validation/shared_ecm_check.json,
per-section ρ_partial and p_partial]`

**The result holds at the patient level.** All seven patients have a positive
median ρ_partial: +0.323 (P4), +0.280 (P2), +0.248 (P6), +0.239 (P10), +0.209
(P9), +0.146 (Patient B, melanoma), +0.141 (P5). No patient contributes a
cohort-level exception.
`[results/validation/paper_stats.json, field patients]`

**And it holds across platform generations.** Within the squamous cohort the four
Visium sections (two patients) give a median ρ_partial of +0.270 and the eleven
first-generation ST sections (four patients) give +0.207 — the same sign, the
same order of magnitude, on assays a platform generation apart.
`[results/validation/shared_ecm_check.json, stratified by data/ledger.csv platform]`

The "dissociation zone" — spots in the lowest quartile of `B_cell` and the
highest quartile of `B_mAb` — occupies **1.1–6.9 % of spots (median 3.7 %)**
after geometric control. Under independence of the two barriers this fraction is
(1 − q)² = 6.25 % by construction at q = 0.75. Observed values are therefore
**0.18–1.10× the chance expectation (median 0.59×), below chance in 18 of 19
sections**: there are *fewer* discordant spots than random, which is a second,
independent expression of the positive association. The single section above
chance (CSCC15, 1.10×) is the same section whose partial correlation is not
positive — the two readouts agree on which section is the exception.
`[results/validation/decoupling.json; chance baseline asserted by
tests/test_statistics.py::test_discordant_fraction_has_a_chance_baseline]`

This result is not a failure of measurement. The framework's own diagnostic
predicted it: the core matrisome score enters both the min-cut edge capacity and
the diffusion edge conductance, so any change in matrix density moves both
barriers in the same direction, and this channel cannot be removed by statistical
control. R3b removes it by construction instead. What the decomposition in R2
adds is that the coupling is **not** merely geometric — after depth is controlled
the association survives — and that the dominant axis of the antibody barrier is
a matrix property (size exclusion) rather than an antigen property.

The conclusion we draw is therefore:

Whether that co-location is a property of the tissue or of our edge-weight
definition is a separate question — see R3b.

---

## R3b. Is the coupling a property of the tissue or of the edge weights?

The core matrisome score enters both the min-cut edge capacity and the diffusion
conductance, so part of the association reported in R3 is guaranteed by
construction. We therefore repeated the analysis with the shared input removed
entirely: the cellular barrier recomputed from the fibroblast signature alone
(`b_ecm` = 0) and the antibody barrier from crosslinking and antigen alone
(`λ` = 0), so that the two operators share no input variable. The geometric
control was computed from the main configuration in both arms, so the same
quantity is being adjusted for.

Across all 19 sections the median partial correlation falls from **+0.217 to
+0.152**; 17/19 remain positive and 15/19 remain significantly positive.

**This section reports a stratified result and deliberately issues no single
verdict.** The decision rule used in earlier drafts required *every* section in a
cohort to remain significant, a threshold whose stringency rises monotonically
with n: adding sections can only ever move the analysis away from the
"tissue-property" label, irrespective of what those sections show. On eight
sections that rule returned SPLIT (on the strength of "CSCC 4/4"); on nineteen it
returns MODEL, with no evidence reversal in between — only a change in n. Rather
than replace one threshold with another after seeing which label it produces, we
report the effect sizes, retention fractions and significance counts by cohort
and by patient, and let the reader apply their own criterion.
`[results/validation/shared_ecm_check.json, field verdict_status = "deprecated_2026-08-27"]`

| Stratum | sections / patients | median ρ_partial (main) | after removal | median retention | significant & positive | patients with a majority significant |
|---|---|---|---|---|---|---|
| cSCC — Visium (P4, P6) | 4 / 2 | +0.270 | **+0.214** | — | **4 / 4** | 2 / 2 |
| cSCC — first-generation ST (P2, P5, P9, P10) | 11 / 4 | +0.207 | **+0.137** | — | **8 / 11** | 3 / 4 |
| cSCC — combined | 15 / 6 | — | **+0.183** | 0.81 | **12 / 15** | **5 / 6** |
| Melanoma (Patient B) | 4 / 1 | — | +0.056 | 0.38 | 3 / 4 | 1 / 1 |

`[results/validation/shared_ecm_check.json]`

**Platform generation is no longer confounded with cohort.** The squamous cohort
is now observed on two platform generations — four Visium sections from two
patients and eleven first-generation ST sections from four further patients — and
the two agree in direction and in magnitude both before removal (median ρ_partial
+0.270 vs +0.207; dissociation-zone enrichment 0.47× vs 0.63× of chance) and
after it (+0.214 vs +0.137). Of the three covariates that moved together in the
eight-section design — tumour type, disease stage and platform generation — this
is the one the extended cohort separates, and it is not the driver.

**Sequencing depth and section size are not the driver either.** Across the 19
sections the partial correlation is essentially independent of both (Spearman
ρ = +0.188 with median UMI per spot, ρ = +0.030 with section size), even though
those two quantities span more than an order of magnitude (median UMI 290–16 686;
462–3 838 spots).
`[results/validation/paper_stats.json + data/ledger.csv]`

**The melanoma side cannot be decided by this design.** All four melanoma
sections are metastatic deposits from a single patient, so the cohort contributes
n = 1 at the level at which we report every other count in this paper. It is also
the only stratum whose classification depends on the geometric control: with the
weighted graph distance 3/4 sections remain significant after removal, with the
hop approximation 1/4 (Table S·x). We therefore state the squamous result and
report the melanoma result as unresolved, rather than presenting the contrast
between them as a finding. The cohort contrast asserted in earlier drafts rested
on "CSCC 4/4 versus MEL 1/4"; the squamous side of that contrast does not survive
the extension to 15 sections (12/15), so the contrast itself was an artefact of
the small design.

**Sensitivity to the geometric control.** `d_vessel_um`, the covariate the partial
correlation adjusts for, is a weighted shortest path along the graph (edge weights
in µm). Re-running the entire analysis with the previous hop-count approximation
(`hops × spacing_um`) changes no stratified conclusion for the squamous cohort —
median ρ_partial after removal +0.176 versus +0.183, 11/15 versus 12/15
significant, the same 5/6 patients — and changes the melanoma classification as
described above. Three sections in total change significance status between the
two conventions (CSCC06, MEL03, MEL04), all marginal in both
(p = 0.116/0.188/0.144 under hops, 0.0073/0.016/0.0057 under the weighted
distance).
`[results/validation/shared_ecm_check_hops.json]`

**The negative part of the conclusion is robust to everything tested.** The
dissociation-zone fraction stays below the chance expectation of (1 − q)² = 6.25 %
in 18/19 sections under the main parameterisation and in **19/19** after the
shared input is removed (median enrichment 0.595× and 0.589× of chance
respectively; 0.626× and 0.647× under the hop approximation, 18/19 in both). The
two barriers do not dissociate, and that statement does not depend on the
edge-weight definition, the geometric control, the platform generation, the
sequencing depth or the section size.

---

## R4. Counterfactual experiments

Three counterfactual experiments probe, respectively, whether the barrier comes
from spatial arrangement (S1), whether continuity of the blockade matters (S2),
and whether the two barriers respond differently to molecular size (S3). All
three are cheap to run and none requires external data.

### S1 — spatial rearrangement

Keeping each spot's scores intact and permuting only their spatial positions
destroys arrangement while leaving composition exactly unchanged. All 19 sections
use 500 permutations, so the empirical p-value resolution (1/501 = 0.002) is
identical across cohorts and platform generations.

Under the `fixed` mode — source and sink held in place, so the test isolates the
arrangement of resistance material — the real `B_cell` exceeds the permutation
null by **0.89× to 4.49×**, and **12 of 19 sections are significant** after
Benjamini–Hochberg correction.
`[results/counterfactual/*.json, n_perm = 500; results/validation/paper_stats.json]`

**The result is organised by patient, not scattered across sections.** Six of the
seven patients are internally consistent: P2 3/3, P10 3/3, P4 2/2 significant;
P6 0/2, P9 0/2 non-significant; melanoma Patient B 3/4. Only P5 is mixed (1/3).
The effective sample size for S1 is therefore closer to seven than to nineteen,
and we report it that way.

| Patient | Cohort / platform | sections significant (`fixed`) | median effect ratio |
|---|---|---|---|
| CSCC_P2 | cSCC / first-gen ST | 3 / 3 | 3.34 |
| CSCC_P4 | cSCC / Visium | 2 / 2 | 2.89 |
| CSCC_P10 | cSCC / first-gen ST | 3 / 3 | 1.57 |
| MEL_PtB | Melanoma / Visium | 3 / 4 | 1.69 |
| CSCC_P5 | cSCC / first-gen ST | 1 / 3 | 1.10 |
| CSCC_P6 | cSCC / Visium | 0 / 2 | 1.10 |
| CSCC_P9 | cSCC / first-gen ST | 0 / 2 | 0.98 |

**S1 is confounded with sequencing depth and with graph connectivity, and we
state this rather than interpret around it.** Across the 19 sections the `fixed`
effect ratio correlates with median UMI per spot at Spearman **ρ = +0.721** and
with the largest-connected-component fraction at **ρ = +0.460**, while showing
essentially no relationship with graph size (ρ = −0.023 against node count,
−0.153 against raw spot count). The two confounds cannot be separated in this
design: depth and connectivity are themselves correlated at ρ = +0.453, because
quality control removes low-count spots and stranded neighbourhoods fragment
(R1).

The three patients with no significant section (P6, P9) or almost none (P5) are
therefore not a tumour-biology subgroup; they sit, with CSCC03/CSCC04, at the
shallow and most fragmented end of the range. Two mechanisms plausibly
contribute and we do not attempt to rank them: signature scores are noisier at
low depth, and a noisier resistance field is closer to its own spatial
permutation by construction; and a fragmented graph both loses source/sink nodes
to stranded components and admits fewer alternative paths, compressing the
difference between the real arrangement and a permuted one.

**We therefore do not use S1 to support the central claim, and we do not present
the S1-negative patients as biologically distinct.** S1 is reported as a positive
control that succeeds where the data are deep enough and the graph connected
enough to support it.
`[Spearman computed over results/validation/paper_stats.json ×
data/interim/*.graph.npz × data/ledger.csv]`

We report effect ratios rather than z-scores as the primary statistic. The
z-scores are large (up to 86) mainly because the permutation null has small
variance; the underlying effect — a real barrier some multiple of the null — is
the interpretable quantity.

The `follow` mode, in which source and sink are re-derived from the permuted
scores, is significant in 10/19 sections and produces much larger ratios in two
sections in particular (CSCC03 8.03×, CSCC04 7.71×). **We do not use it to rescue
the negative `fixed` results.** The two modes test different null hypotheses:
`follow` destroys all spatial structure including the vascular and tumour
compartments, so its null barrier is far lower than the `fixed` null (CSCC03:
0.0137 vs 0.0917), and a section can be significant under `follow` while
genuinely insensitive to the arrangement of resistance material. Note that the
two modes do not even agree on which sections are significant — CSCC03/CSCC04 are
`follow`-only, CSCC10/CSCC14/CSCC15/CSCC16 are `fixed`-only. The modes are
reported separately and never pooled.

Sections insensitive under `fixed` also have a structural correlate — and the
correspondence is stronger on the full 19 sections than it was on the eight
Visium sections where it was first observed. Classifying each section's cut
band by its CAF content and width (continuous CAF band / non-CAF-dominated
diffuse band / signal-sparse narrow band; Suppl. S2):

- **continuous CAF band** (CAF enrichment on the cut 1.39–1.83×; 10 sections,
  both platforms): **10/10 significant under `fixed`** (z = +3.4 to +22.8);
- **non-CAF-dominated diffuse band** (enrichment 1.02–1.24×; 7 sections, all
  but MEL02 first-generation ST): only **2/7 significant** — CSCC15 and
  CSCC16, which come from the same patient (P10);
- **signal-sparse narrow band** (CSCC03, CSCC04; the two narrowest bands at
  9.3 % and 6.7 % of spots, at the two lowest Visium depths): 0/2 significant.

The association is therefore strong but not a clean one-to-one mapping: reading
the taxonomy as a predictor of `fixed`-sensitivity (continuous → significant,
diffuse/narrow → not) calls 17/19 sections correctly, the two misses being
CSCC15 and CSCC16 — diffuse yet significant, and replicates of one patient (P10).
In the other direction the mapping is perfect: no continuous band is
insignificant.
`[results/validation/cut_band_analysis.json — all 19 sections; S1 z read from
results/counterfactual/{sid}.json; Ripley radius scaled to platform pitch
(2.5× spacing)]`

### S2 — blockade continuity

S2 removes resistance material from a contiguous arc of the min-cut and compares
the result with removing the same amount of material from scattered positions
within the same cut set — same material, same amount, only contiguity differs.

The removal size must be defined relative to the blockade, not in absolute spots.
Real cut sets contain 61–422 nodes across the 19 sections; a fixed removal of
3–12 spots (the value appropriate to synthetic sections whose cut sets hold 40–76
nodes) touches only a few per cent of the blockade in the larger sections and
measures the removal fraction rather than any property of the tissue. We
therefore define `k` as a fraction of the cut set and **pre-specify 20 % as the
primary test**, reporting 5 %, 10 % and 30 % as a dose–response.

At the pre-specified 20 % level the contiguous gap leaves a residual barrier
**1.038–1.390× higher** than scattered removal of the same material (median
**1.147×**), and **16 of 19 sections are significant** after per-section BH
correction. The three that are not are marginal (CSCC08 p_FDR = 0.055, CSCC15
0.063, CSCC09 0.119). At the patient level: melanoma Patient B 4/4, P4 2/2,
P6 2/2, P2 3/3, P9 2/2 — **five of seven patients have every section
significant**; P10 2/3 and P5 1/3, and **no patient has zero significant
sections**.

| removal fraction | 5 % | 10 % | **20 %** | 30 % |
|---|---|---|---|---|
| MEL01 | 1.070 | 1.094 | **1.141** | 1.134 |
| MEL02 | 1.067 | 1.080 | **1.136** | 1.033 |
| MEL03 | 1.086 | 1.092 | **1.149** | 1.091 |
| MEL04 | 1.054 | 1.117 | **1.145** | 1.161 |
| CSCC01 | 1.142 | 1.220 | **1.274** | 1.219 |
| CSCC02 | 1.019 | 1.029 | **1.390** | 1.312 |
| CSCC03 | 1.191 | 1.036 | **1.147** | 1.205 |
| CSCC04 | 1.113 | 1.187 | **1.175** | 1.172 |
| CSCC05 | 1.177 | 1.200 | **1.306** | 1.251 |
| CSCC06 | 1.221 | 1.426 | **1.311** | 1.166 |
| CSCC07 | 1.350 | 1.276 | **1.294** | 1.372 |
| CSCC08 | 1.063 | 1.103 | **1.071** | 1.078 |
| CSCC09 | 1.055 | 1.071 | **1.038** | 1.007 |
| CSCC10 | 1.080 | 1.120 | **1.101** | 1.034 |
| CSCC11 | 1.020 | 1.057 | **1.090** | 1.147 |
| CSCC12 | 1.052 | 1.055 | **1.165** | 1.104 |
| CSCC14 | 1.041 | 1.101 | **1.357** | 1.367 |
| CSCC15 | 1.041 | 1.058 | **1.069** | 1.073 |
| CSCC16 | 1.060 | 1.172 | **1.109** | 1.201 |

The residual barrier is larger at 30 % than at 5 % in **15 of 19 sections**. A
dose–response of this shape is the falsifiable prediction of the topological
claim and is not what uncorrelated noise produces; we regard it as stronger
evidence than the significance of any single level.
`[results/counterfactual/*.json, k_mode = frac; run_12_paper_stats.py
--primary-frac 0.20. The empirical p resolution here is 1/101 = 0.010, so a
reported p of 0.010 means "stronger than every one of the 100 scattered
controls", not a smaller true p.]`

Unlike S1, **S2 is not strongly depth-dependent** and it reaches significance in
both platform generations of the squamous cohort (Visium 4/4, first-generation ST
8/11) and in the melanoma cohort (4/4). It is the counterfactual we are willing
to carry.

**Honest framing.** The effect is real and consistent but modest: removing a
fifth of the blockade as a contiguous arc leaves roughly 15 % more barrier than
removing the same fifth at random positions within it. On synthetic sections with
a designed, thin ring the same experiment gives 7.45×. The gap between the two is
itself informative — real blockade bands are thick and redundant, so no single
arc is load-bearing in the way a one-spot-wide synthetic ring is. S2 therefore
supports blockade continuity as a contributing property, and we do not rest the
central argument on it.

### S3 — molecular size scan

Sweeping the hydrodynamic radius from 0.5 nm (small molecule) to 10 nm at fixed
tissue structure raises the mean barrier in the tumour-nest core monotonically in
every section: MEL01 1.48 → 13.20, CSCC03 2.86 → 21.56, CSCC04 2.40 → 19.29
across eight radii, with no saturation at the numerical ceiling.
`[results/counterfactual/*.json]`

S3 is the cleanest evidence in the study that the two barriers are *different
physical problems even though they are driven by the same substrate*: the
identical tissue graph is nearly transparent to a small molecule and strongly
obstructive to an IgG-sized molecule, a distinction no cell-migration model
produces. It is also a deterministic consequence of the operator, not a
statistical test, and we present it as such.

---

## R5. The min-cut is not a rewrite of immune-cell density

A natural objection is that `B_cell` is an elaborate restatement of "how much
CAF / how few T cells are here". Three comparisons argue otherwise, all computed
on the full set of 19 sections.

First, the per-spot `B_cell` field correlates only moderately with the CAF
signature (Spearman ρ = **0.313–0.554**, median 0.458) and essentially not at all
with the T/NK signature (ρ = **−0.275 to +0.140**, median −0.018). The CAF
relationship is consistent across both platform generations, so the field is not
a rescaled CAF map on either.
`[results/validation/benchmark_ext.json, all 19 sections]`

Second, under spatial rearrangement the section-level `B_cell` moves in 12 of 19
sections (z > 2; range −0.68 to +24.9), whereas composition-based measures cannot
move at all: CAF density and T-cell infiltrated fraction are invariant under a
permutation of spot positions **by construction**. Their permutation nulls are
correspondingly degenerate — in 21 of 38 section×measure combinations the null
has exactly zero variance, and in the remaining 17 it takes only two values,
because ties at the quantile threshold let the count wobble by a single spot, so
the "z-score" is exactly ±1.0. We report this as a sanity check, not as a
comparison result: **a metric that cannot move under the null is not evidence
about tissue**, and its z-score is not an effect size.
`[results/validation/benchmark_ext.json: z_dens and z_tinf are NaN in 21/38 and
exactly ±1.0 in 17/38. The degeneracy guard added to `zstat` on 2026-08-27 also
treats a two-valued null as degenerate, so on any subsequent run both cases
return NaN and the ±1.0 values cannot recur.]`

Third, the immune-exclusion proxy most often used in this setting — mean graph
distance from T-cell-high spots to the tumour core — is not stable. Across the 19
sections its permutation z ranges from **−2.4 to +5.9**, is negative in 7 sections
and within ±2 in 12, and the instability is *within* cohorts and within platform
generations, not only between them (negative in 3/8 Visium and 4/11
first-generation ST sections). A summary statistic that changes sign across
replicate sections of the same tumour type is a weak basis for a barrier claim,
which is one motivation for defining the barrier as a transport quantity rather
than a distance summary.

Local spatial statistics do respond to rearrangement, and more strongly than
`B_cell` does: CAF neighbourhood enrichment gives z = 14.2–37.6 on Visium and
z = 1.5–13.8 on the first-generation ST sections. This is definitional — those
statistics measure clustering, and permutation removes clustering. The
contribution of the min cut is not sensitivity to permutation but that it returns
a *semantically interpretable* object: a blockade with a capacity, a source→sink
orientation, and a geometry that can be overlaid on H&E and counterfactually
opened (S2). We state this explicitly rather than presenting a z-score race we
would lose.

**One comparator had to be dropped on the first-generation ST sections.**
Ripley's L peak for CAF-high spots behaves normally on Visium (z = 7.5–64.5,
8/8 sections) but degenerates on the regular 200 µm staggered array: permuting
positions on a regular lattice barely changes the L function, so the null
standard deviation collapses and the z-score either diverges (7/11 sections
returned values of order 10¹⁴–10¹⁵) or is undefined (4/11). We report Ripley's L
for the Visium sections only and have added a degeneracy guard to the statistic
so that a null with effectively zero variance returns NaN rather than a
spuriously enormous z.
`[scripts/run_10_benchmark_ext.py::zstat — relative-variance floor added
2026-08-27; the affected values are recomputed as NaN on the next run]`

---

## R6. Relation to existing spatial-domain and neighbourhood tools

A methodological reviewer will ask what the min cut adds over a standard spatial
domain segmentation. We answer it directly and in both directions.

The question is not idle — cuts have already entered spatial pipelines, but only
as clustering losses. The tissue-cellular-neighbourhood method of Doiron et al.
(2023, *Nat. Methods*) minimises a min-cut objective to assign cells to
neighbourhoods, and Spatial-RGCN (2025) uses a cut objective for domain
identification; in both, the cut optimises an assignment and is discarded once
the labels are produced. Here the cut **is** the quantity of interest — the
cheapest set of edges whose removal separates an immune source from a tumour
sink — and it is paired with a diffusion–absorption operator for the
macromolecular barrier. The output is a transport barrier with a source and a
sink, not a segmentation of the section into regions; the comparison below is
therefore against the segmentation family, which is the one a reviewer would
reach for.

We segmented each section into eight spatial domains using a BANKSY-style
augmented feature representation (each spot's principal components concatenated
with its neighbourhood mean, λ = 0.3, then k-means), and asked how the resulting
domain boundaries relate to the min cut.

| | all 19 sections | Visium (8) | first-generation ST (11) |
|---|---|---|---|
| fraction of edges that are domain boundaries | 0.18 – 0.49 (median 0.36) | 0.18 – 0.37 (median 0.31) | 0.19 – 0.49 (median 0.42) |
| enrichment of min-cut edges on domain boundaries | 0.69 – 1.87× (median 1.15×) | **1.10 – 1.87× (median 1.49×)** | **0.69 – 1.17× (median 1.07×)** |
| recall — cut edges lying on a domain boundary | 22 – 55 % (median 39 %) | 29 – 55 % (median 39 %) | 22 – 52 % (median 40 %) |
| precision — domain-boundary edges lying on the cut | 4.0 – 22.3 % (median 14.3 %) | 4.0 – 17.6 % (median 15.1 %) | 9.9 – 22.3 % (median 14.2 %) |

`[results/validation/benchmark_tools.json, run_13_benchmark_tools.py,
n_domains = 8, λ = 0.3, 20 PCs]`

**The two objects overlap but are not the same, and how much they overlap depends
on the resolution of the assay.** On Visium, domain boundaries are enriched for
cut edges by about 1.5-fold, so a standard segmentation does see part of the
structure and we do not claim to discover it. On the first-generation ST sections
the enrichment collapses to a median of 1.07×, and **4 of 11 sections fall below
1.0×** — the min cut and the domain boundaries are essentially unrelated there.

Part of the mechanism is combinatorial, and we state it rather than
leave the discrepancy hanging: at 200 µm pitch the analysed graphs hold 370–672
nodes, so partitioning them into the same eight k-means domains labels a much
larger share of all edges as boundary (median 42 % versus 31 % on Visium). When
two edges in five are already boundaries, the ceiling on enrichment is about
2.4× and any real signal regresses towards 1. The k = 8 comparison is therefore
resolution-dependent. **It is, however, not the whole story** — the
`n_domains` scan at the end of this section shows that the platform gap
persists even when the boundary fraction is matched to Visium (point 3 below),
so over-partitioning explains the regression towards 1, not the divergence
itself.

What survives both platforms is the number that matters. **Precision is
4–22 % everywhere** (median 14 % on Visium, 14 % on first-generation ST): a
domain segmentation draws a boundary wherever expression changes, and the great
majority of those boundaries carry no transport consequence, on either platform.

What the min cut adds is therefore selection and quantification rather than
detection. Among many places where the tissue changes, it identifies the one
cross-section that actually limits flux from a specified source to a specified
sink, and returns its capacity. Four properties follow that a domain boundary
does not have: a **capacity** (how much flux passes, which is what makes the
molecular-size comparison in R4/S3 possible), an **orientation** (the boundary is
defined relative to a source–sink pair and changes when they change), a
**counterfactual handle** (the blockade can be opened in silico and the barrier
recomputed, as in S2), and **modality transfer** (the same graph under different
edge-weight semantics becomes the antibody problem, which is what R3 and R3b
rest on).

**A neighbourhood-enrichment comparator agrees that the domains are real but
weaker on the coarser platform.** Squidpy's `nhood_enrichment`, run on the *same*
adjacency the min cut uses rather than on a separately rebuilt neighbour graph,
gives maximum off-diagonal z = 7.2–16.5 on the Visium sections (median 12.3) and
2.8–11.7 on the first-generation ST sections (median 5.6); the most negative
off-diagonal z is −36.3 to −20.6 and −27.3 to −9.0 respectively. The domains are
spatially structured on both platforms, and — consistent with the enrichment
figures above — that structure is weaker at 200 µm pitch with eight domains.
`[results/validation/benchmark_tools.json, squidpy_nhood_z_offdiag_*, 19/19
sections]`

Two limitations of this comparison. The number of domains (eight) is a free
choice, and the results above show it is not innocuous, so we scanned
`n_domains` ∈ {4, 6, 8, 10, 12} on all 19 sections, plus a size-scaled arm
(`clip(round(n_nodes/200), 4, 16)`, i.e. one domain per ~200 spots)
`[results/validation/ndomains_sensitivity.json; results/figures/supp_ndomains_scan.png;
run_13b_ndomains_scan.py]`. Three things hold across the entire grid:

1. **Precision — the number this section argues from — stays between 3 and 23 %
   everywhere** (per-k medians 13–15 % on both platforms), so the reading
   "a domain segmentation proposes, the min cut disposes" does not depend on
   the domain count.
2. **Each platform's enrichment is stable in k.** Visium sits at a median of
   1.30–1.48× at every k; first-generation ST sits at 1.02–1.10× at every k.
   No choice of k closes or widens the platform gap materially.
3. **The size-scaled arm rules out over-partitioning as the whole story.**
   Giving the small legacy sections a proportionally coarser segmentation
   brings their boundary-edge fraction down from a median 42 % to 28 % —
   matching Visium's 31 % — yet their enrichment stays at 1.10×.

We retain the combinatorial caveat for the k = 8 comparison in the main table
(42 % boundary edges cap enrichment near 2.4×), but the honest summary of the
scan is stronger than that caveat: **the divergence between the min cut and
domain boundaries on the coarser platform is robust to the segmentation
resolution, including at boundary fractions matched to Visium.** Whether that
divergence reflects the platform's lower molecular resolution (fewer spots,
shallower sequencing) or the tissue cannot be decided from this design.

Second, the comparison is against the domain-segmentation family; we have not
compared against trajectory- or flow-based spatial methods.

---

## Dimension ③ — bulk ICB cohorts (Introduction / Discussion material, not a main result)

No public dataset combines Visium-format spatial transcriptomics with ICB
response labels: the candidate BCC anti-PD-1 Xenium study lacks the ECM and
hypoxia programmes in its 480-gene panel, the melanoma Visium HD ICB study has
nine samples from four patients, and neither of our cohorts carries response
annotation.

As a concept check only, we scored pre-treatment bulk RNA-seq from two anti-PD-1
melanoma cohorts with the same signature gene sets. In GSE78220
(n = 28, 15 responders / 13 non-responders) the composite barrier score
separates non-responders from responders with AUC = 0.759 (p = 0.021), while a
composite immune score does not (AUC = 0.472, p = 0.818). In GSE91061
(n = 49, 10 / 39) the barrier score is not significant and the direction is
reversed (AUC = 0.572, p = 0.495).
`[results/validation/icb_crosscohort_summary.json; the stored value
AUC_R_vs_NR = 0.241 for GSE78220 is the complementary orientation — the reported
0.759 is AUC(NR vs R) and the direction "high barrier → non-response" was fixed
a priori, which must be stated in the text]`

We present this in the Introduction as motivation and in the Discussion as a
limitation, not among the Results. Its instability across cohorts is consistent
with our own conclusion: barrier molecules measured as bulk abundance lose the
spatial arrangement that makes them a barrier.

## Bulk ICB cohorts (Introduction / Discussion material, not a main result)

No public dataset combines Visium-format spatial transcriptomics with ICB
response labels: the candidate BCC anti-PD-1 Xenium study lacks the ECM and
hypoxia programmes in its 480-gene panel, the melanoma Visium HD ICB study has
nine samples from four patients, and neither of our cohorts carries response
annotation.

As a concept check only, we scored pre-treatment bulk RNA-seq from two anti-PD-1
melanoma cohorts with the same signature gene sets. In GSE78220
(n = 28, 15 responders / 13 non-responders) the composite barrier score
separates non-responders from responders with AUC = 0.759 (p = 0.021), while a
composite immune score does not (AUC = 0.472, p = 0.818). In GSE91061
(n = 49, 10 / 39) the barrier score is not significant and the direction is
reversed (AUC = 0.572, p = 0.495).
`[results/validation/icb_crosscohort_summary.json; the stored value
AUC_R_vs_NR = 0.241 for GSE78220 is the complementary orientation — the reported
0.759 is AUC(NR vs R) and the direction "high barrier → non-response" was fixed
a priori, which must be stated in the text]`

We present this in the Introduction as motivation and in the Discussion as a
limitation, not among the Results. Its instability across cohorts is consistent
with our own conclusion: barrier molecules measured as bulk abundance lose the
spatial arrangement that makes them a barrier.

---

# Discussion

## What the coupling means biologically

The finding that the cellular and macromolecular barriers travel together has a
straightforward tissue correlate. A desmoplastic band is not merely dense; it is
crosslinked, and crosslinking simultaneously stiffens the network against cell
migration and contracts the mesh below the size of an IgG. The two obstructions
are mechanistically distinct — one is a matter of a 10 µm cell negotiating a
fibre network, the other of a 5.5 nm molecule negotiating a mesh — but they are
produced by the same remodelling programme, so in tissue they co-localise.

The therapeutic reading follows directly and is the most useful output of this
work. Strategies that normalise the matrix — LOX/LOXL inhibition, TGF-β
blockade, hyaluronidase, or physical modulation of interstitial pressure — are
predicted to improve both T-cell access and antibody exposure in the same
regions. Conversely, interventions aimed at only one modality (for example, a
chemokine-based strategy to increase T-cell recruitment) would not be expected
to improve antibody delivery, and combination designs that assume otherwise may
be over-optimistic about the achievable intratumoural drug concentration.

## Is the coupling a property of the tissue or of the model?

An obvious concern is circularity: the core matrisome score enters both the
min-cut edge capacity and the diffusion conductance, so some coupling is
guaranteed by construction. We tested this by removing the shared input entirely,
recomputing the cellular barrier from the fibroblast signature alone and the
antibody barrier from crosslinking and antigen alone, so that the two operators
share no input variable (Results R3b).

**In the squamous cohort the coupling is tissue-borne.** Across 15 sections from
six patients the partial correlation survives the removal with a median 81 % of
its magnitude and remains significant in 12/15 sections and in five of the six
patients (median ρ_partial after removal = +0.183). It survives on both platform
generations — four Visium sections from two patients (+0.214, 4/4) and eleven
first-generation ST sections from four patients (+0.137, 8/11) — and it survives
whether the geometric control is a weighted graph distance or the older hop
approximation (5/6 patients either way). Neither sequencing depth nor section
size explains it: across all 19 sections the partial correlation is essentially
uncorrelated with both (ρ = +0.19 and +0.03).

We read this as an organised obstruction. A primary squamous carcinoma builds a
desmoplastic front at the tumour margin in which the fibroblast band, the
collagen it deposits and the crosslinking that contracts the mesh are one
physical structure; there, the cellular and macromolecular obstructions coincide
because they *are* the same wall, and removing the shared molecular term does not
dissolve the association. This interpretation is consistent with our cut-band
analysis, in which every squamous section with a continuous CAF band (7
sections, 1.39–1.74× enrichment, both platforms) is also one whose barrier is
sensitive to arrangement, but it is not established by it.

**In the melanoma cohort the question cannot be answered by this design, and we
say so rather than presenting a contrast.** All four melanoma sections are
metastatic deposits from one patient, so the cohort contributes n = 1 at the
level at which we report every other count; and it is the only stratum whose
classification flips with the geometric control (3/4 sections significant after
removal with the weighted distance, 1/4 with the hop approximation). An earlier,
eight-section version of this analysis reported a squamous-versus-melanoma split
on the strength of "4/4 versus 1/4". The squamous side of that contrast does not
survive the extension to 15 sections, so we withdraw the contrast; what remains
is a positive squamous result and an undetermined melanoma one.

What the extended cohort does settle is one of the three covariates that used to
move together. Cohort previously confounded tumour type, disease stage and
platform generation completely; the squamous cohort is now observed on two
platform generations and gives the same answer on both, so platform generation is
not the driver. Tumour type and disease stage remain confounded with each other,
and separating them still requires primary and metastatic sections of the same
tumour type — a concrete and modest follow-up.

What does not depend on any of this is the negative half of the result. The
dissociation-zone fraction stays below the chance expectation in 18/19 sections
before the shared term is removed and in **19/19** after it (median 0.595× and
0.589× of chance; 0.626× and 0.647× under the hop approximation). **That the two
barriers do not dissociate is robust to the edge-weight definition, to the
geometric control, to the platform generation, to sequencing depth and to section
size.**

## Relation to existing tools

SPARTA is not a spatial-domain method and does not compete with one. A
BANKSY-style segmentation recovers part of the same structure — min-cut edges are
enriched on domain boundaries by 1.10–1.87× (median 1.48×) — and we say so rather
than claiming novelty of the band itself.

The informative number is the other direction. Only 4–18 % of domain-boundary
edges lie on the min cut. A segmentation draws a boundary wherever expression
changes; most of those boundaries have no transport consequence. What the
transport formulation contributes is not detection but selection and
quantification: among many places where the tissue changes, it identifies the one
cross-section that limits flux between a specified source and sink, and returns
its capacity.

Four properties follow that a domain boundary does not have: a capacity (which is
what makes the molecular-size comparison possible at all), an orientation
(defined relative to a source–sink pair, and changing when they change), a
counterfactual handle (the blockade can be opened in silico and the consequence
recomputed), and modality transfer (the same graph under different edge-weight
semantics yields the antibody problem). The last is what makes the central
question of this paper askable.

We also note explicitly what we do not claim. Local spatial statistics respond
more strongly to spatial permutation than our min-cut does, which is definitional
— those statistics measure clustering, and permutation destroys clustering. We
report this rather than presenting a comparison we would lose.

## Limitations

**Sample size, and its asymmetry.** This remains the study's principal
limitation and we state it without softening. The 19 sections come from seven
patients: six squamous patients (two contributing two Visium replicates each,
four contributing three first-generation ST replicates each) and **one** melanoma
patient contributing four separate extracranial metastatic deposits. Sections
from one patient are not independent, so we report section-level and
patient-level counts side by side throughout and draw no conclusion that depends
on treating n = 19 as a sample size.

The asymmetry matters more than the total. The squamous side is now six patients
on two platform generations; the melanoma side is one patient, and no addition to
it is available in the source series that would not confound patient with
anatomical compartment (the only other patient in GSE250636 contributes
leptomeningeal deposits exclusively). Every cross-cohort statement in this paper
is therefore a statement about six patients versus one, and we make none that
depends on the comparison. The conclusions that survive at patient level are the
absence of dissociation (7/7 patients with a positive median partial
correlation), the contiguity effect (6/7 patients significant at the
pre-specified removal fraction) and the tissue-borne character of the coupling in
the squamous cohort (5/6 patients). The spatial-rearrangement effect does not:
it is significant in 4/7 patients and tracks sequencing depth (Results R4/S1),
and we present it as a positive control rather than as evidence.

**Sequencing depth spans two platform generations.** Median UMI per spot ranges
from 567 to 16 686 across the admitted sections — a 29-fold span. We used this
deliberately, to test whether conclusions survive a platform generation, and the
central ones do. But depth is not evenly distributed across patients (P9's two
sections are the shallowest and P4's two the deepest), so depth and patient are
partly confounded, and any result that correlates with depth — S1 does, at
ρ = +0.72 — cannot be cleanly attributed to biology in this design. Quality
control also fragments the shallow sections' graphs (largest connected component
90.0–100 % of nodes), which we report per section in Table 1.

**External validation.** No publicly available dataset combines Visium-format
spatial data with ICB response labels; our bulk concept check is therefore
motivational and not validation. Its own instability across two anti-PD-1 cohorts
is informative in the same direction as our main result: barrier molecules
measured as bulk abundance lose the spatial arrangement that makes them a
barrier.

**Treatment status.** The melanoma cohort has no treatment annotation in GEO,
so we cannot exclude sections obtained after therapy, in which barrier structure
may already have been altered. The squamous cohort is treatment-naive, and the
concordance of results across the two cohorts partially mitigates this concern.

**Signature-based scoring.** All inputs are signature scores rather than
deconvolved cell-type proportions. This is adequate for the relative quantities
the operators consume, but it limits resolution — most acutely for the target
antigen, where the ICB target set contains only two genes (`CD274`, `PDCD1LG2`).
The binding-site barrier is weak under every parameter setting we examined, and
we caution that this may reflect the detection floor of a two-gene score rather
than the biology. An ADC-style target with higher and more structured expression
would be the appropriate test, and the operator requires no modification to run
it.

**Scale parameters.** The size-exclusion law is expressed in nanometres but
takes a within-section rank-normalised crosslinking score as input. This fixes
the nominal mesh-size range and the exclusion threshold identically in every
section, so ξ₀ and β govern how sharply the model separates permeable from
impermeable tissue rather than measuring mesh size. Every conclusion that
depends on their values is reported with that dependence, and the full
sensitivity grid is provided.

**The per-spot cellular field.** Section-level topological claims use the
minimum cut. The per-spot cellular field used for the correlation analyses is a
shortest-path accumulated cost, which permits cost-free detours around a
blockade and therefore reports something closer to effective depth than to
blockade. The two quantities are reported for their respective purposes and are
not interchangeable.

**Blockade width.** The min-cut node set spans 6.7–59.1 % of spots depending on
the section (19 sections; the first-generation ST sections run wider, up to
59 %, than the Visium sections, which top out at 37 %), so it is better
described as a blockade *band* than a line; we
report its width distribution rather than presenting it as a thin contour.

## Outlook

Three extensions follow without modification of the operators. First, replacing
the target-antigen set converts the antibody problem into an antibody–drug
conjugate problem (`ERBB2`, `TACSTD2`, `NECTIN4`), where higher target
expression should make the binding-site barrier measurable and thereby provide
the missing test of the antigen channel. Second, cell-type deconvolution would
replace signature scores with proportions and permit a heterogeneous graph.
Third, the same construction applies to any platform that provides coordinates
and expression, including imaging-based assays at single-cell resolution, where
the 10 µm cell and the 5.5 nm molecule are separated by more than two orders of
magnitude in scale and the two barriers might yet prove separable at a
resolution Visium cannot reach.

That last point is the natural successor to this study. We show that the two
barriers are coupled at 55 µm resolution in cutaneous tumours; whether they
remain coupled when the tissue is resolved cell by cell is an open and directly
testable question.

# Methods

### M1. Cohort-level admission thresholds — and what was actually admitted

Sections were screened against seven admission criteria (C1 count matrix,
C2 spatial coordinates, C3 paired H&E, C4 spot number, C5 median UMI,
C6 detectable endothelial signal, C7 known treatment status and site) before any
analysis. The C4 and C5 thresholds are platform-dependent and were relaxed at
the **cohort** level, separately for each of the two GSE144239 platform
generations:

| Cohort key | C4 (spots) | C5 (median UMI) | Justification recorded |
|---|---|---|---|
| `cscc_gse144239` (Visium) | ≥ 500 instead of ≥ 1 000 | ≥ 500 instead of ≥ 1 500 | 2020-generation Visium; section area and depth not comparable with 2024 data |
| `cscc_legacy_gse144239` (1st-gen ST) | ≥ 300 instead of ≥ 1 000 | ≥ 300 instead of ≥ 1 500 | 2016-generation ST chemistry; 1 933-position staggered array, on-tissue coverage 461–1 181 spots |

Both relaxations, with their dates and justifications, are recorded in
`configs/default.yaml` under `admission_overrides`; **no section was admitted by
a command-line `--force` override.** The per-section outcome — which criteria a
section failed and under which thresholds it was nevertheless admitted — is
written to `data/interim/{slide}.admission.json`.

**One section was rejected.** CSCC13 failed C5 (median 289.5 UMI against the
relaxed threshold of 300) and is recorded as `status = rejected` in the ledger;
it appears in Table 1 with its failure reason and is excluded from every
downstream analysis. The batch driver skips non-`ingested` ledger rows by
default, so the exclusion cannot silently reverse itself.

Section-level characteristics, including the per-section admission outcome, are
given in **Table 1** (`docs/table1_sections.md`, regenerated from the artefacts
by `scripts/run_18_table1.py`).

### M2. Parameter identity: what is anchored and what is not

Parameters fall into three classes, and the class determines how the parameter
may be set.

**Physically anchored** parameters take literature values and never participate
in any fit: the IgG hydrodynamic radius `r = 5.5 nm` and the oxygen diffusion
limit `d₀ = 130 µm`.

**Qualitative scale parameters** set the steepness of a monotone relationship and
cannot be calibrated from these data. `ξ₀ = 20 nm` and `β = 3` belong here. We
had previously described them as literature-anchored; that description is
withdrawn, for a reason that must be stated explicitly. The size-exclusion law
`ξ = ξ₀·exp(−β·x)` has nanometre units, but its input `x` is a
**within-section rank-normalised** crosslinking score with no units. Composing
the two fixes the mesh-size range to `[ξ₀·e^(−β), ξ₀] = [1.0, 20] nm` in every
section regardless of that section's actual crosslinking level, and places the
complete-exclusion threshold (`r ≥ ξ`) at the fixed rank
`ln(ξ₀/r)/β = 0.430` — so a comparable fraction of edges is excluded in every
section by construction, not by biology. `ξ₀` and `β` therefore govern how
sharply the model separates permeable from impermeable tissue; they are not
measurements of mesh size. The full sensitivity grid is reported in
Supplementary Figure S_, and every claim that depends on their values is stated
with that dependence (R2).
The predicted constancy has now been verified on all 19 sections: the median
effective mesh size is **4.36–4.68 nm** and the size-excluded edge fraction is
**60.0–65.2 % (median 62.7 %)** — a range of five percentage points across two
tumour types, two platform generations and a 29-fold span of sequencing depth,
which is what "fixed by construction" looks like in the data.
`[results/validation/review_diagnostics.json, run_11_review_diagnostics.py]`

**Sensitivity parameters** must be scanned on a grid and reported: the four
source/sink quantile thresholds, the graph radius, `k_d,eff`, and the three
`B_meta` weights.

No parameter of any class was set using clinical response labels; the two bulk
ICB cohorts were scored only after all spatial parameters were fixed.

### M3. Statistics

Permutation tests use 500 permutations per section and mode, giving an empirical
p-value resolution of 1/501 = 0.002. **All 19 sections use the same number**, so
p-values are comparable across cohorts and platform generations; this was
verified directly from the stored `n_perm` field rather than assumed.
`[results/counterfactual/*.json, field s1.{fixed,follow}.n_perm = 500 in 19/19]`

Multiple testing is controlled by Benjamini–Hochberg **within each test family,
where a family is one test per section**. The `fixed` and `follow` permutation
modes form separate families because they test different null hypotheses; for S2,
the four removal sizes on one section are nested subsets of the same cut set and
are represented by a single test per section rather than four independent ones.

Effect sizes are reported as ratios (observed / null mean for S1; residual
barrier of scattered removal / residual barrier of contiguous removal for S2)
alongside z-scores, because the permutation nulls have small variance and
z-scores alone overstate the practical size of the effects.

Fractions of "discordant" spots are always reported against their chance
expectation `(1 − q)²`, which is 6.25 % at the quartile threshold used
throughout.

### M4. Distance to vessel, and why it is not a hop count

Every partial correlation in R3 and R3b adjusts for `d_vessel`, the distance from
a spot to the nearest vessel node. That distance is computed as the **weighted
shortest path along the graph**, with edge weights equal to the Euclidean
micrometre distance between the spots they join.

An earlier implementation used the hop count times the nominal array pitch. That
approximation is exact only when all edges have equal length, and its bias
**reverses sign between the two platform geometries**: measured against
straight-line distance, `hops × pitch` over-estimates by a median factor of 1.109
on the Visium hexagonal lattice and under-estimates by 0.894 on the
first-generation ST staggered lattice at `radius_um` = 300 — a systematic ~22 %
offset between the two cohorts in the very covariate the cross-cohort comparison
adjusts for. Since `d₀ = 130 µm` is anchored to a straight-line diffusion limit,
we removed the approximation rather than calibrating it.

The change is recorded per section (`d_vessel_mode` in
`data/interim/*.barrier_meta.json`, `"weighted"` in 19/19) and its effect on the
conclusions is reported as a sensitivity analysis in R3b: it changes no
stratified conclusion for the squamous cohort and changes the melanoma
classification, which is one of the reasons the melanoma result is reported as
unresolved.
`[sparta/barrier.py::compute_b_meta, argument D; run_14_shared_ecm_check.py
--d-vessel {weighted,hops}]`

---

# Table 1

**Table 1a. Cohort composition and admission.**

| Section | Patient | Site / replicate | Cohort · platform | GEO | Spots (raw) | Median UMI | Admission |
|---|---|---|---|---|---|---|---|
| MEL01 | MEL_PtB | sternum | Melanoma · Visium | GSM7983359 | 1,667 | 5,514 | admitted |
| MEL02 | MEL_PtB | cecal-nodule | Melanoma · Visium | GSM7983364 | 840 | 3,098 | admitted |
| MEL03 | MEL_PtB | chest-wall | Melanoma · Visium | GSM7983365 | 1,789 | 3,577 | admitted |
| MEL04 | MEL_PtB | ribcage | Melanoma · Visium | GSM7983366 | 3,263 | 1,488 | admitted |
| CSCC01 | CSCC_P4 | rep1 | cSCC · Visium | GSM4565823 | 744 | 15,714 | admitted |
| CSCC02 | CSCC_P4 | rep2 | cSCC · Visium | GSM4565824 | 696 | 16,686 | admitted |
| CSCC03 | CSCC_P6 | rep1 | cSCC · Visium | GSM4565825 | 3,650 | 982 | admitted |
| CSCC04 | CSCC_P6 | rep2 | cSCC · Visium | GSM4565826 | 3,838 | 636 | admitted |
| CSCC05 | CSCC_P2 | rep1 | cSCC · 1st-gen ST | GSM4284316 | 666 | 5,388 | admitted |
| CSCC06 | CSCC_P2 | rep2 | cSCC · 1st-gen ST | GSM4284317 | 646 | 5,737 | admitted |
| CSCC07 | CSCC_P2 | rep3 | cSCC · 1st-gen ST | GSM4284318 | 638 | 7,219 | admitted |
| CSCC08 | CSCC_P5 | rep1 | cSCC · 1st-gen ST | GSM4284319 | 590 | 1,956 | admitted |
| CSCC09 | CSCC_P5 | rep2 | cSCC · 1st-gen ST | GSM4284320 | 521 | 2,846 | admitted |
| CSCC10 | CSCC_P5 | rep3 | cSCC · 1st-gen ST | GSM4284321 | 521 | 1,648 | admitted |
| CSCC11 | CSCC_P9 | rep1 | cSCC · 1st-gen ST | GSM4284322 | 1,145 | 602 | admitted |
| CSCC12 | CSCC_P9 | rep2 | cSCC · 1st-gen ST | GSM4284323 | 1,071 | 567 | admitted |
| CSCC13 | CSCC_P9 | rep3 | cSCC · 1st-gen ST | GSM4284324 | 1,182 | 290 | **rejected — C5** |
| CSCC14 | CSCC_P10 | rep1 | cSCC · 1st-gen ST | GSM4284325 | 608 | 675 | admitted |
| CSCC15 | CSCC_P10 | rep2 | cSCC · 1st-gen ST | GSM4284326 | 621 | 1,817 | admitted |
| CSCC16 | CSCC_P10 | rep3 | cSCC · 1st-gen ST | GSM4284327 | 462 | 1,442 | admitted |

19 sections from 7 patients admitted; 1 rejected at admission.

**Table 1b. Graph properties after quality control.**

| Section | Nodes after QC | Mean degree | Components | Largest component | Source / Sink / Vessel |
|---|---|---|---|---|---|
| MEL01 | 1,660 | 5.74 | 1 | 100.0 % | 133 / 466 / 332 |
| MEL02 | 801 | 5.35 | 8 | 91.0 % | 65 / 232 / 161 |
| MEL03 | 1,724 | 5.40 | 10 | 96.4 % | 138 / 486 / 345 |
| MEL04 | 2,673 | 5.41 | 15 | 98.9 % | 214 / 762 / 535 |
| CSCC01 | 701 | 5.39 | 3 | 99.4 % | 57 / 201 / 141 |
| CSCC02 | 652 | 5.44 | 5 | 98.8 % | 53 / 182 / 131 |
| CSCC03 | 2,645 | 5.29 | 50 | 94.2 % | 212 / 481 / 529 |
| CSCC04 | 2,328 | 4.92 | 54 | 90.0 % | 187 / 358 / 466 |
| CSCC05 | 664 | 7.50 | 1 | 100.0 % | 53 / 197 / 133 |
| CSCC06 | 643 | 7.46 | 1 | 100.0 % | 52 / 190 / 129 |
| CSCC07 | 637 | 7.46 | 1 | 100.0 % | 51 / 187 / 128 |
| CSCC08 | 566 | 7.22 | 2 | 99.6 % | 46 / 156 / 114 |
| CSCC09 | 511 | 7.20 | 1 | 100.0 % | 41 / 142 / 103 |
| CSCC10 | 448 | 6.87 | 1 | 100.0 % | 36 / 128 / 90 |
| CSCC11 | 672 | 6.25 | 15 | 94.5 % | 54 / 185 / 135 |
| CSCC12 | 595 | 6.30 | 16 | 95.1 % | 48 / 161 / 119 |
| CSCC14 | 370 | 5.94 | 3 | 98.6 % | 30 / 106 / 74 |
| CSCC15 | 597 | 7.24 | 1 | 100.0 % | 48 / 166 / 120 |
| CSCC16 | 380 | 6.75 | 2 | 99.7 % | 31 / 103 / 76 |

## Figures

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

![](results/figures/_doc/fig1_framework.png)

**Figure 1. The framework: one graph, two edge-weight semantics, two operators.**
Schematic on a synthetic section (clearly labelled as such; no real data and no
research claim is made from this panel). The same spatial graph is read twice:
T-cell migration as a source–sink minimum cut, which returns a barrier magnitude
and an explicit blockade geometry, and antibody delivery as a screened
Poisson diffusion–absorption field that reduces to the harmonic /
effective-resistance limit when absorption vanishes. Generated by
`scripts/run_15_figure1.py`.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

![](results/figures/_doc/fig2_driver_decomposition.png)

**Figure 2. The antibody barrier decomposes onto one axis, and that axis's
dominance is a modelling choice rather than a measurement.**
(**a**) Share of `B_mAb` variance carried by each of three ablation-defined
channels, per section, on a log axis: crosslinking / size exclusion
95.9–98.9 %, the matrix term shared with the T-cell barrier 0.8–2.8 %, and the
target antigen 0.47–3.06 % over the 15 sections in which the antigen signature
could be scored. In CSCC10, CSCC14, CSCC15 and CSCC16 neither `CD274` nor
`PDCD1LG2` passed the minimum-gene-match threshold; those sections are marked
with a slashed open triangle at the axis edge — the value is **missing, not
near-zero** — and are excluded from the antigen range. Dashed horizontal rule
separates the squamous from the melanoma cohort.
(**b**) Crosslinking share of `B_mAb` variance as a function of the
mesh-contraction steepness β, at λ = 3, for five sections spanning both tumour
types and both platform generations. The five curves are near-identical at every
β, which is the point: the share is a property of the scale parameter, not of the
tissue. Vertical dashed line marks the default β = 3.
`[results/validation/screen_decision.json; results/validation/bmab_sensitivity.json]`

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

![](results/figures/_doc/fig3_coupling.png)

**Figure 3. The two barriers do not dissociate, and in the squamous cohort the
coupling survives removal of every shared input.**
Circles denote primary cSCC sections (6 patients, 15 sections), squares
metastatic melanoma deposits (1 patient, 4 deposits); the dashed vertical rule
separates the cohorts throughout.
(**a**) Partial Spearman correlation between the per-spot `B_cell` field and
`B_mAb`, adjusted for weighted graph distance to the nearest vessel:
positive in 18/19 sections and significant in 17/19 (ρ = −0.000 to +0.385,
median +0.217). Sections not reaching p < 0.05 are marked *n.s.*
(**b**) The same partial correlation before and after removing every input the
two operators share (`b_ecm` = 0 and λ = 0). Solid blue = still significantly
positive; grey = collapses. Squamous 12/15 and melanoma 3/4 survive; median
after removal +0.152. Only the collapsing sections and the extremes of the
retained bundle are labelled; the remaining 13 sections all retain the
association.
(**c**) Dissociation-zone spots — lowest `B_cell` quartile and highest `B_mAb`
quartile — as a percentage of spots (1.1–6.9 %, median 3.7 %) against the 6.25 %
expected under independence at q = 0.75 (orange line). 18/19 sections fall below
the line; the single section above it (CSCC15) is the same section whose partial
correlation is not positive.
`[results/validation/decoupling.json; results/validation/shared_ecm_check.json]`

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

![](results/figures/_doc/fig4_counterfactuals.png)

**Figure 4. Three counterfactuals: arrangement, continuity and molecular size.**
(**a**) S1 — spatial rearrangement with composition held fixed. Observed
`B_cell` divided by the permutation-null mean (500 permutations per section);
blue = significant after Benjamini–Hochberg correction across sections (12/19),
grey = not. Ratios span 0.89–4.49. Direct labels are shown for the significant
sections only.
(**b**) S2 — blockade continuity. Residual barrier after removing a contiguous
arc of the min cut, divided by the residual barrier after removing the same
amount of material at scattered positions within the same cut set, at four
removal fractions. The pre-specified primary test is 20 % (orange dashed line):
16/19 sections significant after BH, ratios 1.038–1.390 (median 1.147); the
effect is larger at 30 % than at 5 % in 15/19 sections.
(**c**) S3 — molecular size scan on fixed tissue. Mean `B_mAb` in the tumour
core rises monotonically with hydrodynamic radius in every section, with no
ceiling in the range examined; the dashed line marks IgG at 5.5 nm.
`[results/counterfactual/*.json; results/validation/paper_stats.json]`

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

![](results/figures/_doc/fig5_vs_domains.png)

**Figure 5. Spatial-domain methods propose candidate boundaries; the min cut
selects the one that limits flux and gives it a capacity.**
(**a**) Fold enrichment of min-cut edges on the boundaries of an
eight-domain BANKSY-style segmentation (λ = 0.3, 20 PCs), per section. Orange
dashed segments give the median within each cohort × platform block; enrichment
is 1.48× on Visium and 1.07× on first-generation ST, and four of the eleven
first-generation ST sections fall below 1.0×. That gap is combinatorial rather
than biological: at 200 µm pitch the analysed graphs hold 370–672 nodes, so the
same eight domains label a median 42 % of edges as boundary versus 31 % on
Visium, which caps attainable enrichment near 2.4×. Only the extremes and the
sub-1.0× sections are labelled.
(**b**) Recall (share of min-cut edges lying on a domain boundary) and precision
(share of domain-boundary edges lying on the min cut). Precision is 4–22 % on
both platforms: a segmentation draws a boundary wherever expression changes, and
the great majority of those boundaries carry no transport consequence.
Squidpy neighbourhood enrichment, computed on the same adjacency, gives maximum
off-diagonal z = 7.2–16.5 on Visium and 2.8–11.7 on first-generation ST
(19/19 sections).
`[results/validation/benchmark_tools.json]`

**Table 1. Sections analysed.**
Per-section provenance, depth, graph size and connectivity, source / sink /
vessel node counts, and the admission outcome under the pre-specified criteria
C1–C7. CSCC13 failed C5 (median 289.5 UMI against a cohort threshold of 300) and
is retained in the table with its failure reason but excluded from every
downstream analysis. Regenerated from the artefacts by
`scripts/run_18_table1.py`.

