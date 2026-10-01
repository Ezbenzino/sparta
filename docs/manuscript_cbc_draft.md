# CBC 投稿初稿（英文版）

> **文档性质**：冲击 *Computational Biology and Chemistry*（Elsevier，订阅路线
> 免发表费）的完整初稿，2026-08-28 成文。Elsevier "Your Paper Your Way"——
> 首投接受 Markdown/任意清晰排版，此处按 CBC 终稿结构组织：Highlights →
> Graphical Abstract（必选项，见编辑注）→ Abstract → Keywords → 1–5 节 →
> CRediT/利益声明/数据可用性 → References（Elsevier 数字引用）→ 图表清单。
>
> **投稿前编辑注（不进入正文）**
> 1. **Graphical abstract 是 CBC 必选项**。建议由 Fig 1（框架图，
>    `results/figures/fig1_framework.png`）裁剪重制：单图呈现"一张空间图、
>    两套边权语义、两个算子、一个共享基质"即可；规格 Elsevier 标准
>    （≥531×1328 px，JPG/PNG，≤10 MB）。
> 2. 文中所有统计数字均取自仓库产物（见各节 `[...]` 留痕标记，投稿时删除）。
> 4. 作者信息、基金号、仓库 URL（GitHub + Zenodo DOI）待 git 首次提交后回填。
> 5. 行号：Elsevier 要求投稿 PDF 带连续行号，LaTeX/Word 排版时打开。
> 6. 中文对照版：`docs/manuscript_cbc_draft_zh.md`。

---

## Title

**Two transport operators, one substrate: a graph model of T-cell migration and
antibody penetration barriers in tumour tissue from spatial transcriptomics**

备选（更收敛、偏方法）：

**A graph-transport model of cell and antibody delivery barriers in tumour
tissue, parameterised from spatial transcriptomics**

---

## Highlights

- Graph operators turn spatial transcriptomics into transport barrier quantities.
- Minimum cut yields the T-cell barrier; a screened Poisson field yields the IgG barrier.
- The antibody barrier is percolation-limited: 60–65% of matrix edges exclude IgG.
- Both barriers share a matrix origin and co-move in 18 of 19 tumour sections.
- Deterministic CPU implementation analyses one section in under 0.2 s.

*(每条 ≤85 字符，符合 Elsevier Highlights 规范)*

---

## Abstract

Monoclonal antibodies and cytotoxic T cells must both cross the tumour
extracellular matrix to reach their targets, yet they differ by three orders of
magnitude in size, and the analytical tools of spatial transcriptomics return
labels rather than transport quantities. We introduce two deterministic graph
operators that share one spatial graph and differ only in edge-weight
semantics. T-cell migration is modelled as a source–sink minimum cut between
immune-entry and tumour-core compartments, returning a barrier strength
together with its blockade geometry. Antibody transport is modelled as a
screened Poisson diffusion–absorption field whose conductances implement size
exclusion against the IgG hydrodynamic radius (5.5 nm) and whose sinks
represent target-antigen binding. Applied to 19 sections from 7 patients
spanning two cutaneous tumour types and two spatial-platform generations, the
model shows that the antibody barrier is percolation-limited — the effective
mesh size falls below 5.5 nm on 60–65% of edges — and that size exclusion, not
antigen availability or bulk matrix density, accounts for a median 97.5% of
barrier variance (at the default scale parameter β = 3; full dependence on β
in Section 3.2). Ablating the matrix channel shared by both operators leaves
their spatial coupling positive and significant in 12 of 15 squamous sections
(median retention 81%), indicating that the two delivery problems are
mathematically distinct but physically co-localised in the tissue. A molecular
size scan on the same tissue confirms the distinction: the identical graph is
nearly transparent to a 0.5 nm solute and strongly obstructive to an IgG. The
full analysis of one section runs in under 0.2 s on CPU with no learned
parameters, making the model usable as a reproducible component of
translational pipelines.

*(241 words；无引用；缩写 IgG 于首现处即 IgG hydrodynamic radius 语境内；满足
"独立成文"要求)*

**Keywords**: spatial transcriptomics; graph algorithms; minimum cut; screened
Poisson equation; antibody transport; tumour microenvironment

---

## 1. Introduction

Solid tumours are difficult places to deliver things to. A therapeutic IgG is a
150 kDa globular protein with a hydrodynamic radius near 5.5 nm; a cytotoxic T
cell is an object of roughly 10 µm. Both must arrive from the vasculature,
through interstitial matrix whose collagen density and crosslinking vary
strongly over distances comparable to a few spot diameters of a spatial assay,
before either antigen binding or tumour-cell killing can occur. The physical
side of this problem has a long modelling tradition — compartmental
pharmacokinetics, Krogh-cylinder penetration models, and reaction–diffusion
descriptions of macromolecule transport through matrix — but those models are
formulated at the tissue or organ level with effective parameters, and they
are not driven by molecular measurements of a specific patient section. The
molecular side now produces such measurements: spatial transcriptomics
resolves, at 10–100 µm pitch, the local abundance of collagen and
crosslinking-associated transcripts, of stromal and immune signatures, and of
the therapeutic target itself.

What spatial transcriptomics has lacked is an operator vocabulary that turns
those measurements into transport quantities. The prevailing analyses return
labels: spatial-domain segmentations and their embeddings, neighbourhood
enrichment scores, cell-type maps. These answer *what is next to what*; they do
not return a conductance, a capacity, or a barrier with a source and a sink.
Even where graph-cut machinery has entered the field, it has entered as a
clustering loss: the tissue-cellular-neighbourhood method of Doiron et al.
[7] and the Spatial-RGCN domain identifier [8] optimise min-cut objectives to
produce assignments, and the cut is discarded once labels are obtained. The
transport question — which cross-section of the matrix limits flux from a
specified source to a specified sink, and by how much — is not asked.

We ask it, and we ask it for two delivery modalities at once, because the two
have a clinically important relationship. If the T-cell barrier and the
antibody barrier were spatially dissociable, matrix-directed intervention
could be timed or dosed to relieve one modality selectively; if they are
co-located because both are expressions of the same matrix substrate, then
matrix normalisation improves both classes of delivery simultaneously and
neither is reachable independently. This is a quantitative question about
transport in measured tissue, which is precisely what a graph-transport model
can answer and a label-based analysis cannot.

The model — we refer to the implementation as SPARTA — is built on one spatial
graph per section and two operators. For cell migration, immune-entry nodes
(endothelial-rich) are the source and tumour-core nodes (malignant-rich,
immune-poor) are the sink; edge capacities decrease with the local matrix
score, and the barrier is the reciprocal of the maximum source–sink flow, with
the minimising cut returned as an explicit blockade geometry. For antibody
transport, the same graph carries a screened Poisson diffusion–absorption
operator: conductances implement size exclusion through an effective mesh size
that contracts with the local crosslinking score, and absorption represents
target-antigen binding. The two operators share no functional form, but they
share one input — the matrix — and the paper's central experiment removes that
shared input entirely and asks whether the coupling survives.

Three design decisions follow from the computational setting. Parameters are
partitioned into physically anchored values (the IgG radius, the oxygen
diffusion limit), qualitative scale parameters that are reported with their
full sensitivity dependence, and sensitivity parameters scanned on grids
(Section 2.4); no parameter was fitted to any outcome, clinical or otherwise.
Every inferential statistic is evaluated with the patient, not the section, as
the independent unit — 19 sections come from 7 patients, and counts are
reported at both levels throughout. And every claim the operators generate is
tested by counterfactual computation on the same tissue: spatial rearrangement
of the resistance field, removal of the blockade as a contiguous arc versus
scattered fragments, and a molecular size sweep from small molecule to large
antibody.

Applied to 19 sections (15 primary cutaneous squamous cell carcinomas from 6
patients, 4 metastatic melanoma deposits from 1 patient; two platform
generations of spatial transcriptomics four years apart), the model yields
four findings. The antibody barrier is percolation-limited: at the default
scale parameters the effective mesh size is below the IgG radius on 60–65% of
edges, and size exclusion accounts for a median 97.5% (at β = 3) of
antibody-barrier variance. The two barriers are mathematically distinct — the same tissue graph
is nearly transparent to a 0.5 nm solute and strongly obstructive to an IgG —
but physically coupled: after controlling for distance from vasculature, the
per-spot barrier fields are positively correlated in 18/19 sections, and the
dissociation of the two barriers is *rarer than chance* in 18/19. Removing the
shared matrix input entirely leaves the coupling significant in 12/15 squamous
sections with a median 81% retention, so the coupling is substantially
tissue-borne rather than an artefact of shared edge weights. And the blockade
behaves as a connected structure rather than a pile of resistant spots:
opening a contiguous gap in the cut band leaves ~15% more residual barrier
than removing the same material at scattered positions, with a dose–response
over removal fraction in 15/19 sections.

Our contributions are:

1. A pair of interpretable, deterministic graph-transport operators for
   spatial molecular data — a source–sink minimum-cut barrier for cell
   migration and a screened Poisson diffusion–absorption barrier for
   macromolecule transport — connected by an explicit parameter-identity
   discipline (anchored / scale / sensitivity);
2. A counterfactual validation framework (rearrangement, blockade continuity,
   molecular size sweep, shared-input ablation) applicable to any operator
   claiming to measure transport in tissue;
3. An application across two platform generations quantifying the coupling of
   the two delivery barriers, with a direct implication for the design of
   matrix-directed combination therapies and for the choice of antibody- vs
   cell-based modalities in matrix-rich tumours;
4. A reference implementation that analyses a full section in 0.02–0.19 s on
   CPU with no learned parameters.

---

## 2. Materials and methods

### 2.1. Data and admission

Nineteen sections from two public cohorts were analysed: fifteen
treatment-naive primary cutaneous squamous cell carcinomas (CSCC01–CSCC16
minus CSCC13) from GSE144239 [2] and four extracranial metastatic deposits of
cutaneous melanoma (MEL01–MEL04) from GSE250636 [3]. The four melanoma
sections are separate deposits (sternum, cecal nodule, chest wall, ribcage) of
a single patient; the second patient in GSE250636 contributes only
leptomeningeal deposits, a different anatomical compartment, and was not
included (Section 4.4). Patient assignment follows the GEO sample metadata and
is recorded with replicate structure in the analysis ledger; **all counts
below are reported at both section and patient level**, and n = 19 is not
treated as a sample size.

Sections passed seven pre-declared admission criteria (tumour-content
fraction, minimum spots, minimum median UMI, detected genes, endothelial
signal, treatment status and site) recorded before analysis. Two cohort-level
threshold relaxations for the 2016- and 2020-generation platforms, with dates
and justifications, were fixed before ingestion; no section was admitted by a
run-time override. One section (CSCC13) failed the relaxed median-UMI
threshold (289.5 vs 300) and was excluded; its admission record is retained.
Per-section outcomes are in Table 1.

The squamous cohort deliberately spans two platform generations: four
2020-generation Visium sections (2 patients) and eleven 2016-generation
first-generation ST sections on a staggered 200 µm array (4 patients, three
technical replicates each). Agreement across that platform shift is used as a
robustness argument throughout.

### 2.2. Spatial graphs

Coordinates were rescaled so that the median nearest-neighbour distance equals
the platform's nominal pitch (100 µm Visium, 200 µm first-generation ST). The
rescaling matters: the first-generation array is staggered (array indices
satisfy x + y ≡ 0 mod 2), so nearest neighbours lie on the diagonal and one
index unit corresponds to 141.4 µm, not 200 µm; assuming otherwise inflates
every micrometre quantity by √2. Graphs were built by radius adjacency
(radius 150 µm Visium, 300 µm first-generation ST; the latter connects both
200 µm and 282.8 µm neighbours; mean degree 4.9–7.5 after quality control).
Spots below 500 UMI were removed, which fragments the shallow sections; we
report the fragmentation rather than assume connectivity: analysed graphs have
1–54 components with 90.0–100% of nodes in the largest, and 8/19 sections have
some source or sink nodes outside the largest component, which biases the
cellular barrier slightly upward there. The operator counts and reports the
affected nodes per section.

### 2.3. Signature scores

Inputs to both operators are within-section signature scores (rank-normalised
within each section): endothelial and T/NK signatures defining source and
sink compartments; a malignant signature defining the tumour core; CAF
(cancer-associated fibroblast), core-matrisome and crosslinking-associated
signatures defining resistance and conductance; a two-gene target-antigen
signature for anti-PD-1 agents (CD274, PDCD1LG2); and hypoxia
(MSigDB v7.1 HALLMARK_HYPOXIA [8]), proliferation and efflux signatures for
the metabolic barrier. Scores are not deconvolved proportions; consequences
are discussed in Section 4.4.

### 2.4. The two operators and parameter identity

**Cell migration — source–sink minimum cut.** Endothelial-rich nodes form the
source (immune entry), malignant-rich immune-poor nodes the sink; quantile
thresholds (scan class below) define membership, with a border fallback when a
compartment is empty. Edge capacities are decreasing functions of the local
matrix resistance (shared core-matrisome and CAF scores). The section-level
barrier `B_cell` is the reciprocal of the maximum source–sink flow, computed
exactly (preflow-push) on the largest component; the minimising cut set is
returned as a blockade band with width distribution (6.7–59.1% of spots across
sections). The per-spot field used in correlation analyses is the
shortest-path accumulated migration cost from immune-entry nodes; the two
quantities are reported for their respective purposes and are not
interchangeable (Section 4.4).

**Antibody transport — screened Poisson diffusion–absorption.** Vessel nodes
hold unit concentration; edge conductance decreases with the local matrix
score and implements size exclusion through an effective mesh size
`ξ = ξ₀·exp(−β·x)`, where `x` is the within-section rank-normalised
crosslinking score and `r = 5.5 nm` the IgG hydrodynamic radius: an edge whose
ξ falls below r contributes only a conductance floor. Target-antigen
expression acts as a distributed absorption sink. The steady state solves a
screened Poisson (Helmholtz-type) system on the graph; the barrier `B_mAb` is
the mean concentration deficit over the tumour-nest core, and the per-spot
deficit field is the correlate used in Section 3.4. With absorption removed
the operator reduces to the harmonic (effective-resistance) limit [11].
Variance of `B_mAb` is decomposed by single-channel ablation into a matrix
channel (shared with the cut), a size-exclusion channel (crosslinking) and an
antigen channel.

**Parameter identity.** Parameters fall into three classes. *Physically
anchored*: IgG hydrodynamic radius r = 5.5 nm; oxygen diffusion limit
d₀ = 130 µm. These take literature values and are never fitted. *Qualitative
scale parameters*: ξ₀ = 20 nm and β = 3 set the steepness of the exclusion
law. They are not mesh-size measurements: because x is rank-normalised, the
model fixes the mesh range to [ξ₀e^−β, ξ₀] and the complete-exclusion rank to
ln(ξ₀/r)/β = 0.430 identically in every section — a constancy verified in the
data (median effective mesh size 4.36–4.68 nm; excluded edge fraction
60.0–65.2% across two tumour types, two platform generations and a 29-fold
depth span). Every claim depending on ξ₀ or β is reported with its dependence
over a β grid (Section 3.2). *Sensitivity parameters*: the four compartment
quantile thresholds, the graph radius, the effective binding rate, and the
metabolic weights, each scanned on a grid. No parameter of any class was set
using clinical response labels.

### 2.5. Statistics

Permutation tests use 500 permutations per section and mode (identical across
sections; empirical resolution 1/501). Multiple testing is controlled by
Benjamini–Hochberg within each test family (one test per section; the
rearrangement modes and removal fractions form separate families). Effect
sizes are reported as ratios alongside z-scores, because the permutation nulls
have small variance and z-scores overstate practical magnitude. Discordant
fractions are compared against their chance expectation (1 − q)² = 6.25% at
the quartile threshold. All correlations are Spearman; partial correlations
residualise both fields on distance to the nearest vessel (quadratic, in rank
space), where that distance is the weighted shortest path in micrometres —
not a hop count: the hop approximation is biased in opposite directions on
the two lattices (median ×1.109 Visium, ×0.894 first-generation ST), and its
effect is reported as a sensitivity analysis wherever it changes anything.

---

## 3. Results

### 3.1. Cohort and graphs

Table 1 summarises the cohort: 19 sections, 7 patients, 2 tumour types, 2
platform generations; 370–2 673 nodes per section after quality control;
median UMI per spot spanning 567–16 686 (a 29-fold span used deliberately as a
robustness axis). Graph construction differs between platforms only in
geometry (Section 2.2), so operator outputs are comparable across the
platform shift.

### 3.2. The antibody barrier is percolation-limited

Single-channel ablation decomposes `B_mAb` variance into the matrix channel
shared with the cut, the size-exclusion channel, and the antigen channel.
Across the 19 sections, size exclusion accounts for **95.9–98.9%** (median
97.5%) and the shared matrix channel for 0.8–2.8% (median 1.8%); the two
inputs are not collinear in any section (|r| = 0.010–0.512). On the 15
sections where the two-gene target-antigen signature is detectable at all,
the antigen channel accounts for 0.47–3.06% (median 0.96%); on the four
shallowest first-generation sections neither gene passes detection, which is
missing data and is reported as such rather than as zero.

The mechanism is measured directly. At default parameters the effective mesh
size falls below the IgG radius on the median edge of every section (median
ξ = 4.36–4.68 nm; 60.0–65.2% of edges size-excluded, median 62.7%), and only
0.0–7.8% of spots (median 1.0%) saturate the numerical conductance floor. The
antibody barrier is therefore not a smooth gradient but a largely
percolation-limited field: transport proceeds through the minority of edges
whose local mesh remains open to a 5.5 nm solute.

This dominance is a property of the scale parameter and we say so. Sweeping β
from 0.25 to 12 at λ = 3 on five sections spanning both tumour types and both
platform generations, the size-exclusion share of `B_mAb` variance moves from
0.2–0.3% (β = 0.25) through 50.3–75.0% (β = 1.5) to 95.9–97.7% (β = 3) and
81.0–96.0% (β = 12), with the five sections agreeing within a few points at
every value. The defensible claim is not an empirical one about skin cancer
but a conditional one about the model: given an exclusion law steep enough
that mid-range crosslinking already contracts the mesh below the IgG radius —
which at β = 3 it does on ~63% of edges — the barrier is dominated by that
axis rather than by bulk matrix density or antigen availability. The antigen
channel is weak under every parameter setting examined (0.3–17.9% over the
full grid); with a two-gene signature this is a detection floor, not evidence
that binding sites are unimportant (Section 4.4).

### 3.3. Molecular size separates the two transport problems

Sweeping the hydrodynamic radius from 0.5 nm to 10 nm at fixed tissue
structure raises the mean barrier in the tumour-nest core monotonically in
every section (MEL01: 1.48 → 13.20; CSCC03: 2.86 → 21.56; CSCC04: 2.40 →
19.29 across eight radii, no saturation at the numerical ceiling). The
identical graph is nearly transparent to a small molecule and strongly
obstructive to an IgG-sized molecule. This is the cleanest demonstration that
the two operators pose different physical problems even though they are driven
by the same substrate — a distinction no cell-migration model produces, and
one that is a deterministic consequence of the operator rather than a
statistical test.

### 3.4. The two barriers do not dissociate

The central question is whether the two barriers are spatially separable in
real tissue. They are not. After residualising both per-spot fields on
distance to the nearest vessel, they remain positively correlated in 18/19
sections and significantly so in 17 (median partial Spearman ρ = +0.217; all
seven patients positive at patient level: +0.141 to +0.323). Within the
squamous cohort the two platform generations agree (Visium median +0.270,
first-generation ST +0.207). The dissociation zone — lowest quartile of the
cell barrier, highest quartile of the antibody barrier — occupies a median
3.7% of spots against a 6.25% chance expectation, below chance in 18/19
sections: discordance is rarer than random, a second independent expression
of the coupling.

Part of this association is guaranteed by construction — the core-matrisome
score enters both edge-capacity and conductance — so the decisive experiment
removes the shared input entirely, recomputing the cellular barrier from the
fibroblast signature alone and the antibody barrier from crosslinking and
antigen alone, such that the two operators share no input variable. Across 19
sections the median partial correlation falls from +0.217 to +0.152; 17/19
remain positive, 15/19 significantly. Stratified:

| Stratum | Sections / patients | ρ after removal | Significant | Patients majority-significant |
|---|---|---|---|---|
| cSCC, Visium | 4 / 2 | +0.214 | 4 / 4 | 2 / 2 |
| cSCC, first-generation ST | 11 / 4 | +0.137 | 8 / 11 | 3 / 4 |
| cSCC combined | 15 / 6 | +0.183 | 12 / 15 | 5 / 6 |
| Melanoma (one patient) | 4 / 1 | +0.056 | 3 / 4 | — |

The squamous coupling survives on both platform generations, under both the
weighted and the hop-count geometric control (5/6 patients either way), and is
essentially independent of sequencing depth (ρ = +0.19) and section size
(ρ = +0.03). We read it as organised obstruction: a desmoplastic front in
which the fibroblast band, its deposited collagen, and the crosslinking that
contracts the mesh are one physical structure, so the cellular and
macromolecular obstructions coincide because they are the same wall. The
melanoma stratum cannot be decided by this design — it is one patient, and
the only stratum whose classification flips with the distance convention —
and we report it as unresolved rather than as a contrast. An earlier
eight-section version of this analysis asserted a squamous-versus-melanoma
split on the strength of 4/4 versus 1/4; the squamous side does not survive
extension to 15 sections (12/15), and we withdraw the contrast.

The negative half of the conclusion is robust to everything tested: the
dissociation-zone fraction stays below chance in 18/19 before the shared-input
removal and 19/19 after it, under both distance conventions, across both
platform generations, independent of depth and size. The two barriers do not
dissociate.

### 3.5. Counterfactual experiments

**Blockade continuity.** Removing resistance material from a contiguous arc
of the min-cut band, compared with removing the same material from scattered
positions within the same cut set (same material, same amount, only
arrangement differs), leaves the residual barrier higher by a factor of
1.038–1.390 (median 1.147) at the pre-specified 20% removal fraction;
16/19 sections are significant after per-section correction, five of seven
patients have every section significant, and no patient has none. The effect
shows a dose–response over removal fraction (residual barrier at 30% > 5% in
15/19 sections) — the falsifiable prediction of a connected blockade and not
the behaviour of uncorrelated noise. We regard the dose–response as the
stronger evidence and the single-level significances as secondary. The effect
is real but modest (real bands are thick and redundant; a synthetic
one-spot-wide ring gives 7.45× under the same experiment), and we do not rest
the central argument on it.

**Spatial rearrangement.** Permuting spot positions while keeping composition
fixed moves the section-level cellular barrier above its permutation null by
0.89×–4.49×, significant in 12/19 sections, organised by patient (three
patients fully significant, two fully non-significant). We report this as a
positive control, not as evidence: the effect ratio correlates with sequencing
depth (ρ = +0.72) and graph connectivity (ρ = +0.46), which are themselves
correlated, and the non-significant patients sit at the shallow, fragmented
end of the range. Composition-based measures cannot move under this null at
all — their permutation nulls are degenerate — which is the point of
difference: a barrier is a statement about arrangement, and only an
arrangement-sensitive quantity can be evidence about it.

### 3.6. The barrier is not a rewrite of density, and its relation to domain tools

The per-spot cellular field correlates only moderately with the CAF signature
(ρ = 0.313–0.554, median 0.458) and essentially not at all with the T/NK
signature (median −0.018), and the commonly used immune-exclusion proxy —
mean distance from T-cell-rich spots to the tumour core — changes sign across
replicate sections of the same tumour type (permutation z from −2.4 to +5.9,
negative in 7 sections), which is one motivation for defining the barrier as
a transport quantity rather than a distance summary.

Relation to domain segmentation follows the same logic. A BANKSY-style
segmentation [4] recovers part of the structure (min-cut edges are enriched
on domain boundaries by 1.10–1.87×, median 1.48×, over the domain-count grid
k = 4–12), but only 4–18% of domain-boundary edges lie on the min cut: a
segmentation draws a boundary wherever expression changes, while the
transport formulation selects the one cross-section that limits flux between
a specified source and sink and returns its capacity. Four properties follow
that a domain boundary does not have: a capacity, an orientation, a
counterfactual handle, and modality transfer. Where cuts have appeared in
spatial pipelines before — as clustering losses in tissue-neighbourhood
assignment [7] and domain identification [8] — the cut optimises an
assignment and is discarded; here it is the quantity of interest.
Computational cost is part of the comparison: on the same graphs, the full
five-stage SPARTA analysis runs in 0.02–0.19 s per section (entire cohort
0.9 s, peak memory < 16 MB, commodity CPU), against 0.16–1.20 s for
BANKSY-style segmentation and 3.8–5.0 s for neighbourhood enrichment [5]
(both reference implementations run within this repository for controlled
comparison; graph-learning domain methods were not run and are not claimed
as comparisons).

---

## 4. Discussion

**What the model says about delivering biologics.** The dominant term in the
antibody barrier is not how much matrix there is, but whether the local mesh
is open to the molecule: transport is percolation-limited, proceeding through
the minority of edges whose effective mesh exceeds the hydrodynamic radius.
Two corollaries follow. First, bulk measurements of matrix abundance are a
poor proxy for penetration — two tissues with identical average matrix but
different mesh topology can differ substantially in delivered fraction, which
is consistent with the instability we observe when scoring bulk RNA-seq
cohorts with the same signature logic (Section 4.4). Second, the model
provides a per-section, parameter-swept way to ask what an intervention buys:
the molecular size scan is the degenerate case, and replacing the target set
(`ERBB2`, `TACSTD2`, `NECTIN4`) converts the antibody problem into an
antibody–drug-conjugate problem with a measurable binding-site channel, at no
change to the operator.

**One substrate, two modalities.** After removing the only input the two
operators share, their coupling retains a median 81% of its magnitude in the
squamous cohort and remains significant in 5/6 patients. The practical
reading for therapy design is that in matrix-rich cutaneous tumours the two
delivery problems are not independently addressable: matrix normalisation is
predicted to relieve both simultaneously, and modality selection (antibody
versus cell) cannot route around the matrix. The corollary for measurement is
that any assay meant to predict delivery must resolve arrangement, not only
abundance — and the dissociation zone, the region where one modality is
delivered and the other is not, occupies less than the chance fraction in
every section examined.

**Relation to continuum transport models.** Macromolecule penetration has a
mature continuum tradition — Krogh cylinders, reaction–diffusion and
pharmacokinetic models — whose parameters are effective and whose geometry is
idealised. The graph formulation is complementary: it inherits the geometry of
a measured section, defines barrier quantities on that geometry directly, and
sacrifices physical continuity (edges are 100–300 µm) for data fidelity. The
size-exclusion channel is where the two views meet, since it carries the only
nanometre-scale quantity (the hydrodynamic radius) in an otherwise
micrometre-scale model.

**Limitations.** The study's principal limitation is sample size and its
asymmetry: six squamous patients across two platform generations against one
melanoma patient; all inferential statistics are patient-level, and no
conclusion is drawn from the cross-cohort comparison. Sequencing depth spans
29-fold and is partly confounded with patient; results that track depth (the
rearrangement control) are presented as controls rather than evidence.
Scale parameters ξ₀ and β are not mesh measurements — the rank-normalised
input fixes the exclusion fraction by construction — so the 97.5% share is a
statement about the model given β, with the full dependence reported. Inputs
are signature scores rather than deconvolved proportions; most acutely, the
two-gene anti-PD-1 target set sits at the detection floor, so the weak
antigen channel is a resolution limit, not a biological finding. The
per-spot cellular field is a shortest-path accumulated cost, which permits
cost-free detours and reads closer to effective depth than to blockade; the
section-level cut is the topological object. No public dataset combines
spatial transcriptomics with checkpoint-response labels, so the bulk
concept check (one cohort significant in the expected direction, one not)
is motivational, and its instability is itself the finding we report: bulk
abundance discards the arrangement that makes matrix a barrier. Finally, the
melanoma cohort lacks treatment annotation, so post-treatment alteration of
barrier structure cannot be excluded there.

**Outlook.** Three extensions require no modification of the operators:
target-set replacement for ADC-style problems; deconvolution-based inputs for
a heterogeneous graph; and application to imaging-based platforms at
single-cell resolution, where the 10 µm cell and the 5.5 nm molecule are
separated by more than two orders of magnitude and the two barriers might yet
prove separable at a resolution spot-based assays cannot reach.

---

## 5. Conclusions

We formulated cell and antibody delivery in tumour tissue as two transport
operators on one spatial-transcriptomic graph: a source–sink minimum cut and
a screened Poisson diffusion–absorption field with size exclusion. The model
is deterministic, interpretable, parameter-honest, and fast enough (sub-second
per section on CPU) to serve as a component of translational analysis
pipelines. On 19 sections from 7 patients across two platform generations it
establishes three quantitative facts: the antibody barrier is
percolation-limited by size exclusion; the two barriers are mathematically
distinct yet physically co-located, with a coupling that survives removal of
every shared input in the squamous cohort; and the blockade behaves as a
connected structure. For matrix-rich tumours the model predicts that
matrix-directed intervention improves both classes of delivery at once — and
that neither quantity is measurable without resolving spatial arrangement.

---

## CRediT authorship contribution statement

**[作者姓名]**: Conceptualisation, Methodology, Software, Validation, Formal
analysis, Investigation, Data curation, Writing — original draft,
Visualisation. *(单作者声明；若导师或其他贡献者列名，按实际 CRediT 角色补)*

## Declaration of competing interest

The authors declare no competing interests.

## Data availability

All analysed data are publicly available: spatial transcriptomic cohorts from
the Gene Expression Omnibus (accessions GSE144239 [2] and GSE250636 [3]); the
hypoxia signature from MSigDB v7.1 [8]. The SPARTA implementation, all
analysis scripts, configuration files and per-section artefacts are available
at **https://github.com/Ezbenzino/sparta** and archived at
**[Zenodo DOI — 待回填]**.

## Acknowledgements

**[基金号/致谢待补]**

---

## References *(Elsevier numbered style)*

[1] P.L. Ståhl, F. Salmén, S. Vickovic, et al., Visualization and analysis of
    gene expression in tissue sections by spatial transcriptomics, Science 353
    (2016) 78–82.

[2] A.L. Ji, D.M. Rubin, K. Tharakan, et al., Multimodal analysis of the
    cellular and molecular landscape of cutaneous squamous cell carcinoma,
    Cell 182 (2020) 497–514.e20.

[3] O.E. Ospina, R. Manjarres-Betancur, G. Gonzalez-Calderon, et al.,
    B.L. Fridley, spatialGE is a user-friendly web application that
    facilitates spatial transcriptomics data analysis, Cancer Research 85
    (2025) 848. (source of GSE250636)

[4] V. Singhal, N. Chou, J. Lee, Y. Yue, J. Liu, W.K. Chock, L. Lin,
    Y.-C. Chang, K.H. Chen, S. Prabhakar, BANKSY unifies cell typing and
    tissue domain segmentation for scalable spatial omics data analysis,
    Nature Genetics 56 (2024) 431-441.

[5] G. Palla, H. Spitzer, M. Klein, et al., Squidpy: a scalable framework for
    spatial omics analysis, Nat. Methods 19 (2022) 171–178.

[6] H. Ren, B.L. Walker, Z. Cang, Q. Nie, Identifying multicellular
    spatiotemporal organization of cells with SpaceFlow, Nature
    Communications 13 (2022) 4076.

[7] K.A. Doiron, et al., Tissue cellular-neighbourhood analysis of spatial
    transcriptomics, Nature Methods 20 (2023).
    https://doi.org/10.1038/s41592-023-02124-2.

[8] A. Liberzon, C. Birger, H. Thorvaldsdóttir, et al., The Molecular
    Signatures Database (MSigDB) hallmark gene set collection, Cell Syst. 1
    (2015) 17–25.

[9] F.A. Wolf, P. Angerer, F.J. Theis, SCANPY: large-scale single-cell gene
    expression data analysis, Genome Biol. 19 (2018) 15.

[10] A.A. Hagberg, D.A. Schult, P.J. Swart, Exploring network structure,
    dynamics, and function using NetworkX, Proc. 7th Python in Science
    Conference (2008) 11–15.

[11] P.G. Doyle, J.L. Snell, Random Walks and Electric Networks, Mathematical
    Association of America, Washington, 1984.

[12] W. Hugo, J.M. Zaretsky, L. Sun, et al., Genomic and transcriptomic
    features of response to anti-PD-1 therapy in metastatic melanoma, Cell
    165 (2016) 35–44.

[13] N. Riaz, J.J. Havel, V. Makarov, A. Desrichard, W.J. Urba, J.S. Sims,
     F.S. Hodi, S. Martin-Algarra, R. Mandal, W.H. Sharfman, T.A. Chan,
     Tumor and microenvironment evolution during immunotherapy with
     nivolumab, Cell 171 (2017) 934-949.e16.

[14] R.K. Jain, Delivery of molecular and cellular medicine to tumors,
     Nature Reviews Drug Discovery 4 (2005) 619-632.

---

## 图表清单（投稿配置）

| 编号 | 内容 | 产物 | 状态 |
|---|---|---|---|
| **Graphical abstract** | **必选**：由 Fig 1 重制（一张图、两套边权、两个算子、一个共享基质） | `results/figures/fig1_framework.png` 改制 | **待制** |
| Table 1 | 队列概况 + 逐张准入结论 + 连通性 | `docs/table1_sections.md`（`run_18_table1.py` 生成） | 已出 |
| Fig 1 | 框架：空间图、两套边权语义、B_cell 与 B_mAb 算子 | `results/figures/fig1_framework.png/.pdf` | 已出 |
| Fig 2 | 三通道份额（对数点图）+ 交联份额对 β 的依赖 | `results/figures/fig2_driver_decomposition.png/.pdf` | 已出 |
| Fig 3 | 尺寸扫描（同一组织的 0.5–10 nm）+ 渗流/网格统计 | 由 `run_17_sensitivity.py` 产物重排 | 待重排 |
| Fig 4 | 耦合与共享输入去除（分层森林图/散点） | 依 R3/R3b 产物 | 待出 |
| Fig 5 | 反事实：S2 剂量反应 + S1 患者分层 | 依 `results/counterfactual/*.json` | 待出 |
| Table 2 | R3b 分层耦合表（正文 3.4 已含，排版为正式表） | `shared_ecm_check.json` | 已有数据 |
| Table 3 | 运行时对比（SPARTA/BANKSY式/Squidpy，19 张） | `runtime_benchmark.json`（run_19） | 已有数据 |

> 正文 3.2/3.5 的 β 表与 S2 剂量反应全表按 Elsevier 惯例移入补充材料
> （Supplementary Tables S1/S2），主文各留摘要行。
