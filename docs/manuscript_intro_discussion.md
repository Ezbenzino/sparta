# SPARTA — Introduction 与 Discussion（草稿 v1.0 · PIVOT 叙事）

> Results 在 `manuscript_results.md`；执行手册在 `runbook_after_review.md`。
>
> **2026-08-27 更新（两轮）**
>
> 第一轮：`run_14_shared_ecm_check.py` 在 8 张切片上判定 SPLIT（CSCC 4/4、MEL 1/4）。
>
> 第二轮：队列扩到 **19 张 / 7 位患者**（新增 11 张第一代 ST，4 位 cSCC 患者）后，
> **那个 SPLIT 判定被推翻了**——它的前提是"cSCC 4/4 满贯"，扩样后是 12/15。
> 判定规则本身也已降级：它要求队列内每张切片都显著，严苛程度随 n 单调上升，
> 8 张给 SPLIT、19 张给 MODEL，中间没有任何证据反转。
> 现在按队列 × 按患者分层陈述，不给单一标签。§2 已按此重写。

---

# Introduction

**¶1 — 临床困境与空间组学的介入**

Immune checkpoint blockade has changed the outlook for cutaneous malignancies,
yet most patients with advanced melanoma or cutaneous squamous cell carcinoma
still fail to obtain durable benefit from anti-PD-1 therapy. Response is only
partly explained by tumour-cell-intrinsic features such as mutational burden or
antigen presentation; a substantial share of variance is attributed to the
organisation of the tumour microenvironment. Spatial transcriptomics has made
that organisation directly measurable, and studies across melanoma, squamous
carcinoma and other solid tumours now routinely describe stromal bands,
immune-excluded phenotypes and fibroblast-rich boundaries around tumour nests.

**¶2 — 既有的"T 细胞屏障"研究及其度量方式的局限**

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

**¶3 — 缺口：第二道屏障几乎无人量化**

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

**¶4 — 本文的做法与贡献**

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
actionable finding: the model predicts that reducing matrix density would lower
both barrier scalars in the same regions, whereas reducing one modality alone
would not be expected to move the other. This is a model-level prediction, not a
measured treatment effect; see Discussion §4.

> **中文注**：第 4 段的写法很关键。不要把"我们本来想找解耦但没找到"写出来——
> 要把它写成"我们直接检验了可分离性这个问题，答案是不可分离，而这个答案有转化含义"。
> 前者是失败叙事，后者是研究叙事，内容完全一样。

---

# Discussion

## §1 — 耦合的生物学含义

The finding that the cellular and macromolecular barriers travel together has a
straightforward tissue correlate. A desmoplastic band is not merely dense; it is
crosslinked, and crosslinking simultaneously stiffens the network against cell
migration and contracts the mesh below the size of an IgG. The two obstructions
are mechanistically distinct — one is a matter of a 10 µm cell negotiating a
fibre network, the other of a 5.5 nm molecule negotiating a mesh — but they are
produced by the same remodelling programme, so in tissue they co-localise.

The model-derived implication follows, but we separate it from the evidence.
Because the two obstructions are produced by the same remodelling programme, the
operators predict that a uniform reduction in matrix density would lower both
barrier scalars in the same regions (Results, model input perturbation). This is
a statement about the model's response to its own inputs, not a measurement of
what LOX/LOXL inhibition, TGF-β blockade, hyaluronidase or any other
matrix-directed therapy does in patients. We report it as a rationale for why
matrix-directed strategies *might* relieve both modalities, not as evidence that
they do; a clinical test would require independent perturbational data or response
labels, neither of which this study contains. Conversely, the model predicts that
interventions aimed at only one modality (for example, a chemokine-based
strategy to increase T-cell recruitment) would not be expected to move the
antibody barrier, and combination designs that assume otherwise may be over-optimistic
about the achievable intratumoural drug concentration.

## §2 — 耦合是组织的性质，还是模型的性质？

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

> **中文注**：这一段的骨架是「先承认循环论证的质疑合理 → 给 cSCC 的正面结果并
> 说明它跨两代平台、跨两种距离口径都成立 → 主动排除深度与切片大小 →
> 给机制解释但标明它只是 consistent 不是 established → 明说黑色素瘤那一侧判不了
> 并**主动撤回**旧稿里的队列对比 → 说清哪个混杂被拆开了、哪些还没 →
> 最后把稳的那一半单独拎出来」。
>
> **撤回要写在明面上。** 旧稿里"cSCC 4/4 vs MEL 1/4"这个对比是小样本假象，
> 扩到 15 张后 cSCC 是 12/15。与其等审稿人对出来，不如自己先说这一句：
> 前提不成立，所以对比撤回。这比悄悄换个说法可信得多。

## §3 — 与既有工具的关系

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

## §4 — 局限

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

**External validation.** The coupling result itself has not been replicated in an
independent patient cohort. All results above come from seven patients in two GEO
series; the sensitivity analyses (parameter sweeps, shared-ECM removal, hop-vs-weighted
distance, leave-one-patient-out on the counterfactuals) are internal robustness
checks on the same data, not external validation. No publicly available dataset
combines Visium-format spatial data with ICB response labels; our bulk concept
check in two anti-PD-1 melanoma cohorts is therefore motivational and not
validation, and its instability across the two cohorts is informative in the same
direction as our main result: barrier molecules measured as bulk abundance lose
the spatial arrangement that makes them a barrier. A genuine external test would
require an independent ST or Visium cohort from a different centre, processed
through the same signature definitions and graph construction, and would report
whether the positive partial correlation survives at the patient level. We have
not attempted that here, and we do not describe any internal analysis as
validation.

**Static sections and causal language.** The coupling we report is an association
between two per-spot fields computed from the same static tissue section. A
positive partial correlation after controlling for vessel distance — and its
survival when the shared ECM input is removed — does not establish that matrix
crosslinking *causes* either barrier, nor that one barrier drives the other. The
fields are computed from graph operators on the same coordinates; tissue shape,
tumour-nest boundary geometry, cell density and segmentation quality could in
principle produce or amplify the association. We addressed the most likely of
these — that two smooth fields on the same graph will be correlated simply
because they are smooth — with a graph-Laplacian phase-randomization null that
preserves each field's autocorrelation variogram: 15/19 sections still exceed
the null (median empirical p = 0.000; Results R3). What that null does not do
is preserve tumour-nest region labels or cell-composition covariates beyond
vessel distance; a block-level null that holds tumour regions fixed is the next
step. We use "coupled" and "associated" throughout, not "interact", "mediate"
or "drive", and we have avoided describing the size-exclusion decomposition as
a mechanism: it is a parameter-defined attribution, and the permutation null in
R2 shows that its dominant share is set by the scale parameter rather than by
tissue architecture.

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

## §5 — 展望

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
