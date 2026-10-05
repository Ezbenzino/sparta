# SPARTA: dual graph operators for modelling immune-cell migration and IgG-sized transport barriers from spatial transcriptomics

**Yize Li**  
Hangzhou Medical College, Hangzhou, China  
Corresponding author: Yize Li (lllyz630031258@gmail.com)

## Abstract

**Background and Objective:** Spatial transcriptomics maps expression onto tissue coordinates, but common analyses identify neighbourhoods or domains rather than quantify transport. We developed SPARTA, a graph-based framework that represents cellular migration and macromolecular transport with separate operators on a shared spatial graph.

**Methods:** We analysed 19 primary-cohort sections from seven patients with cutaneous squamous cell carcinoma or melanoma. A source–sink minimum cut represented a model-defined cellular barrier; a screened diffusion–absorption equation generated an IgG-sized molecule deficit field. Scores were rank-normalised within section. Field association was assessed after adjustment for graph distance to vasculature and tested using 500 graph-spectral sign-randomised surrogates per section, with Benjamini–Hochberg (BH) correction. We repeated the association summary in the 15 sections with a measured `Ag_target` score and assessed the contribution of disconnected nodes by restricting estimates to each largest connected component. Three additional public sections formed a separate exploratory arm.

**Results:** The adjusted field correlation was positive in 18/19 primary sections (median partial ρ = 0.217). Among the 15 sections with a measured `Ag_target` score, all correlations were positive and 13 remained below BH q = 0.05 under the same section-level graph-spectral surrogates. Restricting to largest connected components did not change any correlation's sign; the median changed from 0.217 to 0.201. In a selection-matched in-model counterfactual, the median ratio of residual barrier after scattered versus contiguous removal was 1.053; 10/19 sections had raw p < 0.05 and 8/19 had BH q < 0.05. In the exploratory arm, all three adjusted correlations were positive (median ρ = 0.281); one spatial-null result was at the Monte Carlo boundary (q = 0.0499).

**Conclusions:** SPARTA provides two computationally distinct, model-defined transport summaries on spatial transcriptomic graphs. The observed coupling is not evidence of measured antibody exposure, pharmacological specificity, or clinical response. External replication is limited by three sections from two individuals, and independent functional calibration remains necessary.

**Keywords:** spatial transcriptomics; graph theory; minimum cut; diffusion–absorption; tumour microenvironment; computational biomedicine

## Highlights

- Two graph operators model cell migration and IgG-sized transport fields.
- Model fields correlate positively in 18 of 19 primary sections.
- A graph-spectral null retains BH significance in 14 of 19 sections.
- Selection-matched counterfactual effects are modest and model-defined.

## 1. Introduction

Spatial transcriptomics preserves local gene-expression measurements together with tissue coordinates [1]. In tumour studies, these data support mapping of cell states, neighbourhoods and tissue domains, including cutaneous squamous cell carcinoma (cSCC) and melanoma [2,3]. Most such analyses describe where expression patterns occur. They do not directly quantify a specified source-to-target transport problem.

Transport by a migrating T cell and transport of an antibody-sized molecule differ in scale and mechanism. A graph representation permits both to be posed on the same sampled tissue geometry, while using operators appropriate to each question. For cellular migration, a source–sink minimum cut returns a bottleneck capacity and a corresponding cut set. For a diffusing macromolecule, a conductance-weighted diffusion–absorption equation produces a steady-state field. These graph quantities are computational summaries; their biological meaning depends on how graph edges and node scores are defined.

We introduce SPARTA, which uses a minimum-cut operator for a model-defined cellular barrier and a screened diffusion–absorption operator for an IgG-sized transport deficit. We ask whether the resulting per-spot fields co-vary after adjustment for distance to vasculature, whether the association persists after removing an input shared by both operators, and how a selection-matched spatial counterfactual behaves. We also apply the locked field-analysis pipeline to three external public sections, treating them as exploratory because two are adjacent sections from one breast-cancer patient and the third is a non-tumour lymph node.

The two operators and the cohort structure are summarised in Fig. 1.

The study is deliberately framed as a computational method demonstration. It does not measure extracellular-matrix mesh size, antibody concentration, T-cell passage, or treatment response. In particular, the implemented absorption proxy uses expression of **CD274** and **PDCD1LG2** (PD-L1 and PD-L2 ligands); it is not a receptor-specific or pharmacologically calibrated anti-PD-1 model.

## 2. Materials and methods

### 2.1. Data and analysis cohorts

The primary analysis comprised 19 spatial transcriptomic sections from seven patients. The cSCC arm included 15 treatment-naïve sections from six patients in GSE144239 [2]: four Visium sections from two patients and 11 first-generation spatial-transcriptomics sections from four patients. The melanoma arm included four extracranial metastatic sections from one patient in GSE250636 [3]. Treatment status was not annotated for the melanoma sections; they were retained under a documented cohort-level waiver and are interpreted as exploratory. One cSCC section (CSCC13) was excluded because its median UMI count was below the platform-specific admission threshold.

The separate external arm used two adjacent 10x Genomics breast-cancer sections from one patient and one human lymph-node section. The external sections were not pooled with primary-cohort summaries. The lymph node is non-tumour, and the epithelial region used to define its model sink was not histopathologically verified.

| Analysis arm | Tissue and platform | Sections | Independent individuals | Use in inference |
|---|---|---:|---:|---|
| Primary cSCC | 4 Visium; 11 first-generation ST | 15 | 6 | Primary descriptive analysis |
| Primary melanoma | Visium extracranial metastases | 4 | 1 | Included in section-level primary summary; patient-level inference not performed |
| External | 2 adjacent breast Visium sections; 1 lymph-node Visium section | 3 | 2 | Exploratory, kept separate |

### 2.2. Spatial graph and node compartments

Coordinates were scaled so that the median nearest-neighbour distance matched the nominal platform pitch (100 μm for Visium and 200 μm for first-generation ST). Radius-adjacency graphs used a 150 μm radius for Visium and a 300 μm radius for first-generation ST. Spots below 500 UMI were removed before graph construction. The resulting graphs had 370–2,673 nodes and 1–54 connected components per section; 90.0–100% of nodes belonged to the largest component.

Endothelial-rich nodes were identified using the 80th percentile of the endothelial signature. The immune-entry source comprised nodes in this set whose mean neighbouring T/NK signature exceeded its within-vessel 60th percentile. Malignant nodes were above the 70th percentile of the tumour-type-specific malignant score; the sink comprised malignant nodes at or beyond the median graph distance from non-malignant nodes. A border-based source fallback was available when fewer than five vascular nodes were identified. For the section-level minimum cut, the implementation retained the largest connected component. Eight of 19 primary sections had at least one source or sink outside that component; these nodes cannot contribute to the main component's source–sink flow.

### 2.3. Expression signatures

Signature scores were computed with Scanpy's `score_genes` procedure [6] and rank-normalised to [0,1] within each section. The curated inputs included endothelial and T/NK scores; tumour-type-specific malignant markers; a curated subset of core-matrisome genes; CAF markers; and seven crosslinking-associated genes (LOX, LOXL1–3, PLOD1/2, TGM2). Hypoxia used the MSigDB HALLMARK_HYPOXIA set [12]. These are expression scores, not deconvolved cell proportions or measured matrix concentrations.

The model's target-expression input (`Ag_target`) used CD274 and PDCD1LG2. These genes encode PD-L1 and PD-L2, which are ligands of PD-1; they do not encode PD-1 (PDCD1). The score is therefore described as a ligand-expression proxy in a model absorption term, not as a specific anti-PD-1 binding-site distribution. The score and effective dissociation parameter were not calibrated to protein concentration or binding kinetics. If a signature was unavailable because fewer than three genes were detected, the implementation filled it with 0.5 and recorded the missing key. Four primary sections lacked the two-gene target score; one additional section lacked the efflux score. Results depending on a filled term are treated as degraded or exploratory.

### 2.4. Cellular minimum-cut operator

For each undirected graph edge (u,v), the cellular capacity was defined as

`c_uv = sigmoid[a − b_ECM mean(ECM_u, ECM_v) − c_CAF mean(CAF_u, CAF_v)]`.

The maximum flow between the source and sink sets was computed on the largest connected component, and the section-level barrier score was `B_cell = 1/(maxflow + ε)`. The corresponding minimum-cut edges and incident nodes define a model bottleneck band. For the per-spot field used in association analyses, each edge cost was `1/c_uv`, and Dijkstra's algorithm returned the lowest accumulated cost from the source set. This per-spot shortest-path field is not the same quantity as the section-level minimum-cut score and can reflect low-cost detours around a bottleneck. Nodes disconnected from every source receive the maximum finite field value and are marked unreachable. As a sensitivity analysis, we recomputed the field association after restricting it to the largest connected component in each section.

### 2.5. IgG-sized diffusion–absorption operator

For each edge, crosslink score x was the mean of the two endpoint scores. A qualitative scale law defined an effective mesh parameter `ξ = ξ₀ exp(−βx)`. For molecular radius r, the size factor was `φ_size = (1 − r/ξ)^2` when r < ξ and zero otherwise. Edge conductance was `g_uv = g₀ exp[−λ mean(ECM_u, ECM_v)] φ_size + g_floor`. The default settings were r = 5.5 nm, ξ₀ = 20 nm, β = 3, λ = 3, and `g_floor = 10⁻⁶`.

The vessel nodes had fixed unit concentration. A node absorption proxy was `κ = Ag/(Ag + Kd,eff)`, with default `Kd,eff = 0.5` and weight `κw = 1`. The steady field φ solved a screened graph-Poisson system `(L_g + diag(κwκ + 10⁻⁹))φ = 0` away from vessel nodes, where `L_g` is the conductance-weighted graph Laplacian. The reported field was `B_mAb = −log(φ)`. This is a unitless model score, not a measured concentration or calibrated resistance. If absorption is set to zero, the equation reduces to a harmonic graph field related to effective resistance [8].

The 5.5 nm value is a literature-based IgG Stokes-radius approximation [11]; r is the molecular radius (not diameter) and is compared directly to the effective mesh half-opening ξ, so an edge is excluded when the hydrodynamic radius exceeds half the nominal pore size. By contrast, ξ₀ and β are qualitative scale parameters: x is rank-normalised and dimensionless, so the model fixes the nominal range and complete-exclusion threshold from chosen parameter values. At the defaults, the threshold is `x = ln(ξ₀/r)/β = 0.430`. Thus the observed 60.0–65.2% of edges assigned zero size-exclusion conductance is a model-law output, not a measurement of tissue mesh. We report this parameter dependence explicitly (Fig. 3).

### 2.6. Statistical analysis

The primary association analysis used the per-spot `B_cell` shortest-path field and `B_mAb` deficit field. Both were rank-transformed, residualised on a quadratic function of rank-transformed weighted graph distance to the nearest vessel (edge lengths in μm), and correlated using Spearman's ρ. Nominal point-level p-values treat spots as observations and are reported only as descriptive tests because spatial dependence reduces the effective sample size.

To assess spatial alignment under a graph-based surrogate, the normalised graph Laplacian was eigendecomposed; absolute spectral coefficients of `B_mAb` were retained and assigned independent random signs. The surrogate field was reconstructed and the adjusted correlation recalculated for 500 draws per section. The one-sided empirical p-value used `(k+1)/(N+1)`, giving a minimum attainable p of 1/501. This graph-spectral surrogate follows a graph-signal randomisation approach [9]; it preserves graph-spectral power, but it does not guarantee preservation of every spatial statistic. BH correction was applied across the 19 primary sections and separately across the three external tests [10]. For the input-availability sensitivity analysis, the same per-section surrogate tests were retained for the 15 sections with a measured `Ag_target` score, and BH correction was recalculated within that 15-section family. The `efflux` score is not used in either field or the graph-distance adjustment. Section counts and patient summaries are descriptive. No patient-level inferential model was fitted, and sections from the same individual are not treated as independent patients.

### 2.7. Selection-matched counterfactual

At the selected 20% removal fraction of eligible cut nodes (after excluding protected source/sink nodes), one contiguous arc was selected from up to eight candidate arcs as the candidate yielding the lowest residual `B_cell`. For each matched-control group, eight scattered sets of equal size were generated within the eligible cut set and the lowest residual barrier was retained, applying the same best-of-eight search pressure to both arms. One hundred matched groups were generated per section. The effect ratio was `mean residual barrier after scattered removal / residual barrier after contiguous-gap removal`. A ratio above 1 therefore means the contiguous-gap arm leaves a lower **model** barrier. Empirical p-values use the same lower-tail comparison with finite-sample correction. BH was applied across the 19 primary tests. This experiment alters model scores on a fixed graph; it is not a tissue intervention.

### 2.8. Exploratory bulk response analysis

As an exploratory concept check, a bulk barrier score was compared between responders and non-responders in two public anti-PD-1 cohorts (GSE78220 and GSE91061). We report the direction, group medians, Mann–Whitney p-value, and AUC oriented as the probability that a non-responder has a higher score than a responder. No time-stamped protocol established that the directional comparison was specified before response outcomes were examined. These analyses are not treated as a prediction model or validation of the spatial method.

### 2.9. Software and reproducibility

SPARTA is implemented in Python; the core operators use NumPy, SciPy, and NetworkX [7]. Figures in this submission were generated programmatically from the stored analysis outputs using Matplotlib. The code and analysis instructions are available in the public GitHub repository. A versioned archive containing the exact configuration, per-section outputs, figure inputs, and software manifest must be deposited before submission; the archive DOI will be added to the Data availability statement after deposit. Runtime results are implementation- and machine-specific and are not used as a biological validation claim.

## 3. Results

### 3.1. Cohort and graph coverage

The primary cohort comprised 19 analysed sections from seven patients, spanning two tumour types and two spatial-platform generations. After QC, graph size ranged from 370 to 2,673 spots. The largest component contained 90.0–100% of each section's retained nodes. In eight sections, at least one source or sink node was outside the largest component. The stranded nodes accounted for <5% of combined source–sink nodes in all but one section (MEL02, 9.1%); per-section counts are reported in the repository outputs. Because the section-level minimum cut is calculated on the largest component, the resulting `B_cell` summary can be biased upward where relevant compartments are stranded in smaller components.

### 3.2. Association between model-defined fields

The vessel-distance-adjusted partial correlation was positive in 18/19 primary sections (range −0.0004 to +0.385; median 0.217). Seventeen sections passed the nominal point-level test. Under the graph-spectral surrogate, 15/19 had p < 0.05; 14/19 remained below BH q = 0.05. The exception with near-zero correlation was CSCC15; MEL02 was positive but small. Every patient had a positive median section-level ρ, reported descriptively without patient-level p-values.

Section-level estimates and adjusted q-values are shown in Fig. 2.

| Analysis subset | Sections | Median partial ρ | Positive sections | BH q < 0.05 under graph-spectral null |
|---|---:|---:|---:|---:|
| All primary sections | 19 | 0.217 | 18/19 | 14/19 |
| `Ag_target` score available | 15 | 0.217 | 15/15 | 13/15 |
| Largest connected component only | 19 | 0.201 | 18/19 | Not recalculated |

The first two rows use the same 500-draw section-level surrogate results, with BH correction applied separately to each stated family. The largest-component row is a directional sensitivity analysis only; no component-restricted surrogate tests were run.

In the four Visium cSCC sections the median partial ρ was 0.270; in 11 first-generation ST cSCC sections it was 0.207. This cross-platform agreement is a within-dataset robustness observation, not an estimate from independent platform-level cohorts. The external arm also had positive partial correlations in all three sections (median 0.281). Its two breast sections represent one patient; the third section is a non-tumour lymph node. The BRCA02 spatial-null BH-adjusted value was q = 0.0499 with 500 surrogates, at the Monte Carlo boundary, so it is sensitive to finite permutation resolution.

At the 75th-percentile discordance definition, the fraction of spots in the low-`B_cell`/high-`B_mAb` quadrant was below the 6.25% independence expectation in 18/19 primary sections. This descriptive result is another representation of the field association, not an independent measurement of delivery.

Input availability was evaluated separately by analysis. The field association uses the `Ag_target` term in `B_mAb`, so we excluded the four sections missing that score and retained MEL01, whose only missing input was `efflux`. The remaining 15/15 correlations were positive (median partial ρ = 0.217; range 0.055–0.385); 13/15 remained below BH q = 0.05 under the same graph-spectral surrogate tests with correction recalculated for this subset. By contrast, the matched S2 counterfactual uses ECM/CAF scores and the source–sink graph, and does not use either `Ag_target` or `efflux`; it therefore remains a 19-section result (median ratio 1.053; 8/19 below BH q = 0.05). The earlier five-section `clean14` S2 summary used the legacy un-matched analysis and is not evidence for robustness of the current matched-selection result.

Disconnected nodes also had a small but measurable effect on the field association because they receive finite replacement values. Restricting estimates to the largest connected component left all section-level signs unchanged (18/19 positive) and shifted the median partial ρ from 0.217 to 0.201; the median absolute section-level change was 0.006. MEL02 had 27 stranded source/sink nodes (9.1% of its source/sink set); its partial ρ changed from 0.055 to 0.047 after restriction. The S2 cut and counterfactual are already computed within the largest component, so stranded nodes do not enter that statistic. These checks address graph disconnection within the current dataset and do not establish robustness to tissue segmentation or graph-construction choices.

### 3.3. Shared input and parameter dependence

Both operators use an ECM-related score, creating a direct shared-input route to positive association. An ablation setting the shared ECM contribution to zero reduced the median partial correlation from 0.217 to 0.152 across the 19 sections; 15/19 remained nominally significant in the point-level test. Because the ablated variant was not re-evaluated with the graph-spectral null, this result is treated as exploratory and does not establish that the residual association is independent of spatial autocorrelation.

At the default β = 3, a variance-ablation analysis assigned a median 97.5% of `B_mAb` variance to the crosslink/size channel among sections where the target signature was scoreable. This percentage changes sharply over the β grid. Together with the rank-defined exclusion threshold, this shows that the channel share is highly dependent on a qualitative model parameter; it should not be read as evidence that crosslinking biologically accounts for 97.5% of antibody-delivery resistance.

The selected size-exclusion threshold and its dependence on β are illustrated in Fig. 3.

### 3.4. Selection-matched removal and size scan

In the matched S2 analysis, the median ratio of residual barrier after scattered removal to residual barrier after contiguous removal was 1.053 (range 0.939–1.217); 17/19 ratios exceeded 1. Ten sections had raw p < 0.05 and eight remained below BH q < 0.05. The direction corresponds to a lower residual model barrier after contiguous removal. Its magnitude is modest, and the result depends on candidate search, model edge weights, and the selected 20% removal fraction. The same analysis in the three external sections gave ratios of 1.107, 1.207, and 1.016; two had raw p < 0.05 and the lymph-node section did not. These are graph counterfactuals, not evidence that opening a biological tissue band changes drug penetration.

Across the fixed-input molecular-radius scan, the mean tumour-core `B_mAb` score increased with the specified radius in the analysed sections. This is expected from the programmed size-exclusion law and demonstrates operator sensitivity to r. It does not validate the assumed mesh law or quantify in-vivo antibody exposure.

The selection-matched ratios and fixed-input radius scan are shown in Figs. 4 and 5, respectively.

### 3.5. Relation to domain methods and runtime

Domain-boundary segmentation and transport bottleneck identification address different outputs. We compared the minimum-cut edges against boundaries from a BANKSY-style eight-domain implementation on the same graph: cut edges were enriched among domain boundaries in most sections (median fold enrichment 1.15, range 0.68–1.87), but only a median of 14.3% of domain boundaries coincided with cut edges. This is an implementation-specific descriptive comparison, not the official BANKSY implementation or a head-to-head performance evaluation. The domain-count scan varied the number of clusters from 4 to 12. In an additional sensitivity scan at eight domains, varying the neighbourhood weight λ from 0.1 to 0.7 produced median cut-edge enrichment values between 1.04 and 1.15 across the four settings; median boundary precision ranged from 13.7% to 14.3%. These scans assess selected segmentation choices but do not remove dependence on preprocessing or establish an optimal parameter. Squidpy's neighbourhood enrichment test [5] used the same adjacency graph that contributed to the spatial smoothing of the BANKSY-style labels. We therefore treat it as a spatial-structure check on these labels, not an independent validation. These analyses do not establish superiority over current spatial-domain methods. In a same-machine runtime exercise across the 19 primary sections, the stored core pipeline times were shorter than the selected comparator implementations; these timings are environment- and implementation-dependent and are provided as reproducibility context rather than a general performance claim.

### 3.6. Exploratory bulk response check

In GSE78220, the bulk barrier score was higher in non-responders (median 0.271) than responders (−0.188), with AUC(NR>R) = 0.759 (p = 0.021). In GSE91061, the direction was reversed (AUC = 0.428, p = 0.495); a stricter responder definition gave AUC = 0.426 (p = 0.518). This cross-cohort inconsistency and exploratory design preclude a predictive or comparative claim. The analysis is contextual only and is not evidence that spatial methods outperform bulk scores.

## 4. Discussion

SPARTA demonstrates a way to compute two distinct transport summaries from the same spatial graph: a source–sink cut for a cellular migration proxy and a screened diffusion–absorption score for an IgG-sized molecule. In this dataset, their per-spot fields tended to co-vary after vessel-distance adjustment, including under a graph-spectral surrogate in most sections. The external arm is consistent in direction but is too small and dependent to support independent patient-level validation.

Several design features limit the biological interpretation. First, the source cohort contains only seven patients, including one melanoma patient; 19 sections are not 19 independent subjects. The point-level correlations and multiple-testing correction do not replace a patient-clustered model. Second, the crosslink-to-mesh map uses within-section ranks and qualitative scale parameters. It does not estimate physical mesh size from expression data. Third, `Ag_target` represents PD-L1/PD-L2 ligand expression and the sink strength uses an uncalibrated score and effective dissociation constant. The resulting field is not specific to anti-PD-1 pharmacology and cannot be interpreted as a predicted therapeutic concentration.

The field association is also partly built into the model because matrix expression contributes to both operators. The ECM ablation reduces the association, and its nominal tests have not been repeated with the spatial surrogate. Four primary sections have no measured target signature and receive neutral-value fill; this makes the absorption component constant for those sections. Missing-score diagnostics are preserved in the outputs, but these sections weaken the biological interpretation of the pooled pattern.

External validation is limited to three sections from two individuals, all from the same Visium generation. The breast pair is adjacent tissue from one patient, the lymph node is non-tumour, and the epithelial sink definition has not been verified against pathology. The external results support reproducibility of a model-field calculation across datasets; they do not validate a tumour barrier in independent patients.

The exploratory bulk response comparison was inconsistent between cohorts and is described in Section 3.6. It does not validate clinical discrimination or establish an advantage of spatial over bulk methods.

We also examined exploratory associations between the per-spot model fields and measured expression signatures (T/NK, CD8, proliferation) within each section. Because the immune-entry source points are themselves selected from the vascular T/NK neighbourhood, the correlation between the `B_cell` field and the T/NK signature is not an independent validation of the field; it is partly built into the source definition. Across the 19 sections the median partial Spearman ρ was −0.022 for `B_cell` vs T/NK (6/19 sections p < 0.05) and +0.046 for `B_mAb` vs proliferation (7/19 p < 0.05). One first-generation section (CSCC12) lacked a CD8 field and used T/NK as a proxy. These weak and inconsistent associations are reported as descriptive alignment, not as biological ground truth.

Further evaluation should prioritise patient-level replication, pathology-guided compartments, spatial proteomics or direct mesh measurements, and independent measurement of antibody penetration or concentration. A useful next step is a held-out perturbation experiment that compares the model's predicted changes with measured transport. Until such data exist, the framework is best interpreted as a computational hypothesis generator rather than a clinical decision tool.

## 5. Conclusions

SPARTA combines a minimum-cut operator and a screened diffusion–absorption operator on spatial transcriptomic graphs. The two resulting model-defined fields were positively associated in most primary sections, with 14/19 passing BH correction under a graph-spectral surrogate. A selection-matched counterfactual produced a modest median effect and remains entirely in-model. The method is computationally reproducible, but biological calibration, independent patient-level validation, and functional transport measurements remain necessary before interpreting its scores as tissue permeability or drug delivery.

## CRediT authorship contribution statement

Yize Li: Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Writing – original draft, Visualization.

## Declaration of competing interest

The author declares no competing financial interests or personal relationships that could have appeared to influence the work.

## Funding

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## Ethics statement

Not applicable. This study is a secondary analysis of publicly available, de-identified spatial transcriptomic data and involved no new recruitment or identifiable participant information.

## Data availability

The primary data are available from GEO under accessions GSE144239 and GSE250636. External sections are from the 10x Genomics public spatial-expression datasets listed in the repository documentation. The source code is available at https://github.com/Ezbenzino/sparta. The versioned archive containing the exact analysis code, configuration, per-section outputs, and figure inputs will be deposited before submission; its DOI will be added after the archive has been minted and verified.

## Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During preparation of this submission draft, OpenAI Codex (GPT-6) was used to assist with manuscript organisation, language drafting, and programmatic figure-code authoring. No generative-image model was used to create or edit the figures. The author must review and edit the material, verify every scientific statement and figure against the source data, and take full responsibility for any submitted version.

## References

[1] P.L. Ståhl, F. Salmén, S. Vickovic, et al., Visualization and analysis of gene expression in tissue sections by spatial transcriptomics, *Science* 353 (2016) 78–82. https://doi.org/10.1126/science.aaf2403.

[2] A.L. Ji, A.J. Rubin, K. Thrane, et al., Multimodal analysis of composition and spatial architecture in human squamous cell carcinoma, *Cell* 182 (2020) 497–514.e22. https://doi.org/10.1016/j.cell.2020.05.039.

[3] H. Alhaddad, O.E. Ospina, M.L. Khaled, et al., Spatial transcriptomics analysis identifies a unique tumor-promoting function of the meningeal stroma in melanoma leptomeningeal disease, *Cell Reports Medicine* 5 (2024) 101606. https://doi.org/10.1016/j.xcrm.2024.101606.

[4] V. Singhal, N. Chou, J. Lee, et al., BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis, *Nature Genetics* 56 (2024) 431–441. https://doi.org/10.1038/s41588-023-01576-2.

[5] G. Palla, H. Spitzer, M. Klein, et al., Squidpy: a scalable framework for spatial omics analysis, *Nature Methods* 19 (2022) 171–178. https://doi.org/10.1038/s41592-021-01358-2.

[6] F.A. Wolf, P. Angerer, F.J. Theis, SCANPY: large-scale single-cell gene expression data analysis, *Genome Biology* 19 (2018) 15. https://doi.org/10.1186/s13059-017-1382-0.

[7] A.A. Hagberg, D.A. Schult, P.J. Swart, Exploring network structure, dynamics, and function using NetworkX, in: *Proceedings of the 7th Python in Science Conference*, 2008, pp. 11–15.

[8] P.G. Doyle, J.L. Snell, *Random Walks and Electric Networks*, Mathematical Association of America, Washington, DC, 1984.

[9] D.I. Shuman, S.K. Narang, P. Frossard, A. Ortega, P. Vandergheynst, The emerging field of signal processing on graphs: extending high-dimensional data analysis to networks and other irregular domains, *IEEE Signal Processing Magazine* 30 (2013) 83–98. https://doi.org/10.1109/MSP.2012.2235192.

[10] F. Benjamini, Y. Hochberg, Controlling the false discovery rate: a practical and powerful approach to multiple testing, *Journal of the Royal Statistical Society: Series B* 57 (1995) 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x.

[11] G.D. Thomas, Effect of dose, molecular size, and binding affinity on uptake of antibodies, in: *Methods in Molecular Medicine*, vol. 25, Humana Press, 2000, pp. 115–132. https://doi.org/10.1385/1-59259-075-6:115.

[12] A. Liberzon, C. Birger, H. Thorvaldsdóttir, M. Ghandi, J.P. Mesirov, P. Tamayo, The Molecular Signatures Database (MSigDB) hallmark gene set collection, *Cell Systems* 1 (2015) 417–425. https://doi.org/10.1016/j.cels.2015.12.004.

## Figure legends

[[FIGURE:figure_1_framework]]

**Figure 1. SPARTA's two model-defined transport operators and analysis cohorts.** The schematic shows a shared spatial graph used by a source–sink minimum-cut operator and a screened diffusion–absorption field. The node layout and fields are illustrative and do not depict patient tissue. Primary cohorts and the separate exploratory external arm are shown below.

[[FIGURE:figure_2_coupling]]

**Figure 2. Adjusted association between model fields.** Points show the per-section partial Spearman correlation after quadratic rank-space adjustment for weighted graph distance to vasculature. The BH-adjusted q-value is printed at right. Filled markers indicate BH q < 0.05 and open markers q ≥ 0.05. The external sections are displayed separately and were not pooled with the primary cohort. The BRCA02 value (q = 0.0499) is at the finite-resolution boundary of 500 surrogate draws.

[[FIGURE:figure_3_parameter_sensitivity]]

**Figure 3. The size-exclusion law is parameterised rather than measured.** (A) Model-defined effective mesh parameter versus within-section rank-normalised crosslink score at ξ₀ = 20 nm and β = 3; the vertical threshold follows from the selected parameters. (B) Percentage of graph edges assigned zero size-exclusion conductance over the stored β sensitivity grid at λ = 3 for five example sections. These values are model outputs and should not be interpreted as physical tissue measurements.

[[FIGURE:figure_4_s2_matched]]

**Figure 4. Selection-matched S2 counterfactual.** Points show the residual barrier after scattered removal divided by the residual barrier after contiguous-gap removal of the same fraction of cut nodes. A ratio above 1 means the contiguous-gap arm has a lower model barrier. Filled markers indicate BH q < 0.05 within the primary 19-section family; open markers indicate q ≥ 0.05. External ratios are shown separately; the two breast sections are from one patient and the lymph node is non-tumour.

[[FIGURE:figure_5_size_scan]]

**Figure 5. In-model molecular-radius scan.** Thin lines show section-specific mean tumour-core `B_mAb` scores across molecular radii with tissue inputs held fixed; the heavy line is the median across 19 primary sections. The 5.5 nm reference marks the assumed IgG Stokes radius. The score is not a measured antibody concentration or transport response.
