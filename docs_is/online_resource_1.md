---
title: "Online Resource 1: supplementary methods, tables and figures"
subtitle: "SPARTA: graph transport operators and structural null models for immune-cell and IgG-transport barriers in spatial omics"
---

This document accompanies the main text. Section, equation, table and reference numbers without the prefix S refer to the main text. All numbers are generated from the archived result files (Online Resource 3); per-section values are also given in Online Resource 2.

# S1 Parameter sensitivity (ξ₀/β/λ)

We examined parameter dependence before interpreting the structural fields. Table S1 lists the grids and evaluation metrics; Fig. S1 shows median values across five representative primary sections (MEL01, MEL03, CSCC01, CSCC03 and CSCC05). The default setting is marked on each heatmap. The qualitative pattern is not confined to a narrow parameter ridge: high crosslink contribution persists over the default β range, while dissociation potential changes smoothly with λ, ξ₀ and β. These parameters remain qualitative effective parameters rather than calibrated physical constants in tissue.

**Table S1** Parameter sensitivity grids and evaluation metrics. Crosslink contribution is the percentage of model-derived edge fraction attributed to crosslinking; dissociation potential is the molecular barrier output before rank standardisation. Heatmap values are medians across the five representative sections.

| Parameter | Symbol | Grid | Default | Assessment |
|---|---:|---:|---:|---|
| ECM decay | λ | 0.5, 1, 2, 3, 5, 8, 12, 20 | 3 | Crosslink contribution; log10 potential |
| Mesh contraction | β | 0.25, 0.5, 1, 1.5, 2, 3, 4.5, 6, 9, 12 | 3 | Crosslink contribution; log10 potential |
| Baseline mesh size (nm) | ξ₀ | 8, 12, 16, 20, 30, 40, 60, 80 | 20 | Crosslink contribution; log10 potential |

![**Fig. S1** Parameter sensitivity heatmaps. **a** Crosslink contribution over λ and β; **b** log10 dissociation potential over λ and β; **c** crosslink contribution over ξ₀ and β; **d** log10 dissociation potential over ξ₀ and β. Black rectangles mark the default setting](<<FigS1_parameter_sensitivity>>){width=100%}

# S2 Section admission and quality control

Seven admission criteria were checked for every section: C1, a count matrix; C2, spot coordinates; C3, a paired H&E image; C4, a minimum number of spots; C5, a minimum median UMI count per spot; C6, at least 20 spots with endothelial marker expression; C7, known treatment status and tumour site. Default thresholds were 1,000 spots (Visium), 300 spots (first-generation ST) and a median of 1,500 UMIs per spot. Cohort-level exceptions were written into the configuration with their date and reason before the affected analyses: for the GSE144239 Visium sections (2020 chemistry, small sections) 500 spots and 500 UMIs; for the first-generation ST sections 300 spots and 300 UMIs; and for the 10x Genomics external sections 500 spots and 500 UMIs. The four melanoma sections were retained under a documented cohort-level waiver: their treatment status is not annotated (C7), MEL02 has 801 spots after quality control (C4) and MEL04 has a median of 1,488 UMIs per spot (C5). The third first-generation section of patient P9 (CSCC13; median 289.5 UMIs per spot) failed the predefined threshold and was excluded before any analysis. Table S2 lists the criteria that each analysed section did not meet; no section was admitted silently.

**Table S2 Admission audit. Spots are after quality control; thresholds are those in force for the section's cohort. 10x BC-A1 and BC-A2, 10x Genomics datasets V1_Breast_Cancer_Block_A_Section_1 and V1_Breast_Cancer_Block_A_Section_2; 10x LN, V1_Human_Lymph_Node

<<TABLE_S1>>

# S3 Gene sets

**Table S3 Signature gene sets. Scores were computed with Scanpy's `score_genes` (50 control genes, fixed seed) and rank-normalised within sections. ECM and CAF share five genes (COL1A1, COL1A2, COL3A1, DCN, LUM), which is one reason why their scores are strongly correlated within sections (main text, Fig. 5d)

<<TABLE_S2>>

The ECM set is a curated subset of the core matrisome, not the full core matrisome. Markers labelled "canonical" are widely used lineage markers without a single original source. The ligand score is an expression proxy for PD-L1 and PD-L2 and does not represent PD-1 or any drug's binding sites. Scores are relative abundances within a section; they are not cell proportions, protein concentrations or matrix densities.

# S4 Parameters

**Table S4 Model and analysis parameters. "Fixed a priori" parameters were chosen before the analyses reported here and were not fitted to any outcome or response label

<<TABLE_S3>>

# S5 Simulation benchmark

**Lattice and tissue.** Spots lay on a hexagonal lattice with 100 μm pitch inside a disc of radius 1.5 mm ({{syn_spots}} spots); edges joined spots closer than 150 μm. The tumour nest was a disc of radius 500 μm whose central 300 μm formed the sinks; 14 vessel sources were placed at random in the outer stroma. Matrix-rich spots had true ECM 0.90 and CAF 0.85; the background was 0.12 plus a smooth random field (random Fourier features, length scale 250 μm, amplitude 0.06). The closed capsule occupied the annulus 540–700 μm from the centre; gapped capsules removed a contiguous arc of 5–40% of the circumference starting at a random angle and re-placed the removed spots uniformly in the stroma beyond the capsule; the distant capsule occupied the annulus 840–1000 μm; the band placed the same number of spots uniformly within 530–850 μm; patches grew four compact clusters around random stromal seeds; scatter placed the spots uniformly beyond 620 μm.

**Observation model.** Observed scores were the rank-normalised sum of the true value, a smooth random field with amplitude half the noise level and independent Gaussian noise with standard deviation 0.05, 0.15 or 0.30.

**Ground truth.** 3,000 agents per tissue started at random vessel spots. At each step every agent picked one neighbour uniformly and moved there with probability $c_{uv}$ computed from the true scores (Eq. 1); otherwise it stayed. Agents reaching any core spot were absorbed. Access was the fraction absorbed within 1,500 steps.

**Summaries.** SPARTA's $B_{\mathrm{cell}}$ used the observed scores, vessel sources and core sinks. The SPARTA field summary was the median shortest-path cost over core spots. Neighbourhood enrichment counted tumour–matrix-rich edges, labelling as matrix-rich the top $K$ spots by observed ECM score ($K$, the planted number), and standardised the count against 200 label permutations restricted to non-tumour spots. The BANKSY-style domain summary clustered spots with k-means (three clusters) on own and neighbour-averaged features (observed ECM, CAF and a tumour indicator; neighbour weight 0.3) and reported the fraction of tumour-boundary spots within two edges of the most matrix-rich cluster. Ripley's $L$ used the same $K$ labels at 200 μm. The geometry-only summary, added in this version, was the median distance from the 14 vessel sources to the nest boundary; it uses no score at all. Adding it did not change any other result, because it draws no random numbers.

**Results by noise level.** The AUC of $B_{\mathrm{cell}}$ for separating closed capsules from 5%-gap capsules was {{syn_bcell_auc_gap05_by_noise}} at noise 0.05, 0.15 and 0.30, compared with {{syn_nhood_auc_gap05_by_noise}} for neighbourhood enrichment. Per-tissue values are in the archived file `synthetic_benchmark_tissues.csv`.

**Selection-matched gap experiment on planted capsules.** Under the selection-matched design used for the real sections (best of eight contiguous arcs against best of eight scattered sets, 20% of cut nodes), closed capsules gave a median ratio of {{syn_s2_median}} ({{syn_s2_nsig}}/{{syn_s2_n}} replicates with p < 0.05). {{s2thick_sentence}} Because the matched design gives ratios near 1 even when a continuous capsule is known to exist, the S2 experiment does not test barrier continuity; we report it only as a description of how redundant a barrier band is.

# S6 Calibration of the surrogate test

For each of the {{cal_n_graphs}} real graphs we fixed the observed $B_{\mathrm{mAb}}$ field and the vessel-distance covariate and generated fields that are independent of them by construction. Four scenarios were used: (i) Gaussian fields with the expected graph power spectrum of the observed $B_{\mathrm{mAb}}$ field ({{cal_n_sim}} simulations per graph); (ii) the same fields with the observed marginal distribution of $B_{\mathrm{mAb}}$ imposed by rank matching ({{cal_n_sim}} per graph); (iii) strongly skewed fields, $\exp(2\tilde y)$ of a standardised Gaussian field ({{cal_n_sim_secondary}} per graph); and (iv) heat-kernel Gaussian fields with spectrum $\exp(-\tau\lambda)$, $\tau\in\{2,10,50\}$ ({{cal_n_sim_secondary}} per graph and $\tau$). Each simulated field was tested with {{cal_n_sur}} surrogates drawn from the spectrum of its raw values (the test used for the primary analysis) and with {{cal_n_sur}} surrogates drawn from the spectrum of its normal scores, and with the point-level test.

**Table S5 Pooled false-positive rates at α = 0.05

| Scenario | Raw-spectrum surrogates | Normal-score surrogates | Point-level test |
|:--|--:|--:|--:|
| Spectrum-matched Gaussian | {{cal_pooled_spec}} | {{cal_gauss_nscore}} | {{cal_pooled_naive}} |
| Observed marginal imposed | {{cal_marg_spec}} | {{cal_marg_nscore}} | {{cal_marg_naive}} |
| Strongly skewed | {{cal_skew_spec}} | {{cal_skew_nscore}} | {{cal_skew_naive}} |
| Heat kernel, τ = 2, 10, 50 | {{cal_smooth_spec_min}}–{{cal_smooth_spec_max}} | {{cal_smooth_nscore_min}}–{{cal_smooth_nscore_max}} | {{cal_smooth_naive_min}}–{{cal_smooth_naive_max}} |

The raw-spectrum test is calibrated when the field has a distribution like that of the observed $B_{\mathrm{mAb}}$ fields (median within-section skewness {{bmab_skew_median}}), but it becomes anti-conservative for strongly skewed fields, because sign randomisation of a skewed field's spectrum does not reproduce the spatial structure of its ranks. Drawing surrogates from the spectrum of the normal scores kept the false-positive rate close to nominal in every scenario ({{cal_nscore_min}}–{{cal_nscore_max}}), at the cost of a slightly higher rate than the raw-spectrum test for near-Gaussian fields. Re-testing the primary association with normal-score surrogates gave {{ns_n_q05}}/19 sections with BH q < 0.05 (raw-spectrum test, {{assoc_n_q05}}/19; the difference is {{ns_up_sections}}; Fig. S6c). We report the raw-spectrum test in the main text because it was the primary analysis specified before this calibration study, and we recommend the normal-score variant whenever a field is strongly skewed.

![**Fig. S3** False-positive rate of the surrogate tests on each real graph in the four calibration scenarios. Filled circles, raw-spectrum surrogates; open squares, normal-score surrogates; open triangles, point-level test; the grey band is the 95% binomial range around 0.05 for the scenario's number of simulations](<<FigS3_calibration>>){width=100%}

# S7 Structural nulls and per-section results

The geometry and construction nulls recompute $B_{\mathrm{mAb}}$ from input surrogates while keeping the observed $B_{\mathrm{cell}}$ field, the graph, the vessels and the covariate unchanged. Input surrogates draw independent signs for the graph spectrum of each input score and are rank-normalised to $[0,1]$; a neutral-filled input stays constant. Both nulls used {{dec_n_null}} draws per section. The ablation used $b=0$, $c=12$ in Eq. (1) and $\lambda=0$ in Eq. (2) and was tested with 500 raw-spectrum surrogates.

**Table S6 Per-section results. ρ, partial Spearman correlation; null SD, standard deviation of the 500-draw surrogate null; q, BH-adjusted p-value with raw-spectrum or normal-score surrogates (primary and external sections corrected as separate families); geometry and construction nulls, mean association under each null; p vs constr., one-sided p of the observed association against the construction null; S2 ratio, selection-matched gap experiment

<<TABLE_S5>>

# S8 Graph coverage and disconnection

Sources or sinks outside the largest connected component cannot contribute to the section-level minimum cut, which is computed on that component. Their per-spot field values are finite replacements and are flagged. Table S6 lists these quantities for every section; restricting the association to the largest component changed the median from {{assoc_median}} to {{pll_median}} without changing any sign (main text, Table 3).

**Table S7 Graph coverage. Comp., connected components; LCC %, share of spots in the largest component; Src out and Sink out, sources and sinks outside it; Stranded %, their share of all sources and sinks; Filled inputs, scores replaced by the neutral value

<<TABLE_S4>>

# S9 Additional figures

Figs. S1 and S3–S6 follow; Fig. S3 accompanies Section S6. Map colours follow the main text: grey, ECM score; blue, $B_{\mathrm{cell}}$; orange, $B_{\mathrm{mAb}}$.

![**Fig. S2a** Input and output maps for the primary sections CSCC01–CSCC12. For each section: ECM score (grey), $B_{\mathrm{cell}}$ field with minimum-cut edges (blue; black bars) and $B_{\mathrm{mAb}}$ field (orange), shown as within-section ranks; ρ is the vessel-distance-adjusted field association](<<FigS2a_maps>>){width=100%}

![**Fig. S2b** As Fig. S2a for CSCC14–CSCC16, the four melanoma sections and the three external sections](<<FigS2b_maps>>){width=100%}

![**Fig. S4** Comparison of minimum-cut edges with boundaries between BANKSY-style spatial domains. **a**, **b** Enrichment of cut edges among domain-boundary edges and the share of boundary edges that are cut edges as the number of domains varies (λ = 0.3). **c**, **d** The same as the neighbourhood weight λ varies (eight domains). Grey lines, sections (marker shape, platform); black, medians. This is not the official BANKSY software](<<FigS4_domain_scans>>){width=100%}

![**Fig. S5** Field association as a function of the graph radius; lines, primary sections (marker shape, platform). At 100 μm (Visium) and 200 μm (first-generation ST) the radius equals the spot pitch, so only part of the nearest neighbours are connected and the graph fragments; these configurations are shown for completeness](<<FigS5_radius>>){width=60%}

![**Fig. S6** Reproducibility. **a** Field association recomputed from the archived node tables on a different operating system against the values stored by the original analysis. **b** The same for the standard deviation of the surrogate null. **c** One-sided p-values of the primary association with raw-spectrum and normal-score surrogates; lines mark 0.05](<<FigS6_reproduction>>){width=100%}

![**Fig. S7** Wall time per section for the SPARTA core operators, the BANKSY-style clustering and Squidpy neighbourhood enrichment on one laptop CPU (Intel Core, 20 logical cores, 15.6 GB RAM, Windows 11). Timings are machine- and implementation-specific](<<FigS7_runtime>>){width=60%}

# S10 Exploratory analyses

**Descriptive alignment with immune and proliferation signatures.** Within sections, the $B_{\mathrm{cell}}$ field was weakly and inconsistently associated with the T/NK score (median partial ρ = {{bio_bc_tnk_median}}; significantly negative in {{bio_bc_tnk_nsig}}/19 sections), and the $B_{\mathrm{mAb}}$ field with the proliferation score (median {{bio_bm_prolif_median}}; significantly positive in {{bio_bm_prolif_nsig}}/19). Because immune-entry sources are themselves selected by their T/NK neighbourhood, these associations are partly built into the source definition and are not used as validation.

# S11 Melanoma replication cohort

**Data and plan.** Eight first-generation ST sections of lymph-node metastases from four patients with stage III cutaneous melanoma (two consecutive sections per metastasis; 1,007-spot arrays, 100 μm spots at 200 μm spacing) were obtained from the original study [28 in the main text]. The plan (`docs/replication_melanoma_protocol.md` in Online Resource 3) was written before any graph, barrier or association was computed for these sections. Earlier ingestion had stored gene names as "SYMBOL ENSG…" strings, so no endothelial marker was matched and every section failed criterion C6; nothing downstream had been run. The plan fixed the pipeline (first-generation ST settings of the primary cohort: 200 μm pitch, 300 μm graph radius, unchanged quantile thresholds and operator parameters), the analyses (those of the primary cohort, as a separate family) and the replication criterion (pooled replication-cohort association with a 95% CI excluding zero).

**Admission.** The configuration entry `mel_thrane2018` lowers the minimum spot number to 200 because the array has about half as many spots as the 1,933-spot arrays of the primary first-generation sections; criteria C3 (no paired H&E image in the public release) and C7 (treatment before sampling not reported) are waived at cohort level. All eight sections met the remaining criteria and were analysed ({{rep_spots_min}}–{{rep_spots_max}} spots after quality control).

**Scanpy-free processing.** Quality control, normalisation and signature scoring were run with `sparta/lite.py`, which re-implements `filter_cells`, `filter_genes`, `normalize_total`, `log1p` and `score_genes` (25 expression bins, 50 control genes per bin, legacy global seed) with NumPy and pandas. Re-scoring the primary first-generation sections CSCC05, CSCC08 and CSCC12 from their raw tables reproduced the archived raw signature scores to within {{lite_max_dscore}} (float32 rounding) and their within-section ranks, which are the model inputs, exactly (`results/validation/lite_scoring_equivalence.json`). The ligand score was replaced by the neutral value in {{rep_ag_missing_sections}} because fewer than two of its genes passed filtering.

**Results.** Table S7 lists every section. The surrogate test was calibrated on these graphs (false-positive rate {{cal_rep_spec}}, normal-score variant {{cal_rep_nscore}}, point-level test {{cal_rep_naive}}; `null_calibration_replication.json`). The association was positive in all sections, significant after BH correction in {{rep_n_q05}}/8 and positive in all four patients; the pooled estimate was {{rep_pl_mean}} (95% CI {{rep_pl_lo}}–{{rep_pl_hi}}; leave-one-patient-out estimates ≥ {{rep_lopo_min}}), which meets the pre-specified criterion. The construction null reproduced a median {{rep_share_constr_pct}}% of the observed association, the excess pooled to {{rep_ex_mean}} (95% CI {{rep_ex_lo}}–{{rep_ex_hi}}), and the ECM-ablated association stayed positive in all sections. Spots with low $B_{\mathrm{cell}}$ and high $B_{\mathrm{mAb}}$ were rarer than the 6.25% expected under independence in {{rep_n_discordant_below_chance}}/8 sections.

**Table S8 Replication cohort, per section (LN P1–P4, patients; r1 and r2, consecutive sections). Src and Sink, numbers of sources and sinks; ρ, vessel-distance-adjusted partial Spearman correlation; null SD, standard deviation of the 500-draw surrogate null; q, BH-adjusted p-values within the cohort with raw-spectrum or normal-score surrogates; constr. null, mean association under the construction null (200 draws); ρ ablated, association with the shared ECM term removed; filled inputs, scores replaced by the neutral value

<<TABLE_S7>>

# S12 Validation against measured CD8^+^ T-cell positions (CODEX)

**Protocol.** The protocol (`docs/codex_validation_protocol.md` in Online Resource 3; its SHA-256 is stored with the results) was written before any barrier or outcome was computed; only cell-type counts had been tabulated. Three amendments are recorded in it: before any computation, the sink rule was corrected from a μm-weighted to the hop distance used for the sections; when the results were first summarised, the pre-specified comparator "global matrix mean" turned out to be constant (0.5) in every core because both inputs are within-core ranks, so it is reported as degenerate and the pre-specified partial correlation effectively adjusts for tumour fraction, stromal fraction and peritumoural matrix density; and before the final run, every minimum cut was computed with integer-scaled capacities (Section S14), which leaves every barrier value unchanged and affects only the source and sink sides used by the blockade-line analysis. Analyses added after the comparators had been seen are labelled post hoc: the partial association of $\log B_{\mathrm{rel}}$ given tumour–stroma intermixing; four alternative weightings of Eq. (1); two geometry-only comparators, the median distance from vascular cells to the nearest tumour cell and the median distance from sinks to the nearest non-tumour cell (tumour-core depth); and the partial association of $\log B_{\mathrm{rel}}$ given each simple summary, given both distances (with and without composition) and given matrix density, both distances and intermixing together. Core depth is partly tied to the outcome by construction, because the infiltration ratio is measured in the region assigned to the sinks whose depth it describes.

**Inclusion.** Of {{cx_n_cores_total}} cores, {{cx_n_excluded}} were excluded, all for having fewer than ten sinks, 25 tumour cells or five vascular cells in their scaffold graph (counts by reason in `codex_validation.json`); every patient kept at least one core.

**Results.** Table S8 gives every comparator, the secondary analyses and the sensitivity variants. With all T cells instead of CD8^+^ T cells as the outcome, or with all tumour cells as sinks, the association kept its direction. Measured over the whole tumour region rather than its core, the patient-level correlation was {{cx_tumour_region_rho}}. The blockade-line analysis compared CD8^+^ density beyond the cut with that before it: the SPARTA cut placed a median {{cx_s2_area_sparta_pct}}% of the tissue area on its sink side, against {{cx_s2_area_open_pct}}% for the geometry-only cut, so its sink side contained more stroma, which is where most T cells were; this is one reason why the comparison was negative.

**Table S9 Validation against measured CD8$^+$ T-cell positions: Spearman correlations with the CD8$^+$ infiltration ratio of the tumour core (negative = barrier-like), at core level with patient-cluster bootstrap 95% CIs ({{cx_n_cores}} cores) and at patient level with one-sided permutation p-values ({{cx_n_patients}} patients). The global matrix mean is constant by construction (Section S12)

<<TABLE_S8>>

# S13 Intervention maps

For every spot of the 19 primary, three external and eight replication sections, the ECM, CAF and crosslinking scores of that spot alone were set to the section's 5th percentile, and both barriers were recomputed ({{iv_n_ablations}} single-spot ablations; `run_44_intervention_ranking.py`). Single-spot effects on $B_{\mathrm{cell}}$ were zero for every spot that is not an endpoint of a minimum-cut edge, as max-flow/min-cut duality implies when the cut is unique, so the fraction of spots with an effect equals the fraction of spots on the cut (median {{iv_cut_frac_pct}}% in the primary sections). Joint ablations used the top 1%, 2%, 5% and 10% of each map and budgets of 1, 5, 20 and 50 spots (`run_44`), and a grid of budgets for the half-reduction point with random and density-based selections repeated 20 times (`run_49_intervention_targeting.py`). Table S9 lists every section and Fig. S8 summarises the maps. Small negative single-spot effects on $B_{\mathrm{mAb}}$ occur when opening matrix next to absorbing spots lets more antibody be absorbed there; they are kept as computed.

**Concentration.** In the primary sections, the single most effective spot removed a median {{iv_single_cell_median}}% of a full breach of the cut (up to {{iv_single_cell_max}}% in {{iv_single_cell_max_section}}) but only {{iv_single_mab_median}}% of the antibody reference (at most {{iv_single_mab_max}}%). The top 1% of spots, ablated together, removed a median {{iv_top1_cell_median}}% (range {{iv_top1_cell_min}}–{{iv_top1_cell_max}}%) of the cellular and {{iv_top1_mab_median}}% ({{iv_top1_mab_min}}–{{iv_top1_mab_max}}%) of the antibody reference; the top 5% removed {{iv_top5_cell_median}}% and {{iv_top5_mab_median}}%, and the top 10% {{iv_top10_mab_median}}% of the antibody reference. In the replication sections the top 1% removed {{iv_rep_top1_cell}}% and {{iv_rep_top1_mab}}%. The Jaccard index between the two maps' top 5% of spots had a median of {{iv_jacc_median_pct}}% (range {{iv_jacc_min_pct}}–{{iv_jacc_max_pct}}%) in the primary sections and cohort medians of {{iv_jacc_cohort_min_pct}}–{{iv_jacc_cohort_max_pct}}%.

**Compartments of the high-impact spots.** Spots were labelled tumour (malignant score at or above the 70th percentile, the pipeline's malignant compartment), stromal or immune (the larger of the mean ECM and CAF score and the mean T/NK, myeloid and B/plasma score; `run_53_intervention_domain_enrichment.py`). For each operator and budget, the share of the top spots in each compartment was compared with the share of all spots, pooled over the 29 tumour sections (the lymph node was excluded), with an exact null obtained by convolving the per-section hypergeometric distributions (random placement within sections; this ignores spatial autocorrelation and is descriptive), and across patients with an exact sign-flip test of the mean within-patient difference in shares (12 patients), which is the reported test. Immune-dominated spots were under-represented among the top spots of both operators (top 5%: {{de_c5_imm_top}}% for $B_{\mathrm{cell}}$ and {{de_m5_imm_top}}% for $B_{\mathrm{mAb}}$, against {{de_c5_imm_all}}% of all spots; lower in {{de_c5_imm_nneg}} and {{de_m5_imm_nneg}} of 12 patients; sign-flip p = {{de_c5_imm_p}} and {{de_m5_imm_p}}; top 1%: p = {{de_c1_imm_p}} and {{de_m1_imm_p}}). The cellular top 5% was modestly enriched in stromal spots ({{de_c5_str_top}}% vs {{de_c5_str_all}}%; {{de_c5_str_npos}}/12 patients; p = {{de_c5_str_p}}) and the antibody top 5% in tumour spots ({{de_m5_tum_top}}% vs {{de_m5_tum_all}}%; {{de_m5_tum_npos}}/12 patients; p = {{de_m5_tum_p}}). Because $B_{\mathrm{cell}}$ effects are confined to the cut, its top spots were also compared with the composition of the cut itself ({{de_c5_tum_cut}}% tumour, {{de_c5_str_cut}}% stromal, {{de_c5_imm_cut}}% immune): relative to the cut, the most effective cut spots were more often tumour spots ({{de_c5_tum_top}}%). Compartments are expression-defined and coarse; these shares describe where the model places its bottlenecks, not histology.

![**Fig. S8** In-model intervention maps. **a**, **b** Reduction of $B_{\mathrm{cell}}$ (**a**) and of the mean $B_{\mathrm{mAb}}$ (**b**) when the matrix inputs of a single spot are set to the section's 5th percentile, for every spot of CSCC04 (darker, larger reduction; colour scale clipped at the 97th percentile); circles, top 5% of spots. **c** Joint ablation of the top-ranked spots of each map, as a fraction of a full breach of the minimum cut ($B_{\mathrm{cell}}$) or of ablation of every spot ($B_{\mathrm{mAb}}$); lines, medians over the 19 primary sections; bands, interquartile ranges. **d** Reduction achieved by 20 spots chosen by the map, at random on the cut, by local matrix density, at random or by the other modality's map; dots, primary sections; bars, medians. **e** Compartments of all spots and of the top 1% and 5% of each map, pooled over the 29 tumour sections; p, patient-level sign-flip test for the immune share. All quantities are model-defined](<<FigS8_intervention>>){width=100%}

**Table S10 Intervention maps, all sections. Cut %, share of spots that are endpoints of minimum-cut edges; k50, smallest share of map-ranked spots whose joint ablation recovers half of the reference reduction (grid of budgets); Top 1%, percentage of the reference reduction recovered by jointly ablating the top 1% of each map; 20 spots, percentage recovered by 20 spots chosen by the map, at random on the cut (mean of 20 draws) or by local matrix density; Jaccard, overlap of the two maps' top 5% of spots. In the small replication sections, 20 spots are 5–10% of all spots

<<TABLE_S9>>

# S14 Exact minimum cuts

Up to release 2.1.0, the minimum cut was read off the residual graph of a floating-point preflow-push computation. The flow value, and therefore $B_{\mathrm{cell}}$, was correct (relative difference from an exact integer computation ≤ {{mc_flow_diff_max}}), but the edge set was not always a minimum cut: its total capacity differed from the maximum flow in {{mc_n_old_not_min}} of the {{mc_n_sections}} section graphs analysed here, by up to {{mc_old_gap_max_pct}}%. All cuts in this version are computed with capacities rounded to integer multiples of $2^{-40}$, for which the cut capacity equals the maximum flow (largest deviation {{mc_new_gap_max}}); the median Jaccard index between old and new cut-edge sets was {{mc_median_jaccard}}. Everything that uses the cut geometry (the gap experiments, the domain-boundary comparison, the intervention maps, the CODEX blockade-line analysis and the cut drawings) was recomputed (`run_50_exact_mincut_refresh.py`, `mincut_exactness.json`). The BANKSY-style domains had been computed with Scanpy; because Scanpy was not available where the correction was made, `run_51_domain_comparison_exact_cut.py` re-implements the same preprocessing (Seurat-flavour highly variable genes, scaling with a maximum of 10, ARPACK principal components, k-means with the same seed) with NumPy, pandas and scikit-learn and reads the stored .h5ad files through the HDF5 library. Given the earlier cuts, it reproduced every archived domain-boundary statistic exactly in {{dom_repro_n}} of {{dom_repro_total}} sections over all four neighbourhood weights (`domain_comparison_exact_cut.json`); with the exact cuts, the median enrichment of cut edges on domain boundaries changed from {{dom_enrich_median_old}} to {{dom_enrich_median}} and the median precision did not change. The single-spot intervention effects depend only on flow values and were not recomputed; their cut masks, reference reductions and joint curves were.


# S15 Simple spatial summaries

Three summaries that need no transport model were computed for the 30 sections (`run_52_simple_baselines_spearman.py`, `simple_baselines_spearman.json`). At spot level: stromal density, the mean of the ECM and CAF scores over a spot and its graph neighbours; distance to the tumour, the signed Euclidean distance to the boundary of the malignant compartment (negative inside it); and neighbourhood (niche) enrichment, the binomial z-score of the number of stromal neighbours given the section's stromal share. At section level: the same density over non-malignant spots adjacent to the malignant compartment, the median distance from vessel spots to the nearest malignant spot, and the z-score of tumour–stroma contacts against 1,000 label permutations (as Squidpy's neighbourhood enrichment). Operator summaries were the within-section fields and, per section, $\log B_{\mathrm{rel}}=\log(F_{\mathrm{open}}/F)$ (maximum flow with every capacity at $\sigma(a)$ over the maximum flow with the matrix capacities of Eq. 1) and the median of each field over tumour-core spots. Within sections, the $B_{\mathrm{cell}}$ field correlated with local stromal density (median ρ = {{bl_cell_dens_med}}, interquartile range {{bl_cell_dens_q25}}–{{bl_cell_dens_q75}}, maximum {{bl_cell_dens_max}}), less with niche enrichment ({{bl_cell_niche_med}}) and not with distance to the tumour ({{bl_cell_dist_med}}); the $B_{\mathrm{mAb}}$ field correlated weakly and negatively with all three ({{bl_mab_dens_med}}, {{bl_mab_dist_med}} and {{bl_mab_niche_med}}). Across the 29 tumour sections, $\log B_{\mathrm{rel}}$ was not correlated with any section-level summary (ρ = {{bl_sec_brel_dens_rho}}, {{bl_sec_brel_dist_rho}} and {{bl_sec_brel_niche_rho}}); the core $B_{\mathrm{mAb}}$ correlated with peritumoural stromal density ({{bl_sec_mabcore_dens_rho}}, p = {{bl_sec_mabcore_dens_p}}) and negatively with tumour–stroma enrichment ({{bl_sec_mabcore_niche_rho}}, p = {{bl_sec_mabcore_niche_p}}). Table S10 lists every section.

**Table S11 Operators and simple spatial summaries, all sections. Spot level: within-section Spearman correlation of the $B_{\mathrm{cell}}$ and $B_{\mathrm{mAb}}$ fields with local stromal density (dens.), signed distance to the tumour (dist.) and niche enrichment (niche). Section level: $\log B_{\mathrm{rel}}$, peritumoural stromal density, median vessel-to-tumour distance (μm) and tumour–stroma contact z-score

<<TABLE_S10>>

# S16 2026 extension cohort

This section documents the 50 public sections added as the 2026 extension cohort. QC was performed at the native Visium spot or Slide-seqV2 bead level. Slide-seqV2 counts were subsequently summed into fixed 50-µm grid bins for structural analysis; the h5ad files retain the QC-passed bead-level counts.

**Table S12 Extension cohort metadata. Patient relationship was not provided by GSE289745; those sections are marked as relationship unavailable and were treated as conservative section-level sampling units.

| Project | Section | Patient identifier | Disease / stratum | Platform | Site |
|---|---|---|---|---|---|
| GSE200278 | ECM01_rep1 | MPM01 | melanoma / metastatic | Slide-seqV2 | subcutaneous tissue |
| GSE200278 | ECM01_rep2 | MPM01 | melanoma / metastatic | Slide-seqV2 | subcutaneous tissue |
| GSE200278 | ECM06 | MPM06 | melanoma / metastatic | Slide-seqV2 | subcutaneous tissue |
| GSE200278 | ECM08 | MPM08 | melanoma / metastatic | Slide-seqV2 | subcutaneous tissue |
| GSE200278 | ECM10 | MPM10 | melanoma / metastatic | Slide-seqV2 | subcutaneous tissue |
| GSE200278 | MBM05_rep1 | MBM05 | melanoma / metastatic | Slide-seqV2 | brain |
| GSE200278 | MBM05_rep2 | MBM05 | melanoma / metastatic | Slide-seqV2 | brain |
| GSE200278 | MBM05_rep3 | MBM05 | melanoma / metastatic | Slide-seqV2 | brain |
| GSE200278 | MBM06 | MBM06 | melanoma / metastatic | Slide-seqV2 | brain |
| GSE200278 | MBM07 | MBM07 | melanoma / metastatic | Slide-seqV2 | brain |
| GSE200278 | MBM08 | MBM08 | melanoma / metastatic | Slide-seqV2 | brain |
| GSE200278 | MBM11_rep1 | MBM11 | melanoma / metastatic | Slide-seqV2 | brain |
| GSE200278 | MBM11_rep2 | MBM11 | melanoma / metastatic | Slide-seqV2 | brain |
| GSE200278 | MBM11_rep3 | MBM11 | melanoma / metastatic | Slide-seqV2 | brain |
| GSE200278 | MBM13 | MBM13 | melanoma / metastatic | Slide-seqV2 | brain |
| GSE200278 | MBM18 | MBM18 | melanoma / metastatic | Slide-seqV2 | brain |
| GSE289745 | CSCC289_S1 | relationship unavailable | cscc / cutaneous cSCC (progression not otherwise specified) | Visium | skin |
| GSE289745 | CSCC289_S10 | relationship unavailable | cscc / cutaneous cSCC (progression not otherwise specified) | Visium | skin |
| GSE289745 | CSCC289_S11 | relationship unavailable | cscc / cutaneous cSCC (progression not otherwise specified) | Visium | skin |
| GSE289745 | CSCC289_S15 | relationship unavailable | cscc / cutaneous cSCC (progression not otherwise specified) | Visium | skin |
| GSE289745 | CSCC289_S3 | relationship unavailable | cscc / cutaneous cSCC (progression not otherwise specified) | Visium | skin |
| GSE289745 | CSCC289_S4 | relationship unavailable | cscc / cutaneous cSCC (progression not otherwise specified) | Visium | skin |
| GSE289745 | CSCC289_S5 | relationship unavailable | cscc / cutaneous cSCC (progression not otherwise specified) | Visium | skin |
| GSE289745 | CSCC289_S6 | relationship unavailable | cscc / cutaneous cSCC (progression not otherwise specified) | Visium | skin |
| GSE289745 | CSCC289_S7 | relationship unavailable | cscc / cutaneous cSCC (progression not otherwise specified) | Visium | skin |
| GSE289745 | CSCC289_S8 | relationship unavailable | cscc / cutaneous cSCC (progression not otherwise specified) | Visium | skin |
| GSE289745 | CSCC289_S9 | relationship unavailable | cscc / cutaneous cSCC (progression not otherwise specified) | Visium | skin |
| GSE300445 | MEL300_0019 | MEL300_0019 | melanoma / primary | Visium | skin |
| GSE300445 | MEL300_0022 | MEL300_0022 | melanoma / primary | Visium | skin |
| GSE300445 | MEL300_0113 | MEL300_0113 | melanoma / primary | Visium | skin |
| GSE300445 | MEL300_0133 | MEL300_0133 | melanoma / primary | Visium | skin |
| GSE316760 | MEL316_mel2 | MEL316_mel2 | melanoma / primary | Visium | skin |
| GSE316760 | MEL316_mel3 | MEL316_mel3 | melanoma / primary | Visium | skin |
| GSE320041 | MEL320_WU1340 | MEL320_WU1340 | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE320041 | MEL320_WU1373 | MEL320_WU1373 | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE320041 | MEL320_WU1384_1 | MEL320_WU1384 | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE320041 | MEL320_WU1384_2 | MEL320_WU1384 | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE320041 | MEL320_WU1609 | MEL320_WU1609 | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE320041 | MEL320_WU2130_1 | MEL320_WU2130 | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE320041 | MEL320_WU2130_2 | MEL320_WU2130 | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE320041 | MEL320_WU2415 | MEL320_WU2415 | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE320041 | MEL320_WU3049 | MEL320_WU3049 | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE320041 | MEL320_WU3244 | MEL320_WU3244 | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE320041 | MEL320_YUADD | MEL320_YUADD | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE320041 | MEL320_YUALT | MEL320_YUALT | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE320041 | MEL320_YUBOISE | MEL320_YUBOISE | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE320041 | MEL320_YUMAZO | MEL320_YUMAZO | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE320041 | MEL320_YUSTE | MEL320_YUSTE | melanoma / metastatic | Visium | metastatic melanoma (specific site not annotated) |
| GSE321832 | CSCC321_cut1 | CSCC321_cut1 | cscc / primary | Visium | skin |
| GSE321832 | CSCC321_cut2 | CSCC321_cut2 | cscc / primary | Visium | skin |

**Table S13 Extension cohort QC. Raw and retained locations are Visium spots or Slide-seqV2 beads before the separate Slide-seqV2 50-µm grid aggregation. MT, mitochondrial UMI fraction; NA indicates that the supplied matrix did not contain mitochondrial genes.

| Section | Platform | Raw locations | Retained locations | Retained genes | Median retained UMI | Median MT fraction | QC |
|---|---|---:|---:|---:|---:|---:|---|
| ECM01_rep1 | Slide-seqV2 | 27,325 | 27,325 | 20,782 | 556 | 0.022 | True |
| ECM01_rep2 | Slide-seqV2 | 24,150 | 24,150 | 21,113 | 581 | 0.019 | True |
| ECM06 | Slide-seqV2 | 25,461 | 25,461 | 20,624 | 284 | 0.000 | True |
| ECM08 | Slide-seqV2 | 30,015 | 30,015 | 21,744 | 338 | 0.000 | True |
| ECM10 | Slide-seqV2 | 37,363 | 37,345 | 23,548 | 412 | 0.024 | True |
| MBM05_rep1 | Slide-seqV2 | 29,536 | 29,526 | 22,543 | 344 | 0.025 | True |
| MBM05_rep2 | Slide-seqV2 | 32,228 | 32,146 | 23,365 | 358 | 0.035 | True |
| MBM05_rep3 | Slide-seqV2 | 6,002 | 5,999 | 17,206 | 433 | 0.013 | True |
| MBM06 | Slide-seqV2 | 27,072 | 27,026 | 23,996 | 470 | 0.018 | True |
| MBM07 | Slide-seqV2 | 38,462 | 38,460 | 22,699 | 472 | 0.021 | True |
| MBM08 | Slide-seqV2 | 35,054 | 35,054 | 23,197 | 352 | 0.010 | True |
| MBM11_rep1 | Slide-seqV2 | 27,475 | 27,470 | 23,779 | 280 | 0.031 | True |
| MBM11_rep2 | Slide-seqV2 | 39,049 | 39,044 | 25,886 | 426 | 0.028 | True |
| MBM11_rep3 | Slide-seqV2 | 9,366 | 9,366 | 19,816 | 256 | 0.002 | True |
| MBM13 | Slide-seqV2 | 32,657 | 32,656 | 21,243 | 337 | 0.020 | True |
| MBM18 | Slide-seqV2 | 38,423 | 38,418 | 22,502 | 396 | 0.017 | True |
| CSCC289_S1 | Visium | 2,895 | 2,860 | 15,521 | 5937 | NA | True |
| CSCC289_S10 | Visium | 1,274 | 1,146 | 15,889 | 18764 | NA | True |
| CSCC289_S11 | Visium | 1,514 | 1,483 | 15,581 | 3378 | NA | True |
| CSCC289_S15 | Visium | 2,041 | 2,029 | 17,395 | 7495 | NA | True |
| CSCC289_S3 | Visium | 2,887 | 2,858 | 17,550 | 15935 | NA | True |
| CSCC289_S4 | Visium | 3,195 | 3,195 | 15,928 | 15021 | NA | True |
| CSCC289_S5 | Visium | 2,779 | 2,729 | 15,741 | 7218 | NA | True |
| CSCC289_S6 | Visium | 2,063 | 2,059 | 15,011 | 5758 | NA | True |
| CSCC289_S7 | Visium | 2,748 | 2,746 | 16,674 | 27858 | NA | True |
| CSCC289_S8 | Visium | 2,163 | 2,156 | 15,611 | 17300 | NA | True |
| CSCC289_S9 | Visium | 4,271 | 4,261 | 15,825 | 10234 | NA | True |
| MEL300_0019 | Visium | 2,539 | 2,534 | 18,037 | 34358 | 0.025 | True |
| MEL300_0022 | Visium | 2,783 | 2,762 | 18,052 | 27524 | 0.006 | True |
| MEL300_0113 | Visium | 2,145 | 2,132 | 18,052 | 6802 | 0.025 | True |
| MEL300_0133 | Visium | 4,115 | 4,080 | 18,053 | 37024 | 0.024 | True |
| MEL316_mel2 | Visium | 4,084 | 4,082 | 16,587 | 6694 | NA | True |
| MEL316_mel3 | Visium | 4,949 | 4,933 | 17,195 | 9652 | NA | True |
| MEL320_WU1340 | Visium | 4,991 | 4,989 | 18,052 | 28288 | 0.033 | True |
| MEL320_WU1373 | Visium | 4,708 | 4,651 | 18,047 | 37893 | 0.058 | True |
| MEL320_WU1384_1 | Visium | 4,920 | 4,782 | 18,063 | 52328 | 0.021 | True |
| MEL320_WU1384_2 | Visium | 4,919 | 4,571 | 18,057 | 37535 | 0.021 | True |
| MEL320_WU1609 | Visium | 4,990 | 4,989 | 18,046 | 63900 | 0.096 | True |
| MEL320_WU2130_1 | Visium | 4,910 | 4,251 | 17,685 | 3564 | 0.039 | True |
| MEL320_WU2130_2 | Visium | 3,819 | 3,771 | 18,050 | 42248 | 0.035 | True |
| MEL320_WU2415 | Visium | 4,886 | 4,807 | 18,048 | 60766 | 0.053 | True |
| MEL320_WU3049 | Visium | 3,932 | 3,897 | 18,061 | 51182 | 0.041 | True |
| MEL320_WU3244 | Visium | 4,932 | 4,867 | 18,058 | 51462 | 0.019 | True |
| MEL320_YUADD | Visium | 746 | 482 | 16,826 | 10758 | 0.013 | True |
| MEL320_YUALT | Visium | 2,558 | 2,393 | 18,031 | 6460 | 0.035 | True |
| MEL320_YUBOISE | Visium | 476 | 369 | 17,676 | 56178 | 0.048 | True |
| MEL320_YUMAZO | Visium | 2,597 | 2,355 | 18,017 | 23859 | 0.071 | True |
| MEL320_YUSTE | Visium | 2,785 | 1,458 | 18,011 | 15102 | 0.031 | True |
| CSCC321_cut1 | Visium | 4,454 | 4,258 | 17,513 | 4744 | NA | True |
| CSCC321_cut2 | Visium | 4,197 | 4,195 | 17,636 | 13472 | NA | True |
