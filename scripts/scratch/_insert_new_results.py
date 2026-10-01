# -*- coding: utf-8 -*-
"""把三个新实验结果和定位段落插入 manuscript_cbc_draft.md"""
import io

p = r"D:\sparta\docs\manuscript_cbc_draft.md"
s = io.open(p, encoding="utf-8").read()

# ── 1. Intro 末尾加定位划线（在 "Our contributions are:" 之前）──
anchor1 = "Our contributions are:"
pos1 = s.index(anchor1)
intro_line = (
    "We frame this as a methodological contribution on public data, not as a "
    "clinical biomarker: no parameter was set using clinical outcome labels, and "
    "the biological claims are bounded throughout by the patient counts reported. "
    "The deliverable is a well-defined, deterministic, reproducible transport "
    "operator that turns spatial measurements into barrier quantities, plus a "
    "testable prediction about what stromal modification would do to both "
    "delivery modalities simultaneously.\n\n"
)
s = s[:pos1] + intro_line + s[pos1:]

# ── 2. §3.2 末尾加 null crosslink permutation test（在 "### 3.3." 之前）──
anchor2 = "### 3.3. Molecular size separates"
pos2 = s.index(anchor2)
null_paragraph = (
    "**Permutation control for the rank-normalised input.** Because the "
    "crosslinking score is rank-normalised within each section, the fraction of "
    "edges below the IgG radius is partly fixed by construction. To separate "
    "model artefact from data signal, we permuted the crosslinking score 50 "
    "times within each of six representative sections (spanning both tumour "
    "types and both platform generations) and re-ran the variance decomposition "
    "on each permutation. The size-exclusion share of B_mAb variance is "
    "essentially unchanged: median 97.6% under the real crosslinking field "
    "versus 97.6% under permutation (95% null interval 96.8–98.4%); only 2 of "
    "6 sections show a permutation p < 0.05, and in those the effect size is "
    "less than one percentage point. We therefore report the 97.5% figure "
    "explicitly as a property of the operator given β = 3, not as an empirical "
    "measurement of crosslinking architecture — a framing already stated in "
    "Section 2.4 and confirmed here rather than softened after review.\n\n"
)
s = s[:pos2] + null_paragraph + s[pos2:]

# ── 3. §3.5 末尾加 stromal intervention（在 "### 3.6." 之前）──
anchor3 = "### 3.6. The barrier is not a rewrite"
pos3 = s.index(anchor3)
intervention_paragraph = (
    "**In silico stromal co-targeting.** The coupling implies a directly "
    "testable prediction: if the two barriers share a matrix substrate, "
    "reducing matrix density and crosslinking in silico should lower both "
    "barriers at once. We tested this by scaling the ecm, caf and crosslinking "
    "scores by (1 − reduction) at reduction levels of 20%, 30% and 50%, and "
    "recomputing both operators. At a 30% reduction, the section-level "
    "minimum-cut barrier falls by a median 60.8% (cSCC 61.1%, melanoma 58.1%) "
    "and the mean tumour-core B_mAb by 49.5% (cSCC 49.5%, melanoma 49.6%); "
    "both barriers fall in 19 of 19 sections. At 50% reduction the drops are "
    "80.2% and 71.3% respectively. The two modalities respond in the same "
    "direction and with similar magnitude across both tumour types, which is "
    "the quantitative basis for the paper's clinical claim: stromal-directed "
    "intervention (LOX/TGF-β inhibition, anti-fibrotic combinations) is "
    "predicted to improve cellular and macromolecular delivery simultaneously "
    "rather than relieving one at the expense of the other. This is a "
    "model-based prediction, not a measured treatment response; it is stated "
    "as such.\n\n"
)
s = s[:pos3] + intervention_paragraph + s[pos3:]

# ── 4. §3.6 末尾加 naive baseline 对比（在 "## 4. Discussion" 之前）──
anchor4 = "## 4. Discussion"
pos4 = s.index(anchor4)
naive_paragraph = (
    "**Why a graph at all?** A natural reviewer question is whether a simple "
    "node-level matrix score (0.5·(ECM + CAF)) already captures the barrier, "
    "making the min-cut and screened-Poisson machinery unnecessary. We "
    "compared the per-spot graph B_cell field against this naive score on all "
    "19 sections, using partial Spearman correlation with B_mAb after "
    "residualising on vessel distance. The graph field yields a median "
    "partial ρ of +0.217 versus +0.166 for the naive score; the graph field "
    "is stronger in 14 of 19 sections and reaches p < 0.05 in 17 of 19 "
    "versus 15 of 19 for the naive score. The advantage is modest rather than "
    "transformative (median Δρ = +0.048) and is concentrated in sections "
    "where source–sink geometry is non-trivial (e.g. CSCC11/CSCC12, where "
    "Δρ = +0.14/+0.17). We read this honestly: a matrix score captures most "
    "of the barrier signal because the matrix is the dominant substrate, but "
    "the graph operator adds the spatial arrangement — the capacity of the "
    "narrowest cut and the screened-diffusion field — that a node score "
    "cannot represent.\n\n"
)
s = s[:pos4] + naive_paragraph + s[pos4:]

# ── 5. Limitations 加 single-author（在 "Finally, the melanoma cohort lacks" 之前）──
anchor5 = "melanoma cohort lacks treatment annotation"
pos5 = s.index(anchor5)
single_author = (
    "This is a single-author work: the analysis was not independently "
    "checked by a second analyst, and the code is released for community "
    "verification rather than having passed internal peer review. "
)
s = s[:pos5] + single_author + s[pos5:]

# ── 6. Abstract 加一句 stromal intervention 预测（在最后一句之前）──
anchor6 = "full analysis of one section runs in under 0.2 s"
pos6 = s.index(anchor6)
abstract_add = (
    "In silico reduction of matrix density and crosslinking by 30% lowers "
    "both barriers by a median 61% and 50% respectively, in 19 of 19 sections "
    "across both tumour types. "
)
s = s[:pos6] + abstract_add + s[pos6:]

# ── 7. Conclusions 更新（在最后一句之前加数字）──
anchor7 = "matrix-rich tumours the model predicts"
pos7 = s.index(anchor7)
# 这句已经在，不重复加；只在前面加一句 stromal intervention 数字
conc_add = (
    "Quantitatively, a 30% in silico reduction in matrix density and "
    "crosslinking lowers the minimum-cut barrier by 61% and the antibody "
    "barrier by 50% across all 19 sections. "
)
s = s[:pos7] + conc_add + s[pos7:]

io.open(p, "w", encoding="utf-8", newline="\n").write(s)
print("manuscript updated. New length:", len(s), "chars")
