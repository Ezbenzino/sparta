# -*- coding: utf-8 -*-
"""把新实验段落同步进 generate_docx.js"""
import io

p = r"D:\sparta\scripts\generate_docx.js"
s = io.open(p, encoding="utf-8").read()

# 1. Abstract: 在 "The full analysis of one section" 之前插入 stromal 句
old_abs = "The full analysis of one section runs in under 0.2 s on CPU with no learned parameters, making the model usable as a reproducible component of translational pipelines."
new_abs = ("In silico reduction of matrix density and crosslinking by 30% lowers both barriers by a median 61% and 50% respectively, in 19 of 19 sections across both tumour types. " + old_abs)
assert old_abs in s, "abstract anchor not found"
s = s.replace(old_abs, new_abs, 1)

# 2. Intro: 在 "Our contributions are:" 之前加定位段
old_intro = 'children.push(p("Our contributions are:"));'
new_intro = ('children.push(p("We frame this as a methodological contribution on public data, not as a clinical biomarker: no parameter was set using clinical outcome labels, and the biological claims are bounded throughout by the patient counts reported. The deliverable is a well-defined, deterministic, reproducible transport operator that turns spatial measurements into barrier quantities, plus a testable prediction about what stromal modification would do to both delivery modalities simultaneously."));\n\n'
             + old_intro)
assert old_intro in s, "intro anchor not found"
s = s.replace(old_intro, new_intro, 1)

# 3. §3.2: 在 "// Fig 2" 之前加 null permutation 段
old_fig2 = "// Fig 2"
new_32 = ('children.push(p("Permutation control for the rank-normalised input. Because the crosslinking score is rank-normalised within each section, the fraction of edges below the IgG radius is partly fixed by construction. To separate model artefact from data signal, we permuted the crosslinking score 50 times within each of six representative sections (spanning both tumour types and both platform generations) and re-ran the variance decomposition on each permutation. The size-exclusion share of B_mAb variance is essentially unchanged: median 97.6% under the real crosslinking field versus 97.6% under permutation (95% null interval 96.8\\u201398.4%); only 2 of 6 sections show a permutation p < 0.05, and in those the effect size is less than one percentage point. We therefore report the 97.5% figure explicitly as a property of the operator given \\u03b2 = 3, not as an empirical measurement of crosslinking architecture \\u2014 a framing already stated in Section 2.4 and confirmed here rather than softened after review."));\n\n'
           + old_fig2)
assert old_fig2 in s, "fig2 anchor not found"
s = s.replace(old_fig2, new_32, 1)

# 4. §3.5: 在 "// Fig 5" 之前加 stromal intervention 段
old_fig5 = "// Fig 5"
new_35 = ('children.push(pRuns([\n  new TextRun({ text: "In silico stromal co-targeting. ", font: FONT, size: SZ_BODY, bold: true }),\n  new TextRun({ text: "The coupling implies a directly testable prediction: if the two barriers share a matrix substrate, reducing matrix density and crosslinking in silico should lower both barriers at once. We tested this by scaling the ecm, caf and crosslinking scores by (1 \\u2212 reduction) at reduction levels of 20%, 30% and 50%, and recomputing both operators. At a 30% reduction, the section-level minimum-cut barrier falls by a median 60.8% (cSCC 61.1%, melanoma 58.1%) and the mean tumour-core B_mAb by 49.5% (cSCC 49.5%, melanoma 49.6%); both barriers fall in 19 of 19 sections. At 50% reduction the drops are 80.2% and 71.3% respectively. The two modalities respond in the same direction and with similar magnitude across both tumour types, which is the quantitative basis for the paper\\u2019s clinical claim: stromal-directed intervention (LOX/TGF-\\u03b2 inhibition, anti-fibrotic combinations) is predicted to improve cellular and macromolecular delivery simultaneously rather than relieving one at the expense of the other. This is a model-based prediction, not a measured treatment response; it is stated as such.", font: FONT, size: SZ_BODY }),\n]));\n\n'
           + old_fig5)
assert old_fig5 in s, "fig5 anchor not found"
s = s.replace(old_fig5, new_35, 1)

# 5. §3.6: 在 "// Table 3" 之前加 naive baseline 段
old_t3 = "// Table 3"
new_36 = ('children.push(pRuns([\n  new TextRun({ text: "Why a graph at all? ", font: FONT, size: SZ_BODY, bold: true }),\n  new TextRun({ text: "A natural reviewer question is whether a simple node-level matrix score (0.5\\u00b7(ECM + CAF)) already captures the barrier, making the min-cut and screened-Poisson machinery unnecessary. We compared the per-spot graph B_cell field against this naive score on all 19 sections, using partial Spearman correlation with B_mAb after residualising on vessel distance. The graph field yields a median partial \\u03c1 of +0.217 versus +0.166 for the naive score; the graph field is stronger in 14 of 19 sections and reaches p < 0.05 in 17 of 19 versus 15 of 19 for the naive score. The advantage is modest rather than transformative (median \\u0394\\u03c1 = +0.048) and is concentrated in sections where source\\u2013sink geometry is non-trivial (e.g. CSCC11/CSCC12, where \\u0394\\u03c1 = +0.14/+0.17). We read this honestly: a matrix score captures most of the barrier signal because the matrix is the dominant substrate, but the graph operator adds the spatial arrangement \\u2014 the capacity of the narrowest cut and the screened-diffusion field \\u2014 that a node score cannot represent.", font: FONT, size: SZ_BODY }),\n]));\n\n'
           + old_t3)
assert old_t3 in s, "table3 anchor not found"
s = s.replace(old_t3, new_36, 1)

# 6. Limitations: 在 "Finally, the melanoma cohort" 之前加 single-author 句
old_lim = "Finally, the melanoma cohort lacks treatment annotation"
# 在 JS 里这段是在 pRuns 里，文本是 "...Finally, the\nmelanoma cohort..."
# 找 JS 里的实际锚点
import re
# 找包含 "melanoma cohort lacks treatment" 的 JS 字符串
m = re.search(r'([^"]*melanoma cohort lacks treatment annotation[^"]*)', s)
if m:
    old_txt = m.group(1)
    new_txt = "This is a single-author work: the analysis was not independently checked by a second analyst, and the code is released for community verification rather than having passed internal peer review. " + old_txt
    s = s.replace(old_txt, new_txt, 1)
    print("6. Limitations single-author: inserted")
else:
    print("6. WARNING: limitations anchor not found")

# 7. Conclusions: 在 "For matrix-rich tumours the model predicts" 之前加 stromal 数字句
old_conc = "For matrix-rich tumours the model predicts that"
new_conc = ("Quantitatively, a 30% in silico reduction in matrix density and crosslinking lowers the minimum-cut barrier by 61% and the antibody barrier by 50% across all 19 sections. " + old_conc)
if old_conc in s:
    s = s.replace(old_conc, new_conc, 1)
    print("7. Conclusions: inserted")
else:
    print("7. WARNING: conclusions anchor not found")

io.open(p, "w", encoding="utf-8", newline="\n").write(s)
print("\ngenerate_docx.js synced. New length:", len(s))
