# -*- coding: utf-8 -*-
"""把 A 和 B 段落插入 generate_docx.js"""
import io

p = r"D:\sparta\scripts\generate_docx.js"
s = io.open(p, encoding="utf-8").read()

# B: radius sensitivity —— 在 "// Fig 2" 之前
b_js = (
    'children.push(pRuns([\n'
    '  new TextRun({ text: "Graph-radius sensitivity. ", font: FONT, size: SZ_BODY, bold: true }),\n'
    '  new TextRun({ text: "Because the spatial graph is constructed by a hard radius threshold, we tested whether the main conclusions depend on that choice. We re-graphed all 19 sections at four Visium radii (100, 150, 200, 250 \\u00b5m) and three first-generation-ST radii (200, 300, 400 \\u00b5m), re-defined source/sink/vessel sets at each radius, and recomputed both operators (65 configurations total). The median partial coupling between the B_cell field and B_mAb is +0.225 across configurations (96.9% positive, versus +0.217 at default). The median B_cell reduction after 30% stromal scaling is 59.6% and B_mAb reduction 49.5% (versus 60.8% and 49.5% at default). Both barriers fall simultaneously in 89.2% of configurations. The coupling, the stromal co-targeting prediction, and the direction of the size-exclusion result are therefore not artifacts of a single graph radius.", font: FONT, size: SZ_BODY }),\n'
    ']));\n\n'
)

anchor_b = "// Fig 2"
assert anchor_b in s, "anchor Fig2 not found"
s = s.replace(anchor_b, b_js + anchor_b, 1)

# A: biological validation —— 在 "// Table 3" 之后、"// Discussion" 之前
# 找 §4 Discussion 的 JS 注释
anchor_a = 'children.push(h1("4. Discussion"));'
assert anchor_a in s, "anchor Discussion not found"

a_js = (
    'children.push(pRuns([\n'
    '  new TextRun({ text: "Alignment with measured cell distributions. ", font: FONT, size: SZ_BODY, bold: true }),\n'
    '  new TextRun({ text: "Beyond internal consistency, we asked whether the per-spot barrier fields align with measured cell distributions in the same tissue. After residualising on vessel distance, the B_cell field shows a weak but directionally consistent negative partial correlation with the T/NK signature (median \\u03c1 = \\u22120.022, significant in 6/19 sections) and with the CD8 T-cell signature (median \\u03c1 = \\u22120.013, significant in 5/19): barrier-high spots tend to contain fewer T cells, though the effect is small because T-cell localisation is also driven by antigen availability and inflammatory cues not modelled here. The B_mAb field shows a weak positive partial correlation with proliferation (median \\u03c1 = +0.046, significant in 7/19), consistent with antibody-blocked nests retaining proliferating cells. We do not overstate these alignments: they are weak, they point in the expected direction, and they provide an independent (if modest) check that the operators are not purely mathematical constructs.", font: FONT, size: SZ_BODY }),\n'
    ']));\n\n'
)

s = s.replace(anchor_a, a_js + anchor_a, 1)

io.open(p, "w", encoding="utf-8", newline="\n").write(s)
print("generate_docx.js synced, length:", len(s))
