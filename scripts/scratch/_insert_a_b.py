# -*- coding: utf-8 -*-
"""把 A (biological validation) 和 B (radius sensitivity) 插入 manuscript markdown"""
import io

p = r"D:\sparta\docs\manuscript_cbc_draft.md"
s = io.open(p, encoding="utf-8").read()

# B: radius sensitivity —— 在 "### 3.3." 之前插入
b_text = (
    "**Graph-radius sensitivity.** Because the spatial graph is constructed by a hard radius threshold, "
    "we tested whether the main conclusions depend on that choice. We re-graphed all 19 sections at "
    "four Visium radii (100, 150, 200, 250 \\u00b5m) and three first-generation-ST radii (200, 300, 400 \\u00b5m), "
    "re-defined source/sink/vessel sets at each radius, and recomputed both operators (65 configurations total). "
    "The median partial coupling between the B_cell field and B_mAb is +0.225 across configurations "
    "(96.9% of configurations positive, versus +0.217 at the default radius). The median B_cell reduction after "
    "30% stromal scaling is 59.6% and the median B_mAb reduction is 49.5% (versus 60.8% and 49.5% at default). "
    "Both barriers fall simultaneously in 89.2% of configurations. The coupling, the stromal co-targeting prediction, "
    "and the direction of the size-exclusion result are therefore not artifacts of a single graph radius.\n\n"
)

anchor_b = "### 3.3. Molecular size separates the two transport problems"
assert anchor_b in s, "anchor B not found"
s = s.replace(anchor_b, b_text + anchor_b, 1)

# A: biological validation —— 在 "## 4. Discussion" 之前插入
a_text = (
    "**Alignment with measured cell distributions.** Beyond internal consistency, we asked whether the per-spot "
    "barrier fields align with measured cell distributions in the same tissue. After residualising on vessel distance, "
    "the B_cell field shows a weak but directionally consistent negative partial correlation with the T/NK signature "
    "(median \\u03c1 = \\u22120.022, significant in 6/19 sections) and with the CD8 T-cell signature "
    "(median \\u03c1 = \\u22120.013, significant in 5/19): barrier-high spots tend to contain fewer T cells, "
    "though the effect is small because T-cell localisation is also driven by antigen availability and inflammatory "
    "cues not modelled here. The B_mAb field shows a weak positive partial correlation with proliferation "
    "(median \\u03c1 = +0.046, significant in 7/19), consistent with antibody-blocked nests retaining proliferating cells. "
    "We do not overstate these alignments: they are weak, they point in the expected direction, and they provide an "
    "independent (if modest) check that the operators are not purely mathematical constructs.\n\n"
)

anchor_a = "## 4. Discussion"
assert anchor_a in s, "anchor A not found"
s = s.replace(anchor_a, a_text + anchor_a, 1)

io.open(p, "w", encoding="utf-8", newline="\n").write(s)
print("markdown updated, length:", len(s))
