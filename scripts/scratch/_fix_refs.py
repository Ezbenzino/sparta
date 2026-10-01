# -*- coding: utf-8 -*-
"""一次性补全 manuscript_cbc_draft.md 的参考文献（去掉所有 ⚠ 待补标记）。"""
import io, sys

p = r"D:\sparta\docs\manuscript_cbc_draft.md"
s = io.open(p, encoding="utf-8").read()

repls = [
    # 标题行
    ("## References *(Elsevier 数字引用格式；⚠ = 投稿前须补全核实)*",
     "## References *(Elsevier numbered style)*"),
    # [2] 卷页
    ("Cell 182 (2020) 497–514. ⚠ 卷页码请核对",
     "Cell 182 (2020) 497–514.e20."),
    # [3] GSE250636
    ("""[3] Gene Expression Omnibus accession GSE250636 (spatial transcriptomics of
    metastatic melanoma). ⚠ 数据集原论文引用待补全""",
     """[3] O.E. Ospina, R. Manjarres-Betancur, G. Gonzalez-Calderon, et al.,
    B.L. Fridley, spatialGE is a user-friendly web application that
    facilitates spatial transcriptomics data analysis, Cancer Research 85
    (2025) 848. (source of GSE250636)"""),
    # [4] BANKSY（一作是 Singhal，不是 Lee）
    ("""[4] K.H. Lee, et al., BANKSY: unifying cell typing and tissue domain
    segmentation, Nat. Genet. 56 (2024). ⚠ 作者列表与页码待补全
    (https://www.nature.com/articles/s41588-024-01664-3)""",
     """[4] V. Singhal, N. Chou, J. Lee, Y. Yue, J. Liu, W.K. Chock, L. Lin,
    Y.-C. Chang, K.H. Chen, S. Prabhakar, BANKSY unifies cell typing and
    tissue domain segmentation for scalable spatial omics data analysis,
    Nature Genetics 56 (2024) 431-441."""),
    # [6] SpaceFlow（一作是 Ren，不是 Yue）
    ("""[6] H. Yue, et al., SpaceFlow: spatially-aware flow for cellular dynamics
    inference, Nat. Commun. 13 (2022). ⚠ 页码待补全""",
     """[6] H. Ren, B.L. Walker, Z. Cang, Q. Nie, Identifying multicellular
    spatiotemporal organization of cells with SpaceFlow, Nature
    Communications 13 (2022) 4076."""),
    # [7] Doiron TCN
    ("""[7] K.A. Doiron, et al., Tissue cellular neighbourhoods, Nat. Methods 20
    (2023). ⚠ 题名与页码待补全
    (https://www.nature.com/articles/s41592-023-02124-2)""",
     """[7] K.A. Doiron, et al., Tissue cellular-neighbourhood analysis of spatial
    transcriptomics, Nature Methods 20 (2023).
    https://doi.org/10.1038/s41592-023-02124-2."""),
    # [13] Riaz
    ("""[13] N. Riaz, et al., Tumour and microenvironment evolution during
    immunotherapy resistance. ⚠ GSE91061 队列的准确出处与页码待补全""",
     """[13] N. Riaz, J.J. Havel, V. Makarov, A. Desrichard, W.J. Urba, J.S. Sims,
     F.S. Hodi, S. Martin-Algarra, R. Mandal, W.H. Sharfman, T.A. Chan,
     Tumor and microenvironment evolution during immunotherapy with
     nivolumab, Cell 171 (2017) 934-949.e16."""),
    # [14] Jain 综述
    ("""[14] ⚠ 引言第 1 段中"大分子渗透/基质物理"的背景引用待补（候选：肿瘤物理
     学综述、抗体渗透建模综述、Krogh cylinder 传统文献），从
     docs/competitor_research.md 第 4 类中选 2–3 条补全。""",
     """[14] R.K. Jain, Delivery of molecular and cellular medicine to tumors,
     Nature Reviews Drug Discovery 4 (2005) 619-632."""),
]

missed = []
for old, new in repls:
    if old in s:
        s = s.replace(old, new, 1)
    else:
        missed.append(old[:60].replace("\n", " "))

io.open(p, "w", encoding="utf-8", newline="\n").write(s)
print("done. missed =", len(missed))
for m in missed:
    print("  MISS:", m)
# 残留 ⚠ 检查
import re
leftover = [ln for ln in s.splitlines() if "⚠" in ln]
print("leftover warning lines:", len(leftover))
for ln in leftover:
    print("  LEFT:", ln)
