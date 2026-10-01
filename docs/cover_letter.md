# Cover Letter — Bioinformatics (Original Paper)

> 用法：填掉全部 `<...>` 占位符后随投稿系统提交。正文一页以内。
> 事实性陈述均已按 2026-08-28 产物核对（19 张切片 / 7 位患者）。

---

Dear Editor,

We are pleased to submit our manuscript, "**Coupled, not separable: the
T-cell migration barrier and the antibody mass-transport barrier in cutaneous
tumours share a matrix origin**", for consideration as an Original Paper in
*Bioinformatics*.

Immunotherapy resistance in solid tumours is usually discussed in terms of
cellular immunity, but drug delivery poses two distinct transport problems:
whether a cytotoxic T cell can migrate through the stromal fibre network to
reach the tumour nest, and whether the antibody itself — a 5.5 nm
macromolecule — can diffuse through a matrix whose effective mesh is measured
in tens of nanometres. The second problem is almost never quantified, and
whether the two barriers are separable in real tissue had not been asked.
SPARTA poses both as solvable transport problems on a single
spatial-transcriptomic graph, changing only the edge-weight semantics and the
operator: a source–sink minimum cut for the cellular barrier and a screened
Poisson diffusion–absorption field for the macromolecular barrier. The core
operators depend only on NumPy, SciPy and NetworkX and run entirely on CPU.

We applied the framework to 19 public spatial-transcriptomic sections from
seven patients — six with primary cutaneous squamous carcinoma, profiled
across two platform generations, and one with four melanoma metastases. The
central finding is a clarification rather than a discovery: after controlling
for distance to vasculature, the two barriers are **positively** coupled in
18/19 sections, discordant regions are fewer than chance would predict, and
removing the matrix term shared by the two operators leaves the association
significant in 12/15 squamous sections with a median 81 % of its original
magnitude. The coupling is therefore substantially a property of the tissue,
and its clinical corollary — that matrix-directed intervention could improve
both classes of delivery at once — follows directly from the analysis.

We believe this work suits *Bioinformatics* for three reasons. First, the
methodological contribution is explicit and reusable: two solvable
graph-transport operators with a complete counterfactual framework
(contiguous-arc removal, molecular-radius sweeps) that domain segmentation
cannot support. Second, the validation is deliberately adversarial to our own
framework: every section is reported at both section and patient level,
parameters are pre-designated by role (physically anchored, unsupervised, or
sensitivity-scanned), and the domain-count sensitivity scan, the comparison
against BANKSY-style segmentation and Squidpy neighbourhood enrichment, and
the shared-input removal analysis are all reported with their unfavourable
results. Third, all data are public (GEO: GSE250636, GSE144239) and the
complete pipeline reproduces every figure from raw downloads with one command
per stage.

This manuscript is not under consideration elsewhere, and all authors have
approved the submission. The authors declare no competing interests. All data
used are publicly available; the code is released under the MIT licence at
`<GitHub URL>` and archived at `<Zenodo DOI>`.

Thank you for your consideration.

Sincerely,

`<Corresponding author name>`
`<Affiliation>`
`<Email>`

---

## 提交前自查清单

- [ ] 标题用 results 里的 working title（若编辑偏好 tool-forward，备选
      "SPARTA: a graph-transport framework that separates the operators but
      not the barriers behind immunotherapy resistance in skin cancer"）
- [ ] `12/15 squamous sections / 81 %` 与 Abstract、R3b 一致（2026-08-28 已核对）
- [ ] GitHub URL / Zenodo DOI 填好后再提交（Availability 的硬性条件）
- [ ] 通讯作者信息、日期
- [ ] 建议审稿人：按期刊要求另填（3 名，非合作者、非同机构）
