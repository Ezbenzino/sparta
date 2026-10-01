# SPARTA — Abstract（草稿 v1.0 · 数字全部为 2026-08-27 定稿值）

> 按 Bioinformatics 的 Original Paper 格式（Motivation / Results / Availability）。
> 每个数字后面括号里标了产物来源，投稿前删掉。
> **补入第一代 ST 切片之后必须重写样本量与队列那两句**（见 submission_plan §1d）。

---

## Motivation

Immune checkpoint blockade fails in most patients with advanced cutaneous
malignancy, and spatial transcriptomics has made the microenvironmental
contribution to that failure directly measurable. The field's analytical
vocabulary, however, addresses only one of the two transport problems involved.
Whether a cytotoxic T cell can reach the tumour nest is a cell-migration problem
through a fibre network; whether the anti-PD-1 antibody itself can reach it is a
mass-transport problem for a 5.5 nm macromolecule through a matrix whose
effective mesh is measured in tens of nanometres. The second is almost never
quantified, and whether the two barriers are separable in real tissue has not
been asked.

## Results

We present SPARTA, which poses both questions as transport problems on one
spatial graph, changing only the edge-weight semantics and the operator: T-cell
migration as a source–sink minimum cut, which returns both a barrier magnitude
and an explicit blockade geometry, and antibody delivery as a screened Poisson
diffusion–absorption field that reduces to the harmonic/effective-resistance
limit when absorption vanishes. We applied it to **19 sections from seven
patients** — six primary cutaneous squamous carcinomas profiled across two
platform generations (four Visium, eleven first-generation ST) and one melanoma
patient contributing four extracranial metastatic deposits.

The two barriers prove **not** to be separable. After residualising both on the
weighted graph distance to vasculature the partial correlation is positive in
18/19 sections and significant in 17 (ρ = −0.000 to +0.385, median +0.217), every
one of the seven patients has a positive median, and regions where one barrier is
high while the other is low occupy 1.1–6.9 % of spots against a 6.25 % chance
expectation — fewer discordant spots than independence predicts, in 18/19
sections. The result reproduces across platform generations within the squamous
cohort (median ρ = +0.270 on Visium, +0.207 on first-generation ST), which
separates assay generation from the tumour-type and disease-stage covariates that
were confounded with it in a smaller design.

Removing the matrix term shared by the two operators leaves the association
significant in 15/19 sections — 12/15 squamous sections (a majority in five of
the six squamous patients) with a median 81 % of its original magnitude, and
3/4 melanoma sections with a median 38 % — unchanged whether the geometric control
is a weighted path or a hop count. The single melanoma patient is the one stratum
whose classification depends on that choice, and we report it as unresolved
rather than as a cohort contrast.

Counterfactual experiments support the topological reading: a contiguous gap in
the blockade leaves 1.04–1.39× more residual barrier than removing the same
material scattered within it (16/19 sections, five of seven patients with every
section significant, at a
pre-specified 20 % removal fraction; the effect is larger at 30 % than at 5 % in
15/19), and sweeping molecular radius on fixed tissue raises the tumour-core
barrier monotonically across an order of magnitude. Relative to a spatial-domain
segmentation, just 4–22 % of domain boundaries lie on the cut on either platform
— the transport formulation selects and quantifies the one boundary that limits
flux rather than detecting boundaries.

## Availability and implementation

SPARTA is implemented in Python and runs on CPU; the barrier operators depend
only on NumPy, SciPy and NetworkX. Source code, the analysis pipeline and the
scripts reproducing every figure are available at <GitHub URL> and archived at
<Zenodo DOI>. All data are public (GEO: GSE250636, GSE144239).

---

## 写作说明（投稿前逐条处理）

1. **样本量那句必须先说患者数再说切片数。** 现已写成 "19 sections from seven
   patients"（2026-08-27 扩样后更新）。绝对不要写成 "19 samples"。
   摘要里每一个计数都同时给了切片层与患者层，投稿前逐句核一遍别漏。
2. **交联主导（96–98%）没有进摘要，这是刻意的。** 它几乎只是 β 的单调函数
   （β=0.25 时 0.2%），放进摘要等于把一个建模选择当成发现。正文 R2 里
   带着参数依赖一起讲。
3. **"tissue-borne in one tumour type"这句是全篇最强也最脆的一句。**
   它现在靠 2 位鳞癌患者。补片之后如果 6 位仍然全中，这句可以留；
   如果不是，改成 "in the squamous sections examined here"。
4. Availability 两个 <> 占位符是硬性投稿条件，git/Zenodo 没做完不要投。
5. 字数：Results 段目前约 290 词，Bioinformatics 没有硬性上限但偏好紧凑，
   定稿时可以把 S3 那句压成半句。
