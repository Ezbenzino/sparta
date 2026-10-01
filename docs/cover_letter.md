# Cover Letter — Computational Biology and Chemistry

> 填好后随投稿系统提交，正文一页以内。

---

Dear Editor,

I am pleased to submit my manuscript, "**Two transport operators, one
substrate: a graph model of T-cell migration and antibody penetration
barriers in tumour tissue from spatial transcriptomics**", for consideration
as a research article in *Computational Biology and Chemistry*.

Monoclonal antibodies and cytotoxic T cells must both cross the tumour
extracellular matrix to reach their targets, yet they differ by three orders
of magnitude in size (a 5.5 nm IgG versus a ~10 µm T cell). The classical
modelling tradition for macromolecular delivery — compartmental
pharmacokinetics, Krogh-cylinder penetration models, and reaction–diffusion
descriptions — uses effective parameters on idealised geometry and is not
driven by molecular measurements of a specific patient section. Spatial
transcriptomics now provides those measurements (local collagen, crosslinking,
stromal and immune signatures, and the therapeutic target itself), but
existing analyses return labels rather than transport quantities.

SPARTA closes this gap by formulating the two delivery problems as two
solvable operators on one measured spatial graph, changing only the
edge-weight semantics and the operator. T-cell migration is a source–sink
minimum cut between immune-entry and tumour-core compartments, returning a
barrier strength together with its blockade geometry. Antibody transport is a
screened Poisson diffusion–absorption field whose conductances implement size
exclusion against the IgG hydrodynamic radius and whose sinks represent
target-antigen binding. The core operators depend only on NumPy, SciPy and
NetworkX, run entirely on CPU, and analyse one section in under 0.2 s.

Applied to 19 public sections from 7 patients across two cutaneous tumour
types and two spatial-platform generations, the model establishes three
quantitative facts. First, the antibody barrier is percolation-limited: the
effective mesh falls below the IgG radius on 60–65% of edges. Second, the two
barriers are mathematically distinct (the same graph is nearly transparent
to a 0.5 nm solute and strongly obstructive to an IgG) yet physically
co-localised: they remain positively correlated in 18/19 sections, and
removing the matrix input shared by both operators leaves the coupling
significant in 12/15 squamous sections with a median 81% of its magnitude.
Third — and this is the paper's testable prediction — an in silico 30%
reduction in matrix density and crosslinking lowers the minimum-cut barrier
by a median 61% and the antibody barrier by 50%, in 19 of 19 sections across
both tumour types. The practical corollary is that stromal-directed
intervention is predicted to improve both classes of delivery at once, rather
than relieving one at the expense of the other.

Every analytical decision was fixed before looking at results: seven
pre-declared admission criteria (one section, CSCC13, was rejected on a
pre-set UMI threshold and its record retained), parameters are partitioned
into physically anchored, qualitative scale, and sensitivity-scanned classes,
and no parameter was fitted to any clinical outcome. The paper also reports
its own negative results honestly: the 97.5% size-exclusion share is shown
by a permutation control to be a property of the rank-normalised operator
rather than an empirical measurement; the melanoma stratum is a single patient
and is reported as unresolved rather than as a cross-tumour contrast; and the
weak antigen channel is attributed to the detection floor of a two-gene
signature.

I believe this work fits *Computational Biology and Chemistry* specifically:
its central object is a physicochemical transport model of biologic-drug
delivery — diffusion–absorption, size exclusion, and mesh topology — anchored
to measured tissue geometry rather than a purely statistical or clustering
contribution. It is deterministic, parameter-honest, and reproducible: all
data are public (GEO GSE250636, GSE144239; MSigDB HALLMARK_HYPOXIA), and the
code is released under the MIT licence at https://github.com/Ezbenzino/sparta,
archived at https://doi.org/10.5281/zenodo.23086431.

This manuscript is not under consideration elsewhere, and the author has
approved the submission. The author declares no competing interests and
received no specific funding for this work.

Thank you for your consideration.

Sincerely,

Yize Li
Hangzhou Medical College
lllyz630031258@gmail.com

---

## 提交前自查

- [x] 标题与 manuscript_cbc_draft.md 一致
- [x] GitHub URL / Zenodo DOI 已填
- [x] 通讯作者信息
- [ ] 建议审稿人：投稿系统另填 3 名（非合作者、非同机构）
- [ ] 投稿时在系统选 "Computational Biology and Chemistry"，订阅路线无版面费
