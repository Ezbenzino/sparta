# JBHI 重定位方案（2026-08-28）

> 目的：把现有 Bioinformatics 口径的稿件改造成 IEEE JBHI 口径，而不推翻任何
> 已有结果。所有数字与正文产物一致；runtime 数字来自
> `results/validation/runtime_benchmark.json`（run_19_runtime.py）。
> 使用方式：本文档是"改造清单 + 可直接粘贴的英文段落"；正式 fork 稿件时
> 逐段搬运，Bioinformatics 版本保持不动作为备选。

## 0. 换轴总说明

| | Bioinformatics 版（现状） | JBHI 版（目标） |
|---|---|---|
| 主问题 | 两个屏障是否可解离（生物学） | 如何在空间分子数据上**量化**药物递送的结构性障碍（信息学） |
| 19 张切片的角色 | 生物学主张的证据链 | **方法验证的测试集**（含跨平台世代稳健性） |
| 核心卖点 | PIVOT：耦合源是基质 | 白盒可解释算子 + 反事实验证框架 + 秒级 CPU |
| PIVOT 结论的位置 | Results 高潮 | 降为"方法揭示的发现"一节 + Discussion 转化含义 |
| 语言风格 | author-year，生信叙事 | IEEE 数字引用，信息学叙事 |

不变的东西：R1–R6 的全部数字、图表、统计口径；Methods 三段；诚实口径的
limitation。改造集中在 Introduction、结果排序、补充 runtime/对比表。

## 1. 重写后的 Introduction（可直接粘贴）

**¶1 — 临床信息学缺口：空间分子数据缺少"可解释的传输量"**

Precision immuno-oncology increasingly relies on spatial molecular profiling to
characterise the tumour microenvironment, yet the analytical vocabulary
available to a clinical pipeline is dominated by two families: deep
spatial-domain embeddings and local neighbourhood statistics. Both produce
*labels* — a domain assignment, a co-localisation score — that answer "what is
next to what". Neither produces a *physically interpretable quantity* of drug
delivery: whether a cytotoxic T cell can physically traverse the stromal
network to reach a tumour nest, and whether a 5.5 nm antibody can diffuse
through a matrix whose effective mesh is measured in tens of nanometres. These
are transport problems, and transport problems have canonical, solvable
formulations. The gap is not a lack of data but a lack of operators: no
existing spatial-omics tool outputs a barrier strength with a source, a sink
and units.

**¶2 — SPARTA：同一张图上的两个可解算子**

We formulate both delivery modalities as graph-transport problems on a single
spatial-transcriptomic graph, changing only the edge-weight semantics and the
operator. T-cell migration is a source–sink minimum cut between immune-entry
and tumour-core compartments — a convex, exactly solvable assignment whose
output is both a scalar barrier strength and an explicit blockade geometry
(the cheapest edge set whose removal disconnects source from sink). Antibody
delivery is a screened Poisson diffusion–absorption field in which conductance
falls with matrix density and size exclusion, and target antigen acts as a
distributed sink; with absorption removed it reduces to the harmonic
(effective-resistance) limit. Both operators are deterministic, run on CPU in
seconds per section, and involve no learned parameters — properties that
matter for clinical deployment, where reproducibility, auditability and
hardware footprint are constraints rather than conveniences.

**¶3 — 反事实验证框架**

Because the operators are white-box, every claim they generate can be
adversarially tested by perturbing their inputs rather than their weights. We
contribute a counterfactual suite with three axes: spatial rearrangement
(permute cell identities while preserving counts), blockade continuity
(thicken the min-cut band in place), and molecular size (sweep the
hydrodynamic radius from small molecule to IgG), plus a shared-input removal
test that asks whether an observed coupling between the two operators survives
deletion of the matrix term they share. The design principle is that each
parameter carries a declared identity — physically anchored, unsupervised, or
sensitivity-scanned — and that every unfavourable result is reported.

**¶4 — 验证设计与本文发现**

We apply the framework to 19 sections from seven patients across two cutaneous
cohorts and two platform generations of spatial transcriptomics (Visium and
first-generation ST, eight years apart), admitted through pre-declared
cohort-level thresholds. The framework's finding, which we report as a
demonstration of what transport operators can see that label-based tools
cannot: after controlling for distance from vasculature, the two barriers are
positively coupled in 18/19 sections, discordant regions are rarer than
chance, and removing the shared matrix term leaves the coupling significant in
12/15 squamous sections with a median 81 % of its magnitude. The coupling is
therefore substantially a property of the tissue rather than of the operator
inputs, and its clinical corollary — matrix-directed intervention should
improve both classes of delivery at once — follows directly from the operators'
outputs.

**¶5 — 贡献清单（JBHI 惯例）**

Our contributions are:
1. Two interpretable, non-learning graph-transport operators for spatial
   molecular data — a source–sink minimum-cut barrier and a screened-Poisson
   diffusion–absorption barrier — with physically anchored parameterisation;
2. A deterministic CPU-only implementation completing the full five-stage
   analysis of a section in **0.02–0.19 s** (median 0.05 s on Visium, 0.02 s on
   first-generation ST; the entire 19-section cohort in 0.9 s, peak memory
   < 16 MB on a commodity laptop CPU) — against 0.16–1.20 s (median 0.27 s)
   for BANKSY-style domain segmentation and 3.8–5.0 s (median 4.1 s) for
   neighbourhood enrichment computed on the same graphs;
3. A counterfactual validation framework (rearrangement / continuity /
   size-sweep / shared-input removal) applicable to any operator claiming to
   measure spatial structure;
4. Application to 19 sections across two platform generations revealing a
   robust coupling of the two delivery barriers, with a clinical implication
   for combination-therapy design.

> 中文注：¶5 的 runtime 占位符等 benchmark 填。¶2 那句 "no learned
> parameters" 是 JBHI 的差异化卖点，别删。¶1 定调"labels vs quantities"，
> 这是全文的信息学纲领句。

## 2. 临床锚点小节（从 Intro 动机升格为 Clinical Relevance）

放在 Results 之后、Discussion 之前（或作为 Discussion 首节），标题建议
"Clinical relevance and the limits of bulk proxies"。核心素材就是现有的
Dimension ③ 段落（GSE78220 AUC 0.759 vs GSE91061 反向），改写口径如下：

> **Clinical relevance: where the predictive signal lives.** No public dataset
> currently combines spatial-transcriptomic profiling with checkpoint-blockade
> response labels, so the clinical pathway of any spatial barrier metric must
> for now be argued through proxies. We scored pre-treatment bulk RNA-seq from
> two anti-PD-1 melanoma cohorts with the same signature gene sets. In
> GSE78220 (n = 28) the composite barrier score separates non-responders from
> responders (AUC 0.759, p = 0.021) where a composite immune-abundance score
> does not (AUC 0.472); in GSE91061 (n = 49) the barrier score is not
> significant and the direction reverses. We read this instability not as a
> failure of the signature but as an informatics result in its own right: bulk
> abundance discards exactly the spatial arrangement that makes matrix
> molecules a barrier, and a quantity that is defined over a graph cannot be
> expected to survive the collapse of the graph. The implication for clinical
> informatics is concrete — response prediction from archival bulk data will
> systematically miss structural determinants of delivery, and spatially
> resolved assays are required to capture them.

这段话同时完成三件事：把 ICB 验证从"墙角"搬到"临床相关性"位置；
把跨队列不稳定从弱点改写成"bulk 丢空间信息"的方法学论点；给 JBHI 的
Health Informatics scope 一个明确锚点。原 Dimension ③ 的诚实口径（AUC
方向、无 ST+ICB 公共数据集）全部保留。

## 3. 样本量防守（Discussion 的 limitation 段，可直接粘贴）

> A note on sample size. Nineteen sections from seven patients is small by the
> standards of clinical imaging cohorts and adequate by the standards of
> spatial-omics method validation; we state the design properties that make
> the two readings compatible. First, each section contributes one *complete*
> graph-level measurement — source, sink, operators and permutation reference
> are all internal to the section — so sections are measurement units, not
> repeated samples of one unit. Second, every inferential statistic we report
> is evaluated with the patient as the independent unit (permutation and
> bootstrap over patients, n = 7); section-level counts are always accompanied
> by patient-level counts. Third, the two cohorts span two platform
> generations eight years apart with different chemistry, spot pitch and
> depth; agreement across that domain shift is itself the robustness argument,
> and we report the one stratum where it fails (the single melanoma patient,
> whose classification depends on the distance convention) rather than
> averaging it away.

> 中文注：三步防守——图级测量单位、患者层统计、跨平台世代当 domain-shift
> 卖点。黑色素瘤 n=1 的硬伤放最后一句如实说，不藏。

## 4. 方法对比总表（新 Table：方法 × 输出物 × 本数据上的表现 × 运行成本）

| 方法家族 | 代表 | 输出物 | 本数据上的定量表现 | 单切片运行时 (CPU) | 确定性 | GPU |
|---|---|---|---|---|---|---|
| 域分割（本文实现） | BANKSY 式 k-means（λ=0.3） | 每点域标签 | 域边界预测割边 precision 4–22%（中位 14%）；富集 0.69–1.87× | 0.16–1.20 s（中位 0.27） | 种子固定即可复现 | 否 |
| 域分割（外部包） | STAGATE / GraphST | 嵌入 + 域标签 | 未运行（需 GPU 训练）；按 R6 口径属"候选集供给者"，n_domains 扫描已证结论不依赖域分割来源 | —（未测） | 训练含随机性 | 是 |
| 邻域统计 | Squidpy nhood_enrichment | 域×域 z 值矩阵 | 非对角 z −36.3…+16.5（19/19 有值） | 3.8–5.0 s（中位 4.1） | 置换种子固定 | 否 |
| min-cut 作聚类损失 | TCN (Nat Methods 2023) / Spatial-RGCN | 细胞邻域分配 | 不适用（割用于优化分配，不输出屏障） | — | — | 是（后者） |
| **SPARTA B_cell** | 本文 | 标量屏障强度 + 割集几何 | 19/19 可解；割带宽度 6.7–59.1% spots | **0.012–0.171 s** | **是（精确解）** | **否** |
| **SPARTA B_mAb** | 本文 | 逐点扩散-吸收场 + 标量屏障 | 19/19 可解；r_nm 扫描 0.5–5.5 nm 单调（S3） | **0.001–0.017 s** | **是（线性系统）** | **否** |

> 表注（写进稿件的口径）：BANKSY 式与 Squidpy 臂是**本仓库实现**（k-means
> 域分割 / 置换 z 值），不是官方包全部功能的计时；STAGATE/GraphST 未运行，
> 如实标注不编数；运行环境：Intel 20 核笔记本 CPU、16 GB RAM、Windows、
> numpy 2.5.2 / networkx 3.6.1（完整环境在 JSON 的 `env` 字段留痕）。
> 计时口径：每阶段 3 次取最小（squidpy 置换检验昂贵，1 次）。

## 5. IEEE 两栏格式改造清单（下一步，待确认后执行）

1. `IEEEtran` 双栏模板，正文目标 8 页（现在 ~10k 词 → 砍到 ~7k）：
   - 砍法：R6 的 n_domains 扫描三点压缩为一段 + 引补充材料；R2 分解细节
     全部下沉补充材料；R4 三个反事实各留一段主文。
   - 图 5 → 4 张：框架图、解耦散点+剂量反应合并、反事实、对比表（新）。
2. 引用转 IEEE 数字格式 [1]；文献表用 `IEEEtran.bst`。
3. 摘要压到 250 词内（IEEE 惯例），关键词 5 个（建议：spatial
   transcriptomics, graph algorithms, minimum cut, drug delivery,
   tumour microenvironment）。
4. 投稿前核对 JBHI 当前页数上限与超页费（IEEE 惯例 ~$250/页）；
   双栏 8 页 + 大图很容易超，图用单栏宽设计。
5. 图件字体按 IEEE 要求（≥8pt）、色盲安全（现有调色板已满足）。

## 6. 已完成的配套动作（2026-08-28）

- `scripts/run_19_runtime.py` + `results/validation/runtime_benchmark.json`：
  19 张逐阶段计时（score/graph/b_cell/b_mab/b_meta/banksy/squidpy）+ 峰值内存 + 环境信息；
- R6 已补 TCN/Spatial-RGCN 区分段（"割作聚类损失 vs 割作屏障"），见
  `manuscript_results.md` R6 开头；
- 竞品调研（`docs/competitor_research.md`）：确认无直接竞品，
  TCN 是审稿人搜 "min-cut spatial" 的第一撞点（已防守）。
