# SPARTA v2.2.0 变更说明（2026-10-05）

> 给作者本人看的工作记录。所有数字都由 `scripts/is_figures/facts.py` 从结果文件生成，稿件里没有手填数字。
> 锁定的主队列关联数字（19 张切片 / 7 名患者）**一个都没有变**。

---

## 一、今晚做了什么（按重要性）

### 1. 最小割的精确性 bug（必须知道）

- networkx 的 preflow-push 在**浮点**容量上算出的最大流数值是对的（与整数精确解相差 < 5×10⁻¹²），
  但从残量图读出的割边集合**不一定是最小割**：30 张切片图里 14 张的割容量与最大流相差 0.03%–9.7%。
- 修复：`barrier.exact_min_cut` 把容量取整到 2⁻⁴⁰ 的整数倍再求解，割容量 = 最大流（误差 4×10⁻¹²）。
- 影响：`B_cell` 数值不变；凡是用到"割在哪里"的结果全部重算——S2 缺口实验、BANKSY 式域比较、
  干预图、CODEX 的 blockade-line 次要分析、所有画割的图。
  - S2（主队列）中位比值 1.053 → **1.046**，BH 显著 8/19；患者层面 1.065（1.024–1.108），7/7 患者 > 1。
  - 域边界富集中位 1.15 → **1.14**，precision 不变（14.3%）。
- **更正我昨晚的一句话**：我曾说"干预图里割外 spot 也有非零效应，是因为把基质设为第 5 百分位"——这是错的。
  真正原因是浮点割不精确；改成精确割后，割外 spot 的单点效应严格为 0（最大流-最小割对偶）。第 5 百分位的约定保留不变。

### 2. 你提的四条新要求

| 要求 | 结果 | 位置 |
|---|---|---|
| 干预排序进论文 | 干预图从主图移到 **ESM_1 Fig. S7**（新增 e 面板：区域构成）；结果 3.9 节 + 讨论新增一段。精确割重算后的数字：细胞屏障**高度局灶**——单个最有效 spot 去掉中位 19%（最高 46%，CSCC07），top 1% 联合去掉中位 **62%**（33–86%）；抗体屏障**高度弥散**——top 1% 只去掉 12%，top 5% 31%，top 10% 44%；两图 top 5% 重叠 Jaccard 中位 **5.7%**（2–20%）。你原来记的"46%"对应的是"单个最有效 spot"，不是 top 1%。 | 3.9、讨论第 5 段、ESM_1 S12、Table S9 |
| Thrane 2018 外部验证 | 基因名后缀已修，8 张切片 / 4 名患者走完整流程（QC、打分、两个算子、零模型、患者层面、干预图、简单基线）。作为**独立复现队列**（比"探索性外部臂"更强：有预先写好的复现标准且达标），并在 3.6 节加了一句：主队列以外 6 名患者的 11 张切片关联全部为正。 | 2.1、3.3、3.4、3.6、Table S7 |
| 简单基线对比 | 新 **Fig. 8**（主图）。三个基线：基质密度、到肿瘤边界距离、邻域富集（Squidpy 同算法）。切片内：B_cell 场与局部基质密度中位 ρ = 0.50（共享约 25% 秩方差），与距离 0.07、与邻域富集 0.34；B_mAb 场与三者都只有 −0.14～−0.16。切片间：log B_rel 与三者都不相关。有真值时：模拟里没有任何简单量能分开闭合/5% 缺口包膜（AUC ≤ 0.72，割 0.85）；CODEX 里基质密度与割一样好，**肿瘤核心深度（post hoc）是最强单一相关量（−0.60）**，但给定组成 + 两个距离后割仍保留偏相关 −0.20（−0.35～−0.04），给定全部四个简单量后不再显著（−0.11）。 | 2.14、3.8、讨论第 4 段、Fig. 8、ESM_1 S14、Table S10 |
| top 1%/5% 干预位点的区域富集 | 两个算子的 top 5% 都**避开免疫区**（29%、28% vs 36%；患者层面 sign-flip p = 0.027、0.012）；B_cell top 5% 轻度偏基质（36% vs 34%，9/12 患者，p = 0.024）；B_mAb top 5% 偏肿瘤（37% vs 30%，但只有 8/12 患者，p = 0.115）。区域是表达定义的，不是病理。 | 3.9、Fig. S7e |

### 3. CODEX（你说暂不处理）

在你那条消息之前已经做完（35 名患者、预注册方案、实测 CD8⁺ 为结局），作为独立的 3.7 节 + Fig. 7 保留。
本轮只做了两件事：用精确割重跑（主终点不变：患者层面 ρ = −0.53，p = 0.0006）；为回答"简单方法好在哪"
加了两个 post hoc 几何比较量。**如果你不想要 CODEX，删除 3.7 节、Fig. 7、Fig. 8f、Table S8 和摘要第 4 句即可，
其余部分不依赖它**——但它是全文唯一的"实测结局"验证，我建议保留。

### 4. 其它

- BANKSY 式域比较不再需要回 Windows 重跑：`run_51` 在没有 scanpy 的环境里复刻了同一套预处理
  （HVG seurat / scale / ARPACK PCA / k-means），给旧割时 **19/19 张切片与存档结果逐位一致**，再用精确割重算。
- 模拟基准新增几何量"血管到瘤巢距离"（Table 2 新行），不改变任何其它结果。
- 测试：8 个测试文件共 44 个测试全部通过（新增 BANKSY 预处理、合并超几何零分布两个测试）。
- 版本号 2.2.0（pyproject、`sparta/__init__.py`、CITATION.cff、.zenodo.json）。

---

## 二、需要你本人做的事（约 20 分钟）

1. **Title page**：填 ORCID 和院系（两处方括号），Word 另存 PDF。
2. **Zenodo DOI**：发布 v2.2.0 后把 DOI 填进 Title page 的 Code availability 和 CITATION.cff。
   （git 只能在 Windows 本地做，见下。）
3. **用 Word 打开** `Manuscript_IS_anonymised.docx` 目检公式。
4. **在 Windows 上提交并打 tag**：
   ```powershell
   cd D:\sparta
   .\.venv\Scripts\python.exe tests\test_v22_additions.py
   .\.venv\Scripts\python.exe tests\test_is_revision.py
   git add -A sparta scripts tests configs docs docs_is data/ledger.csv data/interim results/validation results/counterfactual results/intervention results/figures/is README.md CLAUDE.md CITATION.cff .zenodo.json pyproject.toml
   git commit -m "Release v2.2.0: exact minimum cut, CODEX validation, replication cohort, simple baselines, intervention maps"
   git tag -a v2.2.0 -m "SPARTA v2.2.0"
   git push; git push origin v2.2.0
   ```
5. 可选：Thrane 队列的 node tables 含每个 spot 约 100 个模型基因的表达；若担心原始数据的再分发许可，
   公开 Zenodo 版本时可以去掉 `data/interim/MEL_THR*.nodes.npz`（审稿用的 ESM_3 保留不影响）。

## 三、稿件现状

- 摘要约 240 词（上限 250）；正文约 8,600 词；8 图 3 表 47 篇文献；ESM_1 有 S1–S14、Tables S1–S10、Figs. S1–S7。
- 正文比 v2.1 长了不少（新增 CODEX、复现队列、简单基线、干预图）。IS 没有硬性字数上限，但如果编辑嫌长，
  最容易压缩的是 2.7（校准方法）和 3.4（结构零模型）的细节，可移到 ESM_1。
- 新增脚本：`run_51`（域比较）、`run_52`（简单基线）、`run_53`（干预位点区域构成）、`scripts/_h5ad_reader.py`。
