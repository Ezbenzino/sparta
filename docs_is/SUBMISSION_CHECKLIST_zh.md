# SPARTA → *Interdisciplinary Sciences: Computational Life Sciences* 投稿清单（v2.2.1）

生成日期：2026-10-05；2026-10-06 修订重建（新增 Fig. 2e–g / Table 3 分子尺度 ground-truth 实验 run_59；正文压缩至约 9,700 词）。期刊投稿须知：<https://link.springer.com/journal/12539/submission-guidelines>
（双盲评审；摘要 150–250 词；关键词 4–6 个；正文 Word 或 LaTeX；线图 EPS 或 ≥1200 dpi TIFF、组合图 ≥600 dpi；
图宽 84/174 mm；补充材料称 Online Resource，文件名按 ESM_1、ESM_2… 编号；声明放在单独的 Title page；LLM 使用须在 Methods 中说明。）

---

## 1 本文件夹（`D:\sparta\submission_IS\`）里有什么

| 文件 | Editorial Manager 上传类型 | 说明 |
|---|---|---|
| `Cover_letter_IS.pdf`（`.docx` 可编辑） | Cover Letter | 1 页，致两位主编 |
| `Title_page_IS.docx`（`.pdf` 供核对） | Title Page | 作者信息 + 全部 Declarations；**审稿人看不到** |
| `Manuscript_IS_anonymised.docx` | Manuscript（匿名） | 行号、页码、图表嵌入；摘要 234 词，正文约 9,690 词，8 图 4 表 44 篇文献 |
| `Manuscript_IS_anonymised.pdf` | （若系统允许，作为附加的审稿 PDF） | LaTeX 排版的同内容版本（27 页） |
| `Figures/Fig1_overview.tif`、`Fig7_codex.tif`、`Fig8_baselines.tif` | Figure 1、7、8 | 组合图（含切片栅格或栅格化散点），600 dpi，RGB，174 mm |
| `Figures/eps/Fig2–Fig6 *.eps` | Figure 2–6 | 矢量线图（期刊首选格式） |
| `Figures/line_art_1200dpi/*.tif` | 备用 | 若系统不收 EPS，就上传这 5 个 1200 dpi TIFF |
| `Figures/pdf/*.pdf` | 不上传 | 8 张图的矢量 PDF，生产阶段备用 |
| `ESM_1.pdf` | Supplementary Material | Online Resource 1（23 页）：S1–S14，Tables S1–S13，Figs. S1–S9 |
| `ESM_2.xlsx` | Supplementary Material | Online Resource 2：逐切片 / 逐核心结果（含干预图、简单基线、区域构成、CODEX 逐核心） |
| `ESM_3.zip` | Supplementary Material | Online Resource 3：**匿名**代码（含 run_59）+ 30 张切片的 node tables + 全部结果文件（run_59 逐组织数值在 `results/validation/mab_ground_truth_tissues.csv`；含 README_REVIEWERS.md、SHA-256 清单） |
| `suggested_reviewers.md` | Suggested Reviewers（可选） | 需您核实 |
| `_superseded/` | 不上传 | 旧版本文件（v2.1 的 Fig8 等），只为留痕 |

## 2 上传前必须由您完成（约 20 分钟）

1. ~~Title page 两处方括号~~ 已完成：ORCID `0009-0002-9940-5317` 与院系已填入，Title page 已重建（2026-10-06）。
2. ~~Code availability 中的 DOI~~ 已完成：**v2.2.1 已发布，DOI = 10.5281/zenodo.23181026**（GitHub release v2.2.1，Zenodo 自动生成）。
3. **核实声明**：Funding、Competing interests、Author contributions、AI 工具使用说明（正文 2.17 节与 Title page 一致）。
4. **核实推荐审稿人**。
5. **用 Word 打开** `Manuscript_IS_anonymised.docx`，确认公式（Eq. 1–3 和行内符号）正常显示。

## 3 Editorial Manager 填写要点

- Article type：Original Research Article
- Title、Abstract、Keywords：从稿件复制（关键词：Spatial transcriptomics; Multiplexed imaging; Minimum cut; Graph Laplacian; Spatial null model; Tumour microenvironment）
- 作者信息只填在系统表单和 Title page 中；匿名正文、ESM_1/2/3 中没有姓名、单位、邮箱或仓库链接（ESM_3 已自动扫描确认）。
- Data availability / Code availability：与 Title page 一致。

## 4 发布（已完成 2026-10-06）：v2.2.1 → DOI 10.5281/zenodo.23181026

GitHub release v2.2.1: https://github.com/Ezbenzino/sparta/releases/tag/v2.2.1
（历史记录：以下为 v2.2.0/v2.2.1 发布时执行的步骤，仅留档。）

```powershell
cd D:\sparta
.\.venv\Scripts\python.exe tests\test_v22_additions.py
.\.venv\Scripts\python.exe tests\test_is_revision.py
.\.venv\Scripts\python.exe tests\test_barrier.py
git status
git add -A sparta scripts tests configs docs docs_is data/ledger.csv data/interim results/validation results/counterfactual results/intervention results/figures/is README.md CLAUDE.md CITATION.cff .zenodo.json pyproject.toml
git status
git commit -m "Release v2.2.0: exact minimum cut, CODEX validation, replication cohort, simple baselines, intervention maps"
git tag -a v2.2.0 -m "SPARTA v2.2.0 (Interdisciplinary Sciences submission)"
git push
git push origin v2.2.0
```

`git status` 第二次输出里若出现 `.h5ad`、`downloads/` 或 `submission_IS/`，说明 `.gitignore` 没生效，先停下来。
Zenodo 的步骤与 v2.1.0 相同（GitHub → Releases → 选 tag `v2.2.0` → Publish，几分钟后生成 DOI）。

## 5 v2.2 相对 v2.1 的主要变化（详见 `docs/v22_changelog.md`）

| 变化 | 位置 |
|---|---|
| 最小割改为精确整数割（旧浮点割在 14/30 张图上不是最小割）；依赖割几何的结果全部重算 | 2.4 节、ESM_1 S13 |
| CODEX 结直肠癌多重成像验证：35 名患者，割预测实测 CD8⁺ 核心耗竭（患者层面 ρ = −0.53） | 2.12、3.7 节，Fig. 7，Table S8 |
| Thrane 2018 黑色素瘤独立复现队列：8 张切片 / 4 名患者，预设复现标准达标 | 2.1、3.3、3.4 节，Table S7 |
| 两算子 vs 三个简单空间摘要（基质密度、距肿瘤距离、邻域富集） | 2.14、3.8 节，Fig. 8，Table S10 |
| 干预图（细胞屏障局灶、抗体屏障弥散；两图几乎不重叠；高影响位点避开免疫区） | 2.13、3.9 节，Fig. S7，Table S9 |
