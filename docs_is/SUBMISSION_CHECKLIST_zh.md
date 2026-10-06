# SPARTA → *Interdisciplinary Sciences: Computational Life Sciences* 投稿清单（v2.2.2）

生成日期：2026-10-05；2026-10-06 修订重建（新增 Fig. 2e–g / Table 3 分子尺度 ground-truth 实验 run_59；正文压缩至约 9,700 词）；2026-10-06 第二次修订（统一论文标题为匿名稿/cover letter 版本；通讯邮箱统一为学校邮箱 1109240426@hmc.edu.cn；ESM_3 README 版本同步为 v2.2.1；字数统一口径：正文 9,609 词，Title page 显示取整后的 9,610）；2026-10-06 投稿前审查整改（Title page 与匿名稿 Data availability / Ethics 补全 6 个扩展队列 GEO 号；Fig. 4、Fig. 7 由 208/210 mm 重排到 174 mm 版心内；复现协议补记 12 患者合并分析的口径变化；图形摘要按 Springer 规格 13.0×5.2 cm（2.5:1）重绘为图形为主、文字最少）。期刊投稿须知：<https://link.springer.com/journal/12539/submission-guidelines>
（双盲评审；摘要 150–250 词；关键词 4–6 个；正文 Word 或 LaTeX；线图 EPS 或 ≥1200 dpi TIFF、组合图 ≥600 dpi；
图宽 84/174 mm；补充材料称 Online Resource，文件名按 ESM_1、ESM_2… 编号；声明放在单独的 Title page；LLM 使用须在 Methods 中说明。）

---

## 1 本文件夹（`D:\sparta\submission_IS\`）里有什么

| 文件 | Editorial Manager 上传类型 | 说明 |
|---|---|---|
| `Cover_letter_IS.pdf`（`.docx` 可编辑） | Cover Letter | 2 页，致两位主编 |
| `Title_page_IS.docx`（`.pdf` 供核对） | Title Page | 作者信息 + 全部 Declarations；**审稿人看不到** |
| `Manuscript_IS_anonymised.docx` | Manuscript（匿名） | 行号、页码、图表嵌入；摘要 234 词，正文约 9,610 词（与 Title page 同一口径），8 图 4 表 44 篇文献 |
| `Manuscript_IS_anonymised.pdf` | （若系统允许，作为附加的审稿 PDF） | 同内容版本（30 页；由匿名 DOCX 经 LibreOffice 渲染，如本机有 pandoc+pdflatex 可用 `build_is.py` 重建 LaTeX 版） |
| `Figures/Fig1_overview.tif`、`Fig7_codex.tif`、`Fig8_baselines.tif` | Figure 1、7、8 | 组合图（含切片栅格或栅格化散点），600 dpi，RGB，174 mm |
| `Figures/eps/Fig2–Fig6 *.eps` | Figure 2–6 | 矢量线图（期刊首选格式） |
| `Figures/line_art_1200dpi/*.tif` | 备用 | 若系统不收 EPS，就上传这 5 个 1200 dpi TIFF |
| `Figures/pdf/*.pdf` | 不上传 | 8 张图的矢量 PDF，生产阶段备用 |
| `Graphical_Abstract.pdf`（`.tif` 为 500 dpi 位图备份） | Graphical Abstract | 13.0 × 5.2 cm（2.5:1，Springer 图形摘要规格），以图形为主、文字最少 |
| `ESM_1.pdf` | Supplementary Material | Online Resource 1（23 页）：S1–S16，Tables S1–S13，Figs. S1–S9 |
| `ESM_2.xlsx` | Supplementary Material | Online Resource 2：逐切片 / 逐核心结果（含干预图、简单基线、区域构成、CODEX 逐核心） |
| `ESM_3.zip` | Supplementary Material | Online Resource 3：**匿名**代码（含 run_59）+ 30 张切片的 node tables + 全部结果文件（run_59 逐组织数值在 `results/validation/mab_ground_truth_tissues.csv`；含 README_REVIEWERS.md、SHA-256 清单） |
| `suggested_reviewers.md` | Suggested Reviewers（可选） | 需您核实 |
| `_superseded/` | 不上传 | 旧版本文件（v2.1 的 Fig8 等），只为留痕 |

## 2 上传前必须由您完成（约 20 分钟）

1. ~~Title page 两处方括号~~ 已完成：ORCID `0009-0002-9940-5317` 与院系已填入，Title page 已重建（2026-10-06）。
2. ~~Code availability 中的 DOI~~ 已完成：改为引用 **Zenodo concept DOI = `10.5281/zenodo.23086430`**（永远解析到该记录的最新版本，v2.2.2 发布后自动指向 v2.2.2，**以后换版本不用再改稿件**）。
3. **核实声明**：Funding、Competing interests、Author contributions、AI 工具使用说明（正文 2.17 节与 Title page 一致）。
4. **核实推荐审稿人**：期刊要求推荐人完全独立（近 3 年无合作、不同单位、无师生关系），并填写机构邮箱。2026-10-06 已核实 4 人的现任职务与机构邮箱（见更新后的 `suggested_reviewers.md`），仅需作者本人最终确认无合作关系后即可填入系统。
5. **用 Word 打开** `Manuscript_IS_anonymised.docx`，确认公式（Eq. 1–3 和行内符号，尤其是变量、上下标和希腊字母）正常显示。
6. ~~同步 Zenodo 记录元数据~~ 已完成：v2.2.2 记录直接继承 `.zenodo.json`，title/description/version 均已核对。
7. **GitHub release 说明核对**：v2.2.1 release 说明中的正文字数写的是未取整的 9,609，与 Title page 的 9,610（取整到十位）同源不冲突；如希望完全一致，可在 release 页面把该数字改为约 9,610，并顺带确认 release 说明引用的论文标题与新标题一致。
8. ~~发布 v2.2.2~~ 已完成，见第 4 节。

## 3 Editorial Manager 填写要点

- Article type：Original Research Article
- Title、Abstract、Keywords：从稿件复制（关键词：Spatial transcriptomics; Multiplexed imaging; Minimum cut; Graph Laplacian; Spatial null model; Tumour microenvironment）
- 作者信息只填在系统表单和 Title page 中；匿名正文、ESM_1/2/3 中没有姓名、单位、邮箱或仓库链接（ESM_3 已自动扫描确认）。
- Data availability / Code availability：与 Title page 一致。

## 4 发布

**已完成：v2.2.1 → DOI 10.5281/zenodo.23181026**（GitHub release <https://github.com/Ezbenzino/sparta/releases/tag/v2.2.1>）。

**v2.2.2 已发布（2026-10-06，作者授权后由助手完成）**：commit `83617ff` → tag `v2.2.2` → GitHub release <https://github.com/Ezbenzino/sparta/releases/tag/v2.2.2> → Zenodo 自动归档。

- 版本 DOI `10.5281/zenodo.23187800`；**concept DOI `10.5281/zenodo.23086430`（稿件引用的是它，永远指向最新版，换版本无需改稿）**。
- Zenodo 记录的 title/description 直接继承仓库根目录 `.zenodo.json`，已核对为 v2.2.2 与新标题，**网页端无需再手改**（原第 2 节第 6 条作废）。

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
