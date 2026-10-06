# SPARTA 中文说明

SPARTA 是一个用于空间组学数据的图结构方法学框架，用于表征空间组织结构和局部微环境模式。它为每张组织切片构建空间图，并结合结构算子、空间零模型、校准分析和扰动分析描述组织场的空间组织方式。

## 定位声明

SPARTA 不作为治疗响应、预后、治疗疗效或临床效用方法。分子尺度场是计算机内的结构基准，不代表实测药物暴露，也不用于预测治疗结果。

## 目录说明

| 路径 | 内容 |
|---|---|
| `sparta/` | 正式 Python 包 |
| `scripts/` | 正式管线和绘图脚本 |
| `configs/` | 默认参数和准入规则 |
| `docs_is/` | 稿件、补充材料、投稿信和参考文献源文件 |
| `docs/` | 队列修正记录和审稿回应预案 |
| `results/` | 验证结果和稿件图件 |
| `exploratory/` | 探索性脚本，不进入正式投稿代码归档 |
| `legacy/` | 已废弃但保留溯源价值的脚本，禁止作为当前管线运行 |
| `tests/` | 正确性和可复现性测试 |

## 队列概况

| 队列 | 切片数 | 患者/分析单位 | 构成 |
|---|---:|---:|---|
| 主队列 | 19 | 8 位患者 | 15 张 cSCC、4 张黑色素瘤 |
| 复现队列 | 8 | 4 位患者 | 第一代 ST 黑色素瘤淋巴结转移 |
| 2026 扩展队列 | 50 | 32 位可识别患者，另有 11 张患者关系不可得 cSCC | 34 张 Visium、16 张 Slide-seqV2 |
| 合计 | 77 | 44 位可识别患者及 11 个保守切片级单位 | cSCC、原发和转移黑色素瘤 |

扩展队列按疾病分为 13 张 cSCC、6 张原发黑色素瘤、31 张转移黑色素瘤。GSE289745 的 11 张 cSCC 缺少患者 ID，不计为 distinct patients。

## 扩展数据集

| GEO | 平台 | 纳入切片 | 分层 |
|---|---|---:|---|
| GSE289745 | Visium | 11 | cSCC，患者关系不可得 |
| GSE321832 | Visium | 2 | 皮肤 SCC |
| GSE300445 | Visium | 4 | 原发黑色素瘤 |
| GSE316760 | Visium | 2 | 原发黑色素瘤 |
| GSE320041 | Visium | 15 | 转移黑色素瘤 |
| GSE200278 | Slide-seqV2 | 16 | 转移黑色素瘤 |

GSE320041 的 WU1457、WU2109 因属于 Visium HD 16 µm 被排除；GSE321832 的 oral1、oral2、lung1 因不是皮肤 SCC 被排除。

## 核心结果

- 扩展队列 50 张中 48 张为正，中位 partial Spearman ρ = 0.273；患者/分析单位 pooled estimate = 0.28（95% CI 0.23–0.32）。
- 疾病分层 pooled estimates：cSCC 0.37；原发黑色素瘤 0.15；转移黑色素瘤 0.26。
- 平台中位值：Visium 0.320；Slide-seqV2 0.193。
- 扩展校准假阳性率：spectral 0.042；normal score 0.055；point level 0.314。
- 构造零模型和共享 ECM 消融说明部分关联来自图构造和共享基质结构，但并非全部。

以上均为结构验证结果，不是临床预测。

## 运行环境

请使用项目虚拟环境：

```powershell
D:\sparta\.venv\Scripts\python.exe
```

不要直接使用系统 Python；系统 Python 可能缺少项目依赖。

## QC 和平台参数

- Visium spot：UMI ≥ 500。
- Slide-seqV2 bead：UMI ≥ 100。
- 基因：至少在 3 个 retained locations 检出。
- 线粒体比例：若可测则 ≤ 0.20；若矩阵不含线粒体基因，标记为 not assessable。
- 扩展 Visium 切片：≥300 retained spots 且 retained median UMI ≥1500。
- 扩展 Slide-seqV2 切片：≥1000 retained beads 且 retained median UMI ≥100。

Slide-seqV2 先在原生 bead 层面 QC，再将计数汇总到固定 50 µm 网格；正式图半径为 100 µm，不能复用 Visium 的 150 µm。

## 主要复现入口

```powershell
D:\sparta\.venv\Scripts\python.exe scripts\run_54_correct_melanoma_mapping.py
D:\sparta\.venv\Scripts\python.exe scripts\run_55_extension_ingest.py
D:\sparta\.venv\Scripts\python.exe scripts\run_56_extension_cohort.py
D:\sparta\.venv\Scripts\python.exe scripts\run_57_extension_supplement.py
D:\sparta\.venv\Scripts\python.exe scripts\is_figures\facts.py
D:\sparta\.venv\Scripts\python.exe scripts\is_figures\fig4_association.py
D:\sparta\.venv\Scripts\python.exe scripts\is_figures\fig8_baselines.py
D:\sparta\.venv\Scripts\python.exe scripts\is_figures\figS1_parameter_heatmaps.py
```

## 队列变更日志

- 2026-10-05：修正主队列黑色素瘤患者映射。MEL01/MEL04 属于患者 A，MEL02/MEL03 属于患者 B；主队列是 19 张切片、8 位患者。详见 `docs/melanoma_mapping_correction.md`。
- 2026 扩展：新增 6 个 GEO、50 张公共切片。详见 `data/external/extension_2026/README.md`。
- 稿件定位：删除治疗响应、预后、疗效和临床效用表述，并将局限前移到讨论开头。
- 参数敏感性：ξ₀/β/λ 热图位于 Online Resource 1 的 Section S1、Table S1 和 Fig. S1。

## License

MIT。
