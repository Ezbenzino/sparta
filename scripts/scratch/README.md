# scratch/ —— 一次性脚本存档（不属于可复现管线）

这里的脚本是 2026-08 建库与探索阶段的一次性产物。它们：

- **硬编码了 `D:\sparta\...` 绝对路径**，换一台机器就跑不了；
- 不走 `sparta/io_.py` 的路径契约，违反 `CLAUDE.md` 第二节；
- 但其中几个曾经是论文数字的**唯一**来源。

2026-08-26 的投稿前审查把这些统计收编进了正式管线：

| 原脚本 | 现在用什么 |
|---|---|
| `_fdr_correction.py`、`_s2_summary.py`、`_summarize_s2.py`、`_cross_cohort_summary.py` | `scripts/run_12_paper_stats.py` |
| `_sens_bmab.py`、`_sens_bmab_grid.py` | **已收编** → `scripts/run_17_sensitivity.py`（多切片 + 窗口显式化 + 色标改对） |
| `_cut_band_analysis.py` | 仍在此处，待收编 |
| `_icb_*.py` | 维度③的一次性分析，按队列格式而定，保持在此处即可，但须在论文中标明 |
| `_replot_en*.py`、`_screen_dual_plot.py`、`_benchmark_ext_plot.py` | 出图脚本，建议合并为 `run_14_figures.py` |
| `_dbg_*.py`、`_peek_*.py`、`_check_*.py` | 调试用，投稿前可随仓库保留但不必整理 |

**投稿时的处理**：公开仓库里保留本目录并保留这份说明即可——
审稿人看到的是"探索过程留痕 + 正式管线可复现"，这比假装从来没有过临时脚本更可信。
但**论文正文引用的每一个数字，都必须能由 `run_01 … run_13` 复现**。
