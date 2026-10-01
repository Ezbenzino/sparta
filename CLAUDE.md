# CLAUDE.md —— SPARTA 项目规范

> 本文件供所有在本仓库工作的 AI 助手会话遵守。人类协作者也建议先读一遍。
> 修改本文件前请先确认改动确实是项目级的约定，而非一次性的临时决定。

## 项目一句话

SPARTA v2.0：从皮肤肿瘤空间转录组出发，用图论传输量把免疫治疗抵抗的空间屏障
解耦为「T 细胞迁移屏障 `B_cell`」与「抗体传质屏障 `B_mAb`」两个可分离的维度。
全程 CPU，不训练任何深度模型。

---

## 一、硬性约束

**禁止引入的依赖**：`torch` · `torch-geometric` · `transformers` · `geneformer` ·
任何 CUDA / GPU 相关配置。本项目的全部计算是确定性的，不需要也不应该有它们。
如果某个需求看起来需要深度学习，先问：能不能用确定性方法做？大概率能。

**核心模块的依赖边界**：`barrier.py` · `graph.py` · `counterfactual.py` ·
`synthetic.py` 这四个模块**只能依赖 numpy / scipy / networkx**。
不得在其中 import scanpy、anndata、squidpy。理由：
1. 核心数学与数据格式解耦，换平台不用改；
2. 单元测试可在无 scanpy 环境下运行；
3. 零基础读者更容易看懂——输入就是几个数组。

scanpy 相关的代码只能出现在 `signatures.py`、`validate.py` 的部分函数、
`scripts/run_0*.py`，以及各模块末尾明确标注的 `*_from_adata` 适配器中。

---

## 二、目录与文件契约

```
sparta/
├── configs/default.yaml
├── data/{raw,interim,external}/
├── sparta/{io_,loaders,signatures,graph,barrier,counterfactual,validate,viz,synthetic}.py
├── scripts/{setup_check,run_00_demo,run_00b_ingest,run_00c_admission_audit,
│           run_01…run_06,run_07_screen,run_08/09/10_benchmark*,
│           run_11_review_diagnostics,run_12_paper_stats,run_13_benchmark_tools,
│           run_13b_ndomains_scan,run_14_shared_ecm_check,run_15_figure1,
│           run_16_figures,run_17_sensitivity,run_18_table1,run_19_runtime,
│           run_20_pending_figures,run_21_mesh_stats,run_batch}.py
├── scripts/generate_docx.js  Word 文档生成（从正文文本+图表生成 .docx，需 npm install docx）
├── scripts/run_fastpath.py  审查整改的快路径（约 15–25 分钟走到决策岔口）
├── scripts/run_ingest_legacy_cscc.py  批量摄入 GSE144239 第一代 ST（GSM→患者映射写死在代码里）
├── scripts/run_batch_s1s2.py 唯一的长任务：8 张切片重跑 S1+S2（3–5 小时）
├── scripts/fastpath.ps1     纯 ASCII 的一行包装，调 run_fastpath.py
│                            （PowerShell 5.1 按 ANSI 读 .ps1，中文脚本会解析失败）
├── scripts/scratch/         一次性脚本存档（不属于可复现管线，见其 README）
├── tests/{test_barrier,test_counterfactual,test_loaders,test_statistics}.py
├── data/ledger.csv            数据台账（run_00b_ingest 维护）
└── results/{figures,counterfactual,validation}/
```

模块间的数据契约（改动前必须同步更新本文件与 `io_.py`）：

| 模块 | 输入 | 输出 |
|---|---|---|
| M0 `run_00b_ingest` | `data/raw/{任意格式}` | `{sid}.raw.h5ad`、`data/ledger.csv` 一行 |
| M1 `run_01_qc` | `{sid}.raw.h5ad` | `{sid}.qc.h5ad`、`{sid}.admission.json` |
| M2 `run_02_score` | `{sid}.qc.h5ad` | `{sid}.scored.h5ad`（签名在 obs，秩标准化列以 `_n` 结尾） |
| M3 `run_03_graph` | `{sid}.scored.h5ad` | `{sid}.graph.npz`（A、D、source、sink、vessel） |
| M4 `run_04_barrier` | 上两者 | `{sid}.barrier.npz`、`{sid}.mincut.json` |
| M5 `run_05_counterfactual` | 同上 | `results/counterfactual/{sid}.json` |
| M6 `run_06_validate` | 全部切片 | `results/validation/*` |
| M12 `run_12_paper_stats` | 上述全部 JSON | `results/validation/paper_stats.json`（论文统计汇总，只汇总不重算） |
| M13 `run_13_benchmark_tools` | `{sid}.scored.h5ad` + `.mincut.json` | `results/validation/benchmark_tools.json`（与空间域方法的对比） |
| M14 `run_14_shared_ecm_check` | `{sid}.scored.h5ad` + `.graph.npz` | `results/validation/shared_ecm_check.json`（**决策节点**：耦合是组织的性质还是模型自带的） |
| M15 `run_15_figure1` | 无（合成切片） | `results/figures/fig1_framework.png/.pdf` |
| M16 `run_16_figures` | 各汇总 JSON（不重算） | `results/figures/fig{2,3,4,5}_*.png/.pdf` |
| M17 `run_17_sensitivity` | `{sid}.scored.h5ad` + `.graph.npz` | `results/validation/bmab_sensitivity.json`（多切片）+ 逐片热图 |
| M18 `run_18_table1` | `data/ledger.csv` + `admission.json` | `docs/table1_sections.md`（队列概况表） |
| M19 `run_19_runtime` | `{sid}.graph.npz` | `results/validation/runtime_benchmark.json`（SPARTA vs BANKSY vs Squidpy） |
| M20 `run_20_pending_figures` | 各汇总 JSON（不重算） | `results/figures/fig{3,4,5}_*.png/.pdf` + `graphical_abstract.png/.pdf` |
| M21 `run_21_mesh_stats` | `{sid}.scored.h5ad` + `.graph.npz` | `results/validation/mesh_stats.json`（网孔尺寸排除统计，backs 60–65% claim） |

**热图色标的选法**：量级量（如"交联占比 0–100%"）用**单色相 light→dark**；
有阈值的极性量（如"解离潜力，1.0 是 GO 线"）用**双色相 + 中性灰中点**，中点对齐阈值。
两者搞反会让读者在没有极性的地方看出极性——旧的 `_sens_bmab_grid.py` 就是两个都反了。

**出图约定**：色相只承载"哪个屏障"（cell `#2a78d6` / antibody `#eb6834`，已过 CVD 校验），
队列一律用形状 + 分组位置 + 直接标注，不用颜色——颜色被屏障模态占用了，
再拿它编码队列会让读者在两套语义之间来回猜。同一物理量跨面板必须共用色标与轴范围。

**产物文件一律"读回再合并"，禁止整文件覆盖。** `run_06` 跑完 CSCC 就把 MEL 的结果
冲掉过（`decoupling.json` 只剩 MEL、`consistency.json` 只剩 CSCC），
写汇总表时会以为另一个队列没跑过。新写的脚本必须用 `run_06_validate.merge_into` 的模式。

所有磁盘 IO 走 `io_.py`，路径从 config 读，禁止硬编码。

---

## 三、四条工程纪律

1. **逐切片处理，中间产物落盘。** 脚本以单张切片为最小处理单元，处理完保存并释放。
   下游要反复迭代几十次，不落盘就得反复重跑上游。
2. **原始数据只读。** `data/raw/` 下载后不再修改，清洗结果一律写 `data/interim/`。
3. **配置外置。** 参数写 `configs/default.yaml`；出图时把所用 config 复制到
   `results/` 对应目录与图件同存。
4. **随机种子固定并记录。** 用 `io_.set_seed()`；结果元数据用 `io_.stamp_run()` 生成。

---

## 四、参数的三类身份（改参数前必读）

`configs/default.yaml` 中每个参数都标了身份，处理方式完全不同：

- **`[物理锚定]`** —— 有文献值，**直接固定，不参与任何拟合**。
  例：`r_nm: 5.5`（IgG 流体力学半径）、`d0_um: 130`（氧扩散极限）。
  改动它等于改变物理假设，必须有文献依据。
- **`[无监督标定]`** —— 用「同瘤种不同切片间屏障分布一致性」标定，
  **全程不接触任何临床响应标签**。例：`b_ecm`、`lam`、`beta`。
- **`[敏感性]`** —— 必须做网格扫描，在补充材料给出结论随参数变化的热图。
  例：四个分位数阈值、`radius_um`、`kd_eff`、`w1–w3`。

> **信息泄漏红线**：任何用了临床响应标签的参数调整，都会使临床验证失去意义。
> 标定与验证必须彻底分离，且在论文方法部分明确写出「标定使用了什么、
> 没有使用什么」。这是审稿人核查的重点，也是这类研究最常见的隐性造假形式
> （往往是无意的）。

---

## 五、单元测试（不可省略）

```bash
python tests/test_barrier.py          # 6 个测试，无需 pytest
python tests/test_counterfactual.py   # 4 个测试
python tests/test_statistics.py       # 4 个测试（统计口径，2026-08-26 新增）
pytest tests/ -v                       # 有 pytest 时
```

所有测试都写成「pytest 能跑、`python` 直接跑也能跑」的双模式。
新增测试请沿用这个模式（文件末尾的 `if __name__ == "__main__":` 块）。

必须保持全绿的测试及其含义：

| 测试 | 检验什么 | 失败意味着 |
|---|---|---|
| `test_ring_vs_scattered` | 能否区分连续环带与等量散在节点 | 算子测的是组成而非拓扑 |
| `test_ring_with_gap` | 能否感知环带缺口 | 算子已退化为局部统计量 |
| `test_size_dependence` | `B_mAb` 是否随分子尺寸单调 | 尺寸排阻项实现有误 |
| `test_dissociation_is_detectable` | 已知解离场景下能否检出解离 | 无法区分「真的没解离」与「测不出解离」 |
| `test_b_meta_distance_monotonic` | `B_meta` 的几何项 | 图距离实现有误 |
| `test_empty_source_returns_nan` | 空源集的健壮性 | 会静默返回错误值 |
| `test_orientation_always_correct` | 四种原始格式的方向判断 | spot/基因维度弄反，下游全错但不报错 |
| `test_missing_coords_raises` | 缺坐标时显式报错 | 静默产生无意义结果 |

**修改 `barrier.py` 后必须重跑这两个测试文件。**

### 项目的核心风险（改动 barrier.py 前必须知道）

`B_cell` 与 `B_mAb` 共享同一个 ECM 项——它既进了最小割的边容量，也进了扩散的
边电导。因此 ECM 的任何变化会**同时**影响两个屏障，这条耦合通道无法通过统计
校正消除。合成数据上 ECM 解释了 `B_mAb` 方差的 96%，解离潜力比只有 0.04。

含义：本项目的中心主张（两道屏障可解离）只能依靠**只影响一个屏障的因素**
（抗原 / 交联度 / 分子尺寸）来支撑。用 `validate.dissociation_drivers()` 
在建库阶段就量化这一点。详见 README 的「已知的核心风险」。

⚠️ **绝不能为了让解离出现而调 `b_cell.b_ecm` 或 `b_mab.lam`。** 这两个参数的
相对大小直接决定共享通道强度，只能依据文献确定，不能依据结果好看与否确定。
这是信息泄漏的一种变体。

### 统计口径的四条铁律（2026-08-26 投稿前审查后加）

算子写对了，统计口径写错一样会毁掉整篇文章，而且**单元测试抓不到**——
原有的 14 个测试全部检验算子的定性行为，下面四件事它们一件都测不出来。
`tests/test_statistics.py` 就是为守住这四条而写的。

1. **解耦分析必须传 control，且必须报偏相关。**
   `decoupling_stats` 不传 `control` 时 `rho_partial` 是 `NaN`、`controlled` 是 `False`。
   把这种情况下的 `rho` 说成"控制距血管距离后的偏相关"，是事实错误。
2. **逐点 B_cell 只能用 `compute_b_cell_field`。**
   用 `0.5*(ECM_core_n + CAF_n)` 当代理会得到方向相反的结论
   （旧版 `run_06 --dim decoupling` 就是这么写的，产出的负相关是伪结果）。
3. **"解离区占比"必须对照随机期望 `(1-q)^2`。**
   q=0.75 时是 6.25%。实测低于它 = 解离区比随机还少 = 两屏障正相关。
   报一个 5.5% 的解离区占比并称之为"发现"，等于把随机水平当结论。
4. **反事实实验的规模必须相对结构规模定义。**
   S2 的 k 用绝对值时，割集从 76 个节点长到 614 个，同一个 k 只动了 0.5%，
   效应量被稀释到测不出来。用 `counterfactual.s2.k_frac`，并把 `k_over_cut` 记进结果。

同类推广：**任何"我们观察到 X%"的陈述，都要先问 X% 在零假设下是多少。**

### 已被测试抓出来过的三个真 bug（写新代码时留意同类错误）

1. **结合位点屏障只看终点抗原。**
   BSB 的物理本质是抗体在**沿途**被高抗原细胞消耗，深部细胞自己抗原低不低
   不重要——药物在半路就没了。早期实现写成 `kappa[destination]`，
   结果是加不加瘤巢外缘的高抗原环，`B_mAb` 完全不变。
   正确做法：构造有向的"吸收代价图"，用 Dijkstra 求最小累积吸收。
2. **`B_mAb` 用「最少被吸收的最短路径」算透过率。**
   抗体可以**免费绕开**高抗原区域——只要绕路上吸收低，路径再长代价也不增加。
   结果是一个扇区的高抗原完全挡不住深部，与真实扩散物理不符。
   正确做法：解扩散-吸收方程（屏蔽泊松）`(L + diag(k))phi = 0`，
   绕路本身要付扩散阻力，两者自动权衡。k=0 时退化为调和场，
   有效阻抗是其特例，所以"图论传输阻力"的叙事完全保留。
3. **S2 用「随机移除最小割节点」当定向策略。**
   最小割集有几十个节点，随机拿掉 3–5 个只在环带各处打小洞，环带仍闭合。
   必须拿掉空间上**连续的一段弧**，否则会得到与对照无异的结果，
   从而错误地否定拓扑假设。

---

## 六、写代码时的约定

- 每个函数：Python 类型标注 + 中文 docstring，说明输入、输出、单位、
  关键参数含义。`barrier.py` 的三个算子要注释到「零基础者能看懂每一行在算什么」。
- 每个模块文件头部注明：**输入文件 / 输出文件 / 上游模块 / 下游模块**。
- 每个 `run_*.py` 必须有 `--help`，且 `--help` 不应触发重依赖的 import
  （把 `import scanpy` 放在 `parse_args()` 之后）。
- 不要静默失败。缺列、缺基因、源集为空——一律打印告警或返回 `nan`，
  不要用默认值悄悄填过去。下游拿到一堆看似合理其实无意义的数字是最坏的情况。
- 数值边界：`B_mAb` 中与血管不连通的节点在数学上是无穷大，
  实现里赋为有限最大值并在 `reachable` 中标记 `False`，由调用方决定是否纳入统计。

---

## 六之二、git 不能走 Cowork 桥接（2026-08-27 实测）

桥接侧挂载的目录**禁止删除文件**（`rm` 一律 Operation not permitted），
而 git 每次写索引都要先创建 `.git/index.lock`、写完再删掉它。
结果是每条 git 命令都成功创建了锁却删不掉，给下一条留个绊子，仓库直接卡死。

**所有 git 操作（init / add / commit / push）必须在 Windows 本地做。**
桥接侧只能读 `.git`，不要写。发布流程见 `docs/release_checklist.md`。

同类问题：任何"创建临时文件再删除"的工具在挂载目录里都可能卡住，
桥接侧遇到 Operation not permitted 时先想这一条，再想权限。

---

## 七、明确禁止

- 不生成任何"研究结论"或伪造实验数据。合成/占位数据必须标注
  `synthetic` 或 `placeholder`。
- 不自动下载真实数据集，不自动运行耗时计算。
- 不引入深度学习 / GPU 依赖。
- 不改动上述目录结构与文件契约；确需调整先说明理由。
- 不生成超出本项目范围的模块。

---

## 八、当前状态

- 已实现并通过测试（14 个测试）：`barrier` · `graph` · `counterfactual` ·
  `synthetic` · `loaders` · `io_` · `viz` · `validate` · `signatures`，
  `setup_check` · `run_00_demo` · `run_00b_ingest` · `run_01`–`run_06` ·
  `run_07_screen` · `run_batch`。
- `run_07_screen.py --synthetic` 是**早期决策点**，不需要真实数据即可演示。
  它的判定（GO / CAUTION / PIVOT）决定论文叙事走哪条分支，应在拿到
  3–5 张切片后立刻跑，不要等全库建完。
- `run_00_demo.py` 可在无真实数据、无 scanpy 的情况下跑通全流程。
- 2026-08-26 完成投稿前架构审查（`docs/review_for_journal.md`）。**结论已改向 PIVOT**：
  在正确口径下 B_cell 与 B_mAb 控制几何后 8/8 呈正偏相关、解离区占比低于随机期望，
  论文主张改为"两道屏障高度耦合，耦合源是基质而非几何深度"。
  写任何正文之前先读那份报告的 §3 与 §6。
- 已完成的整改：`run_06` 解耦口径、`decoupling_stats` 随机基线、
  `s2` 的 `k_frac`、`run_01` 准入可审计化（`admission_overrides` + 强制放行留痕）、
  `tests/test_statistics.py`、临时脚本收编（`run_12` / `scripts/scratch/`）。
- 待补：`notebooks/` 中的临床队列（维度③）分析示例——不同队列的数据格式
  差异太大，不适合做成统一脚本，因此以 notebook 形式提供更合适。
- 基因集出处已于 2026-08-27 补齐，见 `docs/gene_set_references.md`。
  标注分两级：`[ref]` 有明确文献可直接引用；`[canon]` 是通用谱系标记、
  没有唯一原始出处，只能写成 "canonical lineage markers"。
  **不要给 `[canon]` 硬安一篇文献** —— 一个没核实过的引用比不给更糟。
- 2026-08-28 收尾：run_13 的 squidpy 臂修复（第一代 ST 空图根因是 radius
  小于点间距，现复用 SPARTA 图，19/19 有值）；`run_13b_ndomains_scan.py`
  补 n_domains 敏感性扫描（结论：precision 全网格 3–23%，平台落差对 k
  稳健、边界占比与 Visium 持平时仍在）；cut_band 扩到 19 张（三分法预测
  S1 显著性 17/19 正确）；Supp S3 图重出；Abstract 三处数字按产物校正
  （S2 患者层 5/7、切断共享 ECM 15/19、S2 范围口径）。
- `signatures.py` 的 Hypoxia 目前是占位集合，正式分析前需用
  `--hypoxia-gmt` 传入 MSigDB HALLMARK_HYPOXIA。
