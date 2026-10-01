# SPARTA v2.0

**SPA**tial **R**esistance to **T**herapeutic **A**gents ——
皮肤恶性肿瘤免疫治疗抵抗的双重空间屏障分析框架。

> **核心主张**：免疫治疗抵抗的空间维度不止一层。
> "T 细胞进不去"（细胞迁移屏障 `B_cell`）与"抗 PD-1 抗体自身进不去"
> （大分子传质屏障 `B_mAb`）是两个物理机制不同、在真实肿瘤中可以解离的维度。
> 空间组学界系统性忽视了后者。

全程 **CPU 运行，不需要 GPU**，不训练任何深度模型。核心算子只依赖
`numpy` · `scipy` · `networkx`。

---

## 目录

1. [五分钟上手](#五分钟上手)
2. [这个框架在算什么](#这个框架在算什么)
3. [安装](#安装)
4. [真实数据流程](#真实数据流程)
5. [三个必须做的实验](#三个必须做的实验)
6. [常见问题](#常见问题)
7. [不能跳过的三件事](#不能跳过的三件事)
8. [可选升级](#可选升级)
9. [⚠️ 已知的核心风险](#️-已知的核心风险请在开始真实分析前读完)　← **务必读**
10. [接下来做什么](#接下来做什么按顺序每步都有对应脚本)　← **行动清单**

---

## 五分钟上手

不需要任何真实数据，先确认环境装对了、算子行为正确：

```bash
conda env create -f environment.yml
conda activate sparta
pip install -e .

python scripts/run_00_demo.py --fast     # 在合成数据上跑通全流程
pytest tests/ -v                          # 或：python tests/test_barrier.py
```

`run_00_demo.py` 会在一张人造切片上依次演示八个步骤，并把结果画进
`results/figures/demo_overview.png`。**在动真实数据之前，请把这个脚本
从头到尾读一遍。** 它是理解整个框架最快的路径。

参考输出（合成切片，25×25 网格）：

```
步骤 2｜B_cell：源汇最小割
  环带 B_cell = 155.85   (最大流 0.006)
  散在 B_cell = 0.0955   (最大流 10.475)
  -> 比值 1632x。同样多的阻力细胞，连成一圈就挡得住，散开就挡不住。

步骤 3｜B_mAb：扩散-吸收方程
  瘤巢核心 B_mAb（IgG,  r=5.5nm）=  27.631
  瘤巢核心 B_mAb（小分子, r=0.5nm）=  6.888
  -> 同一张图、同一道基质，IgG 被挡在外面，小分子畅通无阻。

步骤 5｜驱动分解（本项目最重要的诊断）
  共享 ECM 通道 96.2%｜抗原 BSB 3.8%｜解离潜力比 0.040
  -> 见下方「已知的核心风险」
```

---

## 这个框架在算什么

| 分量 | 挡的是什么 | 物理机制 | 图论工具 |
|---|---|---|---|
| `B_cell` | CD8⁺ T 细胞（约 10 μm，主动迁移） | ECM 网孔小于细胞直径、CAF 环带连续性 | **源汇最小割** |
| `B_mAb` | 抗 PD-1 抗体（约 5.5 nm，被动扩散） | 扩散受限 + 尺寸排阻 + 结合位点屏障消耗 | **扩散-吸收方程**（屏蔽泊松；k=0 时退化为有效阻抗） |
| `B_meta` | 小分子（联合化疗场景） | 缺氧静止、外排泵 | **血管图扩散距离** × 代谢状态 |

三者共享同一张空间图，但用不同的边权语义与不同的算子。

**为什么用图论而不是邻域统计**：屏障本质上是拓扑性质。一圈连续闭合的
成纤维细胞构成封锁线，同样数量的散在成纤维细胞不构成封锁线——邻域富集、
Ripley's K、共定位分析对这个区别完全不敏感，因为它们都是局部统计量。
最小割精确刻画"最窄封锁截面的通过容量"；扩散-吸收方程 `(L+diag(k))φ=0`
精确刻画"带吸收的扩散过程"，其中 k=0 的特例就是电阻网络中的有效阻抗。
这不是比喻，是同一个数学对象。

**最小割还额外给了一条几何**：割集本身可以画在 H&E 上，就是一条能与
病理形态直接对照的"封锁线"。这是本项目最有说服力的图件。

---

## 安装

```bash
conda env create -f environment.yml
conda activate sparta
pip install -e .
```

不用 conda 也可以：

```bash
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[all]"
```

**最小安装**（只跑核心算子与单元测试，不处理真实数据）：

```bash
pip install numpy scipy networkx matplotlib pyyaml
python tests/test_barrier.py
```

核心模块（`barrier` / `graph` / `counterfactual` / `synthetic`）不依赖 scanpy，
所以即使 scanpy 装不上，算子验证与合成演示照常可跑。

**WSL 下中文显示为方块**：`sudo apt install fonts-noto-cjk`，
然后 `python -c "import matplotlib; matplotlib.font_manager._load_fontmanager(try_read_cache=False)"`
重建字体缓存。或者调用 `sparta.viz.setup_cjk_font()` 看它报告找到了什么字体。

---

## 真实数据流程

> **当前真实数据状态（2026-08）**：已在两个独立队列跑通 M1–M7 全流程。
> - 队列①转移性黑色素瘤：GSE250636（LMD 研究）4 张颅外转移 Visium（MEL01–04）
> - 队列②原发皮肤鳞癌：GSE144239（Ji et al. 2020 Cell）4 张原发 Visium（CSCC01–04，P4/P6 患者）
> - 双队列决策点均判定 **GO**（中位解离潜力比 58.9 / 41.9），交联主导（96–98%）
> - 详尽的证据汇总见 [`docs/paper_skeleton.md`](docs/paper_skeleton.md)

```
data/raw/{slide_id}/            原始下载（只读）
        │
        ▼  run_01_qc.py         质控 + 准入核查 C1–C7
{slide_id}.qc.h5ad
        │
        ▼  run_02_score.py      签名打分 + 秩标准化
{slide_id}.scored.h5ad
        │
        ▼  run_03_graph.py      空间邻接图 + 源汇定义
{slide_id}.graph.npz
        │
        ▼  run_04_barrier.py    三分量屏障 + 封锁线叠加图
{slide_id}.barrier.npz + .mincut.json
        │
        ├──▶ run_05_counterfactual.py   S1 / S2 / S3
        └──▶ run_06_validate.py         四维验证
```

一张切片的完整流程：

```bash
python scripts/run_01_qc.py  --slide MEL01 --input data/raw/MEL01 \
                             --platform visium --has-he --treatment-known
python scripts/run_02_score.py --slide MEL01 --tumor-type melanoma --plot
python scripts/run_03_graph.py --slide MEL01
python scripts/run_04_barrier.py --slide MEL01
python scripts/run_05_counterfactual.py --slide MEL01
```

每个脚本都有 `--help`。所有参数在 `configs/default.yaml`，代码里不硬编码。

### 准入核查 C1–C7（`run_01` 会逐条打印）

| | 条件 | 不合格的后果 |
|---|---|---|
| C1 | 有 count 矩阵 | 无法分析 |
| C2 | 有空间坐标 | 无法建图，这是命门 |
| C3 | 有配对 H&E | 无法做形态学验证与封锁线叠加图 |
| C4 | spot 数达标（Visium ≥1000，第一代 ST ≥300） | 图过小，最小割不稳定 |
| C5 | 中位 UMI ≥ 1500 | 签名打分不可靠 |
| **C6** | **存在可识别的内皮/血管信号** | **`B_mAb` 与 `B_meta` 都算不出来** |
| C7 | 已知治疗状态与原发/转移 | 治疗后切片的屏障结构已改变 |

**C6 是最容易被忽略、也最致命的一条。** 很多人做到第五个月建图时才发现
切片里内皮信号极弱（第一代 ST 平台尤其常见），此时数据工作已经白做两个月。
`run_01` 会在质控之前就查这一条，不合格的切片直接标记剔除。

降级路径：内皮信号确实不可用时，可以启用"组织最外圈作为血管代理源"
（`configs/default.yaml` 的 `source_sink.use_border_fallback`）。
**一旦启用，必须在论文方法部分写明这一替代及其局限。**

---

## 三个必须做的实验

`run_05_counterfactual.py` 里的三个实验成本极低（CPU 分钟级），
但分别封堵审稿人最可能提出的三个致命质疑。**不要拖到最后做。**

### S1 空间重排对照
> "你的分数是不是只是细胞类型比例的复杂重写？"

保持每个 spot 的全部分数不变、只随机置换空间位置，重算屏障构建零分布。
组成完全没变，只是把细胞随机搬了家——如果屏障消失了，说明它确实来自
空间排布。合成图上 z ≈ 27，p < 0.01。

### S2 环带断裂
> "为什么非要用空间数据？"

在最小割上开一段**连续缺口**，对比移除同样多但**分散**的屏障物质。
材料相同、数量相同，只差"是否连续"。

⚠️ 这里有一个容易写错的地方：随机拿掉几个最小割节点是**不够**的——
最小割集通常有几十个节点，随机拿掉 3–5 个只会在环带各处打几个小洞，
环带整体仍然闭合。必须拿掉空间上连续的一段弧。
（这个坑是在合成图上跑实验时才暴露出来的，按"随机移除"实现会得到与
对照无异的结果，从而错误地否定拓扑假设。）

### S3 分子尺寸扫描
> "两个屏障凭什么说是解耦的？"

把 `r_nm` 从 0.5（小分子）扫到 10（大分子），观察屏障场如何变化。
改一个参数重跑，但它把整篇文章的物理主张可视化成了一条曲线。

---

## 常见问题

**Q：`run_03` 报告"连通分量 5 个"，要紧吗？**
要紧。图碎成多块意味着有效阻抗与图距离在跨分量时是无穷大。增大
`graph.radius_um`（Visium 试 150→180，第一代 ST 试 250→300），
或者检查是不是坐标单位换算错了。

**Q：`run_04` 说"源集与汇集有重叠节点"。**
`source_sink` 的四个分位数阈值设得太宽松，导致同一个 spot 既算免疫入口
又算瘤巢核心。调高 `q_malig` 或 `q_core`。

**Q：`B_mAb` 里有很多 nan / 极大值。**
那些节点与血管不连通。`compute_b_mab` 会把它们赋为有限节点中的最大值，
并在 `reachable` 里标记为 `False`。统计时请按这个掩码筛选，
**不要把它们当正常值**。

**Q：签名分数画出来空间分布很奇怪。**
先怀疑基因名版本（Ensembl ID vs symbol）。`run_02` 会打印每个签名匹配到
几个基因，匹配率低于一半就有问题。

**Q：`pytest` 装不上 / 不想装。**
所有测试文件都能直接跑：`python tests/test_barrier.py`。

---

## 不能跳过的三件事

这个项目能不能成，取决于三件事。其余部分即便打折扣，工作仍然成立。

**1. 准入清单严格执行。** 尤其 C6。在建库阶段花两天核查，能省后面两个月。

**2. 合成图单元测试跑绿。**
```bash
python tests/test_barrier.py          # 6 个测试
python tests/test_counterfactual.py   # 4 个测试
```
屏障算子写错了，在真实数据上几乎不可能被察觉——你会得到一堆看起来合理的
数字，然后基于它们写完整篇论文。合成图是唯一能提前发现这类错误的手段。

其中 `test_dissociation_is_detectable` 会打印**驱动分解诊断**——它告诉你
两个屏障有多少变异来自共享的 ECM、多少来自只影响抗体的因素。
这个数字决定了本课题的中心主张能否成立，务必看懂（见下方「已知的核心风险」）。

**3. S1 与 S2 必须做。** 见上一节。

---

## 可选升级

只有在主线全部跑通、四维验证完成、且仍有时间余量时才考虑。
启动前先问：主线的验证是否已经全部完成？没有就不要启动。

| 升级 | 带来什么 | 成本 | 建议 |
|---|---|---|---|
| 病理基础模型（UNI / CONCH）做形态学特征 | 把维度②从传统特征升级为基础模型特征 | 中 | **性价比最高** |
| 扩展到 ADC 场景 | 把 `Ag_target` 换成 ERBB2 / TACSTD2，公式不用改 | 极低 | 讨论部分的一个小节 |
| 细胞类型反卷积（RCTD / cell2location） | 更严格的细胞比例；可升级为异质图 | 中高 | 审稿人质疑签名打分粗糙时的回应 |
| 扩展到基底细胞癌 | 增加样本量与瘤种覆盖 | 中 | 某瘤种切片不足时的补位首选 |
| Geneformer 节点嵌入 | 若能证明跨切片稳健性提升可作方法学加分 | 高 | **性价比最低**，谨慎 |

> 可选升级的诱惑在于它们看起来能让文章"更高级"。但一篇主线扎实、
> 反事实实验完备、临床验证清晰的文章，远胜过一篇技术堆砌但每处都
> 浅尝辄止的文章。

---

## ⚠️ 已知的核心风险（请在开始真实分析前读完）

这是在构建过程中用合成数据发现的，**不是 bug，而是关于框架本身的结构性发现**。
它直接关系到本课题的中心主张能否成立，建议在投入真实数据之前先看懂。

### 两个屏障之间有两条耦合通道，不是一条

| 通道 | 来源 | 能否校正 |
|---|---|---|
| **共享几何** | `B_cell` 与 `B_mAb` 都是"从血管出发的累积代价"，离血管越远两者都越大 | **能** —— `decoupling_stats(control=d_vessel)` |
| **共享 ECM** | 同一个 core matrisome 分数既进了最小割的边容量（细胞过不去），也进了扩散的边电导（大分子过不去） | **不能** —— 这是框架的结构性质：致密基质本来就同时挡两者 |

第一条通道很强：合成数据上未校正的 Spearman ρ 可达 0.95，校正后降到约 0.35。
**论文里必须报告校正后的偏相关**，否则等于在测深度而不是测机制。

第二条通道是真正的风险。合成实验里，一个 ECM 缺口在放 T 细胞进来的同时，
也让抗体更容易扩散进来——两个屏障一起降低，看不到解离。

### 唯一能产生解离的，是只影响一个屏障的因素

- 只影响 `B_mAb`：**靶抗原**（结合位点屏障，T 细胞不受抗原表达影响）、
  **交联度与分子尺寸**（网孔排阻，对 10 μm 的细胞不是同一个机制）
- 只影响 `B_cell`：趋化因子梯度、细胞尺度的机械阻碍（当前版本未建模）

### 用 `dissociation_drivers()` 先测一测

```python
from sparta.validate import dissociation_drivers
dd = dissociation_drivers(A, ecm, caf, crosslink, ag_target, source, vessel,
                          cfg_cell=..., cfg_mab=...)
print(dd["dissociation_potential"])   # 特异因素 / 共享因素
```

**`dissociation_potential > 1` 才有望在全片尺度观察到解离。**

在当前默认参数下的合成切片上，这个值只有 **0.04**——ECM 解释了 `B_mAb` 变异的
96%，抗原只占 4%。这意味着：

> 如果真实的皮肤肿瘤切片也给出远小于 1 的解离潜力比，那么"两道屏障可解离"
> 这个中心主张在本框架下**不会成立**，而且这与数据质量无关，是模型结构决定的。

### 建议的应对（按顺序）

1. **建库阶段就对每张切片跑 `dissociation_drivers()`**，把结果记进数据台账。
   这只需几秒钟，却能在第二个月就告诉你这条路走不走得通——而不是第八个月。
2. 若普遍远小于 1，先检查 `configs/default.yaml` 里 ECM 在两个算子中的权重
   （`b_cell.b_ecm` 与 `b_mab.lam`）是否有物理依据。这两个参数目前属于
   `[无监督标定]`，其相对大小直接决定了共享通道有多强。
   **注意：不能为了让解离出现而调它们——那是信息泄漏的变体。**
   只能依据文献（ECM 对细胞迁移 vs 对大分子扩散的相对阻碍程度）来定。
3. 若确实解离不了，**把叙事调整为**："两道屏障在皮肤肿瘤中高度耦合，
   且耦合的来源是 ECM 而非几何深度"——这同样是一个有价值、可发表的结论，
   而且有本框架的分解分析作为量化支撑。**但必须在写作之前就知道。**

课题的双向可发表性正体现在这里：解离成立是主线发现，不成立是澄清性结论。
真正致命的只有一种情况——在第八个月才发现这件事。

---

## 接下来做什么（按顺序，每步都有对应脚本）

> **关键改动**：把"这个课题能不能成"的风险检查从第八个月提前到**第三周**。
> 拿到 3–5 张切片就跑 `run_07_screen.py`，不要等全库建完。

### 第 0 步｜环境（1–2 天）

```bash
conda env create -f environment.yml && conda activate sparta && pip install -e .
python scripts/setup_check.py          # 逐项检查依赖，缺什么告诉你装什么
python scripts/run_00_demo.py --fast   # 合成数据跑通全流程
python tests/test_barrier.py && python tests/test_counterfactual.py && python tests/test_loaders.py
```

验收：14 个测试全绿 + demo 出图。**做不到就别往下走。**

### 第 1 步｜先拿 3–5 张，不求全（1 周）

不要一上来就下 20 个数据集。目标是尽快拿到能跑通的少量切片，去做第 3 步的决策。

```bash
# 下载后先扫一遍，只看不加载，很快
python scripts/run_00b_ingest.py --scan data/raw

# 逐个摄入并登记进数据台账
python scripts/run_00b_ingest.py --input data/raw/XXX --slide MEL01 \
    --cancer melanoma --platform legacy_st --source "GEO GSExxxxx" \
    --treatment naive --site primary
```

`--scan` 会告诉你每份数据是什么格式、三要素齐不齐。缺坐标(C2)的直接剔除。
摄入结果记进 `data/ledger.csv`——这是 P1 阶段最重要的产出。

⚠ 摄入后**一定**看一眼 `n_spots` 和 `n_genes`。如果 spot 数大于基因数，
多半是方向判反了，加 `--transpose` 重跑。方向反了下游每一步都能跑通，只是全错。

### 第 2 步｜跑到建图（2–3 天）

```bash
python scripts/run_batch.py --through 03
```

批量跑 QC → 打分 → 建图。单张失败不中断，最后给你一张 `batch_report.csv`。

两个必看：`run_02 --plot` 出的签名空间图（目视确认恶性在瘤巢、内皮呈条索）；
`run_03` 报告的连通分量数（碎成多块说明 `radius_um` 太小）。

### 第 3 步｜⭐ 决策点：解离潜力筛查（半天）

```bash
python scripts/run_07_screen.py --synthetic          # 先看输出长什么样
python scripts/run_07_screen.py --slides MEL01 MEL02 SCC01
```

输出三选一的判定：

| 判定 | 含义 | 怎么办 |
|---|---|---|
| **GO**（潜力比 > 1） | 抗原/交联的贡献超过共享 ECM | 按原计划建全库，主线叙事成立 |
| **CAUTION**（0.3–1） | 解离信号存在但不占优 | 优先纳入抗原异质性大的切片；加大 S3 尺寸扫描的权重 |
| **PIVOT**（< 0.3） | 共享 ECM 主导 | **改叙事**：从"两道屏障可解离"改为"两道屏障高度耦合，耦合源是 ECM 而非几何深度" |

PIVOT 不是失败。它同样是可发表的澄清性结论，而且你已经有全部支撑证据。
**唯一致命的情况是第八个月才发现。**

### 第 4 步｜按判定分叉建全库（3–4 周）

- **GO / CAUTION** → 补到 15–20 张，`run_batch.py --through 05`
- **PIVOT** → 仍然建库（结论需要样本量），但把 S3 尺寸扫描与驱动分解提为主图

### 第 5 步｜反事实与验证（2–3 周）

```bash
python scripts/run_batch.py --through 05                        # S1/S2/S3
python scripts/run_06_validate.py --slides ... --dim consistency
python scripts/run_06_validate.py --slides MEL01 --dim morphology --he path/to/he.png
python scripts/run_06_validate.py --slides MEL01 --dim sensitivity
```

维度③（临床队列）需要按队列格式单独写，`sparta/validate.py` 里的
`project_signature_to_bulk` 与 `incremental_value` 是现成的积木。

### 时间账

| 步骤 | 时间 | 卡住的话 |
|---|---|---|
| 0 环境 | 1–2 天 | `setup_check.py` 会告诉你缺什么 |
| 1 少量数据 | 1 周 | `--scan` 先筛，缺坐标的别浪费时间 |
| 2 建图 | 2–3 天 | 看 `batch_report.csv` 定位 |
| **3 决策点** | **半天** | **这一步决定后面几个月怎么走** |
| 4 全库 | 3–4 周 | 数据获取是瓶颈，不是算力 |
| 5 验证 | 2–3 周 | 反事实实验成本极低，优先做 |

到第 3 步大约第 3 周。**在那之前不要写论文的任何一段。**

---

## 项目结构

```
sparta/
├── configs/default.yaml       所有参数（分物理锚定 / 无监督标定 / 敏感性三类）
├── data/{raw,interim,external}
├── sparta/
│   ├── io_.py                 数据契约与路径管理
│   ├── signatures.py          M2 签名打分 + 基因集清单
│   ├── graph.py               M3 空间图 + 源汇定义
│   ├── barrier.py             M4 三分量算子  ★ 核心
│   ├── counterfactual.py      M5 S1 / S2 / S3
│   ├── validate.py            M6 四维验证
│   ├── viz.py                 出图（封锁线叠加图）
│   ├── loaders.py             原始数据格式适配（Visium/第一代ST/通用矩阵）
│   └── synthetic.py           合成图生成器（测试 + 教学）
├── scripts/
│   ├── setup_check.py         环境自检
│   ├── run_00_demo.py         合成数据端到端演示
│   ├── run_00b_ingest.py      原始数据摄入 + 数据台账
│   ├── run_01 … run_06        主流程
│   ├── run_07_screen.py       ⭐ 早期决策点
│   └── run_batch.py           批量运行
├── tests/  (test_barrier / test_counterfactual / test_loaders)
└── results/
```

详细的开发规范见 [`CLAUDE.md`](CLAUDE.md)。

---

## 引用与免责

本仓库不包含任何真实数据，也不产生任何研究结论。所有合成数据在代码中
明确标记为 `synthetic` / `placeholder`。基因集标注了 `TODO_REF` 的条目，
正式撰写论文时必须补上文献出处。
