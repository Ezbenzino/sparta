# 执行手册：纳入 GSE144239 第一代 ST（12 张）＋ B_meta 距离口径修正

日期 2026-08-27　　适用于 `D:\sparta`　　所有命令在项目根目录下跑

本次改动只有两件事，但第二件会动到已验证的 Visium 数字，所以先说清楚幅度。

---

## 一、今天改了什么

### 1. 加载器：第一代 ST 的选点文件之前**根本没被用上**

`sparta/loaders.py::inspect_raw` 里有一条短路：

```python
if coord_f and not res["has_coords"]:      # ← has_coords 已被 spot 名分支置真
```

第一代 ST 的坐标能直接从 spot 名 `10x20` 解析出来，于是 `has_coords` 提前为真，
`*spot_data-selection-*.tsv.gz` 就永远不会被记进 `files["coords"]`，
`read_table_matrix` 里的选点过滤也就永远不触发。

后果：整块阵列（1933–1934 个点）被当成组织，其中 **39–77% 是背景**。
源汇的分位数阈值、连通性、最小割全部失去意义，**而且一个错都不报**。

同一个函数还有第二个坑：`_IMAGE_PAT` 要求文件名以 `.jpg` 结尾，
而 GEO 上的 H&E 是 `.jpg.gz`，于是准入标准 C3 会对 12 张全部误判为不合格。

两处都已修，并加了 `tests/test_loaders.py::test_legacy_selection_filter` 钉住。

**真实数据验证（CSCC05）**：`1933 -> 666 个 spot`，选点文件里的 666 个点全部能在阵列里对上，一个没漏。

| 切片 | 阵列 | 组织上 | 占比 | | 切片 | 阵列 | 组织上 | 占比 |
|---|---|---|---|---|---|---|---|---|
| CSCC05 | 1933 | 665 | 34% | | CSCC11 | 1933 | 1144 | 59% |
| CSCC06 | 1933 | 648 | 33% | | CSCC12 | 1934 | 1070 | 55% |
| CSCC07 | 1934 | 637 | 32% | | CSCC13 | 1934 | 1181 | 61% |
| CSCC08 | 1934 | 596 | 30% | | CSCC14 | 1934 | 607 | 31% |
| CSCC09 | 1934 | 525 | 27% | | CSCC15 | 1934 | 620 | 32% |
| CSCC10 | 1934 | 525 | 27% | | CSCC16 | 1934 | 461 | 23% |

12 张全部过 `min_spots_legacy_st: 300`。

### 2. `d_vessel_um` 改用加权图距

**起因**：第一代 ST 是**交错阵列**——spot 名里 x+y 恒为偶数，最近邻在对角方向。
`run_01_qc` 把中位最近邻间距归一到 `spacing_um`，所以 1 个阵列索引 = 141.4 μm，
最近邻 200 μm、次近邻 282.8 μm。由此 `radius_um=300` 连的是「最近邻＋次近邻」，均度 7.2–7.6。

`compute_b_meta` 原先用 `d_um = hops × spacing_um`。这个近似只在**所有边等长**时成立，
而它的偏差**随平台阵列几何反号**（实测 `hops×spacing ÷ 直线距离` 的中位数）：

| 平台 | radius_um | 均度 | 比值 | |
|---|---|---|---|---|
| Visium 六方阵列（合成阵列量的几何常数，标注 `synthetic`） | 150 | 5.87 | **1.109** | 高估 11% |
| 第一代 ST 交错阵列（实测 CSCC05/11/16） | 250 | 3.7–3.8 | 1.265 | 高估 27% |
| 第一代 ST 交错阵列（实测 CSCC05/11/16） | 300 | 7.2–7.6 | **0.894** | 低估 11% |

`d0_um = 130` 锚的是**直线**氧扩散极限，于是两个队列的 B_meta 之间凭空多出约 22% 的系统偏移——
而 R3b 比的正是队列差异，这个偏移会被读成生物学差异。

**改法**：`compute_b_meta` 新增 `D=` 参数。传了就用加权图距（`D` 的边权本来就是微米，
`graph.npz` 里一直存着），不传才退回旧的跳数近似（合成图、老单元测试）。
返回值多一个 `d_vessel_mode`（`"weighted"` / `"hops"`），run_04 会打印并写进 `barrier_meta.json`。

**对已验证 Visium 数字的影响**（`B_meta = psi(d) × state`，state 不变，所以变化全在 psi）：

| 切片 | psi 均值（前 → 后） | 逐点 \|Δpsi\| 中位 | 最大 |
|---|---|---|---|
| MEL01 | 0.4217 → 0.4297 | 0.0089 | 0.0385 |
| CSCC01 | 0.4629 → 0.4640 | 0.0009 | 0.0034 |

即 Visium 侧本来就基本无偏（边长几乎都是 100 μm），改动幅度 ≤2%；
差异主要来自跳数是整数、大量并列被打散。第一代 ST 侧则是一次到位，不必先错一遍再改。

**顺带**：`radius_um=300` 的选择因此有了实证依据。之前配置里写的理由是「否则图会碎成很多块」——
这条是错的：250 和 300 都基本连通（250 时 CSCC09/10/13 各裂出 1 个孤点，300 时只剩 CSCC13 裂 1 个，
最大块 ≥99.9%）。真正的理由是上面那个偏差表。已把实测写进 `configs/cscc_legacy_st.yaml`。

---

## 二、跑之前：确认改动没把东西改坏

```
cd D:\sparta
.\.venv\Scripts\python.exe tests\test_loaders.py
.\.venv\Scripts\python.exe tests\test_barrier.py
.\.venv\Scripts\python.exe tests\test_counterfactual.py
.\.venv\Scripts\python.exe tests\test_statistics.py
```

应为 5 / 7 / 4 / 4，共 **20 全绿**（比改动前多 2 条：选点过滤、加权图距）。
不全绿就先停下，别往下跑。

---

## 三、摄入 12 张（分三批，先打通一个患者）

```
.\.venv\Scripts\python.exe scripts\run_ingest_legacy_cscc.py --only P2
```

盯这几行：

1. `[loaders] 选点过滤：1933 -> 666 个 spot` —— 没这行说明过滤没生效
2. `[M1] 坐标换算：中位最近邻间距 1.41 -> 200.0 μm` —— **1.41 是对的**（交错阵列的 √2）；出现 1.00 反而是错的
3. run_03 报的连通分量数应为 1
4. `[M4] B_meta：…（口径 weighted）` —— 出现 `hops` 说明 D 没传进去
5. `n_spots` 落在 596–666，`n_genes` 一万七千多。反过来（一万七千个 spot）是方向判反了，加 `--transpose` 重来

P2 三张都干净之后：

```
.\.venv\Scripts\python.exe scripts\run_ingest_legacy_cscc.py --only P5 P9 P10
```

`--from` 是新加的：某一步单独补跑用，例如只重跑 run_04：

```
.\.venv\Scripts\python.exe scripts\run_ingest_legacy_cscc.py --only P2 P5 P9 P10 --from 04
```

---

## 四、Visium 侧补跑（B_meta 口径变了）

**只有 run_04 及其下游要重跑；run_05（S1 500 次置换，最耗时的一步）不受影响**——
S1/S2/S3 走的是 B_cell 与 B_mAb，跟 B_meta 无关。

先备份将被覆盖的产物：

```
$stamp = Get-Date -Format "yyyyMMdd_HHmm"
Copy-Item data\interim\*.barrier.npz  "data\interim\bak_$stamp\" -Force
Copy-Item results\*.json              "results\bak_$stamp\"      -Force
```

再重跑：

```
foreach ($s in "MEL01","MEL02","MEL03","MEL04","CSCC01","CSCC02","CSCC03","CSCC04") {
  .\.venv\Scripts\python.exe scripts\run_04_barrier.py --slide $s
}
```

然后按顺序补下游（都用到 `d_vessel_um` 或 `b_meta`）：

| 脚本 | 为什么受影响 |
|---|---|
| `run_06_validate.py` | 解耦的 control 直接读 `barrier.npz` 里的 `d_vessel_um` |
| `run_07_screen.py` | 同上，自己调 `compute_b_meta` 取 `d_vessel_um` |
| `run_11_review_diagnostics.py` | 同上 |
| `run_14_shared_ecm_check.py` | **决策节点**，同上 |
| `run_12_paper_stats.py` | 汇总以上 |
| `run_15_figure1.py` / `run_16_figures.py` | 出图 |

12 张第一代 ST 摄入完成后再跑这一轮，避免下游跑两遍。

---

## 五、回滚

三处改动互相独立，可分别回退：

| 改动 | 文件 | 回退办法 |
|---|---|---|
| 选点过滤 | `sparta/loaders.py` | `git checkout` 该文件；或把 `inspect_raw` 里 `if coord_f:` 改回 `if coord_f and not res["has_coords"]:` |
| `.jpg.gz` 识别 | `sparta/loaders.py` | `_IMAGE_PAT` 去掉 `(\.gz)?` |
| 加权图距 | `sparta/barrier.py` ＋ run_04/07/11/14 | 把四处调用的 `D=D,` 删掉即可退回跳数口径，函数本身向后兼容 |

配置只动了注释，没动任何阈值（`yaml.safe_load` 已核对：三个配置的 `b_meta.spacing_um` 与 `graph.radius_um` 数值不变）。

---

## 五之二、P2 跑通后又补的四处（2026-08-27 下午）

摄入 P2 时暴露出来的，都属于「不报错但会静默出错」那一类。

### 1. 缺氧基因集：config 写了却没被读

`run_02_score.py` 原先只认命令行的 `--hypoxia-gmt`，而 `run_batch.py` 与
`run_ingest_legacy_cscc.py` 都不传它。于是 config 里明明写着
`signatures.hypoxia_gmt: data/external/h.all.v7.1.HALLMARK_HYPOXIA.symbols.gmt`
（文件在，200 个基因），实际用的是 `signatures.py` 里那 10 个基因的占位集，
**只在日志里留一行警告**。而且事后查不出来——从前 run_02 不把用了哪个集合写进 `uns`。

两个队列若用了不同的缺氧集，B_meta 会整体偏移，在 R3b 的队列比较里
长得跟生物学差异一模一样。

改法：命令行 > config > 占位集；来源与基因数写进 `uns["sparta_run"]`；
config 指了文件却读不到**直接报错**，绝不静默退回占位集。

追查旧产物用 `scripts/scratch/check_hypoxia_provenance.py`
（拿 `qc.h5ad` 用两个集合各重算一遍，跟 `scored.h5ad` 里存着的 Hypoxia 比 Spearman，
哪个 r≈1.000 当初就是哪个；只读不改，结果落在 `results/validation/hypoxia_provenance.json`）。

**查下来的结论（2026-08-27）**：

| 切片 | r(placeholder) | r(HALLMARK) | 当初用的 |
|---|---|---|---|
| CSCC01–04、MEL01–04（8 张 Visium） | 0.22–0.59 | 1.0000 | HALLMARK |
| CSCC05/06/07 | 1.0000 | 0.37–0.44 | **占位集** |
| CSCC08–16（9 张） | 0.35–0.54 | 1.0000 | HALLMARK |

也就是说 config 里那句「已用于全部 8 张切片」是真的——只有 P2 那三张（在 config 兜底修好之前摄入的）
用了占位集，已用 `--only P2 --from 02` 重刷。CSCC08–16 走的正是修好后的 config 兜底。

顺带一个值得写进 Methods 的数：两个集合给出的 Hypoxia 分数**秩相关只有 0.22–0.59**。
所以这不是无关痛痒的差别——三张切片要是留着占位集，它们的 B_meta 会跟其余 16 张
站在不同基础上，而这三张恰好是 P2 一整个患者。

### 2. 选点文件的反向对账

原先只有「留下的点少于 10 个就报错」这一道闸，拦不住「恰好对上四成」的中间态——
两边坐标约定不一致（1-based/0-based、x/y 互换）时会安静地留下一堆错点。
现在：对不上的点超过两成直接报错，少量则告警并记进 notes。

12 张实测对不上率：CSCC08 最高 7/597 = 1.17%，CSCC09/10 各 5/526，其余 9 张为 0。
没有一张会触发这道闸。

### 3. `--from` 的前置产物检查

从中途某一步开始时，前面几步的产物必须已经在 `data/interim` 下。
不检查的话 run_04 会以一句光秃秃的 `FileNotFoundError` 死掉，12 张一起跑就是
12 条一模一样的报错，看不出是缺文件还是代码有 bug。现在会逐张列出缺哪个文件，
并给出该跑哪条命令。

### 4. Windows 重定向的编码兜底

Python 在 Windows **控制台**走 UTF-16 接口，`⚠ ✓ ✗` 打得出来；
一旦 stdout 被重定向到文件或管道（`*> run.log`、`| tee`、CI 收日志、subprocess 抓输出），
就换成 locale 编码（简中系统 cp936）且 `errors='strict'`，
第一个编不出的字符直接抛 `UnicodeEncodeError` 把脚本打死。
全项目 21 个文件带这类字符。`sparta/__init__.py` 里加了
`stdout/stderr.reconfigure(errors="replace")`，这些字符在 GBK 环境下退化成 `?`，脚本照跑。

收日志时仍建议先设 `$env:PYTHONIOENCODING="utf-8"`，这样日志本身是 UTF-8 而不是退化的 `?`。

---

## 五之三、准入结果与队列账目（2026-08-27 摄入完成后）

### CSCC13 被剔除

准入 C5 未过：中位 UMI **289.5**，队列阈值 300。

12 张第一代 ST 的中位 UMI 排下来是
`289 | 567 602 675 1442 1648 1817 1955 2846 5388 5737 7219`——
289.5 与次低的 567 之间有近 2 倍空档，CSCC13 独自落在空档下面，
不是卡在密集区边缘。阈值是摄入**之前**写进 `admission_overrides` 的（有日期有理由），
所以这条线不是为了卡掉它才画的。

台账已记 `status=rejected`，理由与证据路径写进 notes，
`data/interim/CSCC13.admission.json` 保留作审计凭证，
原台账备份在 `data/ledger.bak_20260827_cscc13.csv`。

代价只是少一张技术重复：P9 还剩 CSCC11/12，**cSCC 患者数不变，仍是 6 位**。

配套修了 `run_batch.py`——它从前**不看台账的 status**，被剔除的切片照样会被批量跑进
下游产物里，而剔除决定只留在 admission.json 与台账里，两边就此对不上且不报错。
现在默认跳过非 `ingested` 的行并打印跳过了谁，要一起跑得显式加 `--include-rejected`。

### 现在的队列账目

| 队列 | 患者 | 切片 |
|---|---|---|
| cSCC | **6 位**（P2 / P4 / P5 / P6 / P9 / P10） | 15 张（4 Visium + 11 第一代 ST） |
| 黑色素瘤 | **1 位**（MEL_PtB） | 4 张 |

### 剩下最大的结构性弱点：黑色素瘤只有 1 位患者

cSCC 从 2 位患者做到了 6 位，黑色素瘤还是 1 位。R3b 比的是队列差异，
而黑色素瘤那一侧是 n=1 位患者的 4 个不同转移灶
（胸骨、盲肠结节、胸壁、肋骨——是四个解剖位点而非技术重复，比技术重复好，但改变不了 n=1）。

查过 GSE250636 的完整样本表：9 个样本、2 位患者。已摄入的 4 张全是患者 B 的颅外转移；
剩下 5 张是软脑膜转移（LMM），其中 GSM7983358（脊髓后部）仍是患者 B，
GSM7983360（右颞）是患者 A。所以这个数据集里确实还有第二位患者，
但他只贡献软脑膜／中枢神经系统组织——血脑屏障、无真皮、基质完全不同。
把它并进来是拿 n=1 换一个「患者与解剖腔室完全共线」的 n=2，对 R3b 不是真改善。

样本页：
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250636>

### 查过但不成立的两个混杂

| 相关 | 全部 19 张 | 第一代 ST 11 张 |
|---|---|---|
| ρ(spot 数, B_cell) | −0.13 | −0.08 |
| ρ(中位 UMI, B_cell) | +0.29 | +0.44 |

图大小基本不影响 B_cell（原以为最小割会随图变大而系统性抬高最大流，实测没有）。
深度与 B_cell 在第一代 ST 里 ρ=+0.44，n=11 时 p≈0.18，不显著——但值得写进 Limitations，
因为 **P9 一个人贡献了 4 张最浅切片里的 3 张**，深度与患者是共变的。

---

## 五之四、切片名单不再硬编码

`run_12_paper_stats.py` 的 `DEFAULT_SLIDES`、`run_16_figures.py` 的 `ORDER`
从前都是手写的切片列表。队列一变就得记得逐个去改，漏改一处的表现是
**图和统计安静地少几张切片**，没有任何报错。CSCC13 被剔除后还多了一面：
剔除只写在台账与 admission.json 里，硬编码名单不会跟着变，
被剔除的切片可能又从某个名单里溜回下游。

现在两处都调 `sparta/io_.py::admitted_slides(P)`——以台账为唯一事实来源，
取 `status == ingested`，cSCC 在前、黑色素瘤在后。两个脚本启动时会打印
「切片 N 张：...」，名单不对一眼就能看见。

`run_batch.py` 也补了同一道闸：默认跳过台账里非 `ingested` 的行，
要一起跑得显式加 `--include-rejected`。

还剩两个**故意**保留的硬编码名单，动之前得先想清楚它们是不是有意的子集：

| 脚本 | 名单 | 说明 |
|---|---|---|
| `run_09_benchmark_real.py` | `SLIDES = MEL01–04` | 与其他空间工具的对比，只在黑色素瘤上做 |
| `run_10_benchmark_ext.py` | 旧的 8 张 | 现在会漏掉 11 张第一代 ST |

---

## 五之五、下游重跑序列

上游（run_00b→05）全部就位后，一趟跑完：

```powershell
cd D:\sparta
$S = "CSCC01","CSCC02","CSCC03","CSCC04","CSCC05","CSCC06","CSCC07","CSCC08",
     "CSCC09","CSCC10","CSCC11","CSCC12","CSCC14","CSCC15","CSCC16",
     "MEL01","MEL02","MEL03","MEL04"
$C = $S | Where-Object { $_ -like "CSCC*" }
$M = $S | Where-Object { $_ -like "MEL*" }

.\.venv\Scripts\python.exe scripts\run_06_validate.py --slides $S --dim decoupling
.\.venv\Scripts\python.exe scripts\run_06_validate.py --slides $C --dim consistency --tag cscc
.\.venv\Scripts\python.exe scripts\run_06_validate.py --slides $M --dim consistency --tag mel
.\.venv\Scripts\python.exe scripts\run_07_screen.py --slides $S
.\.venv\Scripts\python.exe scripts\run_11_review_diagnostics.py --slides $S
.\.venv\Scripts\python.exe scripts\run_12_paper_stats.py
.\.venv\Scripts\python.exe scripts\run_14_shared_ecm_check.py --slides $S
.\.venv\Scripts\python.exe scripts\run_15_figure1.py
.\.venv\Scripts\python.exe scripts\run_16_figures.py
```

片间一致性必须**分队列**跑（`--tag`）：跨队列的一致性没有意义，
而且不分开跑的话第二次会把第一次的结果冲掉。

**run_05 不在这个序列里**——S1/S2/S3 走的是 B_cell 与 B_mAb，跟 B_meta 无关，
所以加权图距那次改动不波及它。19 张的反事实结果已在
`results/counterfactual/`（S1 均 500 次置换，S2 均 `k_mode=frac`）。

---

## 五之六、下游跑完之后做了什么（2026-08-27 下午后半段）

### 1. 缺氧基因集溯源：只有三张切片用错了

`scripts/scratch/check_hypoxia_provenance.py` 拿 `qc.h5ad` 用两个集合各重算一遍，
跟 `scored.h5ad` 里存着的 Hypoxia 比 Spearman。结果：

| 切片 | r(placeholder) | r(HALLMARK) | 当初用的 |
|---|---|---|---|
| CSCC01–04、MEL01–04 | 0.22–0.59 | 1.0000 | HALLMARK |
| CSCC05/06/07 | 1.0000 | 0.37–0.44 | **占位集** |
| CSCC08–16 | 0.35–0.54 | 1.0000 | HALLMARK |

只有 P2 那三张（在 run_02 的 config 兜底修好之前摄入的）用了占位集，已用
`--only P2 --from 02` 重刷。两个集合给出的 Hypoxia 秩相关只有 0.22–0.59，
所以这不是无关痛痒的差别。结果落在 `results/validation/hypoxia_provenance.json`。

### 2. run_14 的单一判定已降级

判定规则要求「队列内每一张切片都显著」，严苛程度随 n 单调上升——8 张时给 SPLIT
（靠「cSCC 4/4 满贯」），19 张时给 MODEL（cSCC 变成 11–12/15），中间没有任何证据反转。
**不是换阈值，是取消「用一个标签概括」这件事**：换阈值就是事后调参，取消才没有旋钮。

`_summarise()` 现在照常算出旧 verdict 但打上 `verdict_status = "deprecated_2026-08-27"`
并附降级理由，输出里新增 `stratified` 结构（按队列 × 按患者的效应量、保留率、
显著性计数）。正文用 `stratified`。

同时给 run_14 加了 `--d-vessel {weighted,hops}` 与 `--tag`，两种距离口径可以分别落盘。
结论：**两种口径都给 MODEL**，所以判定不是那次几何修正造出来的；
cSCC 在两种口径下都是 5/6 位患者，黑色素瘤 1/1 vs 0/1（n=1，且是唯一随口径翻转的分层）。

### 3. 切片名单全部改成台账驱动

`run_12_paper_stats.py` 的 `DEFAULT_SLIDES`、`run_16_figures.py` 的 `ORDER`、
`run_10_benchmark_ext.py` 的 `SLIDES` 都改调 `sparta.io_.admitted_slides(P)`。
`run_09_benchmark_real.py` 保留黑色素瘤子集（那是有意的），但也支持命令行覆盖。

### 4. Table 1 改由脚本生成

`scripts/run_18_table1.py` 从 `ledger.csv` + `graph.npz` + `admission.json` 生成
`docs/table1_sections.md`。手写的 Table 1 在队列扩张后整张都错了且不报错，
凡是「从产物抄进正文」的表都得能一条命令重新生成。

### 5. 连通性：实测比正文原来写的差得多

从 19 份 `graph.npz` 直接算：质控后图有 **1–54 个连通分量**，最大块占 **90.0–100 %**。
最碎的是两张又大又浅的 Visium（CSCC03 50 个分量、CSCC04 54 个）和两张最浅的
第一代 ST（CSCC11 15 个、CSCC12 16 个）。单连通的只有 Visium 1/8、第一代 ST 6/11。

更要紧的是：**8/19 张切片有源或汇落在最大分量之外**（最多 CSCC04 的 23/187 个源、
MEL02 的 17/232 个汇）。最大流因此实际只在主分量（加自足的碎块）上算，
`B_cell = 1/max_flow` 在碎片化切片上被略微抬高。算子本来就统计并报告了这些点
（`barrier_meta.json` 的 `n_unreachable`），但正文原来没写。已写进 R1。

### 6. S1 的混杂查清了（并更正了一个我自己的错误）

先说错误：早先汇总时把 `paper_stats.json` 的 `s1` 列表按切片名建字典，
而每张切片有 `fixed` 和 `follow` 两条，`follow` 把 `fixed` 覆盖掉了。
以下是**正确**的（`fixed` 主模式）：

| | 一度报错的数 | 实际 |
|---|---|---|
| S1 过 FDR | 10/19 | **12/19** |
| vs 切片大小 | +0.584 | **−0.023**（节点数）/ −0.153（原始 spot 数） |
| vs 中位 UMI | +0.316 | **+0.721** |
| vs 最大连通块占比 | — | **+0.460** |
| cSCC 第一代 ST | 3/11 | **7/11** |
| cSCC Visium | 4/4 | **2/4** |

方向整个反了：S1 的混杂**不是切片大小，是测序深度**（并与连通性共线，
深度 ↔ 最大连通块占比 ρ=+0.453）。结论仍然是"S1 不在正文里承重"，但理由要写对。

---

## 五之七、正文重写范围与剩余工作

`docs/manuscript_results.md` 已按 19 张 / 7 位患者重写：
**R1、R2、R3、R3b、R4-S1、R4-S2、R6**，Methods **M1 / M2 / M3**，
新增 Methods **M4**（距血管距离为什么用加权图距而不是跳数），
Table 1 改为脚本生成。三处 `[PENDING]` 全部落实：

| 原 PENDING | 落实结果 |
|---|---|
| 逐张准入结论抄进 Table 1 | `docs/table1_sections.md` 已含准入列与 CSCC13 的剔除理由 |
| 实测排阻边占比 | 60.0–65.2 %（中位 62.7 %），中位网孔 4.36–4.68 nm，19/19 |
| CSCC 只跑了 150 次置换 | 19/19 均为 500 次，直接从 `n_perm` 字段核过 |

### 7. 四张切片根本没有 Ag_target 信号

`screen_decision.json` 里 CSCC10 / CSCC14 / CSCC15 / CSCC16 的 `var_antigen`
**恰好为 0**——`CD274` 与 `PDCD1LG2` 两个基因都没过最小匹配基因数，
`Ag_target` 签名被中性填充。这四张正是全研究检出基因数最少的
（15 383–17 399）。

要紧的是：**那个 0 是缺数据，不是测量值。** 正文原来写"抗原通道占 0.0–3.1%"，
把缺失当成了"测到接近零"。R2 已改成：另外 15 张上是 0.47–3.06%（中位 0.96%），
四张缺失的单独说明。这类"零"在任何汇总里都要先分清是测出来的还是没测到。

---

**剩下三件事：**

1. ~~R5 仍是 8 张的旧数字~~ —— **已重跑并按 19 张重写**（见下面 五之八）。
2. **run_13 的 Squidpy 臂仍只有 8 张 Visium。** 第一次补跑时命令里的省略号被原样
   传了进去（`slides_in_this_run` 记成 `['CSCC05','...','CSCC16']`），要把 11 张显式列全：
   `--slides CSCC05 CSCC06 CSCC07 CSCC08 CSCC09 CSCC10 CSCC11 CSCC12 CSCC14 CSCC15 CSCC16`
   若 squidpy 在第一代 ST 上直接失败，照实写"该臂只在 Visium 上可用"，不必强补。
3. **`n_domains` 敏感性扫描还没做。** R6 里第一代 ST 的富集掉到 1.07×
   （11 张里 4 张低于 1），成因是同样 8 个域切几百节点的图会让边界边占到 42%，
   是组合学而非组织学。正文已主动交代，但补充材料需要一张按 `n_domains` 扫描的图，
   最好把域数按切片大小缩放。

---

## 五之八、R5 重跑后查出的一个统计量退化

`run_10_benchmark_ext.py` 按台账重跑 19 张后，**Ripley's L 的置换 z 在第一代 ST 上
整个失效**：11 张里 7 张给出 10¹⁴–10¹⁵ 量级的 z，另 4 张是 NaN。Visium 侧正常
（7.5–64.5，8/8 张）。

成因是统计量本身退化，不是数据不好：第一代 ST 是规则的 200 μm 交错阵列，
在规则格点上置换位置几乎不改变点过程的二阶结构，零分布的 sd 塌到 1e-14 量级，
`z = (real-mean)/sd` 就爆掉了。而 `zstat` 原来只判 `sd == 0`，拦不住"近零方差"。

组成类指标（CAF 密度、T 浸润比例）有同一个毛病的温和版：它们在置换下**本该完全
不变**，但分位阈值处的并列会让计数抖动一个 spot，零分布只取两个值，z 就恰好是
±1.0。38 个「切片×指标」组合里 21 个是零方差（NaN）、17 个恰好 ±1.0——那不是
效应量，是退化。旧稿写的"z_dens 与 z_tinf 是 NaN 或 −1.0"漏了 +1.0 那一半。

**已在 `zstat` 里加相对判据**：`sd <= 1e-9 × max(|mean|, 1)` 或算出的 `|z| > 1e6`
一律返回 NaN。宁可缺一个数，也不要让一个假的极端值进表——1e15 这种数一旦进了
结果 JSON，看上去还像"极显著"。

R5 已按此重写：Ripley's L 只报 Visium 并说明退化原因；组成类指标的退化写成
"一个动不了的指标不是关于组织的证据，它的 z 也不是效应量"。

**R5 的其余数字（19 张）**：`B_cell` 场与 CAF 的 ρ = 0.313–0.554（中位 0.458），
与 T/NK 的 ρ = −0.275–+0.140（中位 −0.018）；置换下 `B_cell` 动的有 12/19；
免疫排斥代理 z = −2.4–+5.9，**7 张为负、12 张在 ±2 以内，而且同一瘤种的重复切片
之间就变号**——这比 8 张时"两个队列之间不一致"的说法更难反驳。

---

## 五之九、squidpy 臂"不支持第一代 ST"其实是半径传错了（2026-08-28）

补跑 run_13 后第一代 ST 那 11 张的 `squidpy_nhood_z_offdiag_*` 仍然全是 NaN，
而且 **`squidpy_error` 字段是空的**——也就是说 squidpy 没抛异常，是"跑完了但没结果"。

真正的原因在这一行：

```python
sq.gr.spatial_neighbors(adata, coord_type="generic",
                        radius=cfg["graph"]["radius_um"], spatial_key="spatial_um")
```

跑第一代 ST 时若没带 `--config configs/cscc_legacy_st.yaml`，`cfg` 就是 default.yaml，
`radius_um = 150`。而第一代 ST 的点间距是 **200 μm**——150 μm 半径**一个邻居都连不上**。
空图下 `nhood_enrichment` 的置换零分布方差为 0，z 全是 NaN，squidpy 不报错。

这是本项目最典型的一类失败：**结果看上去像"该方法在这个平台上不适用"，
实际是一个配置参数传错了。**如果不去追，它会以"squidpy 臂只在 Visium 上可用"
的形式写进正文，而那是错的。

**改法不是"记得带 --config"，是把这个坑拿掉**：现在直接把 SPARTA 已经建好的图
交给 squidpy（写 `adata.obsp["spatial_connectivities"]` 与 `spatial_distances`，
再补 `adata.uns["spatial_neighbors"]` 的元信息），不再让它按半径重建。
口径上这也更对——域边界与最小割本来就该在**同一个邻接**上比。

同时加了两道留痕：交给 squidpy 的图零边直接抛错；z 全 NaN 也抛错并打印，
不再当成"跑过了但没结果"。

副作用：Visium 那 8 张的 z 也会变（原来是 squidpy 按 radius=150 自己重建的图），
所以要 19 张一起重跑，两个平台同口径。

> **教训写在这里**：凡是"某个外部工具在我们的某个平台/队列上没有输出"的结论，
> 落笔之前必须先确认它是**抛了错**还是**安静地返回了空**。
> 这两者在产物 JSON 里长得一模一样，含义完全相反。

---

## 六、写进 Methods 的话

- 空间图：半径建图，Visium `radius_um=150`（六方阵列，最近邻 100 μm，均度约 5.9）；
  第一代 ST `radius_um=300`（交错阵列，最近邻 200 μm、次近邻 282.8 μm，均度约 7.5）。
- 距最近血管的距离取**沿图的加权最短路**（边权为欧氏微米距），不是跳数 × 名义间距。
  两者在 Visium 上几乎一致（psi 逐点差中位 0.0009–0.0089），但在第一代 ST 的交错阵列上
  跳数口径会低估约 11%，与 Visium 的高估约 11% 反号，直接污染跨队列比较。
- `d0_um = 130 μm` 为文献锚定的氧扩散极限，不参与任何拟合。
