# 审查整改执行手册（2026-08-26）

> 在你自己的机器上按顺序执行。每一步都写了**为什么跑、跑多久、跑完看什么**。
> 全部命令用仓库自带的 venv：`D:\sparta\.venv\Scripts\python.exe`。
> 下文用 `PY` 代指它。

```powershell
cd D:\sparta
```

> **不要用 `$PY` 这类变量简写。** 变量只在当前 PowerShell 会话里存在，
> 换一个窗口就是空的，`& $null scripts\xxx.py` 会报一个和真实原因无关的错误
> 然后什么也不做（2026-08-27 就这么丢过一次）。
> 本文档里所有命令都写成完整路径 `.\.venv\Scripts\python.exe`，
> 复制粘贴即可，不依赖任何前置状态。

---

## 🚀 先跑这一条：快路径（约 15–25 分钟，直接走到决策岔口）

下面第 0–9 步是按逻辑顺序排的，但**执行顺序不该照它来**：
最贵的第 3 步要 3–5 小时，而决定投哪个刊的第 9 步只要 3–5 分钟。
所以把所有便宜的步骤先跑完，再挂长任务。

```powershell
& .\.venv\Scripts\python.exe scripts\run_fastpath.py
```

按顺序做：测试 → 备份 → 准入台账 → 解耦重跑 → 分队列一致性 →
**★ 共享 ECM 通道检查（决策岔口）** → 参数身份核验 → 统计汇总。
跑完会告诉你判定是 TISSUE / MODEL / REVERSED。

中途断了可以接着跑：`& .\.venv\Scripts\python.exe scripts\run_fastpath.py --from 5`（从第 5 步起）。

> **为什么是 .py 不是 .ps1。** Windows PowerShell 5.1 默认按 ANSI/GBK 读 `.ps1`，
> 无 BOM 的 UTF-8 脚本里中文会被拆错字节，进而把引号配对弄断，
> 报一串看不懂的语法错误。换 Python 驱动直接绕开这一整类问题。
> `scripts\fastpath.ps1` 现在只是个纯 ASCII 的一行包装，调的还是这个 .py。

**跑完先把判定结果告诉我**，我按结果调整 Results 的措辞与投稿路线。
下面的分步说明供你想单独重跑某一步时查。

---

## 第 0 步｜确认我的改动没把东西改坏（2 分钟）

```powershell
& .\.venv\Scripts\python.exe tests\test_barrier.py
& .\.venv\Scripts\python.exe tests\test_counterfactual.py
& .\.venv\Scripts\python.exe tests\test_loaders.py
& .\.venv\Scripts\python.exe tests\test_statistics.py      # 新增，4 个统计口径测试
```

**验收**：6/6、4/4、4/4、4/4 全绿（共 18 个）。
我已在隔离环境（numpy 2.4.4 / scipy 1.17.1 / networkx 3.6.1）跑过全绿，
你这边应当一致。任何一个红了先停下来告诉我，不要往下走。

---

## 第 1 步｜先备份将被覆盖的产物（1 分钟，**不要跳过**）

第 3 步会覆盖 `results/counterfactual/*.json`。旧结果是对照基线，必须留着。

```powershell
mkdir results\_archive_2026-08-26 -Force
copy results\counterfactual\*.json results\_archive_2026-08-26\
copy results\validation\*.json      results\_archive_2026-08-26\
```

---

## 第 2 步｜重建准入台账（3–5 分钟）

**为什么**：`run_01_qc.py` 旧版只在准入失败时写记录，而失败分支被 `--force` 跳过，
结果 8 张里只有 MEL01 留下了准入记录。Methods 里"全部通过准入"这句话
在磁盘上没有任何证据。这一步把台账补回来。**它只读 raw.h5ad，不改任何管线产物。**

```powershell
& .\.venv\Scripts\python.exe scripts\run_00c_admission_audit.py
```

**跑完看什么**：终端会打印每张切片的 C1–C7 逐条结论，并写出
`data\admission_audit.csv`。预期结果（按 config 阈值）：

- MEL01 的 C7 不达标（GEO 上治疗状态未知）
- MEL02（840 spot）、CSCC01（744）、CSCC02（696）的 C4 不达标
- MEL04（中位 UMI 1488）、CSCC03（982）、CSCC04（636）的 C5 不达标
- CSCC 四张会显示"已按队列 cscc_gse144239 放宽阈值"（spot≥500、UMI≥500）

把这张表原样做进 **Table 1 的"准入"列**。正文措辞见
`docs/manuscript_results.md` 的 Methods M1。

---

## 第 3 步｜重跑 S1 + S2（最耗时的一步，建议挂着过夜）

**为什么**：两件事一起修。
① S1 的置换次数两队列不一致（MEL 500 次、CSCC 只有 150 次，p 下限不同、不可比）；
② S2 的 k 是绝对值，在 614 节点的割集上只动了 0.5%，效应量被稀释到测不出来。
config 已改为 `s2.k_frac: [0.05, 0.10, 0.20, 0.30]`，`run_05` 会自动使用。

```powershell
& .\.venv\Scripts\python.exe scripts\run_batch_s1s2.py
```

（单张失败不中断，最后给汇总；想先快速试跑用
`& .\.venv\Scripts\python.exe scripts\run_batch_s1s2.py --slides MEL01 --n-perm 200 --n-rand 80`，
但那个结果不能进论文。）

**耗时估计**：每张切片 S1 需要 500×2 次最小割，S2 需要约 1200 次，
按切片大小不同，**单张 20–40 分钟，八张合计 3–5 小时**。建议睡前挂上。

**想先快速看趋势**（约 1/3 时间，结果不能直接进论文）：

```powershell
& .\.venv\Scripts\python.exe scripts\run_05_counterfactual.py --slide MEL01 --n-perm 200 --n-rand 80
```

**跑完看什么**：终端会打印每个 k 占割集的百分比，例如
`k= 92（占割集 30.0%）: 连续缺口 … (1.85x, p=0.010)`。
关键是看 **ratio 是否随移除比例上升而明显增大**：
- 若在 20–30% 移除比例下 ratio 明显 >1（比如 >1.5x）且显著 → S2 救回来了，R4 可以正常写；
- 若仍然只有 1.0–1.2x → 这是一个干净的阴性结果，按 `manuscript_results.md`
  R4 的写法诚实报告，并把"屏障连续性"从论证里拿掉。

**两种结果都可以接受，不要为了好看去调 k_frac。**

---

## 第 4 步｜重跑解耦（口径修正，10–20 分钟）

**为什么**：旧实现用 `0.5*(ECM_core_n + CAF_n)` 当逐点 B_cell 代理、且没传 control，
产出的负相关是伪结果。新实现用 `compute_b_cell_field` + `control=距血管距离`，
并同时报告解离区占比与随机期望 6.25% 的比较。

```powershell
& .\.venv\Scripts\python.exe scripts\run_06_validate.py --slides MEL01 MEL02 MEL03 MEL04 CSCC01 CSCC02 CSCC03 CSCC04 --dim decoupling
```

**跑完看什么**：每张会打印
`ρ=+0.xxx｜偏相关 +0.xxx (p=...)｜解离区 x.x%（随机期望 6.25%，富集 0.xx x）`。
预期与 `screen_decision.json` 一致：**8/8 偏相关为正、富集全部 <1**。
这一步同时修掉了"跑完 CSCC 就冲掉 MEL"的覆盖陷阱（现在是读回再合并）。

顺便把一致性也按队列分别跑一遍（跨队列合并没有意义）：

```powershell
& .\.venv\Scripts\python.exe scripts\run_06_validate.py --slides MEL01 MEL02 MEL03 MEL04 --dim consistency --tag MEL
& .\.venv\Scripts\python.exe scripts\run_06_validate.py --slides CSCC01 CSCC02 CSCC03 CSCC04 --dim consistency --tag CSCC
```

`--tag` 是新加的：不给它两次运行会互相覆盖（旧版的 `consistency.json`
现在只剩 CSCC 四张，就是这么丢的）。

---

## 第 5 步｜论文统计汇总（10 秒）

```powershell
& .\.venv\Scripts\python.exe scripts\run_12_paper_stats.py
```

**为什么**：论文里的 FDR、S2 汇总、跨队列表原本只能由 `scripts\scratch\_*.py`
那批硬编码路径的临时脚本产生，README 的管线复现不出论文数字。现在收编了。
它只汇总不重算，秒级完成，可以随时重跑。

**跑完看什么**：
- S1 表里 `p下限` 一列现在应当全部是 0.0020（第 3 步统一了置换次数）；
- S2 表里 `占割集` 一列应当是 5%/10%/20%/30%，不再是 0.4%–5.1%；
- 解耦表里"来源"一列应当全部是 `run_06(已校正)`，不再有 `未校正!`。

以上三条任何一条不满足，说明前面某一步没跑到位。

---

## 第 6 步｜核验参数身份的两个论断（5–10 分钟）

```powershell
& .\.venv\Scripts\python.exe scripts\run_11_review_diagnostics.py --slides MEL01 MEL02 MEL03 MEL04 CSCC01 CSCC02 CSCC03 CSCC04
```

**跑完看什么**：`排阻边%` 一列。如果 8 张切片的数值高度接近
（预期都在 55–60% 附近），就证实了"排阻比例由秩分布而非生物学决定"，
这是 Methods M2 那段诚实声明的实测依据，直接写进补充材料。

---

## 第 7 步｜与现有空间工具的对比（每张 1–3 分钟）

```powershell
& .\.venv\Scripts\python.exe scripts\run_13_benchmark_tools.py --slides MEL01 MEL02 MEL03 MEL04 CSCC01 CSCC02 CSCC03 CSCC04
```

**为什么**：Bioinformatics / BIB 的方法学审稿人几乎必问
"你的最小割相对 BANKSY 空间域、Squidpy 邻域富集多了什么"。

**可选**：`pip install squidpy` 后会额外记录 `nhood_enrichment` 的 z 值；
不装也能跑（脚本会提示并跳过）。

**跑完看什么**：`富集` 一列（割边落在域边界上的倍数）与 `precision` 一列。
终端最后会直接给出该怎么写进论文的两句话——**两句都要写**，
只写"我们比它强"那半句会被一眼识破。

---

## 第 8 步｜还没做但投稿前必须做的两件事

### 8a｜敏感性网格补跑 ≥3 张切片

现在 `results\validation\bmab_sensitivity.json` 里只有 **MEL01** 一张，
但骨架把它表述成了全局稳健性。脚本已收编为 `scripts\run_17_sensitivity.py`
（走 `io_.Paths`、支持多切片、把"合理邻域"窗口显式化、色标改对）：

```powershell
.\.venv\Scripts\python.exe scripts\run_17_sensitivity.py --slides MEL01 MEL03 CSCC01 CSCC03
```

少于 3 张会告警——正文只能写"以某张为例"，不能称为全局稳健性。

**两组数字都要报。** 脚本同时给"预设窗口内"与"整个网格上"：
只报窗口内的"88–92% 判 GO"而不给全网格范围，等于把结论对 β 的依赖藏起来了。
窗口（lam∈[1,8]、beta∈[1.5,6]、xi0∈[10,40]）现在是命令行参数并会打印出来——
它必须写进 Methods，并说明是**预先定的**，事后按结果划窗口是隐性的挑选。

补第一代 ST 切片之后连新片一起跑，顺带把 Fig 2b 换成多切片曲线版
（`bmab_sensitivity_beta_curves.png`）。

### 8b｜版本控制与代码可获取性（**硬性投稿条件**）

`D:\sparta` 目前完全不在版本控制下。两个期刊都要求代码可公开获取，
Bioinformatics 的 Application Note 还要求可安装、有版本号。

```powershell
git init
git add .
git commit -m "SPARTA v2.0: dual spatial barrier framework"
# 然后在 GitHub 建公开仓库并 push，打 tag v2.0.0，再用 Zenodo 归档拿 DOI
```

`.gitignore` 已经排除了 `data/` 与 `results/`，**不会**把 GEO 原始数据传上去。
请自己再确认一遍 `git status` 里没有大文件。

---

## 第 9 步｜（最值钱的一步）切断共享 ECM 通道

**为什么。** PIVOT 结论"两道屏障高度耦合"会立刻招来一个致命问题：

> "它们正相关，会不会只是因为你把同一个 core matrisome 分数
>  同时放进了最小割的边容量和扩散的边电导？"

这个质疑是**合理的**——按当前公式 ECM 确实是两个算子的共同输入，
所以正相关有一部分是结构上必然的。答不上来，PIVOT 就从
"组织里两道屏障绑在一起"退化成"我们的模型把它们绑在了一起"。

`configs\no_shared_ecm.yaml` 把 `b_cell.b_ecm` 与 `b_mab.lam` 同时置零
（B_cell 只看 CAF，B_mAb 只看交联与抗原），两个算子从此没有任何共同输入。
`run_14` 会在一次运行里同时算两种口径并给出判定，
**它只写一个新文件 `results\validation\shared_ecm_check.json`，
不覆盖任何已有产物**——不需要你手工搬文件。

```powershell
& .\.venv\Scripts\python.exe scripts\run_14_shared_ecm_check.py --slides MEL01 MEL02 MEL03 MEL04 CSCC01 CSCC02 CSCC03 CSCC04
```

想先看输出长什么样（不需要真实数据，10 秒）：

```powershell
& .\.venv\Scripts\python.exe scripts\run_14_shared_ecm_check.py --synthetic
```

**判读**（脚本会直接打印）：

| 判定 | 含义 | 后果 |
|---|---|---|
| **TISSUE** | 切断共同输入后偏相关仍显著为正 | 耦合是组织的性质，是真正的生物学发现 → 值得冲 BIB，本检查作为主图之一 |
| **MODEL** | 切断后偏相关接近 0 | 耦合主要来自共享输入 → 措辞收窄为"在本框架的边权定义下不可分离"，投 Bioinformatics |
| **REVERSED** | 切断后转为负相关 | 共享项掩盖了一个真实的反向关系，比 PIVOT 更有意思但需要重新过一遍论证 → 先找我 |

三种结果都要写进论文。**主动报告远好过被审稿人问出来。**

（合成数据上这个检查是有鉴别力的：为"存在解离"而造的 `SYN_diss` 切片，
切断共享通道后偏相关从 +0.430 掉到 −0.022、解离区从 3.5% 升到 9.1%；
而没有设计解离的 `SYN_ring` 几乎不变。所以它不是一个对什么都给同样答案的检查。）

---

## 跑完之后

把第 3–7、9 步的实际数字填回 `docs\manuscript_results.md` 里所有标 `[PENDING]`
的地方，然后正文就可以定稿了。达标度评估与投稿路线见 `docs\submission_plan.md`。
