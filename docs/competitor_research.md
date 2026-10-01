# 竞品调研纪要（2026-08-28）

> 目的：为 R6（与现有工具的关系）与审稿答辩准备"谁跟我们最像、差在哪"的弹药。
> 结论先行：**没有直接竞品**。最接近的两篇（TCN、Spatial-RGCN）都把 min-cut 用作
> 聚类损失，不是屏障度量；"把迁移屏障与传质屏障写成图上可解算子"这一步没人做过。
> 调研方式：两名并行检索代理 + 人工核实（克隆刊陷阱已排查，见文末）。

## 1. 空间域分割/聚类（R6 已对比的家族，审稿人最可能点名的方法）

| 方法 | 出处 | 与 SPARTA 的差异 |
|---|---|---|
| BANKSY | Nature Genetics 2024, Lee et al. | 输出域标签，不输出屏障强度（正文已对比） |
| GraphST | BIB 2024（自监督图对比学习） | 嵌入+聚类，无传输/割语义 |
| STAGATE | 图注意力自编码器 | 同上；多个 benchmark 前列，审稿人可能点名 |
| SpaGCN / BayesSpace | 图卷积 / MCMC 贝叶斯 | 同上 |
| stKeep | Nature Communications 2024 | 异构图解析 TME，目标是细胞态/基因模块/CCC |

**要点**：这一家族输出的是"域标签或嵌入"，SPARTA 输出的是"物理可解释的屏障量"。
R6 现有的"域分割给候选集、最小割做选择"叙事与该家族定位完全兼容，不用改。

## 2. 图论流/割/传输 × 空间转录组（**直接竞品核查**）

| 方法 | 出处 | 用法 | 与我们的差异 |
|---|---|---|---|
| TCN | Nature Methods 2023, Doiron et al. | min-cut **损失函数**学细胞邻域软分配 | 不是源-汇屏障，不量化迁移阻力 |
| Spatial-RGCN | ACM 2025 | min-cut 目标做空间域识别 | 同上 |
| SpaceFlow | Nature Communications 2022 | 深度图网络 + 伪时空图（名字带 Flow，实无流模型） | 正文已对比 |
| interFLOW | npj Syst Biol Appl 2024 | PPI 网络上 max-flow 找信号汇聚因子 | 非空间转录组、非屏障 |
| DGM | 2025, 单细胞图聚类 + mincut | 无空间坐标 | 背景引用即可 |

**结论**：min-cut 在 ST 里的既有用法全部是聚类/分割目标函数。把割集解释为
**免疫细胞迁移屏障**、再配一个**屏蔽泊松传质场**做双算子对照——没有先例。
新颖性声明站得住。R6 建议补一句与 TCN/Spatial-RGCN 的显式区分（见文末待办）。

## 3. 空间感知细胞通讯 / 免疫排斥生态位

- LIANA+（Nature Methods 2024）、COMMOT（最优传输建 L-R 信号流）、CellPhoneDB、NicheNet：输出通讯强度/路径，不计算"穿过屏障的通过性"。
- 生物学侧（支持前提、不构成方法竞争）：Adv Sci 2025 肝癌 CAF 双屏障生态位（RCTD+免疫荧光）、Comm Biol 2025 TGF-β 成纤维生态型与 ICB 抵抗、黑色素瘤 CAF-M1（PMC 2025）。
- 这些论文证明"CAF/基质屏障"是热点叙事，SPARTA 提供的是它们缺的**量化层**——Introduction/ Discussion 可放心引用来支撑动机。

## 4. 药物渗透 / 大分子传输建模（B_mAb 的邻居学科）

- Theranostics 2024：CRD 对流-反应-扩散 PDE 模拟放射性药物在前列腺癌 ST 的分布——连续场 PDE，非图上可解算子。
- npj Comput Biol 2024：BioFVM/PhysiCell 模拟 GBM 免疫治疗 TME——连续场多尺度模拟。
- ADC 的 PBPK/Krogh cylinder 模型（PMC 2020）、AAPS J 2025 大分子靶结合耗尽——药代动力学传统框架。

**定位话术**：SPARTA 与连续场药代模型**互补而非竞争**——它把传质问题搬到 ST 离散图上，
变成免标定、确定性、逐切片可算的量。Methods/ Discussion 里这段区分值得写清楚。

## 5. 对投稿的启示

1. 新颖性表述模板（供 Abstract/Discussion 用）：
   > We propose two interpretable, non-learning spatial barrier metrics on ST
   > graphs: a source–sink minimum-cut barrier for T-cell migration and a
   > screened-Poisson diffusion–absorption barrier for large-molecule antibody
   > transport — deterministic, physically grounded quantities of access,
   > not domain labels or communication scores.
2. 审稿人可能要求补的对比：STAGATE/GraphST（域分割 SOTA）——应对：R6 的框架
   已把该家族统一为"候选集供给者"，且 n_domains 扫描证明结论不依赖域数；
   可在 rebuttal 指出任何域分割方法都可作为候选集来源，不影响"最小割做选择"。
3. **R6 小补丁（待办）**：加 1–2 句与 TCN（min-cut 作为聚类损失）的显式区分，
   引 Doiron et al. 2023；这是审稿人搜 "min-cut spatial transcriptomics"
   最先撞上的论文。

## 附：调研中排掉的陷阱

- sciexpor.com 与 macwillpublishers.com 上的 "Computational Biology and
  Chemistry"（声称免费到 2028 / $3050 APC）**是克隆刊站点**。正版是 Elsevier
  的混合刊（ISSN 1476-9271，WoS 收录，订阅路线可免费发表）。投稿一律走
  sciencedirect.com 官方入口。

## 关键链接

- BANKSY: https://www.nature.com/articles/s41588-024-01664-3
- TCN: https://www.nature.com/articles/s41592-023-02124-2
- Spatial-RGCN: https://dl.acm.org/doi/10.1145/3820664.3820672
- SpaceFlow: https://pmc.ncbi.nlm.nih.gov/articles/PMC9283532/
- interFLOW: https://www.nature.com/articles/s41540-024-00391-z
- stKeep: https://pmc.ncbi.nlm.nih.gov/articles/PMC11176411/
- LIANA+ 综述: https://pmc.ncbi.nlm.nih.gov/articles/PMC11392821/
- COMMOT: https://pmc.ncbi.nlm.nih.gov/articles/PMC10541142/
- 肝癌 CAF 双屏障（Adv Sci 2025）: https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.202514661
- 放射性药物 CRD 模型: https://www.thno.org/v14p7122.htm
- GBM BioFVM: https://www.nature.com/articles/s41540-024-00419-4
- CBC 正版信息（Elsevier, hybrid）: https://wos-journal.info/journalid/14813
