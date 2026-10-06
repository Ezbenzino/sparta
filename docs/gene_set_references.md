# 基因集出处（Methods 用）

> 2026-08-27 补齐。对应 `sparta/signatures.py` 里的每一个基因集。
> **两级标注，不要混用**：
> `[ref]` = 有明确文献出处，可直接进引用列表；
> `[canon]` = 领域内通用谱系标记，没有唯一原始出处，只能写成
> "canonical lineage markers"。**给一个没核实过的引用，比不给更糟。**

| 基因集 | 级别 | 出处 |
|---|---|---|
| `ECM_core` | `[ref]` | Naba A, et al. The matrisome: in silico definition and in vivo characterization by proteomics of normal and tumor extracellular matrices. *Mol Cell Proteomics* 2012;11:M111.014647.<br>Shao X, et al. MatrisomeDB: the ECM-protein knowledge database. *Nucleic Acids Res* 2020;48:D1136–D1144. |
| `ECM_crosslink` | `[ref]` | Levental KR, et al. Matrix crosslinking forces tumor progression by enhancing integrin signaling. *Cell* 2009;139:891–906.<br>Barker HE, Cox TR, Erler JT. The rationale for targeting the LOX family in cancer. *Nat Rev Cancer* 2012;12:540–552.<br>Cox TR, et al. LOX-mediated collagen crosslinking is responsible for fibrosis-enhanced metastasis. *Cancer Res* 2013;73:1721–1732. |
| `myCAF` / `iCAF` | `[ref]` | Öhlund D, et al. Distinct populations of inflammatory fibroblasts and myofibroblasts in pancreatic cancer. *J Exp Med* 2017;214:579–596.<br>Elyada E, et al. Cross-species single-cell analysis of pancreatic ductal adenocarcinoma reveals antigen-presenting cancer-associated fibroblasts. *Cancer Discov* 2019;9:1102–1123. |
| `Malignant`（melanoma） | `[ref]` | Tirosh I, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science* 2016;352:189–196. |
| `Malignant`（cSCC）、`TSK` | `[ref]` | Ji AL, et al. Multimodal analysis of composition and spatial architecture in human squamous cell carcinoma. *Cell* 2020;182:497–514. |
| `Hypoxia` | `[ref]` | MSigDB HALLMARK_HYPOXIA（v7.1，200 基因）。Liberzon A, et al. The Molecular Signatures Database hallmark gene set collection. *Cell Syst* 2015;1:417–425. |
| `CAF`（泛） | `[canon]` | 通用成纤维/基质标记 |
| `Endothelial` / `T_NK` / `CD8T` / `Myeloid` / `B_Plasma` | `[canon]` | 通用谱系标记。同样的谱系在本项目两个瘤种对应的单细胞研究（Tirosh 2016、Ji 2020）里也用高度重叠的标记定义，Methods 可以这样说明，但**不能写成"取自"那两篇** |
| `Proliferation` | `[canon]` | 与 HALLMARK_G2M_CHECKPOINT / E2F_TARGETS 高度重叠。**建议定稿前直接换成这两个 Hallmark 集合**，这样就升级为 `[ref]`（同 Liberzon 2015） |
| `Efflux` | `[ref]` | ABCB1 = P-gp/MDR1、ABCG2 = BCRP、ABCC1 = MRP1、ABCC2 = MRP2。<br>Gottesman MM, Fojo T, Bates SE. Multidrug resistance in cancer: role of ATP-dependent transporters. *Nat Rev Cancer* 2002;2:48–58.<br>Robey RW, et al. Revisiting the role of ABC transporters in multidrug-resistant cancer. *Nat Rev Cancer* 2018;18:452–464. |
| `Ag_target` | — | `CD274`（PD-L1）/ `PDCD1LG2`（PD-L2），即抗 PD-1/PD-L1 的靶点本身，不需要基因集引用 |

---

## 三条必须写进 Methods 的诚实声明

1. **`ECM_core` 是人工挑出的精简子集，不是完整 core matrisome。**
   正文写 "a curated subset of the core matrisome"，不要写 "the core matrisome"。

2. **myCAF / iCAF 的二分是在胰腺癌里建立的。** 用到皮肤肿瘤上属于跨瘤种借用，
   要写明这一点，并说明本项目只用它做描述性分层，屏障算子并不依赖这个二分。

3. **TSK 与我们的 cSCC 数据同源。** TSK 由 Ji et al. 在**本项目所用的同一批数据**里
   定义。拿它验证封锁带位置既是优势（同源、无需配准）也是局限（不是独立验证），
   两面都要写。

---

## 两个可以顺手升级的地方

- `Proliferation` 换成 HALLMARK_G2M_CHECKPOINT ∪ E2F_TARGETS：出处立刻明确，
  且与已经在用的 HALLMARK_HYPOXIA 同源，Methods 里合成一句话。
  代价：基因数从 7 涨到几百，`score_genes` 的行为会变，需要重跑 M2 之后的全部步骤。
  **建议在补第一代 ST 切片、反正要重跑的时候一起做。**
- `Efflux` 的引用已补齐，但**引用 Robey 2018 时要注意它的立场**：
  ABC 转运体在体外的作用明确，在临床上的贡献一直有争议，那篇正是重新评估这件事的。
  本项目只把它当作 `B_meta` 的一个状态项，不声称它预测耐药——Discussion 里别写过头。
  另外四个基因的均值本来就没比单个 P-gp 多多少信息，若审稿人质疑可以退到 `ABCB1` 单基因。
