"""
M2 签名打分（替代反卷积）
==========================

上游依赖：M1 质控（{slide_id}.qc.h5ad）
下游产出：{slide_id}.scored.h5ad，各签名分数写入 adata.obs，
          秩标准化后的版本以 "_n" 结尾（如 CAF_n），这是 barrier.py 直接消费的列

为什么用签名打分而不是反卷积
----------------------------
反卷积（RCTD、cell2location）能给出更严格的细胞类型比例，但需要配套的单细胞
参考数据、较长的调参过程，且在第一代 ST 平台的大 spot 上常常表现不佳。
本项目的所有屏障算子只需要"这个 spot 的 CAF / 内皮 / T 细胞 / 恶性成分有多强"
这一**相对量**，签名打分完全够用，实现只需几行，且结果完全透明。
反卷积作为可选升级（见 README 的「可选升级」一节）。

基因集来源标注（2026-08-27 补齐，完整表见 docs/gene_set_references.md）
--------------------------------------------------------------------
每个基因集的出处标在下方定义处。标注分两级，**不要把两级混为一谈**：

  [ref]   有明确的文献出处，可以直接写进 Methods 的引用列表。
  [canon] 领域内通用的谱系标记，没有唯一的"原始出处"。这类只能写成
          "canonical lineage markers"，**不能硬安一篇文献上去**——
          给一个没核实过的引用比不给更糟。

Hypoxia 不在此硬编码，用 MSigDB HALLMARK_HYPOXIA（Liberzon et al.,
Cell Syst 2015;1:417-425），config 的 signatures.hypoxia_gmt 指向已下载的 v7.1 文件。
"""

from __future__ import annotations

import numpy as np

__all__ = ["GENE_SETS", "get_gene_sets", "score_all", "rank_normalize", "SCORE_COLUMNS"]


# --------------------------------------------------------------------------
# 跨瘤种通用的签名
# --------------------------------------------------------------------------
_COMMON = {
    # ---- 细胞类型丰度代理 ----
    # [canon] 泛 CAF 标记
    "CAF": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB", "FAP", "THY1", "POSTN"],
    # [ref] myCAF / iCAF 的二分来自 Ohlund et al., J Exp Med 2017;214:579-596；
    #       跨物种单细胞验证与 apCAF 的补充见 Elyada et al., Cancer Discov 2019;9:1102-1123。
    #       注意：该二分是在胰腺癌里建立的，用到皮肤肿瘤上属于跨瘤种借用，Methods 要写明。
    "myCAF": ["ACTA2", "TAGLN", "MYL9", "TPM2", "POSTN", "CTHRC1"],
    "iCAF": ["CXCL12", "C3", "C7", "PDGFRA", "IL6", "CFD"],
    # [canon] 以下五组是通用谱系标记。同样的谱系在本项目两个瘤种对应的单细胞研究里
    #         也是用高度重叠的标记定义的（黑色素瘤 Tirosh et al., Science 2016;352:189-196；
    #         皮肤鳞癌 Ji et al., Cell 2020;182:497-514），可以在 Methods 里这样说明，
    #         但**不要写成"取自"那两篇**——这里的具体基因不是从它们的补充表拷来的。
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "ENG", "EGFL7"],
    "T_NK": ["CD3D", "CD3E", "CD2", "TRAC", "CD8A", "GZMB", "NKG7", "IL7R"],
    "CD8T": ["CD8A", "CD8B", "GZMK", "GZMB", "PRF1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "AIF1", "ITGAX", "CSF1R", "TYROBP", "FCER1G"],
    "B_Plasma": ["MS4A1", "CD79A", "IGHG1", "MZB1", "JCHAIN"],

    # ---- 细胞外基质：屏障边权的核心输入 ----
    # [ref] core matrisome 精简版（胶原 + 糖蛋白 + 蛋白聚糖）。
    #       matrisome 的定义：Naba et al., Mol Cell Proteomics 2012;11:M111.014647；
    #       数据库：Shao et al., MatrisomeDB, Nucleic Acids Res 2020;48:D1136-D1144。
    #       这里是从 core matrisome 中人工挑出的精简子集，不是全集——Methods 必须写明
    #       "a curated subset of the core matrisome"，不能写成"the core matrisome"。
    "ECM_core": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2",
                 "FN1", "LAMB1", "TNC", "THBS2", "FBN1",
                 "VCAN", "BGN", "LUM", "DCN"],
    # [ref] 交联酶：决定局部网孔尺寸，是 B_mAb 尺寸排阻项的输入。
    #       LOX 介导的胶原交联使基质硬化并促进肿瘤进展：
    #         Levental et al., Cell 2009;139:891-906；
    #         Cox et al., Cancer Res 2013;73:1721-1732。
    #       LOX 家族综述（含 LOXL1-3 的分工）：
    #         Barker, Cox & Erler, Nat Rev Cancer 2012;12:540-552。
    #       PLOD1/2（赖氨酰羟化酶）与 TGM2（转谷氨酰胺酶）是另外两类交联途径，
    #       与 LOX 家族并列纳入。
    "ECM_crosslink": ["LOX", "LOXL1", "LOXL2", "LOXL3", "PLOD1", "PLOD2", "TGM2"],

    # ---- 代谢与药理状态：B_meta 的输入 ----
    # [canon] 增殖：与 MSigDB 的 HALLMARK_G2M_CHECKPOINT / E2F_TARGETS 高度重叠，
    #         也可直接换成那两个 Hallmark 集合以获得明确出处（推荐定稿前这么做）。
    "Proliferation": ["MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1", "BIRC5", "TYMS"],
    # [ref] 外排泵：多药耐药中最常被提及的四个 ABC 转运体
    #       ABCB1 = P-gp / MDR1、ABCG2 = BCRP、ABCC1 = MRP1、ABCC2 = MRP2。
    #       经典综述：Gottesman MM, Fojo T, Bates SE. Multidrug resistance in cancer:
    #         role of ATP-dependent transporters. Nat Rev Cancer 2002;2:48-58.
    #       近期再评估（含"临床相关性有多大"的争论）：
    #         Robey RW, et al. Revisiting the role of ABC transporters in
    #         multidrug-resistant cancer. Nat Rev Cancer 2018;18:452-464.
    #       ⚠ 引用后一篇时要注意它的立场：ABC 转运体在体外的作用明确，
    #       在临床上的贡献一直有争议。本项目只把它当作 B_meta 的一个状态项，
    #       不声称它预测耐药——Discussion 里别写过头。
    "Efflux": ["ABCB1", "ABCG2", "ABCC1", "ABCC2"],

    # ---- 靶抗原：B_mAb 的结合位点消耗项 ----
    # 这里是 PD-L1/PD-L2 配体表达代理，不是 PD-1 受体（PDCD1），也不是
    # 抗 PD-1 药物特异性结合位点。仅用于模型中的分布式吸收项；该分数与
    # kd_eff 均未按蛋白浓度/结合动力学标定，不能解释为药物特异性药代参数。
    # 若将来建模 ADC，需重新定义并验证相应靶标表达代理。
    "Ag_target": ["CD274", "PDCD1LG2"],
}

# --------------------------------------------------------------------------
# 瘤种特异的恶性细胞签名
# --------------------------------------------------------------------------
_TUMOR_SPECIFIC = {
    "melanoma": {
        # [ref] 黑色素细胞谱系标记，与 Tirosh et al., Science 2016;352:189-196
        #       用于判定恶性细胞的标记高度一致。
        "Malignant": ["MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10", "S100B", "TYRP1"],
    },
    "cscc": {
        # [ref] 鳞癌角质细胞标记，见 Ji et al., Cell 2020;182:497-514
        #       （本项目 CSCC 队列 GSE144239 即出自该研究）。
        "Malignant": ["KRT5", "KRT14", "KRT6A", "KRT17", "SFN", "S100A2"],
        # 肿瘤特异角质细胞（TSK）：位于 cSCC 肿瘤前沿，与纤维血管龛共同构成天然边界结构。
        # 检验 TSK 富集区是否与最小割屏障线在空间上重合，是一个来自已发表生物学的
        # 独立验证点，建议在分析中专门做。
        # [ref] TSK（tumour-specific keratinocyte）由 Ji et al., Cell 2020;182:497-514
        #       在本项目所用的同一批 cSCC 数据中定义，位于肿瘤前沿的纤维血管龛。
        #       用它做验证等于用同一篇论文自己的发现来核对我们的封锁带位置——
        #       这既是优势（同源、无需配准）也是局限（不是独立数据），Methods 要说清。
        "TSK": ["MMP10", "PTHLH", "FEZ1", "IL24", "KCNMA1", "INHBA",
                "MAGEA4", "NT5E", "LAMC2"],
    },
    "bcc": {
        "Malignant": ["KRT15", "KRT17", "BNC1", "PTCH1", "GLI1", "EPCAM"],
    },
    # 乳腺癌（外部验证队列，2026-10-03 补）。
    # 用泛上皮角蛋白 + EPCAM + MUC1 作为恶性上皮代理——不区分腔面/HER2/三阴性，
    # 因为外部验证只需要"哪里是肿瘤上皮"来定义 sink，不需要分子分型。
    # 注意：不要把黑色素瘤标记（MLANA/PMEL/TYR...）用在乳腺癌上——那些基因在
    # 乳腺上皮里不表达，会把 sink 选在随机位置（2026-10-03 实测踩过这个坑）。
    "brca": {
        "Malignant": ["EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"],
    },
}

GENE_SETS = {k: {**_COMMON, **v} for k, v in _TUMOR_SPECIFIC.items()}

# barrier.py 需要的秩标准化列名
SCORE_COLUMNS = [
    "Malignant", "CAF", "myCAF", "iCAF", "Endothelial", "T_NK", "CD8T", "Myeloid",
    "ECM_core", "ECM_crosslink", "Hypoxia", "Proliferation", "Efflux", "Ag_target",
]


def get_gene_sets(tumor_type: str, hypoxia_genes: list[str] | None = None) -> dict:
    """取出某瘤种的全部基因集。

    参数
    ----
    tumor_type    : "melanoma" / "cscc" / "bcc"
    hypoxia_genes : 缺氧基因集。建议从 MSigDB 读取 HALLMARK_HYPOXIA 后传入。
                    不传则使用一个精简的占位集合（标记为 placeholder，
                    正式分析前务必替换）。
    """
    if tumor_type not in GENE_SETS:
        raise KeyError(f"未知瘤种 '{tumor_type}'，可选：{list(GENE_SETS)}")
    sets = dict(GENE_SETS[tumor_type])
    if hypoxia_genes:
        sets["Hypoxia"] = list(hypoxia_genes)
    else:
        # placeholder —— 正式分析前请替换为 MSigDB HALLMARK_HYPOXIA 全集
        sets["Hypoxia"] = ["VEGFA", "SLC2A1", "CA9", "LDHA", "PGK1", "ADM",
                           "NDRG1", "P4HA1", "BNIP3", "ANKRD37"]
    return sets


# --------------------------------------------------------------------------
def score_all(adata, tumor_type: str, hypoxia_genes=None, min_genes: int = 3,
              ctrl_size: int = 50, seed: int = 0):
    """为一张切片计算全部签名分数，写入 adata.obs。

    匹配到的基因少于 min_genes 的签名会被跳过并告警 —— 不要静默跳过，
    否则你会在下游拿到一列全 0.5 的填充值却毫无察觉。
    """
    import scanpy as sc

    sets = get_gene_sets(tumor_type, hypoxia_genes)
    skipped = []
    for name, genes in sets.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) < min_genes:
            skipped.append((name, len(present), len(genes)))
            continue
        sc.tl.score_genes(adata, gene_list=present, score_name=name,
                          ctrl_size=ctrl_size, random_state=seed)
    if skipped:
        print("[warn] 以下签名因匹配基因过少被跳过：")
        for name, k, tot in skipped:
            print(f"       {name}: 只匹配到 {k}/{tot} 个基因")
        print("       常见原因：基因名版本不一致（如用了 Ensembl ID 而非 symbol），"
              "或该切片测序深度过浅。")
    return adata


def rank_normalize(adata, keys: list[str] | None = None, suffix: str = "_n"):
    """把 adata.obs 中的签名分数秩标准化到 [0, 1]，写入 <key>_n 列。

    barrier.py 消费的是这些 _n 列。跳过这一步会导致边权公式被量纲隐性加权。
    """
    from .graph import rank_normalize_array

    keys = keys or [k for k in SCORE_COLUMNS if k in adata.obs]
    for k in keys:
        if k not in adata.obs:
            continue
        adata.obs[k + suffix] = rank_normalize_array(
            np.asarray(adata.obs[k].values, dtype=float)
        )
    return adata
