"""
M6 四维验证矩阵
================

上游依赖：barrier.py 的输出、原始 H&E 图像、外部 bulk 队列
下游产出：results/validation/*

四个维度及其回应的质疑
----------------------
  ① 跨切片一致性     -> "是不是单张切片的偶然结构？"
  ② 同片 H&E 形态学  -> "转录组分数有没有独立证据支持？"
  ③ 临床队列与头对头 -> "和已有签名比增量在哪？"
  ④ 生物学一致性     -> "生物学上讲得通吗？"

维度③的核心提醒
----------------
临床验证的关键不是"能预测"，而是"比现有方法多提供了什么"。三件事缺一不可：
多变量校正后仍显著、与已有签名头对头、增量价值检验（似然比）。
本模块的 incremental_value 函数就是为第三件事准备的。

本课题在维度③上有一个天然优势：B_mAb 是已有方法完全不覆盖的维度。
即使 B_cell 与 TIDE 等高度相关（这很可能），只要 B_mAb 能提供独立的预测信息，
增量价值就成立。所以分析时应重点报告 B_mAb 的独立贡献，而非合并分数。
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "cross_slide_consistency", "decoupling_stats", "dissociation_drivers",
    "he_patch_features", "correlate_with_morphology",
    "project_signature_to_bulk", "incremental_value",
    "parameter_sensitivity",
]


# --------------------------------------------------------------------------
# 维度① 跨切片一致性
# --------------------------------------------------------------------------
def cross_slide_consistency(fields: dict[str, np.ndarray], n_bins: int = 30) -> dict:
    """同瘤种不同切片的屏障分数分布一致性。

    做法：把每张切片的屏障场标准化后做直方图，两两计算分布距离
    （这里用 1 - Jensen-Shannon 距离作为相似度，取值 0–1，越大越一致）。

    参数
    ----
    fields : {slide_id: (n,) 屏障值数组}
    """
    from scipy.spatial.distance import jensenshannon

    sids = sorted(fields)
    hists = {}
    for s in sids:
        v = np.asarray(fields[s], float)
        v = v[np.isfinite(v)]
        if len(v) < 10:
            continue
        v = (v - v.mean()) / (v.std() + 1e-12)
        h, _ = np.histogram(v, bins=n_bins, range=(-3, 3), density=True)
        hists[s] = h / (h.sum() + 1e-12)

    keys = sorted(hists)
    m = np.eye(len(keys))
    for i, a in enumerate(keys):
        for j, b in enumerate(keys):
            if j <= i:
                continue
            d = jensenshannon(hists[a], hists[b], base=2)
            m[i, j] = m[j, i] = 1.0 - (0.0 if np.isnan(d) else d)

    off = m[~np.eye(len(keys), dtype=bool)] if len(keys) > 1 else np.array([np.nan])
    return dict(slides=keys, similarity=m,
                mean_similarity=float(np.nanmean(off)),
                min_similarity=float(np.nanmin(off)) if len(keys) > 1 else np.nan)


# --------------------------------------------------------------------------
# 解耦性统计（本课题的中心结论）
# --------------------------------------------------------------------------
def _rankz(v):
    from scipy.stats import rankdata
    r = rankdata(v).astype(float)
    return (r - r.mean()) / (r.std() + 1e-12)


def _residualize(v, control):
    """把 v 中能被 control 解释的部分去掉，返回残差（均在秩空间做，抗离群值）。"""
    import numpy as np
    zv, zc = _rankz(v), _rankz(control)
    # 二次项足以吸收单调的几何趋势，又不至于过拟合掉真实信号
    X = np.column_stack([np.ones_like(zc), zc, zc ** 2])
    beta, *_ = np.linalg.lstsq(X, zv, rcond=None)
    return zv - X @ beta


def decoupling_stats(b_cell_like: np.ndarray, b_mab: np.ndarray,
                     *, q: float = 0.75, control: np.ndarray | None = None) -> dict:
    """量化两个屏障的解耦程度，并定位关键的解离区域。

    关键输出是 frac_discordant —— "T 细胞屏障低但抗体屏障高"的 spot 占比。
    这类区域在现有框架下会被判为"热而无效"，本框架给出了可测量的机制解释。
    它是本课题最有临床意义的一个数字。

    ⚠ 必须校正的共同混杂：距血管的几何距离
    ----------------------------------------
    B_cell 与 B_mAb 都是"从血管/免疫入口出发的累积代价"，因此**天然共享一个
    很强的几何成分**——离血管越远，两者都越大。在合成数据上，未校正的
    Spearman ρ 可达 0.95，几乎全部来自这个共同成分。

    如果不校正就报告"两个屏障高度相关，所以不解耦"，那是在测深度而不是测机制；
    反过来，如果不校正就报告"存在解离区"，审稿人会立刻问"是不是只是深浅不同"。

    **正确做法**：传入 control（通常是 barrier.compute_b_meta 返回的
    d_vessel_um，或图上到血管的跳数），本函数会在秩空间对 control 做二次回归
    并取残差，额外报告：
      · rho_partial       —— 去掉几何成分后的偏相关，这才是"机制是否不同"的证据
      · frac_discordant_r —— 基于残差定义的解离区占比，这才是可写进论文的数字

    参数
    ----
    b_cell_like : (n,) 逐 spot 的细胞迁移屏障。**用 barrier.compute_b_cell_field()**，
                  不要用局部 ECM/CAF 分数——那几乎是二值的，分位数切分会退化
    b_mab       : (n,) 抗体屏障
    q           : 判定"高/低"的分位阈值（high >= q，low <= 1-q，对称）
    control     : (n,) 需要校正的共同混杂，通常是距血管的图距离

    返回
    ----
    dict。未传 control 时 rho_partial / frac_discordant_r 为 nan，
    并在返回值中把 controlled 标为 False —— 论文里请报告校正后的版本。
    """
    import numpy as np
    from scipy.stats import spearmanr

    x, y = np.asarray(b_cell_like, float), np.asarray(b_mab, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if control is not None:
        c = np.asarray(control, float)
        ok &= np.isfinite(c)
    x, y = x[ok], y[ok]
    if len(x) < 10:
        return dict(rho=np.nan, p=np.nan, n=int(len(x)), controlled=False)

    rho, p = spearmanr(x, y)

    def _quad(u, v):
        """给定两个量，返回 (低u&高v 的占比, 高u&高v, 低u&低v, 是否退化)。"""
        lo_u, hi_u = np.quantile(u, 1 - q), np.quantile(u, q)
        lo_v, hi_v = np.quantile(v, 1 - q), np.quantile(v, q)
        il, ih = u <= lo_u, u >= hi_u
        jl, jh = v <= lo_v, v >= hi_v
        return (float(np.mean(il & jh)), float(np.mean(ih & jh)),
                float(np.mean(il & jl)), bool(lo_u == hi_u or lo_v == hi_v),
                il & jh)

    fd, fbh, fbl, degen, mask_raw = _quad(x, y)
    if degen:
        print("[warn] decoupling_stats：输入近似二值或分位数退化，解离比例不可解释。"
              "请改用连续量（barrier.compute_b_cell_field）。")

    # 解离区占比必须对照随机期望，否则数字没有意义。
    # 两个屏障相互独立时，"低 u 且高 v"的期望占比就是 (1-q)*(1-q)。
    # q=0.75 时是 6.25% —— 实测低于它，说明解离区比随机还少，
    # 也就是两个屏障是正相关的（2026-08-26 审查发现，见 review §3 B1）。
    chance = float((1.0 - q) ** 2)

    out = dict(
        rho=float(rho), p=float(p), n=int(len(x)), degenerate=degen,
        frac_discordant=fd, frac_both_high=fbh, frac_both_low=fbl,
        chance_discordant=chance,
        enrichment_vs_chance=float(fd / chance) if chance > 0 else np.nan,
        controlled=False, rho_partial=np.nan, p_partial=np.nan,
        frac_discordant_r=np.nan, enrichment_vs_chance_r=np.nan,
        idx_discordant=np.where(ok)[0][mask_raw],
    )

    if control is not None:
        c = np.asarray(control, float)[ok]
        rx, ry = _residualize(x, c), _residualize(y, c)
        rp, pp = spearmanr(rx, ry)
        fdr, fbhr, fblr, _, mask_r = _quad(rx, ry)
        out.update(
            controlled=True,
            rho_partial=float(rp), p_partial=float(pp),
            frac_discordant_r=fdr,
            enrichment_vs_chance_r=float(fdr / chance) if chance > 0 else np.nan,
            frac_both_high_r=fbhr, frac_both_low_r=fblr,
            idx_discordant_r=np.where(ok)[0][mask_r],
        )
    return out


# --------------------------------------------------------------------------
# 解离的驱动因素分解（本课题中心风险的诊断工具）
# --------------------------------------------------------------------------
def dissociation_drivers(A, ecm, caf, crosslink, ag_target, source, vessel,
                         *, cfg_cell: dict, cfg_mab: dict, r_nm: float = 5.5,
                         collinearity_thresh: float = 0.80) -> dict:
    """分解 B_mAb 的变异来源，判断这张切片有没有观察到解离的潜力。

    为什么需要这个诊断
    ------------------
    B_cell 与 B_mAb 之间有**两条**耦合通道，不是一条：
      1. 共享几何 —— 都是"从血管出发的累积代价"，离血管越远两者都越大。
         这条可以用 decoupling_stats(control=d_vessel) 校正掉。
      2. **共享 ECM** —— 同一个 core matrisome 分数既进了最小割的边容量
         （细胞过不去），也进了扩散的边电导（大分子过不去）。
         这条**校正不掉**，它是框架的结构性质：致密基质本来就同时挡两者。

    因此一个 ECM 缺口会同时降低两个屏障。要在全片尺度上观察到解离，
    必须依靠**只影响 B_mAb 而不影响 B_cell** 的因素：
      · 靶抗原（结合位点屏障）—— T 细胞不受抗原表达影响
      · 交联度与分子尺寸    —— 决定网孔排阻，对 10 μm 的细胞不是同一个机制

    本函数通过"消融各个因素后 B_mAb 变化多少"来估计三者的相对贡献。

    返回
    ----
    dict:
      frac_shared_ecm       共享 ECM 通道解释的 B_mAb 变异比例
      frac_antigen          抗原（BSB）特异贡献
      frac_crosslink        交联/尺寸排阻特异贡献
      ecm_crosslink_corr      ECM 与交联度的相关性
      crosslink_collinear     交联度是否与 ECM 共线（超过 collinearity_thresh）。
                              为 True 时交联贡献被保守地并入共享桶——分不开
                              的时候就不声称分得开
      dissociation_potential  (特异 / 共享)。**大于 1 才有望在全片尺度
                              观察到解离**；远小于 1 说明这张切片上两个屏障
                              主要由同一个 ECM 驱动，解耦分析会得到阴性结果

    使用建议
    --------
    在建库阶段就对每张切片跑一遍。若全部切片的 dissociation_potential 都远小于 1，
    应当尽早把论文叙事从"两道屏障可解离"调整为"两道屏障在皮肤肿瘤中高度耦合，
    且耦合的来源是 ECM"——后者同样是一个有价值、可发表的结论，
    但需要在写作前就知道，而不是在第八个月才发现。
    """
    import numpy as np
    from .barrier import compute_b_mab

    def _bm(e, x, a):
        return compute_b_mab(A, e, x, a, vessel, r_nm=r_nm, **cfg_mab)["b_mab"]

    full = _bm(ecm, crosslink, ag_target)
    flat = np.full_like(np.asarray(ecm, float), 0.5)

    # 逐一把某个因素拉平，看 B_mAb 变了多少（变得越多说明该因素贡献越大）
    no_ecm = _bm(flat, crosslink, ag_target)
    no_xl = _bm(ecm, np.full_like(flat, float(np.mean(crosslink))), ag_target)
    no_ag = _bm(ecm, crosslink, np.full_like(flat, float(np.mean(ag_target))))

    def _delta(v):
        d = np.asarray(full, float) - np.asarray(v, float)
        d = d[np.isfinite(d)]
        return float(np.var(d)) if len(d) else 0.0

    v_ecm, v_xl, v_ag = _delta(no_ecm), _delta(no_xl), _delta(no_ag)

    # ---- 共线性校正（关键）----
    # 消融法有一个陷阱：如果交联度与 ECM 高度相关（真实组织里它们确实相关——
    # 胶原多的地方 LOX 家族通常也高），那么"拉平交联度"和"拉平 ECM"改变的是
    # 同一件事，会把共享通道的贡献重复计入"交联特异"，虚高解离潜力。
    #
    # 保守处理：两者相关性超过阈值时，**把交联贡献并入共享桶**。
    # 分不开的时候就不声称分得开——这个方向的错误代价小得多。
    e = np.asarray(ecm, float); x = np.asarray(crosslink, float)
    ok = np.isfinite(e) & np.isfinite(x)
    if ok.sum() > 3 and np.std(e[ok]) > 1e-9 and np.std(x[ok]) > 1e-9:
        r_ex = float(abs(np.corrcoef(e[ok], x[ok])[0, 1]))
    else:
        r_ex = 1.0                      # 无变异时按最保守处理
    collinear = r_ex >= collinearity_thresh

    tot = v_ecm + v_xl + v_ag + 1e-12
    if collinear:
        shared = v_ecm + v_xl
        specific = v_ag
    else:
        shared = v_ecm
        specific = v_xl + v_ag

    return dict(
        var_ecm=v_ecm, var_crosslink=v_xl, var_antigen=v_ag,
        frac_shared_ecm=(v_ecm + (v_xl if collinear else 0.0)) / tot,
        frac_crosslink=(0.0 if collinear else v_xl / tot),
        frac_crosslink_raw=v_xl / tot,
        frac_antigen=v_ag / tot,
        ecm_crosslink_corr=r_ex,
        crosslink_collinear=bool(collinear),
        dissociation_potential=float(specific / (shared + 1e-12)),
    )


# --------------------------------------------------------------------------
# 维度② 同片 H&E 形态学
# --------------------------------------------------------------------------
def he_patch_features(image: np.ndarray, coords_px: np.ndarray,
                      patch_px: int = 64) -> dict:
    """按 spot 坐标裁 H&E 图像块，计算传统形态学特征（全程 CPU）。

    每张空间切片自带配对 H&E，所以这一步不需要额外找数据、不需要配准，
    是本项目性价比最高的正交验证。

    计算三类特征：
      nuclei_density : 核密度（Otsu 阈值分割后的连通域计数 / 面积）
      stroma_frac    : 基质面积占比（用苏木素/伊红的颜色分离近似）
      texture_*      : 灰度共生矩阵纹理（对比度、同质性）

    参数
    ----
    image     : (H, W, 3) uint8 的 H&E 图像
    coords_px : (n, 2) 每个 spot 在该图像上的像素坐标 (x, y)
    patch_px  : 图像块边长。Visium 的 55 μm spot 在全片图上通常对应
                数十到上百像素，请按 scalefactors 换算后再定

    返回
    ----
    dict of (n,) 数组；越界的 spot 记为 nan
    """
    from skimage.color import rgb2gray, rgb2hed
    from skimage.feature import graycomatrix, graycoprops
    from skimage.filters import threshold_otsu
    from skimage.measure import label

    img = np.asarray(image)
    H, W = img.shape[:2]
    n = len(coords_px)
    half = patch_px // 2
    out = {k: np.full(n, np.nan) for k in
           ("nuclei_density", "stroma_frac", "texture_contrast", "texture_homogeneity")}

    for i, (x, y) in enumerate(np.asarray(coords_px, int)):
        if not (half <= x < W - half and half <= y < H - half):
            continue
        p = img[y - half:y + half, x - half:x + half]
        if p.size == 0:
            continue
        try:
            hed = rgb2hed(p)
            h_ch = hed[..., 0]                     # 苏木素通道 ~ 细胞核
            e_ch = hed[..., 1]                     # 伊红通道   ~ 胞质/基质
            thr = threshold_otsu(h_ch)
            mask = h_ch > thr
            lab = label(mask)
            out["nuclei_density"][i] = lab.max() / (patch_px ** 2) * 1e4
            out["stroma_frac"][i] = float(np.mean(e_ch > np.quantile(e_ch, 0.5)))

            g = (rgb2gray(p) * 255).astype(np.uint8)
            glcm = graycomatrix(g, distances=[1], angles=[0], levels=256,
                                symmetric=True, normed=True)
            out["texture_contrast"][i] = float(graycoprops(glcm, "contrast")[0, 0])
            out["texture_homogeneity"][i] = float(graycoprops(glcm, "homogeneity")[0, 0])
        except Exception:  # noqa: BLE001 —— 单个 patch 失败不应中断全片
            continue
    return out


def correlate_with_morphology(barrier: np.ndarray, morph: dict) -> dict:
    """检验屏障分数与形态学特征的相关性（Spearman）。

    预期方向：高屏障区应对应更高的基质占比、更低的核密度（致密纤维基质）。
    方向不符时不要硬解释——那可能是坐标未对齐，先检查 patch 是否裁对了位置。
    """
    from scipy.stats import spearmanr

    b = np.asarray(barrier, float)
    res = {}
    for k, v in morph.items():
        ok = np.isfinite(b) & np.isfinite(v)
        if ok.sum() < 20:
            res[k] = dict(rho=np.nan, p=np.nan, n=int(ok.sum()))
            continue
        rho, p = spearmanr(b[ok], v[ok])
        res[k] = dict(rho=float(rho), p=float(p), n=int(ok.sum()))
    return res


# --------------------------------------------------------------------------
# 维度③ 临床队列
# --------------------------------------------------------------------------
def project_signature_to_bulk(bulk_expr, gene_weights: dict, *, z_score: bool = True):
    """把空间屏障对应的基因签名投射到 bulk 队列。

    参数
    ----
    bulk_expr    : DataFrame，行=基因 symbol，列=样本
    gene_weights : {gene: weight}。权重可来自"屏障高分区 vs 低分区"的
                   差异表达统计量（如 log2FC 或 t 值）

    返回
    ----
    (n_sample,) 的签名分数 Series

    局限（必须写进论文）
    --------------------
    bulk 投射不等价于空间计算：它丢失了全部空间信息，只保留了"屏障区域的
    分子特征在整个样本中有多强"。应在有配对数据处报告二者的相关性，
    并明确声明这一局限。
    """
    import pandas as pd

    genes = [g for g in gene_weights if g in bulk_expr.index]
    if len(genes) < 5:
        raise ValueError(f"投射失败：只有 {len(genes)} 个签名基因出现在 bulk 数据中")
    X = bulk_expr.loc[genes]
    if z_score:
        X = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1) + 1e-12, axis=0)
    w = pd.Series({g: gene_weights[g] for g in genes})
    return (X.T * w).T.sum() / (np.abs(w).sum() + 1e-12)


def incremental_value(y, base_covariates, new_score, *, family: str = "binomial") -> dict:
    """增量价值检验：在已有签名/协变量之上，加入新分数是否显著提升模型。

    这是维度③三件事中最关键的一件。做法是嵌套模型的似然比检验：
        基础模型 : y ~ 已有协变量（TMB、PD-L1、CD8、TIDE ...）
        完整模型 : y ~ 已有协变量 + 新的屏障分数
    若似然比检验显著，说明新分数带来了独立信息。

    参数
    ----
    y               : (n,) 响应标签（0/1）或生存时间（family="cox" 时另见说明）
    base_covariates : DataFrame，已有协变量
    new_score       : (n,) 待检验的新分数
    family          : "binomial"（响应）；生存分析请用 lifelines 的
                      CoxPHFitter 自行做嵌套比较，本函数只处理二分类

    返回
    ----
    dict: llr_stat, p_value, auc_base, auc_full, delta_auc, coef, coef_p
    """
    import pandas as pd
    import statsmodels.api as sm
    from sklearn.metrics import roc_auc_score

    if family != "binomial":
        raise NotImplementedError("目前只实现二分类；生存分析请用 lifelines 做嵌套 Cox 比较")

    X0 = sm.add_constant(pd.DataFrame(base_covariates).astype(float), has_constant="add")
    X1 = X0.copy()
    X1["sparta_score"] = np.asarray(new_score, float)

    y = np.asarray(y, float)
    m0 = sm.Logit(y, X0).fit(disp=0)
    m1 = sm.Logit(y, X1).fit(disp=0)

    llr = 2.0 * (m1.llf - m0.llf)
    from scipy.stats import chi2
    p = float(chi2.sf(llr, df=1))

    auc0 = roc_auc_score(y, m0.predict(X0))
    auc1 = roc_auc_score(y, m1.predict(X1))
    return dict(
        llr_stat=float(llr), p_value=p,
        auc_base=float(auc0), auc_full=float(auc1), delta_auc=float(auc1 - auc0),
        coef=float(m1.params["sparta_score"]), coef_p=float(m1.pvalues["sparta_score"]),
    )


# --------------------------------------------------------------------------
# 参数敏感性（补充材料必备）
# --------------------------------------------------------------------------
def parameter_sensitivity(fn, grid: dict, *, base: dict | None = None) -> list[dict]:
    """在参数网格上重跑某个函数，收集结果，用于画敏感性热图。

    参数
    ----
    fn   : 接受 **kwargs、返回标量或 dict 的函数
    grid : {参数名: [取值列表]}。会做全组合（笛卡尔积），注意组合数
    base : 固定不变的其余参数

    用法示例
    --------
        parameter_sensitivity(
            lambda **kw: compute_b_cell(A, ecm, caf, S, T, **kw)["b_cell"],
            grid={"b_ecm": [4, 6, 8, 10], "c_caf": [2, 4, 6]},
            base={"a": 3.0},
        )
    """
    import itertools

    base = base or {}
    keys = sorted(grid)
    rows = []
    for vals in itertools.product(*[grid[k] for k in keys]):
        kw = {**base, **dict(zip(keys, vals))}
        try:
            out = fn(**kw)
        except Exception as e:  # noqa: BLE001
            out = {"error": f"{type(e).__name__}: {e}"}
        row = {k: v for k, v in zip(keys, vals)}
        row["result"] = out
        rows.append(row)
    return rows
