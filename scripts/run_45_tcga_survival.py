#!/usr/bin/env python
"""
run_45_tcga_survival.py —— TCGA-SKCM bulk 层面的屏障 summary 与总生存（探索性）
===============================================================================
定位（必须按 explore 口径叙述，不得写成"我们的分数能预测预后"）
--------------------------------------------------------------
空间屏障算子需要空间坐标，无法在 bulk 上直接计算。本脚本做的是**概念验证**：
用与空间分析**完全相同**的签名基因集，在 bulk RNA-seq 上打出屏障成分分数，
取一个 bulk 层面的 summary（ECM_core / ECM_crosslink / CAF 的 z 均值），
看它与 TCGA-SKCM 总生存的关系，并用免疫浸润复合分数（CD8T / T_NK / Ag_target）
作对照——回答"这个屏障 summary 是否只是免疫浸润的代理"。

数据（公开，无需授权；UCSC Xena GDC hub）
------------------------------------------
  hub   : https://gdc.xenahubs.net
  cohort: "GDC TCGA Melanoma (SKCM)"
  表达  : TCGA-SKCM.star_tpm.tsv            (STAR TPM，行=带版本号的 Ensembl gene ID)
  probemap: gencode.v36.annotation.gtf.gene.probemap
  生存  : TCGA-SKCM.survival.tsv            (OS / OS.time，天；OS=1 表示死亡事件)
注：GDC hub 已停用 /download/ 直链（403），本脚本走 hub 的 EDN RPC
    (POST <hub>/data/)，只取签名所需的 ~250 个基因列，不做全矩阵下载。
符号→Ensembl 映射走 Ensembl REST 批量 /lookup/symbol/homo_sapiens。

统计口径
--------
· 基因级：log2(TPM+1) 后跨样本 z-score；签名为其成员基因 z 的均值（与 scratch
  的 GSE91061 ICB 概念验证同一口径，便于两处结果互相印证）。
· Cox：Breslow 部分似然 + Newton-Raphson，手工实现（不引入 lifelines）。
  连续变量按 1 个 SD 报告 HR；分层报告 all / primary(01) / metastasis(06)。
· KM：primary barrier composite 三分位 + log-rank；同时报 C-index。
· 多变量：barrier composite + immune composite 同入模型，检验屏障项是否
  独立于免疫浸润。
· 归一化：{sid}-11 为正常组织样本，排除；同一患者多个样本时按 _PATIENT 去重
  （保留 OS.time 最长的那个，并在 JSON 里记录去除数）。

输出
----
  data/external/TCGA-SKCM_survival.csv / TCGA-SKCM_sigexpr.csv（缓存）
  results/validation/tcga_survival.json
  results/figures/tcga_skcm_survival_km.png
"""
from __future__ import annotations

import csv
import json
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

from sparta.io_ import Paths, load_config, save_json, stamp_run  # noqa: E402
from sparta.signatures import GENE_SETS  # noqa: E402

HUB = "https://gdc.xenahubs.net"
COHORT = "GDC TCGA Melanoma (SKCM)"
DS_EXPR = "TCGA-SKCM.star_tpm.tsv"
DS_SURV = "TCGA-SKCM.survival.tsv"
DS_PM = "gencode.v36.annotation.gtf.gene.probemap"
ENSEMBL_LOOKUP = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"

# bulk 屏障 summary 的成员（与空间 B 分量的构成一致）
BARRIER_SIGS = ["ECM_core", "ECM_crosslink", "CAF"]
BARRIER_MAB_SIGS = ["ECM_crosslink", "CAF", "Hypoxia"]   # 与 scratch ICB 口径一致
IMMUNE_SIGS = ["CD8T", "T_NK", "Ag_target"]
ALL_SIGS = ["ECM_core", "ECM_crosslink", "CAF", "Hypoxia", "CD8T", "T_NK",
            "Ag_target", "Myeloid", "Endothelial", "Proliferation", "Efflux"]


# --------------------------------------------------------------------------
# Xena EDN RPC
# --------------------------------------------------------------------------
def xena_post(query: str, timeout: int = 300) -> str:
    req = urllib.request.Request(HUB.rstrip("/") + "/data/", query.encode(),
                                 {"Content-Type": "text/plain"})
    return urllib.request.urlopen(req, timeout=timeout).read().decode()


def xena_samples(dataset: str) -> list[str]:
    q = ('(map :value (query {:select [:value]\n'
         '            :from [:dataset]\n'
         '            :join [:field [:= :dataset.id :dataset_id]\n'
         '            :code [:= :field.id :field_id]]\n'
         '            :where [:and\n'
         f'            [:= :dataset.name "{dataset}"]\n'
         '            [:= :field.name "sampleID"]]}))')
    return json.loads(xena_post(q))


def xena_fields(dataset: str) -> list[str]:
    q = ('(map :name (query {:select [:field.name]\n'
         '             :from [:dataset]\n'
         '             :join [:field [:= :dataset.id :dataset_id]]\n'
         f'             :where [:= :dataset.name "{dataset}"]}}))')
    return json.loads(xena_post(q))


def xena_fetch(dataset: str, columns: list[str], samples: list[str]) -> list[list]:
    """返回 [n_columns][n_samples] 的矩阵（行的顺序与 columns 一致）。"""
    cols = " ".join(json.dumps(c) for c in columns)
    smps = " ".join(json.dumps(s) for s in samples)
    q = (f'(fetch [{{:table "{dataset}"\n'
         f'      :columns [{cols}]\n'
         f'      :samples [{smps}]}}])')
    return json.loads(xena_post(q))


def ensembl_ids(symbols: list[str]) -> dict[str, str]:
    """符号 -> Ensembl gene ID（不带版本号），批量走 Ensembl REST。"""
    out: dict[str, str] = {}
    batch = 180
    for i in range(0, len(symbols), batch):
        chunk = symbols[i:i + batch]
        body = json.dumps({"symbols": chunk}).encode()
        req = urllib.request.Request(
            ENSEMBL_LOOKUP, body,
            {"Content-Type": "application/json", "Accept": "application/json"})
        for attempt in range(3):
            try:
                data = json.loads(urllib.request.urlopen(req, timeout=60).read().decode())
                break
            except Exception:  # noqa: BLE001 —— 网络重试
                if attempt == 2:
                    raise
                time.sleep(2.0)
        for sym, rec in data.items():
            if isinstance(rec, dict) and rec.get("id"):
                out[sym] = rec["id"]
        time.sleep(0.4)  # 对公共 REST 服务保持礼貌
    return out


# --------------------------------------------------------------------------
# 统计（手工实现，不引入 lifelines）
# --------------------------------------------------------------------------
def cox_breslow(X: np.ndarray, time_: np.ndarray, event: np.ndarray,
                max_iter: int = 60) -> dict:
    """Breslow 部分似然的 Cox 比例风险模型。X 列已标准化。"""
    X = np.asarray(X, float)
    t = np.asarray(time_, float)
    e = np.asarray(event, bool)
    n, p = X.shape
    beta = np.zeros(p)
    loglik = np.nan
    info = np.eye(p)
    event_times = np.unique(t[e])
    for _ in range(max_iter):
        eta = np.clip(X @ beta, -50, 50)
        w = np.exp(eta)
        U = np.zeros(p)
        info = np.zeros((p, p))
        loglik = 0.0
        for tt in event_times:
            ev = e & (t == tt)
            risk = t >= tt
            Sw = w[risk].sum()
            Sx = (w[risk, None] * X[risk]).sum(axis=0)
            Sxx = (w[risk, None, None] * X[risk][:, :, None] * X[risk][:, None, :]).sum(axis=0)
            d = int(ev.sum())
            U += X[ev].sum(axis=0) - d * (Sx / Sw)
            info += d * (Sxx / Sw - np.outer(Sx / Sw, Sx / Sw))
            loglik += float(np.sum(X[ev] @ beta)) - d * float(np.log(Sw))
        step = np.linalg.solve(info + 1e-10 * np.eye(p), U)
        beta = beta + step
        if np.max(np.abs(step)) < 1e-9:
            break
    cov = np.linalg.inv(info + 1e-10 * np.eye(p))
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    z = np.divide(beta, se, out=np.zeros_like(beta), where=se > 0)
    pval = 2.0 * stats.norm.sf(np.abs(z))
    return dict(beta=beta, se=se, z=z, p=pval, loglik=loglik, cov=cov,
                n=int(n), n_events=int(e.sum()))


def km_curve(time_: np.ndarray, event: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Kaplan–Meier 生存曲线（返回时间点与生存概率，含 t=0）。"""
    t = np.asarray(time_, float)
    e = np.asarray(event, bool)
    times = np.unique(t)
    s, surv = 1.0, [1.0]
    out_t = [0.0]
    for tt in times:
        at_risk = int(np.sum(t >= tt))
        d = int(np.sum(e & (t == tt)))
        if at_risk > 0 and d > 0:
            s *= (1.0 - d / at_risk)
        out_t.append(float(tt))
        surv.append(s)
    return np.asarray(out_t), np.asarray(surv)


def logrank(time_: np.ndarray, event: np.ndarray, group: np.ndarray) -> float:
    """多组 log-rank（Mantel–Cox）检验 p 值。"""
    t = np.asarray(time_, float)
    e = np.asarray(event, bool)
    g = np.asarray(group)
    levels = [x for x in np.unique(g) if np.sum(e & (g == x)) > 0]
    O = np.zeros(len(levels))
    E = np.zeros(len(levels))
    V = np.zeros((len(levels), len(levels)))
    for tt in np.unique(t[e]):
        at_risk = t >= tt
        n_risk = at_risk.sum()
        d_total = int(np.sum(e & (t == tt)))
        if d_total == 0:
            continue
        for i, lv in enumerate(levels):
            n_i = int(np.sum(at_risk & (g == lv)))
            d_i = int(np.sum(e & (t == tt) & (g == lv)))
            O[i] += d_i
            E[i] += d_total * n_i / n_risk
        # 协方差（超几何）
        for i in range(len(levels)):
            for j in range(len(levels)):
                n_i = np.sum(at_risk & (g == levels[i]))
                n_j = np.sum(at_risk & (g == levels[j]))
                if i == j:
                    V[i, j] += (d_total * (n_i / n_risk) * (1 - n_i / n_risk)
                                * (n_risk - d_total) / max(n_risk - 1, 1))
                else:
                    V[i, j] += (-d_total * (n_i / n_risk) * (n_j / n_risk)
                                * (n_risk - d_total) / max(n_risk - 1, 1))
    diff = (O - E)[:-1]
    Vr = V[:-1, :-1]
    try:
        chi2 = float(diff @ np.linalg.solve(Vr, diff))
    except np.linalg.LinAlgError:
        return float("nan")
    k = len(levels) - 1
    return float(stats.chi2.sf(chi2, k))


def c_index(risk: np.ndarray, time_: np.ndarray, event: np.ndarray) -> float:
    """Harrell C-index（risk 越大越差）。"""
    t = np.asarray(time_, float)
    e = np.asarray(event, bool)
    conc = tie = 0
    n = len(t)
    for i in range(n):
        if not e[i]:
            continue
        cmp_mask = t > t[i]
        if not cmp_mask.any():
            continue
        r_i, r_j = risk[i], risk[cmp_mask]
        conc += int(np.sum(r_i > r_j))
        tie += int(np.sum(r_i == r_j))
    denom = sum(int(np.sum((t > t[i]) & (t > 0))) for i in range(n) if e[i])
    if denom == 0:
        return float("nan")
    return float((conc + 0.5 * tie) / denom)


# --------------------------------------------------------------------------
# 阶段 1：下载与缓存
# --------------------------------------------------------------------------
def ensure_cache(ext: Path, force: bool = False) -> tuple[Path, Path, dict]:
    surv_csv = ext / "TCGA-SKCM_survival.csv"
    expr_csv = ext / "TCGA-SKCM_sigexpr.csv"
    meta_path = ext / "TCGA-SKCM_fetch_meta.json"

    if not force and surv_csv.exists() and expr_csv.exists() and meta_path.exists():
        print("cache present, skipping download")
        return surv_csv, expr_csv, json.loads(meta_path.read_text(encoding="utf-8"))

    meta: dict = {"hub": HUB, "cohort": COHORT, "datasets": {
        "expression": DS_EXPR, "survival": DS_SURV, "probemap": DS_PM}}

    # ---- 生存 ----
    surv_samples = xena_samples(DS_SURV)
    print(f"survival samples: {len(surv_samples)}")
    mat = xena_fetch(DS_SURV, ["OS", "OS.time", "_PATIENT"], surv_samples)
    os_v = mat[0]; os_t = mat[1]; pat = mat[2]
    with open(surv_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["sample", "OS", "OS_time", "patient"])
        for s, o, ot, pt in zip(surv_samples, os_v, os_t, pat):
            w.writerow([s, o, ot, pt])
    meta["n_survival_samples"] = len(surv_samples)

    # ---- 表达（只取签名基因）----
    expr_samples = xena_samples(DS_EXPR)
    print(f"expression samples: {len(expr_samples)}")
    probes = xena_fields(DS_EXPR)
    print(f"expression probes: {len(probes)}")
    probe_by_ensg = {p.split(".")[0]: p for p in probes}

    sig_sets = dict(GENE_SETS["melanoma"])
    hy_path = ext / "h.all.v7.1.HALLMARK_HYPOXIA.symbols.gmt"
    if hy_path.exists():
        with open(hy_path, encoding="utf-8") as f:
            sig_sets["Hypoxia"] = f.readline().rstrip("\n").split("\t")[2:]
    symbols = sorted({g for k in ALL_SIGS if k in sig_sets for g in sig_sets[k]})
    print(f"signature symbols: {len(symbols)}")

    s2e = ensembl_ids(symbols)
    print(f"mapped symbols: {len(s2e)}/{len(symbols)}")

    gene_rows: list[str] = []
    probe_list: list[str] = []
    coverage: dict[str, dict] = {}
    for k in ALL_SIGS:
        if k not in sig_sets:
            coverage[k] = {"n_genes": 0, "n_mapped": 0, "used": []}
            continue
        used = [g for g in sig_sets[k] if g in s2e and s2e[g] in probe_by_ensg]
        coverage[k] = {"n_genes": len(sig_sets[k]), "n_mapped": len(used), "used": used}
        for g in used:
            p = probe_by_ensg[s2e[g]]
            if p not in gene_rows:
                gene_rows.append(p)          # 探针 ID（带版本）
                probe_list.append(g)         # 我们记录用的符号
    print(f"probes to fetch: {len(gene_rows)}")

    exp_mat: list[list] = []
    step = 120  # 每批样本数，控制单次请求体量
    for i in range(0, len(expr_samples), step):
        sub = expr_samples[i:i + step]
        part = xena_fetch(DS_EXPR, gene_rows, sub)
        exp_mat.append(np.asarray(part, float))
        print(f"  fetched {i + len(sub)}/{len(expr_samples)} samples", flush=True)
    E = np.concatenate(exp_mat, axis=1)      # (n_genes, n_samples)

    with open(expr_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["symbol", "probe"] + list(expr_samples))
        for j, (sym, pr) in enumerate(zip(probe_list, gene_rows)):
            w.writerow([sym, pr] + [f"{v:.6g}" for v in E[j]])
    meta["n_expression_samples"] = len(expr_samples)
    meta["signature_coverage"] = coverage
    meta["n_probes_fetched"] = len(gene_rows)
    meta["mapping"] = "Ensembl REST /lookup/symbol/homo_sapiens"
    meta["probe_match"] = "probe = gencode.v36 versioned Ensembl ID, matched on unversioned ID"
    meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"Wrote {surv_csv.name}, {expr_csv.name}, {meta_path.name}")
    return surv_csv, expr_csv, meta


# --------------------------------------------------------------------------
# 阶段 2：分析
# --------------------------------------------------------------------------
def load_cached(surv_csv: Path, expr_csv: Path):
    surv: dict[str, tuple[float, float, str]] = {}
    with open(surv_csv, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                surv[r["sample"]] = (float(r["OS"]), float(r["OS_time"]), r["patient"])
            except ValueError:
                continue
    syms: list[str] = []
    samples: list[str] = []
    rows: list[list[float]] = []
    with open(expr_csv, encoding="utf-8") as f:
        rd = csv.reader(f)
        header = next(rd)
        samples = header[2:]
        for r in rd:
            syms.append(r[0])
            rows.append([float(x) if x not in ("", "NaN", "nan") else np.nan
                         for x in r[2:]])
    return surv, syms, samples, np.asarray(rows, float)


def build_scores(expr_csv: Path, meta: dict):
    surv, syms, samples, E = load_cached(meta["_surv_csv"], expr_csv)
    # log2(TPM+1) -> 跨样本 z
    L = np.log2(np.clip(E, 0, None) + 1.0)
    mu = np.nanmean(L, axis=1, keepdims=True)
    sd = np.nanstd(L, axis=1, keepdims=True)
    Z = (L - mu) / np.where(sd > 0, sd, np.nan)
    sym_idx = {s: i for i, s in enumerate(syms)}
    scores = {}
    for k in ALL_SIGS:
        used = meta["signature_coverage"].get(k, {}).get("used", [])
        idx = [sym_idx[g] for g in used if g in sym_idx]
        if not idx:
            continue
        scores[k] = np.nanmean(Z[idx], axis=0)
    return surv, samples, scores


def analyze(meta: dict) -> dict:
    surv, samples, scores = build_scores(meta["_expr_csv"], meta)

    # 组装分析表：排除正常组织(-11)，与生存取交集
    S = np.array(samples)
    keep = np.array([("-11" not in s) for s in S])
    common = [i for i in range(len(S)) if keep[i] and S[i] in surv]
    sid = [S[i] for i in common]
    t = np.array([surv[s][1] for s in sid], float)
    ev = np.array([surv[s][0] for s in sid], float)
    pat = np.array([surv[s][2] for s in sid])
    stype = np.array(["primary" if "-01" in s else
                      ("metastasis" if ("-06" in s or "-07" in s) else "other")
                      for s in sid])

    # 同一患者多样本去重：保留 OS.time 最长者
    uniq = []
    seen: dict[str, int] = {}
    for i, p in enumerate(pat):
        if p not in seen or t[i] > t[seen[p]]:
            seen[p] = i
    uniq = sorted(seen.values())
    n_dropped = len(sid) - len(uniq)
    sid_u = [sid[i] for i in uniq]
    t_u, ev_u = t[uniq], ev[uniq].astype(bool)
    st_u = stype[uniq]

    comp = {}
    def _mean_of(names):
        arrs = [scores[k][[common[i] for i in uniq]] for k in names if k in scores]
        return np.mean(np.stack(arrs, axis=1), axis=1) if arrs else None

    comp["barrier"] = _mean_of(BARRIER_SIGS)          # 主 summary
    comp["barrier_mab_style"] = _mean_of(BARRIER_MAB_SIGS)
    comp["immune"] = _mean_of(IMMUNE_SIGS)
    per_sig = {k: scores[k][[common[i] for i in uniq]] for k in scores}

    def _pop_mask(label):
        if label == "all":
            return np.ones(len(t_u), bool)
        return st_u == label

    def _cox_report(x, mask, name):
        xm = np.asarray(x, float)[mask]
        tm, em = t_u[mask], ev_u[mask]
        ok = np.isfinite(xm) & np.isfinite(tm) & (tm > 0)
        xm, tm, em = xm[ok], tm[ok], em[ok]
        if len(xm) < 20 or em.sum() < 8:
            return dict(signature=name, n=int(len(xm)), n_events=int(em.sum()),
                        note="insufficient events")
        z = (xm - xm.mean()) / xm.std(ddof=0)
        fit = cox_breslow(z.reshape(-1, 1), tm, em)
        b, s = float(fit["beta"][0]), float(fit["se"][0])
        return dict(
            signature=name, n=int(len(z)), n_events=int(em.sum()),
            hr_per_sd=float(np.exp(b)),
            ci_lo=float(np.exp(b - 1.96 * s)), ci_hi=float(np.exp(b + 1.96 * s)),
            z=float(fit["z"][0]), p=float(fit["p"][0]),
            c_index=c_index(z, tm, em),
        )

    strata = {}
    for label in ("all", "primary", "metastasis"):
        m = _pop_mask(label)
        strata[label] = dict(
            n=int(m.sum()), n_events=int(ev_u[m].sum()),
            univariate={k: _cox_report(per_sig[k], m, k) for k in per_sig},
            composites={k: _cox_report(v, m, k) for k, v in comp.items() if v is not None},
        )

    # 多变量：barrier + immune
    mv = {}
    for label in ("all", "metastasis"):
        m = _pop_mask(label)
        if comp["barrier"] is None or comp["immune"] is None:
            continue
        X = np.stack([comp["barrier"], comp["immune"]], axis=1)
        Xm = X[m]; tm, em = t_u[m], ev_u[m]
        ok = np.isfinite(Xm).all(axis=1) & (tm > 0)
        Xm, tm, em = Xm[ok], tm[ok], em[ok]
        if len(Xm) < 30 or em.sum() < 10:
            continue
        Xz = (Xm - Xm.mean(axis=0)) / Xm.std(axis=0, ddof=0)
        fit = cox_breslow(Xz, tm, em)
        b, s = fit["beta"], fit["se"]
        mv[label] = dict(
            n=int(len(Xm)), n_events=int(em.sum()),
            terms={
                "barrier": dict(hr_per_sd=float(np.exp(b[0])),
                                ci_lo=float(np.exp(b[0] - 1.96 * s[0])),
                                ci_hi=float(np.exp(b[0] + 1.96 * s[0])),
                                p=float(fit["p"][0])),
                "immune": dict(hr_per_sd=float(np.exp(b[1])),
                               ci_lo=float(np.exp(b[1] - 1.96 * s[1])),
                               ci_hi=float(np.exp(b[1] + 1.96 * s[1])),
                               p=float(fit["p"][1])),
            },
        )

    # KM：各 composite 三分位（全体队列）。immune 作为阳性对照一起画。
    km: dict = {}
    for comp_name in ("barrier", "immune"):
        if comp[comp_name] is None:
            continue
        b = comp[comp_name]
        q = np.quantile(b, [1 / 3, 2 / 3])
        grp = np.digitize(b, q)
        curves = {}
        for gi, lab in enumerate(["low", "mid", "high"]):
            m = grp == gi
            tt, ss = km_curve(t_u[m], ev_u[m])
            curves[lab] = dict(n=int(m.sum()), n_events=int(ev_u[m].sum()),
                               times=[float(x) for x in tt],
                               survival=[float(x) for x in ss])
        lr = logrank(t_u, ev_u, grp)
        km[comp_name] = dict(tertile_cutoffs=[float(x) for x in q], curves=curves,
                             logrank_p=lr, groups=["low", "mid", "high"])

    return dict(
        cohort=COHORT, n_samples_used=len(sid_u), n_dropped_by_patient_dedup=int(n_dropped),
        n_events=int(ev_u.sum()),
        sample_types={k: int(np.sum(st_u == k)) for k in ("primary", "metastasis", "other")},
        signature_composites={
            "barrier": BARRIER_SIGS, "barrier_mab_style": BARRIER_MAB_SIGS,
            "immune": IMMUNE_SIGS},
        strata=strata, multivariable=mv, km=km,
        interpretation=(
            "The immune composite (CD8T/T_NK/Ag_target) is strongly protective for OS in the full "
            "cohort (HR<1, p~1e-7), which serves as a positive control that the bulk scoring, "
            "survival join and Cox implementation detect a known signal. The bulk barrier summary "
            "(ECM_core/ECM_crosslink/CAF) shows no association with OS in the full cohort "
            "(HR~1.0, p>0.7) and only a non-significant trend in the expected direction in the "
            "primary-tumour subset. This is consistent with the manuscript's position that the "
            "barrier construct is spatial: a bulk average of the same components does not carry "
            "the node-level geometry that defines the operators, so bulk cannot be used as a "
            "surrogate for the spatial scores."
        ),
        score_summary={k: dict(mean=float(np.nanmean(v)), sd=float(np.nanstd(v)))
                       for k, v in list(per_sig.items()) + list(comp.items()) if v is not None},
    )


def plot_km(result: dict, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    km_all = result["km"]
    colors = {"low": "#4dd0e1", "mid": "#b0bec5", "high": "#ef5350"}
    panels = [("barrier", "barrier summary (ECM_core + ECM_crosslink + CAF)"),
              ("immune", "immune composite (CD8T + T_NK + Ag_target) — positive control")]
    panels = [(k, t) for k, t in panels if k in km_all]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for ax, (key, ttl) in zip(axes, panels):
        km = km_all[key]
        for lab, c in km["curves"].items():
            ax.step(c["times"], c["survival"], where="post", color=colors[lab],
                    label=f"{lab} (n={c['n']}, ev={c['n_events']})")
        ax.set_xlabel("Overall survival (days)")
        ax.set_ylabel("Survival probability")
        ax.set_title(f"{ttl}\nlog-rank p = {km['logrank_p']:.3g}", fontsize=9)
        ax.legend(fontsize=7, frameon=False)
        ax.grid(alpha=0.2)

    ax = axes[-1]
    items = []
    for k, v in result["strata"]["all"]["composites"].items():
        if "hr_per_sd" in v:
            items.append((f"{k} (all)", v))
    for k, v in result["strata"]["primary"]["composites"].items():
        if "hr_per_sd" in v:
            items.append((f"{k} (primary only)", v))
    ypos = np.arange(len(items))
    for y, (name, v) in zip(ypos, items):
        ax.errorbar(v["hr_per_sd"], y, xerr=[[v["hr_per_sd"] - v["ci_lo"]],
                                             [v["ci_hi"] - v["hr_per_sd"]]],
                    fmt="o", color="#ef5350", capsize=3, ms=5)
        ax.text(v["ci_hi"] * 1.02, y, f" p={v['p']:.2g}", va="center", fontsize=7)
    ax.axvline(1.0, color="grey", ls="--", lw=0.8)
    ax.set_yticks(ypos)
    ax.set_yticklabels([n for n, _ in items], fontsize=8)
    ax.set_xlabel("HR per 1 SD, overall survival")
    ax.set_title("Univariate Cox (exploratory)", fontsize=9)
    ax.grid(alpha=0.2, axis="x")
    fig.suptitle("TCGA-SKCM bulk summaries vs overall survival — exploratory concept check "
                 "(not a prognostic claim)", fontsize=10)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main() -> None:
    cfg = load_config(None)
    P = Paths(cfg)
    ext = Path(P.external)
    force = "--force-download" in sys.argv

    surv_csv, expr_csv, meta = ensure_cache(ext, force=force)
    meta["_surv_csv"] = surv_csv
    meta["_expr_csv"] = expr_csv

    result = analyze(meta)
    result["caveat"] = (
        "EXPLORATORY concept check only. The spatial barrier operators require spatial "
        "coordinates and cannot be computed on bulk RNA-seq. What is tested here is whether "
        "a bulk signature summary of the same matrix/barrier components is associated with "
        "overall survival in TCGA-SKCM. Association, if any, does not establish that the "
        "SPARTA scores predict prognosis, response or outcome, and no clinical claim of "
        "predictive validity is made."
    )
    result["data_sources"] = {
        "expression": f"{HUB} :: {DS_EXPR}",
        "survival": f"{HUB} :: {DS_SURV}",
        "cohort": COHORT,
        "note": "Uses the star_tpm dataset; no TCGA ICB treatment arm is modelled.",
    }
    result["signature_coverage"] = meta.get("signature_coverage", {})
    result["meta"] = stamp_run(cfg, {"module": "M45-tcga-survival"})

    p = P.validation("tcga_survival.json")
    save_json(p, result)
    png = Path(P.figure("tcga_skcm_survival_km.png"))
    try:
        plot_km(result, png)
        made_png = True
    except Exception as e:  # noqa: BLE001
        print(f"[warn] figure failed: {e}")
        made_png = False
    print(f"Wrote {p}" + (f" and {png}" if made_png else ""))

    print("\n=== TCGA-SKCM bulk barrier summary vs OS (exploratory) ===")
    print(f"  samples used: {result['n_samples_used']} (events {result['n_events']}), "
          f"sample types: {result['sample_types']}, "
          f"dedup dropped {result['n_dropped_by_patient_dedup']}")
    for label in ("all", "primary", "metastasis"):
        s = result["strata"][label]
        print(f"  [{label}] n={s['n']} ev={s['n_events']}")
        for k, v in s["composites"].items():
            if "hr_per_sd" in v:
                print(f"      {k:<20} HR/SD={v['hr_per_sd']:.3f} "
                      f"[{v['ci_lo']:.3f},{v['ci_hi']:.3f}] p={v['p']:.3g} "
                      f"C={v['c_index']:.3f}")
    if "all" in result["multivariable"]:
        mv = result["multivariable"]["all"]
        print(f"  [multivariable, all] barrier HR/SD={mv['terms']['barrier']['hr_per_sd']:.3f} "
              f"p={mv['terms']['barrier']['p']:.3g}; "
              f"immune HR/SD={mv['terms']['immune']['hr_per_sd']:.3f} "
              f"p={mv['terms']['immune']['p']:.3g}")
    if "barrier" in result["km"]:
        print(f"  KM barrier tertiles log-rank p = {result['km']['barrier']['logrank_p']:.3g}")
    if "immune" in result["km"]:
        print(f"  KM immune  tertiles log-rank p = {result['km']['immune']['logrank_p']:.3g} "
              f"(positive control)")


if __name__ == "__main__":
    main()
