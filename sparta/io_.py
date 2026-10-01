"""
数据契约的统一读写与路径管理
==============================

所有磁盘 IO 都走这里。理由不是洁癖，而是：
  1. 路径统一从 config 读，换机器、换盘符只改一个 yaml；
  2. 中间产物的命名规则集中在一处，改名不会漏掉某个脚本；
  3. 每次写盘同时记录随机种子与所用配置，结果可追溯。

文件契约（模块间的接口，改动前请先更新 CLAUDE.md）
----------------------------------------------------
  M1  ->  interim/{slide_id}.qc.h5ad
  M2  ->  interim/{slide_id}.scored.h5ad
  M3  ->  interim/{slide_id}.graph.npz      (A, D, source, sink, vessel)
  M4  ->  interim/{slide_id}.barrier.npz    (b_cell, b_mab, b_meta, ...)
          interim/{slide_id}.mincut.json    (最小割边集，供 viz)
  M5  ->  results/counterfactual/{slide_id}.json
  M6  ->  results/validation/*
"""

from __future__ import annotations

import json
import os
import random
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import yaml

__all__ = [
    "load_config", "set_seed", "Paths",
    "save_graph", "load_graph", "save_json", "load_json", "stamp_run",
]

DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "configs" / "default.yaml"


# --------------------------------------------------------------------------
def load_config(path: str | os.PathLike | None = None) -> dict:
    """读取 YAML 配置。不传路径则用 configs/default.yaml。"""
    path = Path(path) if path else DEFAULT_CONFIG
    if not path.exists():
        raise FileNotFoundError(f"找不到配置文件 {path}")
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    cfg["_config_path"] = str(path)
    return cfg


def set_seed(seed: int = 0) -> int:
    """统一设置随机种子。反事实实验涉及随机重排，不固定种子将无法复现。"""
    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    return seed


class Paths:
    """按 config 组织的路径管理器。所有目录不存在时自动创建。"""

    def __init__(self, cfg: dict):
        root = Path(cfg["paths"]["root"]).expanduser().resolve()
        self.root = root
        self.raw = root / cfg["paths"]["raw"]
        self.interim = root / cfg["paths"]["interim"]
        self.external = root / cfg["paths"]["external"]
        self.results = root / cfg["paths"]["results"]
        for p in (self.raw, self.interim, self.external, self.results):
            p.mkdir(parents=True, exist_ok=True)

    # --- 各模块的产物路径 ---
    def qc(self, sid):        return self.interim / f"{sid}.qc.h5ad"
    def scored(self, sid):    return self.interim / f"{sid}.scored.h5ad"
    def graph(self, sid):     return self.interim / f"{sid}.graph.npz"
    def barrier(self, sid):   return self.interim / f"{sid}.barrier.npz"
    def mincut(self, sid):    return self.interim / f"{sid}.mincut.json"
    def counterfactual(self, sid):
        d = self.results / "counterfactual"; d.mkdir(parents=True, exist_ok=True)
        return d / f"{sid}.json"
    def validation(self, name):
        d = self.results / "validation"; d.mkdir(parents=True, exist_ok=True)
        return d / name
    def figure(self, name):
        d = self.results / "figures"; d.mkdir(parents=True, exist_ok=True)
        return d / name


# --------------------------------------------------------------------------
def save_graph(path, A: sp.spmatrix, D: sp.spmatrix, source, sink, vessel, meta: dict | None = None):
    """把邻接矩阵、距离矩阵与源汇索引存成一个 npz。"""
    A, D = A.tocoo(), D.tocoo()
    np.savez_compressed(
        path,
        A_row=A.row, A_col=A.col, A_data=A.data, A_shape=np.array(A.shape),
        D_row=D.row, D_col=D.col, D_data=D.data,
        source=np.asarray(source, int), sink=np.asarray(sink, int),
        vessel=np.asarray(vessel, int),
        meta=np.array([json.dumps(meta or {}, ensure_ascii=False)], dtype=object),
    )


def load_graph(path):
    """读回 save_graph 存的图。返回 (A, D, source, sink, vessel, meta)。"""
    z = np.load(path, allow_pickle=True)
    shape = tuple(z["A_shape"])
    A = sp.coo_matrix((z["A_data"], (z["A_row"], z["A_col"])), shape=shape).tocsr()
    D = sp.coo_matrix((z["D_data"], (z["D_row"], z["D_col"])), shape=shape).tocsr()
    meta = json.loads(str(z["meta"][0])) if "meta" in z else {}
    return A, D, z["source"], z["sink"], z["vessel"], meta


def save_json(path, obj):
    """存 JSON，自动把 numpy 类型转成 Python 原生类型。"""
    def _default(o):
        if isinstance(o, (np.integer,)):   return int(o)
        if isinstance(o, (np.floating,)):  return float(o)
        if isinstance(o, np.ndarray):      return o.tolist()
        raise TypeError(f"无法序列化 {type(o)}")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1, default=_default)


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def stamp_run(cfg: dict, extra: dict | None = None) -> dict:
    """生成一份运行元数据，随每次结果一起归档。

    审稿人要求复现时，这份记录是唯一的救命稻草。它至少要能回答：
    用了哪份配置、什么种子、什么时间、什么版本的依赖。
    """
    import platform
    import scipy
    meta = dict(
        timestamp=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        config_path=cfg.get("_config_path"),
        seed=cfg.get("seed"),
        python=platform.python_version(),
        numpy=np.__version__,
        scipy=scipy.__version__,
    )
    try:
        import networkx as nx
        meta["networkx"] = nx.__version__
    except Exception:  # noqa: BLE001
        pass
    if extra:
        meta.update(extra)
    return meta


# --------------------------------------------------------------------------
# 数据台账：切片 -> 患者
# --------------------------------------------------------------------------
def load_ledger(P) -> list[dict]:
    """读 data/ledger.csv，返回逐行的 dict 列表。台账不存在时返回空列表。

    用 stdlib csv 而不是 pandas，免得核心 IO 模块多一个依赖。
    """
    import csv as _csv

    p = P.root / "data" / "ledger.csv"
    if not p.exists():
        return []
    with open(p, encoding="utf-8") as f:
        return list(_csv.DictReader(f))


def admitted_slides(P, cancer: str | None = None) -> list[str]:
    """台账里 status == "ingested" 的切片 ID，cSCC 在前、黑色素瘤在后，各自按名字排。

    为什么要有这个函数（2026-08-27）
    --------------------------------
    切片名单从前散在各脚本里硬编码（run_12 的 DEFAULT_SLIDES、run_16 的 ORDER
    等等）。队列一变就得逐个记得去改，漏改一处的表现是**图和统计安静地少几张切片**，
    没有任何报错。CSCC13 被准入剔除之后，这个隐患还多了一面：
    剔除决定只写在台账与 admission.json 里，硬编码名单不会跟着变，
    于是被剔除的切片可能又从某个名单里溜回下游。

    以台账为唯一事实来源，两个方向的漏就都堵住了。
    """
    rows = load_ledger(P)
    sids = [r["slide_id"] for r in rows
            if str(r.get("status", "")).strip() == "ingested"
            and (cancer is None or str(r.get("cancer_type", "")) == cancer)]
    cscc = sorted(s for s in sids if s.upper().startswith("CSCC"))
    mel = sorted(s for s in sids if s.upper().startswith("MEL"))
    rest = sorted(s for s in sids if s not in set(cscc) | set(mel))
    return cscc + mel + rest


def patient_map(P) -> dict:
    """切片 ID -> 患者 ID。

    **为什么这个函数很重要。** 同一个患者的多张切片（技术重复、相邻切面、
    同一肿瘤的不同区域）不是独立样本。把它们当 n 个独立观测做统计，
    是把伪重复（pseudo-replication）当样本量，审稿人一查 GEO 元数据就会发现。

    2026-08-27 查 GEO 才发现：本项目的 8 张切片实际来自 **3 个患者**——
    CSCC01/02 = GSE144239 的 P4 两个重复，CSCC03/04 = P6 两个重复，
    MEL01–04 = GSE250636 同一位患者（Patient B）的四处颅外转移灶。
    所有"n=8"的陈述都必须同时给出患者数。

    台账里没有 patient 列时，退化为"每张切片自成一个患者"并打印告警——
    不要静默假装它们独立。
    """
    rows = load_ledger(P)
    if not rows:
        return {}
    if "patient" not in rows[0] or not any(r.get("patient") for r in rows):
        print("[warn] data/ledger.csv 里没有 patient 列。"
              "统计将按切片计数，这会把同一患者的重复切片当成独立样本。"
              "请用 run_00b_ingest.py --patient 补上，或直接编辑台账。")
        return {r["slide_id"]: r["slide_id"] for r in rows}
    return {r["slide_id"]: (r.get("patient") or r["slide_id"]) for r in rows}
