"""
原始数据摄入与格式适配
========================

上游依赖：无（这是流水线的最前端）
下游模块：run_00b_ingest.py → run_01_qc.py

为什么需要这个模块
------------------
公开空间转录组数据的格式是**混乱的**。`sc.read_visium()` 只能读 spaceranger
的标准输出，而你要用的皮肤肿瘤数据里有相当一部分不是这个格式：

  · 第一代 ST 平台（Ji et al. cSCC、Thrane et al. melanoma）通常给的是
    一个 TSV 计数矩阵，spot 名字形如 "10x20"（坐标编码在名字里），
    有的另给一个坐标文件，有的没有；
  · GEO 的补充文件经常是 genes×spots 而不是 spots×genes，需要转置；
  · 有的只给处理后的 h5ad，有的给 mtx 三件套；
  · 基因名有的是 symbol，有的是 Ensembl ID。

如果不先解决摄入，你会在第二周就卡住，而且是以"每个数据集单独写一段代码"
的方式卡住——三个月后没人（包括你自己）看得懂那些脚本。

本模块的策略
------------
`inspect_raw()`  只看目录结构与文件头，**不加载数据**，快速判断这份数据
                 是什么格式、三要素齐不齐。用于建数据台账。
`load_raw()`     真正加载并转成统一的 AnnData（含 obsm["spatial"]）。

两者分开，是因为筛选阶段你要过几十个候选数据集，不该为了知道"能不能用"
就把每个都完整加载一遍。
"""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np

__all__ = ["inspect_raw", "load_raw", "read_table_matrix", "FORMATS", "parse_spot_coords"]

FORMATS = ("visium", "legacy_st", "matrix_coords", "h5ad", "mtx", "unknown")

# 常见的文件名模式
_COUNT_PAT = re.compile(
    r"(count|matrix|expression|expr|umi|gene[_-]?bc|filtered|raw)", re.I)
_COORD_PAT = re.compile(
    r"(tissue_position|position|coord|spot|barcode.*loc|location)", re.I)
_IMAGE_PAT = re.compile(r"\.(png|jpg|jpeg|tif|tiff)(\.gz)?$", re.I)
_HE_PAT = re.compile(r"(hires|lowres|tissue|he|h&e|histo|image)", re.I)

# 形如 "10x20" / "10.5x20.3" / "10_20" 的 spot 名（第一代 ST 常见）
_SPOTNAME_PAT = re.compile(r"^(\d+(?:\.\d+)?)[x_×](\d+(?:\.\d+)?)$")


# --------------------------------------------------------------------------
def parse_spot_coords(names) -> np.ndarray | None:
    """从形如 "10x20" 的 spot 名中解析坐标。

    第一代 ST 平台（以及不少 GEO 上的补充矩阵）把坐标编码在行名/列名里。
    解析成功返回 (n, 2) 数组，失败返回 None。

    注意：这里解析出来的是**阵列坐标**（第几行第几列），不是微米也不是像素。
    run_01_qc.py 会用最近邻间距把它标定到微米。
    """
    names = [str(x).strip() for x in names]
    out = []
    for s in names:
        m = _SPOTNAME_PAT.match(s)
        if not m:
            return None
        out.append((float(m.group(1)), float(m.group(2))))
    return np.asarray(out, dtype=float)


def _open_maybe_gz(p: Path, mode="rt"):
    return gzip.open(p, mode) if p.suffix == ".gz" else open(p, mode, encoding="utf-8")


def _peek_table(p: Path, n_lines: int = 3):
    """只读表格前几行，判断分隔符、方向、维度量级。不加载全表。"""
    try:
        with _open_maybe_gz(p) as f:
            lines = [f.readline().rstrip("\n") for _ in range(n_lines)]
    except Exception:  # noqa: BLE001
        return None
    lines = [x for x in lines if x]
    if not lines:
        return None
    sep = "\t" if lines[0].count("\t") > lines[0].count(",") else ","
    header = lines[0].split(sep)
    first_col = [ln.split(sep)[0] for ln in lines[1:]]
    return dict(sep=sep, n_header=len(header), header=header[:8],
                first_col=first_col, path=p)


# --------------------------------------------------------------------------
def inspect_raw(path: str | Path) -> dict:
    """快速检查一份原始数据，**不加载矩阵**。用于建数据台账。

    参数
    ----
    path : 数据集目录，或单个 .h5ad 文件

    返回
    ----
    dict，含：
      format        : FORMATS 之一
      has_counts    : 是否找到计数矩阵
      has_coords    : 是否能拿到空间坐标（含"编码在 spot 名里"的情况）
      has_image     : 是否找到组织学图像
      coord_source  : 坐标从哪来（"file" / "spot_names" / "obsm" / None）
      files         : 识别到的关键文件路径
      notes         : 人类可读的说明与警告

    典型用法
    --------
        for d in Path("data/raw").iterdir():
            info = inspect_raw(d)
            print(d.name, info["format"], info["has_counts"],
                  info["has_coords"], info["has_image"])
    """
    p = Path(path)
    res = dict(format="unknown", has_counts=False, has_coords=False,
               has_image=False, coord_source=None, files={}, notes=[])

    if not p.exists():
        res["notes"].append(f"路径不存在：{p}")
        return res

    # ---- 单个 h5ad ----
    if p.is_file() and p.suffix == ".h5ad":
        res["format"] = "h5ad"
        res["has_counts"] = True
        res["files"]["h5ad"] = str(p)
        try:
            import anndata as ad
            a = ad.read_h5ad(p, backed="r")
            if "spatial" in a.obsm:
                res["has_coords"] = True
                res["coord_source"] = "obsm"
            res["notes"].append(f"{a.n_obs} spot × {a.n_vars} 基因")
            if a.uns.get("spatial"):
                res["has_image"] = True
            a.file.close()
        except Exception as e:  # noqa: BLE001
            res["notes"].append(f"无法读取 h5ad 头部：{e}")
        return res

    if not p.is_dir():
        res["notes"].append(f"既不是目录也不是 .h5ad：{p}")
        return res

    files = sorted([f for f in p.rglob("*") if f.is_file()])
    names = {f.name.lower(): f for f in files}

    # ---- Visium (spaceranger) ----
    has_spatial_dir = any("spatial" in str(f.parent).lower() for f in files)
    has_h5 = any(f.suffix == ".h5" and _COUNT_PAT.search(f.name) for f in files)
    has_tp = any(_COORD_PAT.search(f.name) and f.suffix in (".csv", ".parquet", ".gz")
                 for f in files)
    if has_spatial_dir and (has_h5 or "matrix.mtx.gz" in names or "matrix.mtx" in names):
        res.update(format="visium", has_counts=True)
        if has_tp:
            res["has_coords"] = True
            res["coord_source"] = "file"
        res["has_image"] = any(_IMAGE_PAT.search(f.name) and _HE_PAT.search(f.name)
                               for f in files)
        res["files"] = {"dir": str(p)}
        res["notes"].append("识别为 spaceranger 标准输出，可直接用 sc.read_visium")
        return res

    # ---- mtx 三件套（无 spatial 目录）----
    if ("matrix.mtx" in names or "matrix.mtx.gz" in names):
        res.update(format="mtx", has_counts=True)
        res["files"] = {"dir": str(p)}
        coord_f = [f for f in files if _COORD_PAT.search(f.name)
                   and f.suffix in (".csv", ".tsv", ".txt", ".gz")]
        if coord_f:
            res["has_coords"] = True
            res["coord_source"] = "file"
            res["files"]["coords"] = str(coord_f[0])
        res["has_image"] = any(_IMAGE_PAT.search(f.name) for f in files)
        return res

    # ---- 表格矩阵（第一代 ST / GEO 补充文件）----
    tabs = [f for f in files
            if f.suffix in (".tsv", ".csv", ".txt") or
            (f.suffix == ".gz" and f.stem.endswith((".tsv", ".csv", ".txt")))]
    count_tabs = [f for f in tabs if _COUNT_PAT.search(f.name)] or tabs
    if count_tabs:
        # 取最大的那个表当计数矩阵
        cf = max(count_tabs, key=lambda f: f.stat().st_size)
        pk = _peek_table(cf)
        res["files"]["counts"] = str(cf)
        res["has_counts"] = True
        if pk:
            sep_name = "TAB" if pk["sep"] == "\t" else "COMMA"
            res["notes"].append(f"计数表 {cf.name}：{pk['n_header']} 列，分隔符 {sep_name}")
            # 坐标是否编码在表头里（spots 在列）或首列里（spots 在行）
            if parse_spot_coords(pk["header"][1:]) is not None:
                res.update(has_coords=True, coord_source="spot_names", format="legacy_st")
                res["notes"].append("坐标编码在**列名**中（spots 在列，需转置）")
            elif pk["first_col"] and parse_spot_coords(pk["first_col"]) is not None:
                res.update(has_coords=True, coord_source="spot_names", format="legacy_st")
                res["notes"].append("坐标编码在**行名**中（spots 在行）")

        coord_f = [f for f in tabs if _COORD_PAT.search(f.name) and f != cf]
        if coord_f:
            # 即使坐标已经能从 spot 名里解析出来，这个文件也必须记下来。
            # 第一代 ST（GSE144239）的 *spot_data-selection-* 是"哪些阵列点真的
            # 落在组织上"的唯一依据；早期版本这里有个 `and not has_coords` 的短路，
            # 于是 1933 个阵列点全被当成组织，图建在一大片背景上而且一个错都不报。
            res["files"]["coords"] = str(coord_f[0])
            if not res["has_coords"]:
                res.update(has_coords=True, coord_source="file", format="matrix_coords")
                res["notes"].append(f"坐标文件：{coord_f[0].name}")
            else:
                res["notes"].append(
                    f"点位文件：{coord_f[0].name}（坐标仍取自 spot 名，此文件用于裁掉组织外的点）")

        if res["format"] == "unknown":
            res["format"] = "matrix_coords" if res["has_coords"] else "unknown"

    res["has_image"] = any(_IMAGE_PAT.search(f.name) for f in files)
    if res["has_image"]:
        img = [f for f in files if _IMAGE_PAT.search(f.name)]
        res["files"]["image"] = str(max(img, key=lambda f: f.stat().st_size))

    if not res["has_counts"]:
        res["notes"].append("未找到可识别的计数矩阵。请检查目录内容，"
                            "或手工整理成 spots×genes 的 csv/tsv 后重试。")
    if not res["has_coords"]:
        res["notes"].append("⚠ 未找到空间坐标 —— 这是准入标准 C2，不合格则整份数据不可用。")
    if not res["has_image"]:
        res["notes"].append("⚠ 未找到组织学图像 —— 准入标准 C3，"
                            "该切片无法做形态学验证与封锁线叠加图。")
    return res


# --------------------------------------------------------------------------
def _attach_visium_spatial(a, p):
    """为 mtx 目录形态的 Visium 数据补上空间坐标、缩放因子与组织学图像。

    兼容 spaceranger 目录（filtered_feature_bc_matrix/ + spatial/）：
      · tissue_positions_list.csv —— 每个捕获位点的行列坐标与像素坐标
      · scalefactors_json.json   —— 图像缩放因子
      · tissue_hires_image.png   —— 高分辨率组织学图像
    """
    import json
    spatial_dir = p / "spatial"
    tp = spatial_dir / "tissue_positions_list.csv"
    if tp.exists():
        import pandas as pd
        # 旧版 spaceranger 输出无表头：barcode,in_tissue,array_row,array_col,pxl_col,pxl_row
        raw = tp.read_text().splitlines()
        if raw and raw[0].startswith("barcode"):
            pos = pd.read_csv(tp, index_col=0)
            xcol, ycol = "pxl_col_in_fullres", "pxl_row_in_fullres"
        else:
            pos = pd.read_csv(tp, header=None,
                              names=["barcode", "in_tissue", "array_row",
                                     "array_col", "pxl_col", "pxl_row"])
            pos = pos.set_index("barcode")
            xcol, ycol = "pxl_col", "pxl_row"
        pos.index = pos.index.astype(str)
        common = [b for b in a.obs_names if b in pos.index]
        if common:
            sub = pos.loc[common]
            a.obsm["spatial"] = sub[[xcol, ycol]].to_numpy(float)
    # 缩放因子与图像（供 run_01 标定与 run_06 形态学验证）
    sf = spatial_dir / "scalefactors_json.json"
    if sf.exists():
        try:
            a.uns["scalefactors"] = json.loads(sf.read_text())
        except Exception:  # noqa: BLE001
            pass
    img = spatial_dir / "tissue_hires_image.png"
    if img.exists():
        try:
            from PIL import Image
            import numpy as np
            a.uns["hires_image"] = np.asarray(Image.open(img))
        except Exception:  # noqa: BLE001
            pass


def load_raw(path: str | Path, *, coords_file: str | None = None,
             transpose: bool | None = None, gene_col: str | None = None):
    """把一份原始数据加载为统一的 AnnData。

    返回的 AnnData 保证：
      · adata.X 是计数（spots × genes）
      · adata.obsm["spatial"] 存在（阵列坐标或像素坐标，单位由 run_01 标定）
      · adata.var_names 是基因名

    参数
    ----
    path        : 数据集目录或 .h5ad
    coords_file : 显式指定坐标文件。不给则按 inspect_raw 的判断自动找
    transpose   : 强制转置。None 时自动判断（看行名像基因还是像 spot）
    gene_col    : 若计数表有多列注释，指定哪一列是基因名

    注意
    ----
    自动判断有出错的可能。加载后**务必**检查：
      · adata.n_obs 是 spot 数（几百到几千），adata.n_vars 是基因数（上万）
      · 若两者反了，说明转置判断错了，手工传 transpose=True/False
    """
    import anndata as ad
    import pandas as pd

    p = Path(path)
    info = inspect_raw(p)

    # ---- h5ad ----
    if info["format"] == "h5ad":
        return ad.read_h5ad(p)

    # ---- Visium ----
    if info["format"] == "visium":
        import scanpy as sc
        # spaceranger 输出有两种形态：filtered_feature_bc_matrix.h5（h5）
        # 或 filtered_feature_bc_matrix/ 目录（mtx 三件套）。
        # 老版本 scanpy 的 read_visium 只认 h5，这里手动兼容 mtx 目录形态。
        h5 = list(p.glob("filtered_feature_bc_matrix.h5")) or \
             list(p.glob("filtered_feature_bc_matrix/*.h5"))
        if not h5:
            mtx_dir = p / "filtered_feature_bc_matrix"
            if (mtx_dir / "matrix.mtx").exists() or (mtx_dir / "matrix.mtx.gz").exists():
                a = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols", make_unique=True)
                _attach_visium_spatial(a, p)
                return a
        return sc.read_visium(p)

    # ---- mtx ----
    if info["format"] == "mtx":
        import scanpy as sc
        a = sc.read_10x_mtx(p, var_names="gene_symbols", make_unique=True)
        cf = coords_file or info["files"].get("coords")
        if cf:
            a.obsm["spatial"] = _read_coords(cf, a.obs_names)
        return a

    # ---- 表格矩阵 ----
    vals, obs, var, xy = read_table_matrix(
        info, coords_file=coords_file, transpose=transpose, gene_col=gene_col)
    a = ad.AnnData(vals)
    a.obs_names = obs
    a.var_names = var
    a.var_names_make_unique()
    a.obsm["spatial"] = xy
    return a



_GENE_PAT = re.compile(r"^[A-Z][A-Z0-9._-]{1,14}$")


def _looks_like_genes(labels) -> float:
    """一组标签里有多大比例像基因 symbol。用于方向判断的兜底信号。"""
    labels = [str(x) for x in labels][:200]
    if not labels:
        return 0.0
    return sum(bool(_GENE_PAT.match(x)) for x in labels) / len(labels)


def _decide_orientation(df, coords_file):
    """判断表格是 spots×genes 还是 genes×spots，返回 (need_transpose, 判断依据)。

    按可靠性排序使用四种信号：
      1. 坐标文件的 ID 与哪一维对得上 —— **确定性信号，最可靠**
      2. 行名/列名是否形如 "10x20"（第一代 ST 把坐标编码在名字里）
      3. 哪一维更像基因 symbol
      4. 维度大小（真实数据里基因数远多于 spot 数）—— 最不可靠，只作兜底

    早期版本只用信号 4，在"基因被过滤到很少"的数据上会判反。
    """
    idx, col = list(df.index), list(df.columns)

    # 信号 1：坐标文件
    # 选点文件（x/y 是阵列坐标，第一列是纯数字）没有 ID 列，
    # 拿它的第一列去跟行名/列名比对没有意义，直接跳过。
    if coords_file and _read_spot_selection(coords_file) is None:
        try:
            import pandas as pd
            cp = Path(coords_file)
            sep = "\t" if cp.name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")) else ","
            cdf = pd.read_csv(cp, sep=sep, nrows=500)
            ids = set(str(x) for x in cdf.iloc[:, 0])
            hit_idx = sum(str(x) in ids for x in idx[:200]) / max(min(len(idx), 200), 1)
            hit_col = sum(str(x) in ids for x in col[:200]) / max(min(len(col), 200), 1)
            if max(hit_idx, hit_col) > 0.5:
                return (hit_col > hit_idx), "坐标文件的 ID 与列名匹配" if hit_col > hit_idx \
                    else "坐标文件的 ID 与行名匹配"
        except Exception:  # noqa: BLE001
            pass

    # 信号 2：坐标编码在名字里
    row_is_spot = parse_spot_coords(idx[:20]) is not None
    col_is_spot = parse_spot_coords(col[:20]) is not None
    if col_is_spot and not row_is_spot:
        return True, "列名形如 10x20，spots 在列"
    if row_is_spot and not col_is_spot:
        return False, "行名形如 10x20，spots 在行"

    # 信号 3：哪一维更像基因名
    g_idx, g_col = _looks_like_genes(idx), _looks_like_genes(col)
    if abs(g_idx - g_col) > 0.3:
        return (g_idx > g_col), ("行名更像基因 symbol" if g_idx > g_col
                                 else "列名更像基因 symbol")

    # 信号 4：兜底
    return (df.shape[0] > df.shape[1]), "维度启发式（最不可靠，请核对 n_obs/n_vars）"


def read_table_matrix(info: dict, *, coords_file=None, transpose=None, gene_col=None):
    """把表格式的计数矩阵读成 (values, obs_names, var_names, coords)。

    从 load_raw 中拆出来，是为了让方向判断与坐标解析这两处**最容易出错**的
    逻辑可以在没有 anndata 的环境下单独测试（见 tests/test_loaders.py）。

    方向判断的规则（按优先级）：
      1. 行名像 "10x20" -> spots 在行，不转置
      2. 列名像 "10x20" -> spots 在列，转置
      3. 都不像 -> 用维度启发式：基因数通常远多于 spot 数

    规则 3 会出错（比如一张只有 300 spot 的第一代 ST 切片、基因也被过滤到
    只剩几百个）。加载后**必须**检查 n_obs 是 spot 数、n_vars 是基因数，
    反了就手工传 transpose=。
    """
    import pandas as pd

    cf = info["files"].get("counts")
    if cf is None:
        raise FileNotFoundError(
            f"找不到可识别的计数矩阵。inspect_raw 的说明：{info['notes']}")

    sep = "\t" if Path(cf).name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")) else ","
    df = pd.read_csv(cf, sep=sep, index_col=0)
    if gene_col and gene_col in df.columns:
        df = df.set_index(gene_col)

    cfile = coords_file or info["files"].get("coords")

    if transpose is None:
        transpose, why = _decide_orientation(df, cfile)
        info.setdefault("notes", []).append(f"方向判断：{'转置' if transpose else '不转置'}（依据：{why}）")

    if transpose:
        df = df.T

    obs = [str(x) for x in df.index]
    var = [str(x) for x in df.columns]
    vals = df.values.astype(np.float32)

    xy = parse_spot_coords(obs)
    if xy is None:
        if not cfile:
            raise ValueError(
                "没有可用的空间坐标（既不在 spot 名里，也没有坐标文件）。"
                "这是准入标准 C2，该切片不可用。")
        xy = _read_coords(cfile, obs)

    # ---- 第一代 ST：按选点文件裁掉组织外的点 ----
    # 坐标能从 spot 名里解析出来时，上面那个分支不会去读坐标文件，
    # 于是选点信息就被跳过了。必须在这里单独处理一次。
    sel = _read_spot_selection(cfile)
    if sel is not None:
        keep = np.array([o in sel for o in obs], dtype=bool)
        n_drop = int((~keep).sum())

        # 反向对账：选点文件里有多少个点在计数矩阵里根本找不到。
        # 两边坐标约定不一致（比如一边 1-based 一边 0-based、或 x/y 反了）时，
        # 下面那个 keep.sum() < 10 的门槛拦不住"恰好对上四成"这种中间态——
        # 那会安静地留下一堆错点，图照样建得出来，只是全错。
        orphan = sel - set(obs)
        if len(orphan) > 0.2 * len(sel):
            raise ValueError(
                f"选点文件里有 {len(orphan)}/{len(sel)} 个点在计数矩阵里找不到"
                f"（超过两成），两边的坐标约定多半对不上，不能就这么往下走。"
                f"矩阵前几个名字：{obs[:3]}；对不上的前几个：{sorted(orphan)[:3]}。")
        if orphan:
            info.setdefault("notes", []).append(
                f"选点文件有 {len(orphan)}/{len(sel)} 个点不在计数矩阵里，已忽略")
            print(f"[loaders] 注意：选点文件有 {len(orphan)}/{len(sel)} 个点"
                  f"在计数矩阵里找不到，已忽略")

        if keep.sum() < 10:
            raise ValueError(
                f"选点文件与计数矩阵的 spot 名几乎对不上（只剩 {int(keep.sum())} 个）。"
                f"矩阵前几个名字：{obs[:3]}；选点文件里前几个：{sorted(sel)[:3]}。")
        if n_drop:
            info.setdefault("notes", []).append(
                f"按选点文件裁掉组织外的 {n_drop} 个 spot，保留 {int(keep.sum())} 个")
            print(f"[loaders] 选点过滤：{len(obs)} -> {int(keep.sum())} 个 spot"
                  f"（裁掉组织外 {n_drop} 个）")
        vals, obs, xy = vals[keep], [o for o, k in zip(obs, keep) if k], xy[keep]

    return vals, obs, var, xy


def _read_spot_selection(path):
    """第一代 ST 的选点文件 -> 组织上 spot 名的集合（形如 "10x20"）。

    为什么必须有这一步
    ------------------
    GSE144239（Ji et al. 2020）的每张第一代 ST 切片有两个文件：
      · `*_stdata.tsv`               **整块阵列**的计数，行名形如 "10x10"
      · `*spot_data-selection-*.tsv`  真正落在组织上的点，列是 x/y/new_x/new_y/pixel_x/pixel_y
    两者数量差很多——P2_rep1 是 1933 个阵列点里只有 665 个在组织上。

    而 read_table_matrix 里"坐标能从 spot 名解析出来就不读坐标文件"这条捷径，
    正好会把选点信息整个跳过。后果是图建在一大片背景上：
    源汇的分位数阈值、连通性、最小割全部失去意义，
    **而且一个错都不报**——只会安静地给出一堆看着挺合理的数字。
    （2026-08-27 纳入 P2/P5/P9/P10 时发现。）

    返回 None 表示"这不是一个第一代 ST 的选点文件"（比如 Visium 的
    tissue_positions），此时调用方不做任何过滤。
    """
    if not path:
        return None
    import pandas as pd

    p = Path(path)
    if "selection" not in p.name.lower():
        return None
    try:
        sep = "\t" if p.name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")) else ","
        df = pd.read_csv(p, sep=sep)
    except Exception:  # noqa: BLE001
        return None

    cols = {str(c).lower(): c for c in df.columns}
    if "x" not in cols or "y" not in cols:
        return None
    # 有 barcode/spot 这类 ID 列的是别的格式，交给 _read_coords 处理
    if any(k in cols for k in ("barcode", "spot", "id", "name", "cell")):
        return None

    out = set()
    for xv, yv in zip(df[cols["x"]], df[cols["y"]]):
        # 阵列坐标可能写成 7 或 7.0，spot 名里是整数
        out.add(f"{int(round(float(xv)))}x{int(round(float(yv)))}")
    return out or None


def _read_coords(path, obs_names) -> np.ndarray:
    """读坐标文件并按 obs_names 对齐。

    支持常见的几种列名组合：
      (barcode, x, y) / (spot, pxl_col, pxl_row) / Visium 的 tissue_positions
    """
    import pandas as pd

    p = Path(path)
    sep = "\t" if p.name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")) else ","
    df = pd.read_csv(p, sep=sep)

    # Visium 的 tissue_positions_list.csv 无表头
    if df.shape[1] >= 6 and not any(isinstance(c, str) and c.isalpha()
                                    for c in df.columns[:1]):
        df = pd.read_csv(p, sep=sep, header=None)
        df.columns = ["barcode", "in_tissue", "array_row", "array_col",
                      "pxl_row_in_fullres", "pxl_col_in_fullres"][:df.shape[1]]

    cols = {c.lower(): c for c in df.columns}
    id_col = next((cols[k] for k in ("barcode", "spot", "id", "name", "cell")
                   if k in cols), df.columns[0])
    x_col = next((cols[k] for k in ("x", "pxl_col_in_fullres", "imagecol",
                                    "array_col", "col") if k in cols), None)
    y_col = next((cols[k] for k in ("y", "pxl_row_in_fullres", "imagerow",
                                    "array_row", "row") if k in cols), None)
    if x_col is None or y_col is None:
        raise ValueError(f"坐标文件 {p} 中找不到 x/y 列。实际列名：{list(df.columns)}")

    df = df.set_index(id_col)
    idx = [str(x) for x in obs_names]
    missing = [i for i in idx if i not in df.index.astype(str)]
    if len(missing) > len(idx) * 0.5:
        raise ValueError(
            f"坐标文件与计数矩阵的 spot 名对不上（{len(missing)}/{len(idx)} 缺失）。"
            f"坐标文件前几个名字：{list(df.index[:3])}；"
            f"矩阵前几个：{idx[:3]}。请检查是否需要去掉后缀或前缀。"
        )
    df.index = df.index.astype(str)
    sub = df.reindex(idx)
    return sub[[x_col, y_col]].values.astype(float)
