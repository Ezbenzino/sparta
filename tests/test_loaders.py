"""
loaders.py 的格式适配测试
===========================

运行方式：
    pytest tests/test_loaders.py -v
    python tests/test_loaders.py

公开空间数据的格式是混乱的，摄入环节出错会让后面**所有**结果静默失效
（比如把 genes×spots 当成 spots×genes，你会得到 20000 个"spot"和 300 个"基因"，
而下游每一步都能跑通，只是全错）。这里用四种模拟格式检验自动识别。
"""
import os
import gzip
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from sparta.loaders import inspect_raw, parse_spot_coords, read_table_matrix  # noqa: E402

GENES = ["MLANA", "PMEL", "TYR", "COL1A1", "COL1A2", "DCN", "PECAM1", "VWF",
         "CD3D", "CD8A", "LOX", "LOXL2", "CD274", "MKI67", "ABCB1", "VEGFA"]
N_SPOT = 400
N_ON_TISSUE = 120     # legacy_sel 里真正落在组织上的点数


def _mock_root():
    """造四种格式的模拟原始数据。返回临时目录（调用方负责清理）。"""
    root = Path(tempfile.mkdtemp(prefix="sparta_mock_"))
    rng = np.random.default_rng(0)
    spots = [f"{i}x{j}" for i in range(5, 25) for j in range(5, 25)]
    M = rng.poisson(4, size=(len(GENES), len(spots)))
    bc = [f"BC{i:04d}" for i in range(len(spots))]
    xy = np.array([[float(s.split("x")[0]), float(s.split("x")[1])] for s in spots])

    d = root / "legacy_cols"; d.mkdir()          # spots 在列，坐标在列名
    pd.DataFrame(M, index=GENES, columns=spots).to_csv(d / "counts.tsv", sep="\t")
    (d / "HE_image.jpg").write_bytes(b"\xff\xd8\xff" + b"0" * 2000)

    d = root / "legacy_rows"; d.mkdir()          # spots 在行，坐标在行名
    pd.DataFrame(M.T, index=spots, columns=GENES).to_csv(d / "expression.tsv", sep="\t")

    d = root / "matrix_coords"; d.mkdir()        # 独立坐标文件
    pd.DataFrame(M.T, index=bc, columns=GENES).to_csv(d / "counts.csv")
    pd.DataFrame({"barcode": bc, "x": xy[:, 0], "y": xy[:, 1]}).to_csv(
        d / "spot_coordinates.csv", index=False)
    (d / "tissue_hires_image.png").write_bytes(b"\x89PNG" + b"0" * 3000)

    d = root / "no_coords"; d.mkdir()            # 缺坐标，应判为不可用
    pd.DataFrame(M.T, index=bc, columns=GENES).to_csv(d / "counts.csv")

    # 第一代 ST 的真实形态（GSE144239）：整块阵列的计数 + 只含组织上点的选点文件
    # + .gz 包着的 H&E。三个文件名都不含 count/matrix 之类的词。
    d = root / "legacy_sel"; d.mkdir()
    pd.DataFrame(M.T, index=spots, columns=GENES).to_csv(d / "P9_ST_rep1_stdata.tsv", sep="\t")
    on_tissue = spots[:N_ON_TISSUE]
    sxy = np.array([[float(t.split("x")[0]), float(t.split("x")[1])] for t in on_tissue])
    pd.DataFrame({"x": sxy[:, 0].astype(int), "y": sxy[:, 1].astype(int),
                  "new_x": sxy[:, 0], "new_y": sxy[:, 1],
                  "pixel_x": sxy[:, 0] * 194.0, "pixel_y": sxy[:, 1] * 194.0}).to_csv(
        d / "spot_data-selection-P9_ST_rep1.tsv", sep="\t", index=False)
    with gzip.open(d / "P9_ST_rep1.jpg.gz", "wb") as fh:
        fh.write(b"\xff\xd8\xff" + b"0" * 2000)
    return root


def test_parse_spot_coords():
    """坐标能从 "10x20" 形式的 spot 名中解析出来，非坐标名返回 None。"""
    xy = parse_spot_coords(["10x20", "11.5x20", "12_21"])
    assert xy is not None and xy.shape == (3, 2)
    assert np.allclose(xy[1], [11.5, 20.0])
    print(f"  解析 3 个 spot 名 -> {xy.tolist()}")
    assert parse_spot_coords(["AAACCTGAGC-1", "AAACGGGTCA-1"]) is None
    print("  10x barcode 正确返回 None")


def test_inspect_detects_formats():
    """四种格式都能被识别，缺坐标的被正确判为不可用。"""
    root = _mock_root()
    try:
        exp = {"legacy_cols": True, "legacy_rows": True,
               "matrix_coords": True, "legacy_sel": True, "no_coords": False}
        for name, should_have_coords in exp.items():
            info = inspect_raw(root / name)
            print(f"  {name:<16} format={info['format']:<14} "
                  f"coords={info['has_coords']} src={info['coord_source']}")
            assert info["has_counts"], f"{name} 未识别出计数矩阵"
            assert info["has_coords"] == should_have_coords, (
                f"{name} 的坐标识别错误：期望 {should_have_coords}，"
                f"实际 {info['has_coords']}"
            )
        assert inspect_raw(root / "legacy_cols")["has_image"]
        assert not inspect_raw(root / "legacy_rows")["has_image"]
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_orientation_always_correct():
    """三种可用格式都必须加载成 (400 spot, 16 基因)，方向不能弄反。

    这是本文件最重要的一条：方向弄反后下游每一步都能跑通，只是全错。
    """
    root = _mock_root()
    try:
        for name in ("legacy_cols", "legacy_rows", "matrix_coords"):
            info = inspect_raw(root / name)
            vals, obs, var, xy = read_table_matrix(info)
            print(f"  {name:<16} {vals.shape}  依据：{info['notes'][-1].split('：')[-1]}")
            assert vals.shape == (N_SPOT, len(GENES)), (
                f"{name} 方向错误：得到 {vals.shape}，应为 ({N_SPOT}, {len(GENES)})"
            )
            assert len(obs) == N_SPOT and len(var) == len(GENES)
            assert xy.shape == (N_SPOT, 2)
            assert set(var) == set(GENES), "基因名没落在 var 上，说明方向反了"
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_legacy_selection_filter():
    """第一代 ST：必须按选点文件裁掉组织外的阵列点，且 .jpg.gz 也算图像。

    为什么单独测这一条
    ------------------
    GSE144239 的 stdata 是**整块阵列**（P2_rep1 有 1933 个点），
    只有选点文件里那 665 个真的落在组织上。漏掉这一步不会报任何错，
    只会把空间图建在一大片背景上——源汇分位数、连通性、最小割全部失去意义。
    """
    root = _mock_root()
    try:
        info = inspect_raw(root / "legacy_sel")
        print(f"  format={info['format']} src={info['coord_source']} "
              f"coords={Path(info['files'].get('coords', '')).name}")
        assert info["format"] == "legacy_st"
        assert info["coord_source"] == "spot_names", "坐标应仍取自 spot 名"
        assert info["files"].get("coords"), "选点文件没有被记录，过滤不会生效"
        assert info["has_image"], ".jpg.gz 应算作组织学图像（准入标准 C3）"

        vals, obs, var, xy = read_table_matrix(info)
        print(f"  裁剪后 {vals.shape}（整块阵列是 {N_SPOT} 个点）")
        assert vals.shape == (N_ON_TISSUE, len(GENES)), (
            f"选点过滤没生效：得到 {vals.shape[0]} 个 spot，应为 {N_ON_TISSUE}")
        assert len(obs) == N_ON_TISSUE and xy.shape == (N_ON_TISSUE, 2)
        assert set(var) == set(GENES), "方向判断被选点文件带偏了"

        # 坐标约定对不上时必须报错，不能安静地留下"恰好对上的那部分"
        d = root / "legacy_sel"
        sel = pd.read_csv(d / "spot_data-selection-P9_ST_rep1.tsv", sep="\t")
        sel["x"] = sel["x"] + 100          # 整体平移，模拟约定不一致
        sel.to_csv(d / "spot_data-selection-P9_ST_rep1.tsv", sep="\t", index=False)
        try:
            read_table_matrix(inspect_raw(d))
        except ValueError as e:
            print(f"  坐标约定不一致时正确报错：{str(e)[:34]}...")
        else:
            raise AssertionError("选点与矩阵完全对不上时应抛 ValueError，却成功返回了")
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_missing_coords_raises():
    """缺坐标时必须显式报错，不能静默返回。"""
    root = _mock_root()
    try:
        info = inspect_raw(root / "no_coords")
        try:
            read_table_matrix(info)
        except ValueError as e:
            print(f"  正确报错：{str(e)[:44]}...")
            return
        raise AssertionError("缺坐标时应抛 ValueError，却成功返回了")
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        print(f"\n[RUN ] {t.__name__}")
        print(f"       {(t.__doc__ or '').strip().splitlines()[0]}")
        try:
            t(); print(f"[PASS] {t.__name__}")
        except AssertionError as e:
            failed += 1; print(f"[FAIL] {t.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            failed += 1; print(f"[ERROR] {t.__name__}: {type(e).__name__}: {e}")
    print("\n" + "=" * 60)
    print(f"{len(tests) - failed} / {len(tests)} passed")
    sys.exit(1 if failed else 0)
