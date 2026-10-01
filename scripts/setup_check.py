#!/usr/bin/env python
"""
setup_check.py —— 环境自检（第 0 步，装完就跑这个）
=====================================================

输入：无      输出：终端报告
上游模块：无  下游模块：run_00_demo.py

逐项检查依赖是否可用，缺什么就告诉你装什么。分三档：
  [必需]   缺了核心算子都跑不了
  [真实数据] 缺了只能跑合成演示，处理不了真实切片
  [验证]   缺了做不了 M6 的形态学与统计验证

用法
----
    python scripts/setup_check.py
"""
import importlib
import platform
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

CHECKS = [
    ("必需", "numpy", "pip install numpy", "数组运算"),
    ("必需", "scipy", "pip install scipy", "稀疏矩阵、图算法、线性求解"),
    ("必需", "networkx", "pip install networkx", "最大流/最小割"),
    ("必需", "pandas", "pip install pandas", "表格处理"),
    ("必需", "yaml", "pip install pyyaml", "读配置文件"),
    ("必需", "matplotlib", "pip install matplotlib", "出图"),
    ("真实数据", "anndata", "pip install anndata", "空间数据容器"),
    ("真实数据", "scanpy", "pip install scanpy", "质控、签名打分"),
    ("真实数据", "squidpy", "pip install squidpy", "空间邻接图（可选，有替代实现）"),
    ("验证", "skimage", "pip install scikit-image", "H&E 形态学特征"),
    ("验证", "sklearn", "pip install scikit-learn", "AUC 等指标"),
    ("验证", "statsmodels", "pip install statsmodels", "多变量回归、似然比检验"),
    ("验证", "lifelines", "pip install lifelines", "生存分析"),
    ("开发", "pytest", "pip install pytest", "跑测试（不装也能用 python 直接跑）"),
]


def main():
    print("=" * 68)
    print(f"SPARTA 环境自检   Python {platform.python_version()}  {platform.system()}")
    print("=" * 68)

    if sys.version_info < (3, 9):
        print(f"\n✗ Python 版本过低（{platform.python_version()}），需要 >= 3.9")
        return 1

    missing = {"必需": [], "真实数据": [], "验证": [], "开发": []}
    cur = None
    for tier, mod, fix, why in CHECKS:
        if tier != cur:
            print(f"\n[{tier}]")
            cur = tier
        try:
            m = importlib.import_module(mod)
            v = getattr(m, "__version__", "?")
            print(f"  ✓ {mod:<14} {v:<12} {why}")
        except ImportError:
            print(f"  ✗ {mod:<14} {'缺失':<12} {why}")
            missing[tier].append((mod, fix))

    # 中文字体
    print("\n[出图]")
    try:
        from sparta.viz import setup_cjk_font
        f = setup_cjk_font()
        print(f"  {'✓' if f else '✗'} 中文字体      {f or '未找到（图中中文会显示为方块）'}")
    except Exception as e:  # noqa: BLE001
        print(f"  ? 中文字体      检查失败：{e}")

    print("\n" + "=" * 68)
    if missing["必需"]:
        print("✗ 核心依赖缺失，先装这些：")
        for mod, fix in missing["必需"]:
            print(f"    {fix}")
        return 1

    print("✓ 核心依赖齐全。现在可以跑：")
    print("    python scripts/run_00_demo.py --fast")
    print("    python tests/test_barrier.py")

    if missing["真实数据"]:
        print("\n处理真实切片还需要：")
        for mod, fix in missing["真实数据"]:
            print(f"    {fix}")
    if missing["验证"]:
        print("\n做 M6 验证还需要：")
        for mod, fix in missing["验证"]:
            print(f"    {fix}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
