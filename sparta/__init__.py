"""
SPARTA —— SPAtial Resistance to Therapeutic Agents
====================================================

皮肤恶性肿瘤免疫治疗抵抗的双重空间屏障分析框架。

核心主张
--------
免疫治疗抵抗的空间维度不止一层：
  · "T 细胞进不去"（细胞迁移屏障，B_cell）
  · "抗 PD-1 抗体自身进不去"（大分子传质屏障，B_mAb）
这两者物理机制不同，在真实肿瘤中可以解离，而空间组学界系统性忽视了后者。

模块地图
--------
  io_          数据契约的统一读写与路径管理
  signatures   M2 签名打分（替代反卷积）
  graph        M3 空间图构建与源汇定义
  barrier      M4 三分量屏障算子  <- 核心，改动前请先看懂 tests/test_barrier.py
  counterfactual M5 反事实实验（S1 空间重排 / S2 环带断裂 / S3 尺寸扫描）
  validate     M6 四维验证矩阵
  viz          出图（屏障线叠加图是本课题最有说服力的图件）
  synthetic    合成图生成器，用于单元测试与教学

快速上手
--------
    python scripts/run_00_demo.py       # 在合成数据上跑通全流程，看看每步长什么样
    python tests/test_barrier.py        # 验证核心算子（不需要 pytest 也能跑）

注意：barrier / graph / counterfactual / synthetic 只依赖 numpy·scipy·networkx，
      不需要 scanpy。scanpy 只在处理真实 AnnData 时用到。
"""

__version__ = "2.2.3"


# --------------------------------------------------------------------------
# Windows 控制台编码兜底
# --------------------------------------------------------------------------
# 为什么需要这几行（2026-08-27 踩到）
# ----------------------------------
# Python 3.6+ 在 Windows **控制台**上走 UTF-16 接口，打印 ⚠ ✓ ✗ 都没问题；
# 可一旦 stdout 被**重定向到文件或管道**（`... *> run.log`、`| tee`、CI 收日志、
# subprocess 抓输出），就改用 locale 编码（简中系统是 cp936）且 errors='strict'，
# 于是第一个 GBK 编不出的字符直接抛 UnicodeEncodeError 把脚本打死。
# 全项目有 20 多个文件带这类字符，日志里连报错都看不全——最坏的一种失败方式：
# 在控制台好好的，一收日志就挂，而且挂在跟数据毫无关系的地方。
#
# errors="replace" 只改「编不出时怎么办」，不改编码本身，
# 因此对已经能正常输出的环境完全没有影响。
import sys as _sys

for _stream in (_sys.stdout, _sys.stderr):
    try:
        _stream.reconfigure(errors="replace")
    except Exception:  # noqa: BLE001  —— 比如 stdout 被换成了不支持 reconfigure 的对象
        pass
del _stream, _sys

from . import barrier, counterfactual, graph, synthetic  # noqa: F401,E402

__all__ = ["barrier", "counterfactual", "graph", "synthetic", "__version__"]
