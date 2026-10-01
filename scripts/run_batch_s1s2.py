#!/usr/bin/env python
"""
run_batch_s1s2.py —— 唯一的长任务：8 张切片重跑 S1 + S2（3–5 小时）
=====================================================================

输入：{sid}.scored.h5ad + {sid}.graph.npz
输出：results/counterfactual/{sid}.json（**会覆盖**，请先跑 run_fastpath.py 的备份步骤）
上游模块：run_03_graph.py
下游模块：run_12_paper_stats.py

为什么要重跑
------------
① S1 的置换次数两队列不一致（MEL 500 次、CSCC 只有 150 次，p 下限不同、不可比）；
② S2 的 k 是绝对值，在 614 节点的割集上只动了 0.5%，效应量被稀释到测不出来。
config 已改为 counterfactual.s2.k_frac，run_05 会自动使用。

单张失败不中断，最后给汇总。建议睡前挂上。

用法
----
    python scripts/run_batch_s1s2.py
    python scripts/run_batch_s1s2.py --slides MEL01 MEL02      # 只跑两张
    python scripts/run_batch_s1s2.py --n-perm 200 --n-rand 80  # 快速试跑，结果不能进论文
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors="replace")
    except Exception:  # noqa: BLE001
        pass

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable
SLIDES = ["MEL01", "MEL02", "MEL03", "MEL04", "CSCC01", "CSCC02", "CSCC03", "CSCC04"]


def main():
    ap = argparse.ArgumentParser(description="批量重跑 S1 + S2",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slides", nargs="+", default=SLIDES)
    ap.add_argument("--n-perm", type=int, default=500,
                    help="S1 置换次数。两个队列必须用同一个值，否则 p 不可比")
    ap.add_argument("--n-rand", type=int, default=None, help="S2 随机对照次数，默认取 config")
    args = ap.parse_args()

    print(f"重跑 S1 + S2　{len(args.slides)} 张切片　n_perm={args.n_perm}")
    print(f"仓库 {ROOT}\n解释器 {PY}")
    print("提示：这会覆盖 results/counterfactual/*.json。"
          "没备份的话先 Ctrl+C，跑 run_fastpath.py 会自动备份。\n")

    t0 = time.time()
    ok, bad = [], []
    for i, sid in enumerate(args.slides, 1):
        print("=" * 78)
        print(f"[{i}/{len(args.slides)}] {sid}　（已用时 {(time.time()-t0)/60:.0f} 分钟）")
        print("=" * 78, flush=True)
        cmd = [PY, "scripts/run_05_counterfactual.py", "--slide", sid,
               "--n-perm", str(args.n_perm)]
        if args.n_rand:
            cmd += ["--n-rand", str(args.n_rand)]
        code = subprocess.call(cmd, cwd=str(ROOT))
        (ok if code == 0 else bad).append(sid)
        if code:
            print(f"  [!] {sid} 返回退出码 {code}，继续下一张。")

    print()
    print("=" * 78)
    print(f"完成　成功 {len(ok)}/{len(args.slides)}　总用时 {(time.time()-t0)/60:.0f} 分钟")
    if bad:
        print(f"失败：{', '.join(bad)}")
    print("=" * 78)
    print("\n下一步：python scripts/run_12_paper_stats.py")
    print("  看两件事：S1 表的 'p下限' 是否全部 0.0020；"
          "S2 表的 '占割集' 是否变成 5%/10%/20%/30%。")


if __name__ == "__main__":
    main()
