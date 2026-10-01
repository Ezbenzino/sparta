#!/usr/bin/env python
"""
run_fastpath.py —— 审查整改快路径（约 15–25 分钟，一路跑到决策岔口）
=======================================================================

输入：已有的 data/interim/* 与 results/*
输出：见各步骤；本脚本自己只写 results/_archive_<日期>/ 的备份
上游模块：run_03_graph.py（各步骤的实际上游见其各自文件头）
下游模块：无

为什么是 Python 而不是 .ps1
---------------------------
Windows PowerShell 5.1 默认按 ANSI/GBK 读 .ps1，无 BOM 的 UTF-8 脚本里
中文会被拆错字节，进而把引号配对弄断，报一串莫名其妙的语法错误。
换成 Python 驱动就完全绕开了这件事：源码编码由 Python 规定（UTF-8），
用哪个解释器跑就用哪个解释器调子进程（sys.executable），不依赖任何 shell。

它做什么
--------
把所有**不耗时**的整改按"先便宜后昂贵"跑完，最后落在决策节点上。
唯一的长任务（S1/S2 重跑，3–5 小时）不在这里，跑完本脚本再挂。

用法
----
    python scripts/run_fastpath.py
    python scripts/run_fastpath.py --from 5      # 从第 5 步接着跑
    python scripts/run_fastpath.py --skip-tests  # 跳过第 0 步
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# 控制台编码兜底：GBK 控制台遇到不可映射的字符（如 ⚠）会抛 UnicodeEncodeError，
# 换成 replace 之后最多显示成 '?'，不会把整个脚本弄崩。
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors="replace")
    except Exception:  # noqa: BLE001 —— 老版本 Python 没有 reconfigure
        pass

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable
SLIDES = ["MEL01", "MEL02", "MEL03", "MEL04", "CSCC01", "CSCC02", "CSCC03", "CSCC04"]
TESTS = ["test_barrier", "test_counterfactual", "test_loaders", "test_statistics"]


def banner(n, title, why=None):
    print()
    print("=" * 78)
    print(f"[{n}] {title}")
    if why:
        for line in why.splitlines():
            print(f"    {line}")
    print("=" * 78, flush=True)


def run(cmd):
    """跑一个子进程，输出直接流到当前控制台（不捕获，编码问题最少）。"""
    print(f"  $ {' '.join(str(c) for c in cmd)}", flush=True)
    return subprocess.call([str(c) for c in cmd], cwd=str(ROOT))


def _capture(cmd):
    """跑一个子进程并收集输出。

    强制子进程用 UTF-8 写管道（PYTHONIOENCODING），这边也按 UTF-8 解，
    避免在 GBK 控制台上管道两端编码不一致导致乱码或抛异常。
    """
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    r = subprocess.run([str(c) for c in cmd], cwd=str(ROOT), env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return r.returncode, r.stdout.decode("utf-8", errors="replace")


def step_tests():
    banner(0, "确认整改没把东西改坏（约 2 分钟）")
    bad = []
    for t in TESTS:
        # 测试输出很长（每个断言都打印诊断），这里只收不显；失败了才把尾巴摊开
        code, out = _capture([PY, ROOT / "tests" / f"{t}.py"])
        tail = [ln for ln in out.strip().splitlines() if ln.strip()]
        summary = tail[-1].strip() if tail else "(无输出)"
        print(f"  {t:<22} {'OK  ' if code == 0 else '失败'}  {summary}")
        if code != 0:
            bad.append(t)
            print("  " + "-" * 60)
            for ln in tail[-25:]:
                print("  | " + ln)
            print("  " + "-" * 60)
    if bad:
        print(f"\n  {', '.join(bad)} 没过。先停下来，把上面的失败输出发给我，不要继续。")
        return False
    print("  -> 18 个测试全绿。")
    return True


def step_backup():
    banner(1, "备份将被覆盖的产物（10 秒）")
    arch = ROOT / "results" / f"_archive_{datetime.now():%Y-%m-%d}"
    if arch.exists():
        print(f"  {arch.name} 已存在，跳过（不覆盖已有备份）。")
        return True
    arch.mkdir(parents=True, exist_ok=True)
    n = 0
    for sub in ("counterfactual", "validation"):
        d = ROOT / "results" / sub
        if not d.exists():
            continue
        for f in d.glob("*.json"):
            shutil.copy2(f, arch / f"{sub}__{f.name}")
            n += 1
    print(f"  -> 已备份 {n} 个 JSON 到 results/{arch.name}/")
    return True


STEPS = [
    (2, "重建准入台账（3–5 分钟）",
     "旧版 run_01 只在准入失败时写记录，而失败分支被 --force 跳过，\n"
     "8 张里只有 MEL01 留下了记录。Methods 里'全部通过准入'没有证据。",
     lambda: run([PY, "scripts/run_00c_admission_audit.py"])),

    (3, "重跑解耦（口径修正，3–5 分钟）",
     "旧版用 0.5*(ECM+CAF) 当逐点 B_cell 代理且没传 control，\n"
     "产出的负相关是伪结果。",
     lambda: run([PY, "scripts/run_06_validate.py", "--slides", *SLIDES,
                  "--dim", "decoupling"])),

    (4, "分队列跑一致性（1 分钟）",
     "跨队列合并没有意义；--tag 让两次运行分别落盘，不再互相覆盖。",
     lambda: (run([PY, "scripts/run_06_validate.py", "--slides", *SLIDES[:4],
                   "--dim", "consistency", "--tag", "MEL"])
              or run([PY, "scripts/run_06_validate.py", "--slides", *SLIDES[4:],
                      "--dim", "consistency", "--tag", "CSCC"]))),

    (5, "★ 决策岔口：切断共享 ECM 通道后两屏障还相关吗（3–5 分钟）",
     "这一步的结果决定投 Bioinformatics 还是冲 BIB。\n"
     "它只写一个新文件 results/validation/shared_ecm_check.json，不覆盖任何产物。",
     lambda: run([PY, "scripts/run_14_shared_ecm_check.py", "--slides", *SLIDES])),

    (6, "核验参数身份的论断（5–10 分钟）",
     "看'排阻边%'一列：若 8 张高度接近，就证实排阻比例由秩分布而非生物学决定。",
     lambda: run([PY, "scripts/run_11_review_diagnostics.py", "--slides", *SLIDES])),

    (7, "论文统计汇总（10 秒）",
     "FDR 按切片成族、报 p 值分辨率下限、效应量改报比值。",
     lambda: run([PY, "scripts/run_12_paper_stats.py"])),
]


def main():
    ap = argparse.ArgumentParser(
        description="审查整改快路径（跑到决策岔口）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--from", dest="start", type=int, default=0,
                    help="从第几步开始（0=测试, 1=备份, 2..7=各整改步骤）")
    ap.add_argument("--skip-tests", action="store_true")
    args = ap.parse_args()

    print(f"SPARTA 快路径　仓库 {ROOT}")
    print(f"解释器 {PY}")

    if args.start <= 0 and not args.skip_tests:
        if not step_tests():
            sys.exit(1)
    if args.start <= 1:
        step_backup()

    failed = []
    for n, title, why, fn in STEPS:
        if n < args.start:
            continue
        banner(n, title, why)
        code = fn()
        if code:
            print(f"\n  [!] 第 {n} 步返回非零退出码 {code}，继续往下跑，"
                  f"但请把这一段输出发给我。")
            failed.append(n)

    print()
    print("=" * 78)
    print("快路径跑完了。接下来两件事：")
    print("=" * 78)
    if failed:
        print(f"  [!] 第 {', '.join(map(str, failed))} 步有非零退出码，先看那几段输出。")
        print()
    print("  1) 把第 5 步的判定（TISSUE / MODEL / REVERSED）告诉我 —— 它决定投哪个刊。")
    print("     完整结果在 results/validation/shared_ecm_check.json")
    print()
    print("  2) 睡前挂上唯一的长任务（S1/S2 重跑，3–5 小时）：")
    print()
    print("     python scripts/run_batch_s1s2.py")
    print()
    print("     跑完第二天早上再跑一次 python scripts/run_12_paper_stats.py 看新数字。")
    print()


if __name__ == "__main__":
    main()
