#!/usr/bin/env python
"""
run_batch.py —— 批量跑全流程
==============================

输入：data/ledger.csv 中登记的切片
输出：各切片的中间产物 + results/batch_report.csv
上游模块：run_00b_ingest.py
下游模块：run_07_screen.py / run_06_validate.py

按台账逐张切片依次跑 run_01 → run_02 → run_03 →（可选）run_04 → run_05。
单张失败不中断整批，最后汇总成一张报告表，让你一眼看出哪几张卡在哪一步。

用法
----
    python scripts/run_batch.py --through 03          # 只跑到建图（筛查前够用）
    python scripts/run_batch.py --through 05          # 跑完全套
    python scripts/run_batch.py --slides MEL01 MEL02  # 只跑指定几张
    python scripts/run_batch.py --through 03 --dry-run
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402

from sparta.io_ import Paths, load_config  # noqa: E402

SCRIPTS = Path(__file__).resolve().parent
STEPS = ["01", "02", "03", "04", "05"]
NAMES = {"01": "run_01_qc.py", "02": "run_02_score.py", "03": "run_03_graph.py",
         "04": "run_04_barrier.py", "05": "run_05_counterfactual.py"}


def build_cmd(step, row, P):
    sid = row["slide_id"]
    py = sys.executable
    if step == "01":
        cmd = [py, str(SCRIPTS / NAMES[step]), "--slide", sid,
               "--input", str(P.interim / f"{sid}.raw.h5ad"),
               "--platform", "legacy_st" if row.get("platform") == "legacy_st" else "visium"]
        if row.get("has_image"):
            cmd.append("--has-he")
        if str(row.get("treatment", "unknown")) != "unknown":
            cmd.append("--treatment-known")
        return cmd
    if step == "02":
        ct = str(row.get("cancer_type", "melanoma"))
        return [py, str(SCRIPTS / NAMES[step]), "--slide", sid,
                "--tumor-type", ct if ct in ("melanoma", "cscc", "bcc") else "melanoma",
                "--plot"]
    return [py, str(SCRIPTS / NAMES[step]), "--slide", sid]


def main():
    ap = argparse.ArgumentParser(description="批量跑全流程",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--through", default="03", choices=STEPS,
                    help="跑到第几步为止。筛查阶段跑到 03 就够")
    ap.add_argument("--slides", nargs="+", default=None, help="只跑这几张")
    ap.add_argument("--config", default=None)
    ap.add_argument("--dry-run", action="store_true", help="只打印命令不执行")
    ap.add_argument("--timeout", type=int, default=3600, help="单步超时（秒）")
    ap.add_argument("--include-rejected", action="store_true",
                    help="连台账里 status 不是 ingested 的切片也一起跑（默认跳过）")
    args = ap.parse_args()

    cfg = load_config(args.config)
    P = Paths(cfg)
    led_p = P.root / "data" / "ledger.csv"
    if not led_p.exists():
        sys.exit(f"找不到数据台账 {led_p}。请先用 run_00b_ingest.py 摄入数据。")
    led = pd.read_csv(led_p)
    if args.slides:
        led = led[led["slide_id"].isin(args.slides)]

    # 台账里 status != ingested 的是准入没过、已经明确剔除的切片。
    # 从前这里不过滤，被剔除的切片照样会被批量跑进下游产物里，
    # 而剔除决定只留在 admission.json 与台账里——两边就此对不上，且不报错。
    if not args.include_rejected and "status" in led.columns:
        bad = led[led["status"].astype(str) != "ingested"]
        if len(bad):
            print("已跳过台账中未准入的切片（要一起跑请加 --include-rejected）：")
            for _, r in bad.iterrows():
                print(f"  {r['slide_id']:<8} status={r['status']}  {str(r.get('notes', ''))[:60]}")
            print()
            led = led[led["status"].astype(str) == "ingested"]

    if led.empty:
        sys.exit("台账里没有匹配的切片。")

    steps = STEPS[:STEPS.index(args.through) + 1]
    print(f"批量运行：{len(led)} 张切片 × {len(steps)} 步（{'、'.join(steps)}）\n")

    report = []
    for _, row in led.iterrows():
        sid = row["slide_id"]
        rec = {"slide_id": sid, "cancer_type": row.get("cancer_type")}
        for step in steps:
            cmd = build_cmd(step, row, P)
            if args.dry_run:
                print("  " + " ".join(cmd))
                rec[step] = "dry"
                continue
            print(f"[{sid}] 第 {step} 步 ...", end=" ", flush=True)
            try:
                r = subprocess.run(cmd, capture_output=True, text=True,
                                   timeout=args.timeout)
                if r.returncode == 0:
                    print("OK")
                    rec[step] = "ok"
                else:
                    tail = (r.stderr or r.stdout).strip().splitlines()
                    msg = tail[-1][:110] if tail else "未知错误"
                    print(f"失败 — {msg}")
                    rec[step] = f"fail: {msg}"
                    break               # 上游失败就不跑下游了
            except subprocess.TimeoutExpired:
                print("超时")
                rec[step] = "timeout"
                break
        report.append(rec)

    if args.dry_run:
        return

    df = pd.DataFrame(report)
    out = P.results / "batch_report.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    print("\n" + "=" * 60)
    done = sum(1 for r in report if r.get(steps[-1]) == "ok")
    print(f"完成全部 {len(steps)} 步的切片：{done} / {len(report)}")
    fails = [r for r in report if any(str(v).startswith(("fail", "timeout"))
                                      for v in r.values())]
    if fails:
        print(f"\n有问题的 {len(fails)} 张：")
        for r in fails:
            bad = [f"{k}={v}" for k, v in r.items()
                   if str(v).startswith(("fail", "timeout"))]
            print(f"  {r['slide_id']}: {'; '.join(bad)}")
    print(f"\n完整报告：{out}")
    if done >= 3:
        ok_ids = [r["slide_id"] for r in report if r.get(steps[-1]) == "ok"]
        print(f"\n下一步（重要）——跑早期决策点筛查：")
        print(f"  python scripts/run_07_screen.py --slides {' '.join(ok_ids[:8])}")


if __name__ == "__main__":
    main()
