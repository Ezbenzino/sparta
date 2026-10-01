#!/usr/bin/env python
"""
run_ingest_legacy_cscc.py —— 批量摄入 GSE144239 的第一代 ST 切片（P2/P5/P9/P10）
==================================================================================

输入：data/raw/ 下已下载好的 12 个样本目录（本脚本**不下载任何数据**）
输出：data/interim/{sid}.{raw,qc,scored,graph,barrier}.* + 台账 + 准入记录
上游模块：无（建库入口）
下游模块：run_05_counterfactual / run_12 / run_14 / run_16

为什么要有这个脚本
------------------
把 12 张切片手打一遍就是 60 条命令，每条都要带对 --patient / --replicate /
--cohort / --config 四个参数。打错一个 --patient，下游的患者层统计就全错，
而且**不会报错**——只会安静地把两个患者算成一个。
所以这里把 GSM → 患者 / 重复 / 切片 ID 的映射写死在代码里，人不再手填。

映射依据（2026-08-27 从 GEO 元数据核对，不是推测）
--------------------------------------------------
GSE144239 的 Visium 只有 4 个样本（P4/P6，已用作 CSCC01–04）。
另有 12 个第一代 ST 样本，来自**另外 4 位患者**，每人 3 张重复。
纳入它们一次解决两件事：
  1. cSCC 患者数 2 -> 6；
  2. 第一代 ST 是第三个平台世代——若其结果与 Visium 的 cSCC 一致，
     就拆掉了 R3b 里"队列差异可能是平台世代造成的"这条混杂，
     而那是瘤种/分期/平台三个共变量里唯一能拆开的一个。

第一代 ST 最容易翻车的两处（本脚本会逐张检查并汇总）
----------------------------------------------------
  · **方向**：stdata.tsv 是 spot × gene，行名形如 "10x20"。若摄入后
    n_spots > n_genes，多半是判反了，要加 --transpose 重来。
  · **连通性**：点间距约 200 μm，radius_um 必须放大到 300。
    建图后连通分量数 > 1 说明还不够大，往 350–400 调。

用法
----
    python scripts/run_ingest_legacy_cscc.py --dry-run       # 只看要做什么
    python scripts/run_ingest_legacy_cscc.py --only P2       # 先跑通一个患者
    python scripts/run_ingest_legacy_cscc.py                 # 全部 12 张
    python scripts/run_ingest_legacy_cscc.py --through 03    # 只跑到建图
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors="replace")
    except Exception:  # noqa: BLE001
        pass

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable
CFG = "configs/cscc_legacy_st.yaml"
COHORT = "cscc_legacy_gse144239"

# (slide_id, GSM, patient, replicate)  —— 顺序即 GEO 上的样本顺序
SAMPLES = [
    ("CSCC05", "GSM4284316", "CSCC_P2", "rep1"),
    ("CSCC06", "GSM4284317", "CSCC_P2", "rep2"),
    ("CSCC07", "GSM4284318", "CSCC_P2", "rep3"),
    ("CSCC08", "GSM4284319", "CSCC_P5", "rep1"),
    ("CSCC09", "GSM4284320", "CSCC_P5", "rep2"),
    ("CSCC10", "GSM4284321", "CSCC_P5", "rep3"),
    ("CSCC11", "GSM4284322", "CSCC_P9", "rep1"),
    ("CSCC12", "GSM4284323", "CSCC_P9", "rep2"),
    ("CSCC13", "GSM4284324", "CSCC_P9", "rep3"),
    ("CSCC14", "GSM4284325", "CSCC_P10", "rep1"),
    ("CSCC15", "GSM4284326", "CSCC_P10", "rep2"),
    ("CSCC16", "GSM4284327", "CSCC_P10", "rep3"),
]

STEPS = ["00b", "01", "02", "03", "04"]


def _find_raw(sid: str, gsm: str) -> Path | None:
    """在 data/raw 下找这个样本的目录。

    容忍三种命名：直接叫切片 ID、叫 GSM 号、或目录名里包含 GSM 号。
    找不到就返回 None，由调用方报告缺哪些——不要猜。
    """
    raw = ROOT / "data" / "raw"
    if not raw.exists():
        return None
    for cand in (raw / sid, raw / gsm):
        if cand.exists():
            return cand
    hits = [d for d in raw.iterdir() if d.is_dir() and gsm in d.name]
    return hits[0] if len(hits) == 1 else None


def run(cmd) -> int:
    print(f"    $ {' '.join(str(c) for c in cmd[1:])}", flush=True)
    return subprocess.call([str(c) for c in cmd], cwd=str(ROOT))


def main():
    ap = argparse.ArgumentParser(
        description="批量摄入 GSE144239 第一代 ST 切片",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--only", nargs="+", default=None,
                    help="只处理这些患者，如 --only P2（建议先跑通一个）")
    ap.add_argument("--through", default="04", choices=STEPS,
                    help="跑到哪一步为止")
    ap.add_argument("--from", dest="from_step", default=STEPS[0], choices=STEPS,
                    help="从哪一步开始。前面几步已经跑过时用，比如只补 run_04："
                         "--only P2 --from 04")
    ap.add_argument("--dry-run", action="store_true", help="只打印计划，不执行")
    ap.add_argument("--transpose", action="store_true",
                    help="强制转置（摄入后发现 n_spots > n_genes 时用）")
    args = ap.parse_args()

    todo = [s for s in SAMPLES
            if not args.only or any(p in s[2] for p in args.only)]
    if not todo:
        sys.exit(f"--only {args.only} 没匹配到任何样本。可选：P2 P5 P9 P10")

    stop = STEPS.index(args.through)
    start = STEPS.index(args.from_step)
    if start > stop:
        sys.exit(f"--from {args.from_step} 排在 --through {args.through} 后面，跑不了。")

    print("=" * 78)
    print(f"GSE144239 第一代 ST 摄入　{len(todo)} 张切片／"
          f"{len({s[2] for s in todo})} 位患者　"
          f"run_{args.from_step} → run_{args.through}")
    print(f"配置 {CFG}　准入队列 {COHORT}")
    print("=" * 78)

    missing = []
    for sid, gsm, pat, rep in todo:
        src = _find_raw(sid, gsm)
        flag = "" if src else "  ← 找不到原始数据"
        print(f"  {sid:<8} {gsm:<12} {pat:<9} {rep:<6} "
              f"{src.name if src else '—':<28}{flag}")
        if not src:
            missing.append((sid, gsm))

    if missing:
        print()
        print(f"  有 {len(missing)} 个样本在 data/raw 下找不到。请先从 GEO 下载，")
        print("  每个样本解压到 data\\raw\\<切片ID>\\ 或 data\\raw\\<GSM号>\\，")
        print("  目录里应当有这三个文件（解压后）：")
        print("    *_stdata.tsv.gz                 spot × gene 计数，行名形如 10x20")
        print("    *_spot_data-selection-*.tsv.gz  点选与像素坐标")
        print("    *.jpg.gz                        组织图像（准入 C3 需要）")
        print("  样本页：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<GSM号>")
        if not args.dry_run:
            sys.exit("\n  缺数据，先不跑。补齐后重来。")

    # ---- --from 的前置产物检查 ----
    # 直接从中途某一步开始时，前面几步的产物必须已经在 data/interim 下。
    # 不检查的话，run_04 会以一句光秃秃的 FileNotFoundError 死掉，
    # 12 张一起跑就是 12 条一模一样的报错，看不出到底是缺文件还是代码有 bug。
    NEEDS = {
        "00b": [],
        "01":  ["{sid}.raw.h5ad"],
        "02":  ["{sid}.qc.h5ad"],
        "03":  ["{sid}.scored.h5ad"],
        "04":  ["{sid}.scored.h5ad", "{sid}.graph.npz"],
    }
    if start > 0:
        lacking = []
        for sid, *_ in todo:
            gone = [f.format(sid=sid) for f in NEEDS[args.from_step]
                    if not (ROOT / "data" / "interim" / f.format(sid=sid)).exists()]
            if gone:
                lacking.append((sid, gone))
        if lacking:
            print()
            print(f"  --from {args.from_step} 要求前面几步的产物已经存在，但这些切片缺文件：")
            for sid, gone in lacking:
                print(f"    {sid:<8} 缺 " + "、".join(gone))
            prev = STEPS[start - 1]
            print()
            print("  也就是说前面的步骤没跑成（或者根本没跑）。去掉 --from 从头跑：")
            only = " ".join(sorted({t[2].replace("CSCC_", "") for t in todo}))
            print(f"    python scripts/run_ingest_legacy_cscc.py --only {only} --through {prev}")
            if not args.dry_run:
                sys.exit("\n  前置产物不全，先不跑。")

    if args.dry_run:
        print("\n  --dry-run：以上是计划，没有执行任何命令。")
        return

    ok, bad = [], []
    for i, (sid, gsm, pat, rep) in enumerate(todo, 1):
        src = _find_raw(sid, gsm)
        print(f"\n{'=' * 78}\n[{i}/{len(todo)}] {sid}　{pat} {rep}　{gsm}\n{'=' * 78}", flush=True)

        cmds = [
            [PY, "scripts/run_00b_ingest.py", "--input", src, "--slide", sid,
             "--patient", pat, "--replicate", rep,
             "--cancer", "cscc", "--platform", "legacy_st",
             "--source", "GEO GSE144239", "--accession", gsm,
             "--treatment", "naive", "--site", "primary", "--config", CFG]
            + (["--transpose"] if args.transpose else []),
            [PY, "scripts/run_01_qc.py", "--slide", sid,
             "--input", f"data/interim/{sid}.raw.h5ad",
             "--platform", "legacy_st", "--has-he", "--treatment-known",
             "--cohort", COHORT, "--config", CFG],
            [PY, "scripts/run_02_score.py", "--slide", sid, "--tumor-type", "cscc",
             "--config", CFG],
            [PY, "scripts/run_03_graph.py", "--slide", sid, "--config", CFG],
            [PY, "scripts/run_04_barrier.py", "--slide", sid, "--config", CFG],
        ]
        failed_at = None
        for j, cmd in enumerate(cmds[start:stop + 1], start=start):
            if run(cmd) != 0:
                failed_at = STEPS[j]
                break
        if failed_at:
            print(f"  [!] {sid} 在 run_{failed_at} 失败，跳过它继续下一张。")
            bad.append((sid, failed_at))
        else:
            ok.append(sid)

    print(f"\n{'=' * 78}\n完成　成功 {len(ok)}/{len(todo)}\n{'=' * 78}")
    if bad:
        for sid, st in bad:
            print(f"  失败：{sid}（run_{st}）")
    if ok:
        print("\n跑完必看两件事：")
        print("  1) 台账里的 n_spots 与 n_genes —— 若 spot 数大于基因数，方向判反了，")
        print("     加 --transpose 重跑这几张（下游一切照常，只是全错）。")
        print("  2) run_03 打印的连通分量数 —— 大于 1 说明 radius_um=300 还不够，")
        print(f"     把 {CFG} 里的 graph.radius_um 调到 350–400 再来一次。")
        print("\n两件都没问题之后：")
        print("  .\\.venv\\Scripts\\python.exe scripts\\run_05_counterfactual.py --slide <每一张> --config " + CFG)
        print("  .\\.venv\\Scripts\\python.exe scripts\\run_06_validate.py --slides <全部> --dim decoupling")
        print("  .\\.venv\\Scripts\\python.exe scripts\\run_14_shared_ecm_check.py --slides <全部>   ← 决策")
        print("  .\\.venv\\Scripts\\python.exe scripts\\run_16_figures.py                            ← 更新图")


if __name__ == "__main__":
    main()
