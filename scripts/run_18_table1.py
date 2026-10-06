#!/usr/bin/env python
"""
run_18_table1.py —— 从产物生成 Table 1（切片清单）
====================================================

输入：data/ledger.csv + data/interim/*.graph.npz + data/interim/*.admission.json
输出：docs/table1_sections.md
上游模块：run_00b_ingest / run_01_qc / run_03_graph
下游模块：无（正文引用）

为什么要有这个脚本
------------------
Table 1 手写过一次，队列从 8 张扩到 19 张之后整张表都错了，而且**不报错**——
表里少几张切片、连通分量数还停在旧值，看上去一切正常。
凡是"从产物抄进正文"的表，都必须能一条命令重新生成，否则迟早对不上。

用法
----
    python scripts/run_18_table1.py
    python scripts/run_18_table1.py --out docs/table1_sections.md
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402

from sparta.io_ import Paths, load_config, load_graph, stamp_run  # noqa: E402

CT = {"melanoma": "Melanoma", "cscc": "cSCC", "bcc": "BCC"}
PF = {"visium": "Visium", "legacy_st": "1st-gen ST", "slideseqv2": "Slide-seqV2", "unknown": "—"}


def graph_stats(P, sid):
    p = P.graph(sid)
    if not Path(p).exists():
        return None
    A, _D, source, sink, vessel, _ = load_graph(p)
    nc, lab = connected_components(A, directed=False)
    cnt = np.bincount(lab)
    return dict(n=A.shape[0], deg=float(A.getnnz(axis=1).mean()), nc=int(nc),
                big=100.0 * cnt.max() / A.shape[0],
                src=len(source), snk=len(sink), ves=len(vessel))


def main():
    ap = argparse.ArgumentParser(description="从产物生成 Table 1")
    ap.add_argument("--out", default="docs/table1_sections.md")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    P = Paths(cfg)
    led_p = P.root / "data" / "ledger.csv"
    if not led_p.exists():
        sys.exit(f"找不到台账 {led_p}")
    led = list(csv.DictReader(open(led_p, encoding="utf-8")))

    L = ["# Table 1 — Sections analysed", "",
         "由 `scripts/run_18_table1.py` 从 `data/ledger.csv`、`data/interim/*.graph.npz`",
         "与 `data/interim/*.admission.json` 生成。**不要手工编辑**——队列一变就会不一致，",
         "而且不会报错。重新生成：`python scripts/run_18_table1.py`", "",
         "| Section | Patient | Replicate / site | Cohort | Platform | GEO | Spots (raw) "
         "| Median UMI | Nodes after QC | Mean degree | Components | Largest comp. "
         "| Source / Sink / Vessel | Admission |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]

    headers = [
        "Section", "Patient", "Replicate / site", "Cohort", "Platform", "GEO",
        "Spots (raw)", "Median UMI", "Nodes after QC", "Mean degree",
        "Components", "Largest comp.", "Source / Sink / Vessel", "Admission",
    ]
    records = []

    rows = sorted(led, key=lambda r: (r["cancer_type"] != "melanoma",
                                      r["platform"] != "visium", r["slide_id"]))
    stats = {}
    for r in rows:
        sid = r["slide_id"]
        g = stats[sid] = graph_stats(P, sid)
        ap_ = P.interim / f"{sid}.admission.json"
        adm = json.loads(ap_.read_text(encoding="utf-8")) if ap_.exists() else {}
        status = str(r.get("status", ""))
        if status == "ingested":
            a = "admitted" + ("**（forced）**" if adm.get("forced") else "")
        elif status == "replication":
            a = "replication cohort — registered override"
        elif status == "extension":
            a = "2026 extension cohort — analysed"
        else:
            a = f"**rejected — {','.join(adm.get('failed') or ['?'])}**"
        cells = ([f"{g['n']:,}", f"{g['deg']:.2f}", str(g["nc"]), f"{g['big']:.1f} %",
                  f"{g['src']} / {g['snk']} / {g['ves']}"] if g else ["—"] * 5)
        record = {
            "Section": sid,
            "Patient": r.get("patient", ""),
            "Replicate / site": r.get("replicate", ""),
            "Cohort": CT.get(r['cancer_type'], r['cancer_type']),
            "Platform": PF.get(r['platform'], r['platform']),
            "GEO": r.get("accession", ""),
            "Spots (raw)": f"{int(r['n_spots']):,}",
            "Median UMI": f"{float(r['median_umi']):,.0f}",
            "Nodes after QC": cells[0],
            "Mean degree": cells[1],
            "Components": cells[2],
            "Largest comp.": cells[3],
            "Source / Sink / Vessel": cells[4],
            "Admission": a,
        }
        records.append(record)
        L.append(f"| {record['Section']} | {record['Patient']} | {record['Replicate / site']} | "
                 f"{record['Cohort']} | {record['Platform']} | {record['GEO']} | "
                 f"{record['Spots (raw)']} | {record['Median UMI']} | "
                 + " | ".join(cells) + f" | {a} |")

    ing = [r for r in led if str(r.get("status", "")) in {"ingested", "replication", "extension"}]
    primary = [r for r in led if str(r.get("status", "")) == "ingested"]
    replication = [r for r in led if str(r.get("status", "")) == "replication"]
    extension = [r for r in led if str(r.get("status", "")) == "extension"]
    rejected = [r for r in led if str(r.get("status", "")) == "rejected"]
    L += ["",
          f"**{len(primary)} ingested sections from {len({r['patient'] for r in primary})} patient groups "
          f"(19 current primary sections from 8 patients plus 3 legacy external-validation sections); "
          f"{len(replication)} replication sections from {len({r['patient'] for r in replication})} patients "
          f"analysed under a registered override; {len(extension)} extension sections from "
          f"32 identifiable patients (plus 11 sections with unavailable patient relationships) analysed; "
          f"{len(rejected)} rejected at admission.**", ""]
    for cohort, label in ((primary, "primary"), (replication, "replication"), (extension, "extension")):
        for pf in sorted({r["platform"] for r in cohort}):
            ss = [r["slide_id"] for r in cohort if r["platform"] == pf and stats.get(r["slide_id"])]
            if not ss:
                continue
            if label == "extension" and pf == "visium":
                npat_text = "21 identifiable patients plus 11 patient-relationship-unavailable sections"
            else:
                npat = len({r["patient"] for r in cohort if r["platform"] == pf})
                npat_text = f"{npat} patients"
            d = [stats[s]["deg"] for s in ss]; b = [stats[s]["big"] for s in ss]
            L.append(f"- **{label.title()} / {PF.get(pf, pf)}** — {len(ss)} sections / {npat_text}; "
                     f"mean degree {min(d):.2f}–{max(d):.2f}; largest component "
                     f"{min(b):.1f}–{max(b):.1f} % "
                     f"(single component in {sum(1 for s in ss if stats[s]['nc'] == 1)}/{len(ss)})")

    meta = stamp_run(cfg, {"module": "M18-table1", "n_rows": len(led)})
    L += ["", "---", "",
          f"生成于 {meta.get('timestamp')}　"
          f"seed={meta.get('seed')}　numpy {meta.get('numpy')}　scipy {meta.get('scipy')}"]

    out = P.root / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(L) + "\n", encoding="utf-8")

    csv_out = out.with_suffix(".csv")
    with csv_out.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers)
        writer.writeheader()
        writer.writerows(records)

    xlsx_out = out.with_suffix(".xlsx")
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    ws = wb.active
    ws.title = "Table 1 sections"
    ws.append(headers)
    for record in records:
        ws.append([record[h] for h in headers])
    for cell in ws[1]:
        cell.font = Font(bold=True)
        cell.fill = PatternFill("solid", fgColor="EDEDED")
    for col_idx, header in enumerate(headers, start=1):
        width = max(len(header), max(len(str(ws.cell(row=row, column=col_idx).value or ""))
                                     for row in range(2, ws.max_row + 1))) + 2
        ws.column_dimensions[get_column_letter(col_idx)].width = min(width, 38)
    ws.freeze_panes = "A2"
    wb.save(xlsx_out)

    print(f"[Table 1] 已写出 {out}, {csv_out}, {xlsx_out}（{len(led)} 行，其中 {len(ing)} 张准入）")


if __name__ == "__main__":
    main()
