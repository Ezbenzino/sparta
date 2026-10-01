# 批量管线：对每个 slide 依次跑 摄入→QC(force)→打分→建图→决策
import subprocess, sys, os

PY = r"D:\sparta\.venv\Scripts\python.exe"
SCRIPTS = r"D:\sparta\scripts"

SLIDES = {
    "MEL02": "GSM7983364",
    "MEL03": "GSM7983365",
    "MEL04": "GSM7983366",
}

def run(cmd, check_file=None):
    print(f"\n>>> {' '.join(cmd)}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    # 只打印告警/关键行（去重、截断）
    lines = [l for l in (r.stdout + r.stderr).splitlines() if l.strip()]
    for l in lines:
        if any(k in l for k in ("[M", "[warn]", "警告", "✓", "✗", "已写出", "判定", "中位", "解离", "连通", "节点", "跳", "失败", "Error", "Traceback")):
            print(l[:200])
    if check_file and not os.path.exists(check_file):
        print(f"!! 产物缺失: {check_file}")
        return False
    return True

for sid, acc in SLIDES.items():
    print(f"\n{'='*60}\n处理 {sid} ({acc})\n{'='*60}")
    ok = run([PY, os.path.join(SCRIPTS, "run_00b_ingest.py"),
              "--input", rf"D:\sparta\data\raw\{sid}", "--slide", sid,
              "--cancer", "melanoma", "--platform", "visium",
              "--source", "GEO GSE250636", "--accession", acc,
              "--treatment", "unknown", "--site", "metastasis"],
             rf"D:\sparta\data\interim\{sid}.raw.h5ad")
    if not ok: continue

    ok = run([PY, os.path.join(SCRIPTS, "run_01_qc.py"),
              "--slide", sid, "--input", rf"D:\sparta\data\interim\{sid}.raw.h5ad",
              "--platform", "visium", "--has-he", "--force"],
             rf"D:\sparta\data\interim\{sid}.qc.h5ad")
    if not ok: continue

    ok = run([PY, os.path.join(SCRIPTS, "run_02_score.py"),
              "--slide", sid, "--tumor-type", "melanoma"],
             rf"D:\sparta\data\interim\{sid}.scored.h5ad")
    if not ok: continue

    ok = run([PY, os.path.join(SCRIPTS, "run_03_graph.py"), "--slide", sid],
             rf"D:\sparta\data\interim\{sid}.graph.npz")
    if not ok: continue

    ok = run([PY, os.path.join(SCRIPTS, "run_07_screen.py"), "--slides", sid],
             rf"D:\sparta\results\validation\screen_decision.json")

print("\n全部完成")
