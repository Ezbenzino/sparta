# 对 MEL02/03/04 跑 run_04 barrier + run_05 counterfactual
import subprocess, os

PY = r"D:\sparta\.venv\Scripts\python.exe"
SCRIPTS = r"D:\sparta\scripts"
SLIDES = ["MEL02", "MEL03", "MEL04"]

def run(cmd, check_file=None):
    print(f"\n>>> {' '.join(cmd)}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    for l in (r.stdout + r.stderr).splitlines():
        l = l.strip()
        if any(k in l for k in ("[M", "[warn]", "警告", "已写出", "已保存",
                                "S1", "S2", "S3", "mode=", "k=", "r=",
                                "失败", "Error", "Traceback")):
            print(l[:200])
    if check_file and not os.path.exists(check_file):
        print(f"!! 产物缺失: {check_file}")
        return False
    return True

for sid in SLIDES:
    run([PY, os.path.join(SCRIPTS, "run_04_barrier.py"), "--slide", sid],
        rf"D:\sparta\data\interim\{sid}.barrier.npz")
    run([PY, os.path.join(SCRIPTS, "run_05_counterfactual.py"), "--slide", sid],
        rf"D:\sparta\results\counterfactual\{sid}.json")
print("\n全部完成")
