# 批量下载 GSE250636 其余样本（并行）
import subprocess, os, sys, time

SAMPLES = [
    # (slide_id, gsm, sample_prefix, desc)
    ("MEL02", "GSM7983364", "sample14", "cecal nodule"),
    ("MEL03", "GSM7983365", "sample15", "right upper chest wall"),
    ("MEL04", "GSM7983366", "sample16", "ribcage"),
]

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7983nnn/{gsm}/suppl/{gsm}_{pre}_"
FILES = [
    ("filtered_feature_bc_matrix.h5", "filtered_feature_bc_matrix.h5"),
    ("tissue_positions_list.csv.gz", "spatial/tissue_positions_list.csv.gz"),
    ("tissue_hires_image.png.gz", "spatial/tissue_hires_image.png.gz"),
    ("tissue_lowres_image.png.gz", "spatial/tissue_lowres_image.png.gz"),
    ("scalefactors_json.json.gz", "spatial/scalefactors_json.json.gz"),
]


def download_one(slide, gsm, pre):
    raw_dir = os.path.join(r"D:\sparta\data\raw", slide)
    os.makedirs(os.path.join(raw_dir, "spatial"), exist_ok=True)
    for remote, local in FILES:
        url = BASE.format(gsm=gsm, pre=pre) + remote
        dst = os.path.join(raw_dir, local)
        if os.path.exists(dst) and os.path.getsize(dst) > 1000:
            print(f"[{slide}] 已存在，跳过: {local}")
            continue
        # 反复续传直到成功
        for attempt in range(5):
            r = subprocess.run(
                ["curl.exe", "-s", "-L", "-C", "-", "--connect-timeout", "30",
                 "--max-time", "1800", "-o", dst, url],
                capture_output=True, text=True)
            if r.returncode == 0 and os.path.exists(dst) and os.path.getsize(dst) > 1000:
                print(f"[{slide}] OK {local} ({os.path.getsize(dst)} bytes)")
                break
            print(f"[{slide}] 重试 {attempt+1}: {local}")
            time.sleep(2)
        else:
            print(f"[{slide}] 失败: {local}")
    print(f"[{slide}] {desc} 完成")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        slides = sys.argv[1:]
        tasks = [s for s in SAMPLES if s[0] in slides]
    else:
        tasks = SAMPLES
    import concurrent.futures as cf
    with cf.ThreadPoolExecutor(max_workers=len(tasks)) as ex:
        futs = [ex.submit(download_one, sid, gsm, pre) for sid, gsm, pre, desc in tasks]
        for f in cf.as_completed(futs):
            f.result()
    print("ALL DONE")
