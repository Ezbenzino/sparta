# 解压 MEL02/03/04 的 gz 文件
import gzip, shutil, os

SLIDES = ["MEL02", "MEL03", "MEL04"]
GZ_MAP = {
    "tissue_positions_list.csv.gz": "tissue_positions_list.csv",
    "tissue_hires_image.png.gz": "tissue_hires_image.png",
    "tissue_lowres_image.png.gz": "tissue_lowres_image.png",
    "scalefactors_json.json.gz": "scalefactors_json.json",
}

for sid in SLIDES:
    spatial = os.path.join(r"D:\sparta\data\raw", sid, "spatial")
    for gz, out in GZ_MAP.items():
        src = os.path.join(spatial, gz)
        dst = os.path.join(spatial, out)
        if os.path.exists(dst) and os.path.getsize(dst) > 0:
            print(f"[{sid}] 已存在 {out}")
            continue
        if not os.path.exists(src):
            print(f"[{sid}] 缺少 {gz}")
            continue
        with gzip.open(src, "rb") as f_in, open(dst, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(src)
        print(f"[{sid}] 解压 {gz} -> {out} ({os.path.getsize(dst)} bytes)")
print("DONE")
