#!/usr/bin/env python3
"""Build nbsdataset index files for the multi-screen frames, one per device.

Rows are frame pairs (reference frame, upscaled rendition frame) carrying the
MOS of their clip on that device; group = source content, so that
--folds=3 in nbstune/nbsjoint is leave-one-content-out.
"""
import os, sys, json, glob, csv
import openpyxl

# usage: msds_build_index.py WORKDIR   (WORKDIR holds frames/ from msds_pipeline.py and Results.xlsx)
ROOT = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("MSDS_WORK", ".")
FRAMES = os.path.join(ROOT, "frames")
DEV = {"mobile": ("M", 97.9), "tablet": ("T", 57.4), "tv": ("TV", 56.55)}

wb = openpyxl.load_workbook(os.path.join(ROOT, "Results.xlsx"), data_only=True)
ws = wb["MOS and VQA Scores"]
mos = {}
for r in ws.iter_rows(min_row=3, values_only=True):
    if r[0] and str(r[0]).endswith(".mp4"):
        name = r[0][:-4]
        mos[name] = {"M": float(r[6]), "T": float(r[7]), "TV": float(r[8]),
                     "psnr": float(r[2]), "ssim_xlsx": float(r[3]), "vif": float(r[4]), "vmaf": float(r[5])}
json.dump(mos, open(os.path.join(ROOT, "mos.json"), "w"), indent=1)
manifest = json.load(open(os.path.join(FRAMES, "manifest.json")))
for dev, (col, ppd) in DEV.items():
    path = os.path.join(ROOT, "msds_%s.csv" % dev)
    n = 0
    with open(path, "w") as f:
        f.write("# nbsdataset 1\n# name: MultiScreen-%s\n# root: %s\n# score: mos\n# lower_better: 0\n# ppd: %g\n" % (dev, FRAMES, ppd))
        for m in manifest:
            pvs = m["pvs"]
            refs = sorted(glob.glob(os.path.join(FRAMES, pvs, "ref", "*.bmp")))
            dists = sorted(glob.glob(os.path.join(FRAMES, pvs, "dist", "*.bmp")))
            for r, d in zip(refs, dists):
                f.write("%s,%s,%.6f,%s\n" % (os.path.relpath(r, FRAMES), os.path.relpath(d, FRAMES), mos[pvs][col], pvs.split("_")[0]))
                n += 1
    print("wrote", path, n, "pairs")
