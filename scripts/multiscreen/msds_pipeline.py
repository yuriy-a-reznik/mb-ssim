#!/usr/bin/env python3
"""Frame-wise SSIM / MS-SSIM / MB-SSIM on the BC-KU Multi-Screen dataset.

Pipeline (per processed video sequence, PVS):
  1. decode the HEVC rendition and upscale it to 1920x1080 (bicubic, as a
     player would), decode the 1080p source; dump every STEP-th frame of
     both as 8-bit BMP (luma extraction happens inside nbsmetric, BT.601).
  2. run nbsmetric on each frame pair once per device, at that device's
     *effective* angular sampling density on the 1080p grid:
        TV     55" 3840x2160 at 3H      -> ppd 113.1, video x2.0   -> 56.5
        Tablet 12.4" 2800x1752 at 18in  -> ppd  83.7, video x1.458 -> 57.4
        Mobile 6.1" 2400x1080 at 13in   -> ppd  97.9, video x1.0   -> 97.9
     SSIM and MS-SSIM are viewing-blind and computed once.
  3. average over frames -> one score per (PVS, metric, device).

Usage:
  msds_pipeline.py frames  SRC_DIR TEST_DIR OUT_DIR [--step=10]
  msds_pipeline.py score   OUT_DIR NBSMETRIC [--weights=W.txt] [--jobs=2]
  msds_pipeline.py collect OUT_DIR RESULTS.xlsx
"""
import os, sys, re, glob, subprocess, json, csv, math
from concurrent.futures import ThreadPoolExecutor

DEV = {"tv": 56.55, "tablet": 57.4, "mobile": 97.9}   # effective ppd on the 1080p grid


def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(cmd + "\n" + r.stderr[-800:])
    return r.stdout


def nframes(path):
    out = run('ffprobe -v error -count_frames -select_streams v:0 -show_entries stream=nb_read_frames -of csv=p=0 "%s"' % path)
    return int(out.strip())


def dump(path, outdir, step, scale=True):
    os.makedirs(outdir, exist_ok=True)
    vf = "select='not(mod(n\\,%d))'" % step
    if scale:
        vf += ",scale=1920:1080:flags=bicubic"
    run('ffmpeg -v error -y -i "%s" -vf "%s" -vsync vfr -pix_fmt bgr24 "%s/f%%04d.bmp"' % (path, vf, outdir))
    return sorted(glob.glob(outdir + "/f*.bmp"))


def cmd_frames(src_dir, test_dir, out, step):
    os.makedirs(out, exist_ok=True)
    srcs = {os.path.splitext(os.path.basename(p))[0]: p for p in glob.glob(src_dir + "/*")}
    tests = sorted(glob.glob(test_dir + "/*"))
    manifest = []
    for t in tests:
        name = os.path.splitext(os.path.basename(t))[0]
        # match the source by the longest source name that prefixes the PVS name
        key = name.split("_")[0].lower()           # s01_H_540p -> s01 ; source s01_1080p
        src = next((s for s in srcs if s.lower().split("_")[0] == key), None)
        if src is None:
            print("no source for", name); continue
        d = os.path.join(out, name)
        if not os.path.exists(d + "/done"):
            dump(t, d + "/dist", step, scale=True)
            dump(srcs[src], d + "/ref", step, scale=True)
            open(d + "/done", "w").write("ok")
        nd, nr = len(glob.glob(d + "/dist/*.bmp")), len(glob.glob(d + "/ref/*.bmp"))
        manifest.append(dict(pvs=name, src=src, frames=min(nd, nr)))
        print("%-40s src=%-20s frames %d/%d" % (name, src, nd, nr), flush=True)
    json.dump(manifest, open(out + "/manifest.json", "w"), indent=1)


def cmd_score(out, nbsmetric, weights, jobs):
    manifest = json.load(open(out + "/manifest.json"))
    base = "--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --quiet"
    if weights:
        base += " --weight-file=%s" % weights
    tasks = []
    for m in manifest:
        d = os.path.join(out, m["pvs"])
        refs = sorted(glob.glob(d + "/ref/*.bmp")); dists = sorted(glob.glob(d + "/dist/*.bmp"))
        for r, x in zip(refs, dists):
            for dev, ppd in DEV.items():  # per device
                resf = x.replace(".bmp", ".%s.txt" % dev)
                if not os.path.exists(resf):
                    tasks.append((m["pvs"], dev, ppd, r, x, resf))
    print(len(tasks), "frame x device evaluations to run", flush=True)

    def work(t):
        pvs, dev, ppd, r, x, resf = t
        o = run('%s %s --ppd=%s "%s" "%s"' % (nbsmetric, base, ppd, r, x))
        open(resf, "w").write(o)
        return pvs, dev
    done = 0
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        for _ in ex.map(work, tasks):
            done += 1
            if done % 50 == 0: print("  %d / %d" % (done, len(tasks)), flush=True)


def cmd_collect(out, xlsx):
    manifest = json.load(open(out + "/manifest.json"))
    rows = []
    for m in manifest:
        d = os.path.join(out, m["pvs"])
        for dev in DEV:
            vals = []
            for f in sorted(glob.glob(d + "/dist/*.%s.txt" % dev)):
                s = open(f).read().split()
                if len(s) >= 3: vals.append([float(v) for v in s[:3]])
            if not vals: continue
            n = len(vals)
            ssim, ms, mb = (sum(v[i] for v in vals) / n for i in range(3))
            rows.append(dict(pvs=m["pvs"], device=dev, frames=n, ssim=ssim, msssim=ms, mbssim=mb, ppd=DEV[dev]))
    with open(out + "/scores.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote", out + "/scores.csv", len(rows), "rows")


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a: print(__doc__); sys.exit(1)
    opts = {k.lstrip("-").split("=")[0]: k.split("=", 1)[1] for k in a if k.startswith("--") and "=" in k}
    pos = [k for k in a if not k.startswith("--")]
    if pos[0] == "frames":
        cmd_frames(pos[1], pos[2], pos[3], int(opts.get("step", 10)))
    elif pos[0] == "score":
        cmd_score(pos[1], pos[2], opts.get("weights"), int(opts.get("jobs", 2)))
    elif pos[0] == "collect":
        cmd_collect(pos[1], pos[2] if len(pos) > 2 else None)
