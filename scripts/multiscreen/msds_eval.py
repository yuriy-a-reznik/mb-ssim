#!/usr/bin/env python3
"""Evaluate SSIM, MS-SSIM and MB-SSIM variants on the BC-KU Multi-Screen dataset
at the clip level, with leave-one-content-out (LOCO) refits of MB-SSIM.

MB-SSIM variants:
  transfer   CID:IQ joint weights, no fitting on this data
  joint      one weight set refit across the three devices (nbsjoint, LOCO)
  perdev     weights refit separately per device (LOCO)
  uniform    w_ij = 2/N^2 (no fitting)

Scores for MB-SSIM come from the cached log-xi tensors (nbsjoint --dump-scores);
SSIM / MS-SSIM from nbsbench.  Everything is averaged over a clip's frames
before any correlation is computed.
"""
import os, sys, json, csv, subprocess, math, warnings
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import pearsonr, spearmanr

ROOT = next((a for a in sys.argv[1:] if not a.startswith("--")), None) or os.environ.get("MSDS_WORK", ".")   # work dir with cache_*.bin, bench_ssim.csv, mos.json, frames/manifest.json
NBS = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DEVS = ["mobile", "tablet", "tv"]
COL = {"mobile": "M", "tablet": "T", "tv": "TV"}
CONTENTS = ["s01", "s02", "s03"]
FIT = "--combine=mean --objective=plcc --folds=0 --restarts=2 --basin=4 --maxit=3000 --eval-stride=4"
SYM = "--preset=wilson --weight-sym=off"
mos = json.load(open(ROOT + "/mos.json"))
manifest = json.load(open(ROOT + "/frames/manifest.json"))
PVS = [m["pvs"] for m in manifest]


def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(cmd + "\n" + r.stdout[-500:] + r.stderr[-500:])
    return r.stdout


def cache(dev): return "%s/cache_%s.bin" % (ROOT, dev)


def subset(dev, keep=None, drop=None):
    tag = "keep_" + keep if keep else "drop_" + drop
    out = "%s/cache_%s_%s.bin" % (ROOT, dev, tag)
    if not os.path.exists(out):
        gid = {"s01": 0, "s02": 1, "s03": 2}
        arg = "--keep=%d" % gid[keep] if keep else "--drop=%d" % gid[drop]
        sh("%s/bin/nbscachesel %s %s %s" % (NBS, cache(dev), out, arg))
    return out


def dump(cachefile, weights, out):
    """per-pair MB-SSIM scores from a cache with given weights -> dict index->score"""
    sh("%s/bin/nbsjoint %s --cond=%s --folds=0 --load-weights=%s --dump-scores=%s > /dev/null"
       % (NBS, SYM, cachefile, weights, out))
    d = {}
    for row in csv.DictReader(l for l in open(out) if not l.startswith("#")):
        d[int(row["index"])] = float(row["score"])
    return d


def clip_means(scores_by_index, pvs_order):
    """index -> score, 15 frames per clip in manifest order -> clip means"""
    per = {}
    n = len(scores_by_index)
    k = n // len(pvs_order)
    for i, p in enumerate(pvs_order):
        vals = [scores_by_index[j] for j in range(i * k, (i + 1) * k) if j in scores_by_index]
        per[p] = float(np.mean(vals))
    return per


def logistic(s, b1, b2, b3, b4, b5):
    return b1 * (0.5 - 1.0 / (1.0 + np.exp(np.clip(b2 * (s - b3), -50, 50)))) + b4 * s + b5


def fit_eval(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    best = None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for b2 in (1.0, 10.0, 50.0, 200.0):
            for sign in (1, -1):
                p0 = [sign * (y.max() - y.min()), b2, x.mean(), 0.0, y.mean()]
                try:
                    p, _ = curve_fit(logistic, x, y, p0=p0, maxfev=20000)
                    q = logistic(x, *p)
                    rmse = math.sqrt(np.mean((q - y) ** 2))
                    if best is None or rmse < best[0]: best = (rmse, q)
                except Exception:
                    continue
    if best is None:   # linear fallback
        A = np.vstack([x, np.ones_like(x)]).T
        c = np.linalg.lstsq(A, y, rcond=None)[0]; q = A @ c
        best = (math.sqrt(np.mean((q - y) ** 2)), q)
    rmse, q = best
    return pearsonr(q, y)[0], spearmanr(x, y)[0], rmse


def main():
    # ---------- per-clip scores for every metric and device ----------
    S = {}   # S[metric][dev][pvs]
    # SSIM / MS-SSIM (viewing-blind): one bench run
    ss, ms = {}, {}
    for row in csv.DictReader(l for l in open(ROOT + "/bench_ssim.csv") if not l.startswith("#")):
        p = row["dist"].split("/")[0]
        ss.setdefault(p, []).append(float(row["ssim"])); ms.setdefault(p, []).append(float(row["msssim"]))
    S["SSIM"] = {d: {p: float(np.mean(ss[p])) for p in PVS} for d in DEVS}
    S["MS-SSIM"] = {d: {p: float(np.mean(ms[p])) for p in PVS} for d in DEVS}

    # uniform and transferred weights
    wu = ROOT + "/w_uniform.txt"
    open(wu, "w").write("# nbssim exponents, 6 bands, asymmetric, full matrix\n" + "\n".join(" ".join("%.6f" % (2.0 / 36) for _ in range(6)) for _ in range(6)) + "\n")
    S["MB uniform"] = {d: clip_means(dump(cache(d), wu, ROOT + "/dump_unif_%s.csv" % d), PVS) for d in DEVS}
    S["MB transfer (CID:IQ w)"] = {d: clip_means(dump(cache(d), NBS + "/results/cidiq/w_cidiq_joint.txt",
                                                      ROOT + "/dump_xfer_%s.csv" % d), PVS) for d in DEVS}

    if "--refit" not in sys.argv:          # the paper reports only the transferred weights
        json.dump(S, open(ROOT + "/clip_scores.json", "w"), indent=1)
        return report(S)
    # per-device LOCO refit
    S["MB per-device refit"] = {d: {} for d in DEVS}
    for d in DEVS:
        for g in CONTENTS:
            w = "%s/w_perdev_%s_drop_%s.txt" % (ROOT, d, g)
            if not os.path.exists(w):
                sh("%s/bin/nbsjoint %s --label=%s --cond=%s %s --out=%s > %s.log" % (NBS, SYM, d, subset(d, drop=g), FIT, w, w))
            held = dump(subset(d, keep=g), w, ROOT + "/dump_perdev_%s_%s.csv" % (d, g))
            order = [p for p in PVS if p.startswith(g)]
            S["MB per-device refit"][d].update(clip_means(held, order))
    # joint LOCO refit across devices
    S["MB joint refit"] = {d: {} for d in DEVS}
    for g in CONTENTS:
        w = "%s/w_joint_drop_%s.txt" % (ROOT, g)
        if not os.path.exists(w):
            conds = " ".join("--label=%s --cond=%s" % (d, subset(d, drop=g)) for d in DEVS)
            sh("%s/bin/nbsjoint %s %s %s --out=%s > %s.log" % (NBS, SYM, conds, FIT, w, w))
        for d in DEVS:
            held = dump(subset(d, keep=g), w, ROOT + "/dump_joint_%s_%s.csv" % (d, g))
            order = [p for p in PVS if p.startswith(g)]
            S["MB joint refit"][d].update(clip_means(held, order))
    json.dump(S, open(ROOT + "/clip_scores.json", "w"), indent=1)
    report(S)


def report(S):
    # ---------- evaluation ----------
    print("\nper-device logistic (27 clips each); PLCC / SRCC / RMSE")
    hdr = "%-26s" % "metric" + "".join("%26s" % d for d in DEVS) + "%26s" % "mean of 3"
    print(hdr)
    table = {}
    for m in S:
        cells = []; accs = []
        for d in DEVS:
            x = [S[m][d][p] for p in PVS]; y = [mos[p][COL[d]] for p in PVS]
            plcc, srcc, rmse = fit_eval(x, y); accs.append((plcc, srcc, rmse))
            cells.append("%7.3f %7.3f %7.3f" % (plcc, srcc, rmse))
        mp = np.mean([a[0] for a in accs]); msr = np.mean([a[1] for a in accs]); mr = np.mean([a[2] for a in accs])
        table[m] = dict(perdev=accs, mean=(mp, msr, mr))
        print("%-26s" % m + "".join("%26s" % c for c in cells) + "%26s" % ("%7.3f %7.3f %7.3f" % (mp, msr, mr)))

    print("\npooled: one logistic over all 81 (clip, device) points; PLCC / SRCC / RMSE")
    for m in S:
        x = [S[m][d][p] for d in DEVS for p in PVS]; y = [mos[p][COL[d]] for d in DEVS for p in PVS]
        plcc, srcc, rmse = fit_eval(x, y); table[m]["pooled"] = (plcc, srcc, rmse)
        print("%-26s %7.3f %7.3f %7.3f" % (m, plcc, srcc, rmse))

    print("\ndevice separation: mean score difference mobile - tv, and sign agreement with MOS (mobile > tv in all 27 clips)")
    for m in S:
        dm = np.array([S[m]["mobile"][p] - S[m]["tv"][p] for p in PVS])
        print("%-26s mean diff %+.4f   clips with mobile>tv: %2d/27" % (m, dm.mean(), int((dm > 0).sum())))
    json.dump(table, open(ROOT + "/eval_table.json", "w"), indent=1)


if __name__ == "__main__":
    main()
