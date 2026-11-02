#!/usr/bin/env python3
"""SSIM / MS-SSIM on exactly the pairs and folds the LIVE sweep used.

nbstune --subset=2 keeps every second entry of the index; nbsjoint numbers
the surviving reference images in order of first appearance and splits the
folds by that number.  Reproducing both here is what makes the baseline
columns of the LIVE table comparable to the MB-SSIM column.
"""
import sys, os
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import spearmanr, kendalltau

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUBSET = int(sys.argv[1]) if len(sys.argv) > 1 else 2
FOLDS  = int(sys.argv[2]) if len(sys.argv) > 2 else 3

def vqeg(x, b1, b2, b3, b4, b5):
    return b1 * (0.5 - 1.0 / (1.0 + np.exp(b2 * (x - b3)))) + b4 * x + b5

def stats(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    best, bp = -1, None
    for sc in (0.02, 0.2, 1.0, 5.0, 20.0, 100.0):
        for sgn in (1, -1):
            try:
                p, _ = curve_fit(vqeg, x, y, maxfev=200000,
                                 p0=[y.max()-y.min(), sgn*sc, x.mean(), 0.0, y.mean()])
                pr = vqeg(x, *p)
                if np.all(np.isfinite(pr)):
                    r = abs(np.corrcoef(pr, y)[0, 1])
                    if r > best: best, bp = r, p
            except Exception:
                pass
    pred = vqeg(x, *bp)
    return (best, abs(spearmanr(x, y).statistic), abs(kendalltau(x, y).statistic),
            float(np.sqrt(np.mean((pred - y) ** 2))))

rows = []
for line in open(f"{ROOT}/runs/live_base.csv"):
    if line.startswith('#') or line.startswith('ref,'):
        continue
    f = line.rstrip('\n').split(',')
    rows.append((f[2], float(f[3]), float(f[4]), float(f[5])))
rows = rows[::SUBSET]
gid, g = {}, []
for r in rows:
    gid.setdefault(r[0], len(gid)); g.append(gid[r[0]])
print("%d pairs, %d references, %d folds" % (len(rows), len(gid), FOLDS))

for name, col in (("SSIM", 2), ("MS-SSIM", 3)):
    acc = []
    for f in range(FOLDS):
        sel = [k for k, x in enumerate(g) if x % FOLDS == f]
        acc.append(stats([rows[k][col] for k in sel], [rows[k][1] for k in sel]))
    m = np.mean(acc, axis=0)
    print("%-8s held-out  PLCC %.4f  SRCC %.4f  KRCC %.4f  RMSE %.4f" % (name, *m))
