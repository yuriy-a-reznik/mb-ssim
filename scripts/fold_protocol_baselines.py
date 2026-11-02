#!/usr/bin/env python3
"""SSIM / MS-SSIM on CID:IQ under both protocols (all-data logistic, and the
paper's five content folds with a logistic per fold), from results/cidiq/full_*.csv."""
import numpy as np, csv, sys
sys.path.insert(0, "scripts")
from csf_weights import fit_plcc

def load(path):
    rows = [r for r in csv.DictReader(l for l in open(path) if not l.startswith("#"))]
    refs = sorted(set(r["refid"] for r in rows), key=lambda s: int(s.replace("final", "")))
    gid = {r: i for i, r in enumerate(refs)}
    g = np.array([gid[r["refid"]] for r in rows])
    y = np.array([float(r["subj"]) for r in rows])
    return g, y, {k: np.array([float(r[k]) for r in rows]) for k in ("ssim", "msssim", "nbssim")}

R = "results/cidiq/"
g50, y50, s50 = load(R + "full_50cm.csv")
g100, y100, s100 = load(R + "full_100cm.csv")
print("groups:", g50.max() + 1, "n:", len(y50), len(y100))
for m in ("ssim", "msssim", "nbssim"):
    a = fit_plcc(s50[m], y50)[0]; b = fit_plcc(s100[m], y100)[0]
    p = fit_plcc(np.concatenate([s50[m], s100[m]]), np.concatenate([y50, y100]))[0]
    ho = {"50": [], "100": [], "pool": []}
    for f in range(5):
        m50 = (g50 % 5) == f; m100 = (g100 % 5) == f
        ho["50"].append(fit_plcc(s50[m][m50], y50[m50])[0])
        ho["100"].append(fit_plcc(s100[m][m100], y100[m100])[0])
        ho["pool"].append(fit_plcc(np.concatenate([s50[m][m50], s100[m][m100]]),
                                   np.concatenate([y50[m50], y100[m100]]))[0])
    print("%-7s all data: 50cm %.4f 100cm %.4f pooled %.4f | fold-mean: 50cm %.4f 100cm %.4f pooled %.4f"
          % (m, a, b, p, np.mean(ho["50"]), np.mean(ho["100"]), np.mean(ho["pool"])))
