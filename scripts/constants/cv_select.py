#!/usr/bin/env python3
"""Cross-validated choice of the stabilising-constant level K and the lowest
sub-band weight beta for the CSF + lowest-sub-band form of MB-SSIM.
For each of the paper's five content folds, (K, beta) is chosen by nbsjoint's
objective (mean per-condition PLCC after the logistic) on the training groups
only, and the held-out fold is scored with that choice."""
import sys, numpy as np
sys.path.insert(0, "scripts")
from csf_weights import read_cache, scores, fit_plcc, csf_weights

LEVELS = {  # scale of (K1, K2) relative to SSIM's (0.01, 0.03) -> cache dir
    0.1: "runs/constants/tenth", 0.2: "runs/constants/fifth", 1/3: "runs/constants/third",
    0.5: "runs/constants/half", 1.0: None, 3.0: "runs/constants/triple"}
betas = np.arange(0, 2.01, 0.05)

def load(d):
    if d is None:
        return [read_cache("runs/csf/cache_50cm.bin"), read_cache("runs/csf/cache_100cm.bin")]
    return [read_cache(d + "/c50.bin"), read_cache(d + "/c100.bin")]

import os, pickle
S = {}   # (level, beta) -> [scores50, scores100]
if os.path.exists("cv_scores.pkl"):
    S, groups, subj = pickle.load(open("cv_scores.pkl", "rb"))
else:
    conds = None
    for lev, d in LEVELS.items():
        c = load(d); conds = c; fb = c[1]["f"]
        for b in betas:
            W = csf_weights(fb, L=40, Xo=20, total=2.0 - b); W[0, 0] += b
            S[(lev, b)] = [scores(cc, W) for cc in c]
        print("scored level %.3f" % lev, flush=True)
    groups = [c["group"] for c in conds]; subj = [c["subj"] for c in conds]
    pickle.dump((S, groups, subj), open("cv_scores.pkl", "wb"))

def train_obj(key, f):
    acc = 0.0
    for i in range(2):
        m = (groups[i] % 5) != f
        acc += fit_plcc(S[key][i][m], subj[i][m])[0]
    return acc / 2

ho = {"50cm": [], "100cm": [], "pooled": []}; chosen = []
for f in range(5):
    best = max(S.keys(), key=lambda k: train_obj(k, f))
    chosen.append((round(best[0], 3), round(float(best[1]), 2)))
    xs, ys = [], []
    for i, lab in enumerate(("50cm", "100cm")):
        m = (groups[i] % 5) == f
        ho[lab].append(fit_plcc(S[best][i][m], subj[i][m])[0]); xs.append(S[best][i][m]); ys.append(subj[i][m])
    ho["pooled"].append(fit_plcc(np.concatenate(xs), np.concatenate(ys))[0])
    print("fold %d: chosen (K scale, beta) = %s  held-out 50cm %.4f 100cm %.4f pooled %.4f"
          % (f, chosen[-1], ho["50cm"][-1], ho["100cm"][-1], ho["pooled"][-1]), flush=True)
print("K and beta chosen on training folds: mean held-out 50cm %.4f  100cm %.4f  pooled %.4f"
      % (np.mean(ho["50cm"]), np.mean(ho["100cm"]), np.mean(ho["pooled"])))

# the same with K fixed at each level (beta chosen on training folds), for the table
for lev in LEVELS:
    ho2 = {"50cm": [], "100cm": [], "pooled": []}
    for f in range(5):
        best = max([k for k in S if k[0] == lev], key=lambda k: train_obj(k, f))
        xs, ys = [], []
        for i, lab in enumerate(("50cm", "100cm")):
            m = (groups[i] % 5) == f
            ho2[lab].append(fit_plcc(S[best][i][m], subj[i][m])[0]); xs.append(S[best][i][m]); ys.append(subj[i][m])
        ho2["pooled"].append(fit_plcc(np.concatenate(xs), np.concatenate(ys))[0])
    print("K scale %.3f, beta chosen per fold: 50cm %.4f  100cm %.4f  pooled %.4f"
          % (lev, np.mean(ho2["50cm"]), np.mean(ho2["100cm"]), np.mean(ho2["pooled"])))
