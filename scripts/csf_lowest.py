#!/usr/bin/env python3
"""Table 7, row "CSF + lowest sub-band": Barten CSF weights scaled to sum to
2 - beta, plus beta on the lowest sub-band (0.8 x 0.8 c/deg).  beta is the one
fitted number; it is chosen on the training groups of each of the paper's five
content folds by nbsjoint's objective (mean per-condition PLCC after the
logistic) and the fold is scored held out.  Also prints the single-ingredient
rows quoted in the text (lowest sub-band alone; CSF without it).

usage: python3 scripts/csf_lowest.py [cache_dir]
"""
import sys, numpy as np
sys.path.insert(0, "scripts")
from csf_weights import read_cache, scores, fit_plcc, csf_weights, fmt, evaluate

root = sys.argv[1] if len(sys.argv) > 1 else "runs/csf"
c50 = read_cache(root + "/cache_50cm.bin"); c100 = read_cache(root + "/cache_100cm.bin")
conds = [c50, c100]; fb = c100["f"]
betas = np.arange(0, 1.21, 0.05)

W = np.zeros((6, 6)); W[0, 0] = 1.0
print("lowest sub-band alone      ", fmt(evaluate(conds, W)))
W = csf_weights(fb, L=40, Xo=20, total=1.4); W[0, 0] = 0
print("CSF without lowest sub-band", fmt(evaluate(conds, W)))

for kind, p in (("uniform", 1.0), ("csf", 1.0), ("csf", 2.0)):
    S = {}
    for b in betas:
        W = csf_weights(fb, power=p, L=40, Xo=20, total=2.0 - b) if kind == "csf" else np.full((6, 6), (2.0 - b) / 36)
        W[0, 0] += b
        S[b] = [scores(c, W) for c in conds]
    ho = {"50cm": [], "100cm": [], "pooled": []}; chosen = []
    for f in range(5):
        best = None
        for b in betas:
            acc = 0.0
            for c, s in zip(conds, S[b]):
                m = (c["group"] % 5) != f
                acc += fit_plcc(s[m], c["subj"][m])[0]
            if best is None or acc > best[0]: best = (acc, b)
        b = best[1]; chosen.append(round(float(b), 2))
        xs, ys = [], []
        for c, s, lab in zip(conds, S[b], ("50cm", "100cm")):
            m = (c["group"] % 5) == f
            ho[lab].append(fit_plcc(s[m], c["subj"][m])[0]); xs.append(s[m]); ys.append(c["subj"][m])
        ho["pooled"].append(fit_plcc(np.concatenate(xs), np.concatenate(ys))[0])
    print("%s^%.0f + beta on lowest, beta per fold %s: held-out 50cm %.4f 100cm %.4f pooled %.4f"
          % (kind, p, chosen, np.mean(ho["50cm"]), np.mean(ho["100cm"]), np.mean(ho["pooled"])))
