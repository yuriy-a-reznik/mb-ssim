#!/usr/bin/env python3
"""Refit the 36 sub-band weights under a per-band compression exponent gamma,
    score = mean_p exp( - sum_ij w_ij (-log xi_ij)^gamma ),
gamma = 1 is MB-SSIM as published, gamma = 1/2 the SQRI-style square root.
Same five content folds as the paper, objective = mean of per-condition PLCC
after the logistic on the training groups (nbsjoint --combine=mean
--objective=plcc), search on every 8th cached pixel, held-out on all pixels.
A deliberately modest optimiser (bounded Powell from the uniform start) applied
identically to every gamma, so the comparison is like for like."""
import sys, time, numpy as np
sys.path.insert(0, "scripts")
from csf_weights import read_cache, fit_plcc
from scipy.optimize import minimize

gammas = [float(g) for g in sys.argv[1:]] or [1.0, 0.5]   # usage: python3 scripts/csf_refit.py 1.0 0.5
c50 = read_cache("runs/csf/cache_50cm.bin")
c100 = read_cache("runs/csf/cache_100cm.bin")
conds = [c50, c100]

def prep(c, gamma, sub):
    """per-image list of (rows) of (-L)^gamma, subsampled by `sub`."""
    T = np.power(np.maximum(-c["L"], 0.0), gamma).astype(np.float32)
    segs = []
    for i in range(c["n"]):
        a = c["off"][i] // c["nw"]
        segs.append(T[a:a + c["npx"][i]:sub])
    lens = np.array([len(s) for s in segs])
    return np.concatenate(segs), np.repeat(np.arange(c["n"]), lens), lens

def score_all(P, idx, lens, w, nb, n):
    s = -(P @ w[:nb, :nb].reshape(-1))
    e = np.exp(s)
    return np.bincount(idx, weights=e, minlength=n) / lens

def run(gamma):
    t0 = time.time()
    search = [prep(c, gamma, 8) for c in conds]
    full = [prep(c, gamma, 1) for c in conds]
    ho = {"50cm": [], "100cm": [], "pooled": []}
    W = []
    for f in range(5):
        tr = [((c["group"] % 5) != f) for c in conds]
        te = [~m for m in tr]
        def obj(v):
            w = v.reshape(6, 6)
            acc = 0.0
            for c, (P, idx, lens), m in zip(conds, search, tr):
                s = score_all(P, idx, lens, w, c["nbands"], c["n"])
                acc += fit_plcc(s[m], c["subj"][m])[0]
            return -acc / len(conds)
        x0 = np.full(36, 2.0 / 36.0)
        r = minimize(obj, x0, method="Powell", bounds=[(0.0, 3.0)] * 36,
                     options={"maxfev": 6000, "xtol": 1e-3, "ftol": 1e-5})
        w = r.x.reshape(6, 6)
        W.append(w)
        xs, ys = [], []
        for c, (P, idx, lens), m, lab in zip(conds, full, te, ("50cm", "100cm")):
            s = score_all(P, idx, lens, w, c["nbands"], c["n"])
            ho[lab].append(fit_plcc(s[m], c["subj"][m])[0]); xs.append(s[m]); ys.append(c["subj"][m])
        ho["pooled"].append(fit_plcc(np.concatenate(xs), np.concatenate(ys))[0])
        print("  gamma %.2f fold %d: train obj %.4f  held-out 50cm %.4f 100cm %.4f pooled %.4f  (%d evals, %.0fs)"
              % (gamma, f, -r.fun, ho["50cm"][-1], ho["100cm"][-1], ho["pooled"][-1], r.nfev, time.time() - t0), flush=True)
    print("gamma %.2f  mean held-out: 50cm %.4f  100cm %.4f  pooled %.4f" %
          (gamma, np.mean(ho["50cm"]), np.mean(ho["100cm"]), np.mean(ho["pooled"])), flush=True)
    np.save("runs/csf/w_gamma_%.2f.npy" % gamma, np.array(W))
    np.set_printoptions(precision=3, suppress=True)
    print("mean fold weights:\n", np.mean(W, axis=0), flush=True)

for g in gammas:
    run(g)
