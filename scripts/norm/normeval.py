import sys, numpy as np
sys.path.insert(0, "scripts")
from csf_weights import read_cache, scores, fit_plcc, csf_weights, evaluate
d = sys.argv[1]
c50 = read_cache(d + "/c50.bin"); c100 = read_cache(d + "/c100.bin"); conds = [c50, c100]; fb = c100["f"]
U = np.full((6, 6), 2.0 / 36); r = evaluate(conds, U)["heldout"]
print("uniform      %.4f %.4f %.4f" % (r["50cm"], r["100cm"], r["pooled"]))
r = evaluate(conds, csf_weights(fb, L=40, Xo=20))["heldout"]
print("CSF          %.4f %.4f %.4f" % (r["50cm"], r["100cm"], r["pooled"]))
betas = np.arange(0, 1.21, 0.05); S = {}
for b in betas:
    W = csf_weights(fb, L=40, Xo=20, total=2.0 - b); W[0, 0] += b; S[b] = [scores(c, W) for c in conds]
ho = {"50cm": [], "100cm": [], "pooled": []}; chosen = []
for f in range(5):
    best = None
    for b in betas:
        acc = sum(fit_plcc(s[(c["group"] % 5) != f], c["subj"][(c["group"] % 5) != f])[0] for c, s in zip(conds, S[b]))
        if best is None or acc > best[0]: best = (acc, b)
    b = best[1]; chosen.append(round(float(b), 2)); xs = []; ys = []
    for c, s, lab in zip(conds, S[b], ("50cm", "100cm")):
        m = (c["group"] % 5) == f; ho[lab].append(fit_plcc(s[m], c["subj"][m])[0]); xs.append(s[m]); ys.append(c["subj"][m])
    ho["pooled"].append(fit_plcc(np.concatenate(xs), np.concatenate(ys))[0])
print("CSF+lowest   %.4f %.4f %.4f  beta %s" % (np.mean(ho["50cm"]), np.mean(ho["100cm"]), np.mean(ho["pooled"]), chosen))
