#!/usr/bin/env python3
"""CSF-weights / square-root-compression experiment on the CID:IQ log-xi caches
(Table 7 of the paper, rows "CSF" and "CSF + lowest sub-band" in csf_lowest.py).

usage: python3 scripts/csf_weights.py [cache_dir] [weights_file]
  cache_dir holds cache_50cm.bin and cache_100cm.bin from REPRODUCE.md section 7c.

Reads the NBSCACHE3 tensors written by nbstune, scores every image under a
given master 6x6 weight matrix and a given per-band compression, and reports
PLCC after the VQEG logistic (a) on all data per condition and pooled, as the
paper does for SSIM, MS-SSIM and the uniform row, and (b) as the mean over
the paper's five content folds (group % 5), which is how the fitted rows of
the paper are reported.
"""
import sys, struct, numpy as np
from scipy.optimize import curve_fit
from scipy.stats import pearsonr, spearmanr

MAXB = 12

def read_cache(path):
    with open(path, "rb") as f:
        magic = f.read(9)
        assert magic == b"NBSCACHE3", magic
        n, nw, stride = struct.unpack("iii", f.read(12))
        total, = struct.unpack("q", f.read(8))
        nbands, w_sym = struct.unpack("ii", f.read(8))
        ppd, = struct.unpack("d", f.read(8))
        fb = np.frombuffer(f.read(8 * MAXB), dtype=np.float64)[:nbands].copy()
        rho = np.frombuffer(f.read(8 * MAXB), dtype=np.float64)[:nbands].copy()
        kc = np.frombuffer(f.read(8 * MAXB), dtype=np.float64)[:nbands].copy()
        off = np.frombuffer(f.read(8 * n), dtype=np.int64).copy()
        npx = np.frombuffer(f.read(4 * n), dtype=np.int32).copy()
        subj = np.frombuffer(f.read(8 * n), dtype=np.float64).copy()
        group = np.frombuffer(f.read(4 * n), dtype=np.int32).copy()
        L = np.frombuffer(f.read(4 * total), dtype=np.float32)
    assert w_sym == 0 and nw == nbands * nbands
    return dict(n=n, nw=nw, nbands=nbands, ppd=ppd, f=fb, off=off, npx=npx,
                subj=subj, group=group, L=L.reshape(-1, nw))

def cond_weights(wmaster, nb):
    """master 6x6 -> flat row-major vector over the condition's nb bands."""
    return wmaster[:nb, :nb].reshape(-1)

def scores(c, wmaster, comp="log"):
    """MB-SSIM per image.  comp='log': mean_p exp(sum w L);
    comp='sqrt': mean_p exp(-sum w sqrt(-L))."""
    w = cond_weights(wmaster, c["nbands"])
    L = c["L"]
    if comp == "log":
        s = L @ w
    elif comp == "sqrt":
        s = -(np.sqrt(np.maximum(-L, 0.0)) @ w)
    else:
        raise ValueError(comp)
    e = np.exp(s)
    out = np.empty(c["n"])
    for i in range(c["n"]):
        a = c["off"][i] // c["nw"]
        out[i] = e[a:a + c["npx"][i]].mean()
    return out

def logistic(x, b1, b2, b3, b4, b5):
    return b1 * (0.5 - 1.0 / (1.0 + np.exp(np.clip(b2 * (x - b3), -500, 500)))) + b4 * x + b5

def fit_plcc(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    best = None
    for b2sign in (1.0, -1.0):
        p0 = [np.max(y) - np.min(y), b2sign * 4.0 / (np.std(x) + 1e-12), np.mean(x), 0.0, np.mean(y)]
        try:
            p, _ = curve_fit(logistic, x, y, p0=p0, maxfev=20000)
            pred = logistic(x, *p)
            r = pearsonr(pred, y)[0]
            rmse = np.sqrt(np.mean((pred - y) ** 2))
            if best is None or rmse < best[1]:
                best = (abs(r), rmse)
        except Exception:
            pass
    if best is None:
        best = (abs(pearsonr(x, y)[0]), float("nan"))
    return best

def evaluate(conds, wmaster, comp="log", folds=5, verbose=False):
    """Returns dict with all-data per-condition and pooled PLCC/SRCC/RMSE and
    mean held-out (per fold) PLCC for each condition and pooled."""
    sc = [scores(c, wmaster, comp) for c in conds]
    res = {}
    for c, s, lab in zip(conds, sc, ("50cm", "100cm")):
        r, e = fit_plcc(s, c["subj"])
        res[lab] = (r, abs(spearmanr(s, c["subj"])[0]), e)
    allx = np.concatenate(sc); ally = np.concatenate([c["subj"] for c in conds])
    r, e = fit_plcc(allx, ally)
    res["pooled"] = (r, abs(spearmanr(allx, ally)[0]), e)
    # held-out protocol: fold f = groups with g % folds == f, logistic fitted on the fold
    acc = {"50cm": [], "100cm": [], "pooled": []}
    for f in range(folds):
        xs, ys = [], []
        for c, s, lab in zip(conds, sc, ("50cm", "100cm")):
            m = (c["group"] % folds) == f
            acc[lab].append(fit_plcc(s[m], c["subj"][m])[0])
            xs.append(s[m]); ys.append(c["subj"][m])
        acc["pooled"].append(fit_plcc(np.concatenate(xs), np.concatenate(ys))[0])
    res["heldout"] = {k: float(np.mean(v)) for k, v in acc.items()}
    return res

def fmt(res):
    a = res
    return ("all data: 50cm %.4f  100cm %.4f  pooled %.4f (SRCC %.4f, RMSE %.4f)   |  "
            "held-out mean: 50cm %.4f  100cm %.4f  pooled %.4f" %
            (a["50cm"][0], a["100cm"][0], a["pooled"][0], a["pooled"][1], a["pooled"][2],
             a["heldout"]["50cm"], a["heldout"]["100cm"], a["heldout"]["pooled"]))

# ---------------------------------------------------------------- Barten CSF
def barten_S(u, L=40.0, Xo=20.0, k=3.0, T=0.1, Xmax=12.0, Nmax=15.0, eta=0.03,
             Phi0=3e-8, u0=7.0, sigma0=0.5 / 60.0, Cab=0.08 / 60.0, p=1.2e6):
    """Barten (1999) eq. 3.26.  u in c/deg, L in cd/m2, Xo field in deg."""
    u = np.asarray(u, float)
    d = 5.0 - 3.0 * np.tanh(0.4 * np.log10(L * Xo ** 2 / 40.0 ** 2))
    sigma = np.sqrt(sigma0 ** 2 + (Cab * d) ** 2)
    Mopt = np.exp(-2.0 * np.pi ** 2 * sigma ** 2 * u ** 2)
    E = (np.pi * d ** 2 / 4.0) * L * (1.0 - (d / 9.7) ** 2 + (d / 12.4) ** 4)
    A = 1.0 / Xo ** 2 + 1.0 / Xmax ** 2 + u ** 2 / Nmax ** 2
    uu = np.maximum(u, 1e-6)
    N = 1.0 / (eta * p * E) + Phi0 / (1.0 - np.exp(-(uu / u0) ** 2))
    return (Mopt / k) / np.sqrt((2.0 / T) * A * N)

def read_weights(path):
    rows = [l for l in open(path) if not l.startswith("#") and l.strip()]
    return np.array([[float(v) for v in l.split()] for l in rows])

def csf_weights(fb, power=1.0, radial=True, total=2.0, **kw):
    """w_ij proportional to S(f_ij)^power, normalised to sum `total` over 6x6."""
    n = len(fb)
    W = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            fij = np.hypot(fb[i], fb[j]) if radial else max(fb[i], fb[j])
            W[i, j] = barten_S(fij, **kw) ** power
    return W * (total / W.sum())

if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) > 1 else "runs/csf"
    wfile = sys.argv[2] if len(sys.argv) > 2 else "results/cidiq/w_cidiq_joint.txt"
    c50 = read_cache(root + "/cache_50cm.bin")
    c100 = read_cache(root + "/cache_100cm.bin")
    conds = [c50, c100]
    fb = c100["f"]
    print("bands:", fb, " n:", c50["n"], c100["n"], " ppd:", c50["ppd"], c100["ppd"])
    print("stride-16 pixels per image:", c50["npx"][0], c100["npx"][0])

    U = np.full((6, 6), 2.0 / 36.0)
    print("\n[A] uniform 2/36            ", fmt(evaluate(conds, U)))
    Wfit = read_weights(wfile)
    print("[G] fitted (all-data weights)", fmt(evaluate(conds, Wfit)))

    for L in (20.0, 40.0, 80.0):
        for power, tag in ((1.0, "S  "), (0.5, "S^½"), (2.0, "S² ")):
            W = csf_weights(fb, power=power, L=L, Xo=20.0)
            print("[B] w∝%s L=%3.0f Xo=20   " % (tag, L), fmt(evaluate(conds, W)))
    W = csf_weights(fb, power=1.0, L=40.0, Xo=20.0, radial=False)
    print("[B'] w∝S, max(fi,fj) L=40 ", fmt(evaluate(conds, W)))
    # per-condition field: 24.9 deg at 50 cm, 12.6 deg at 100 cm
    W50 = csf_weights(fb, L=40.0, Xo=24.9); W100 = csf_weights(fb, L=40.0, Xo=12.6)
    sc = [scores(c50, W50), scores(c100, W100)]
    r50 = fit_plcc(sc[0], c50["subj"]); r100 = fit_plcc(sc[1], c100["subj"])
    rp = fit_plcc(np.concatenate(sc), np.concatenate([c50["subj"], c100["subj"]]))
    print("[E] w∝S per-field 24.9/12.6  all data: 50cm %.4f  100cm %.4f  pooled %.4f" % (r50[0], r100[0], rp[0]))
    print("    S(f) at 24.9 vs 12.6 deg:", np.round(barten_S(fb, 40, 24.9) / barten_S(fb, 40, 12.6), 3))

    print("\nsquare-root compression  exp(-sum w sqrt(-log xi)):")
    print("[F1] uniform, sqrt           ", fmt(evaluate(conds, U, comp="sqrt")))
    print("[F2] w∝S L=40, sqrt          ", fmt(evaluate(conds, csf_weights(fb, L=40, Xo=20), comp="sqrt")))
    print("[F3] w∝S^½ L=40, sqrt        ", fmt(evaluate(conds, csf_weights(fb, power=0.5, L=40, Xo=20), comp="sqrt")))
    print("[F4] fitted weights, sqrt    ", fmt(evaluate(conds, Wfit, comp="sqrt")))
    print("\nCSF weights (w∝S, L=40, Xo=20), rows i = horizontal band, cols j = vertical band:")
    np.set_printoptions(precision=4, suppress=True)
    print(csf_weights(fb, L=40, Xo=20))
    print("S(f) at band peaks, L=40, Xo=20:", np.round(barten_S(fb, 40, 20), 1))
