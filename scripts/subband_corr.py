#!/usr/bin/env python3
"""How independent are the N^2 sub-bands, as the fit sees them?

The fitter only ever sees each image through its per-sub-band mean log-xi.
If two sub-bands give nearly the same value on every image, no amount of
data can decide how to split an exponent between them -- the objective is
flat along that direction, and the specific sparse pattern a search returns
says more about the search than about vision.  This reads the cache and
reports the correlation matrix of those per-image values.
"""
import struct, sys
import numpy as np

PATH = sys.argv[1] if len(sys.argv) > 1 else "runs/cache_100cm.bin"
PEAKS = ["0.8", "1.7", "2.8", "4", "8", "16"]

with open(PATH, "rb") as f:
    assert f.read(9) == b"NBSCACHE3"
    n, nw, stride = struct.unpack("<iii", f.read(12))
    total, = struct.unpack("<q", f.read(8))
    nbands, w_sym = struct.unpack("<ii", f.read(8))
    ppd, = struct.unpack("<d", f.read(8))
    f.read(8 * 12 * 3)                      # f[], rho[], kc[]
    off = np.frombuffer(f.read(8 * n), dtype="<i8")
    npx = np.frombuffer(f.read(4 * n), dtype="<i4")
    subj = np.frombuffer(f.read(8 * n), dtype="<f8")
    f.read(4 * n)                           # group
    L = np.frombuffer(f.read(4 * total), dtype="<f4")

print("%s: %d pairs, %d bands, %d slots%s, ppd %.3f"
      % (PATH, n, nbands, nw, "" if w_sym else " (free)", ppd))

# per-image mean log-xi for each slot: exactly what a unit exponent buys
M = np.empty((n, nw))
for e in range(n):
    blk = L[off[e]:off[e] + npx[e] * nw].reshape(npx[e], nw)
    M[e] = blk.mean(axis=0)

C = np.corrcoef(M.T)
lab = ["%s x %s" % (PEAKS[i], PEAKS[j]) for i in range(nbands) for j in range(nbands)]

print("\nmost collinear sub-band pairs (|r| over the 690 images):")
pairs = [(abs(C[a, b]), a, b) for a in range(nw) for b in range(a + 1, nw)]
pairs.sort(reverse=True)
for r, a, b in pairs[:10]:
    print("  r = %.4f   (%s)  vs  (%s)" % (r, lab[a], lab[b]))

print("\nmirror pairs (i,j) vs (j,i) -- the ones the symmetry constraint ties:")
for i in range(nbands):
    for j in range(i + 1, nbands):
        a, b = i * nbands + j, j * nbands + i
        print("  r = %.4f   (%s)  vs  (%s)" % (C[a, b], lab[a], lab[b]))

print("\ncorrelation with the subjective score (single-sub-band predictors):")
order = sorted(range(nw), key=lambda k: -abs(np.corrcoef(M[:, k], subj)[0, 1]))
for k in order[:8]:
    print("  r = %+.4f   %s" % (np.corrcoef(M[:, k], subj)[0, 1], lab[k]))

ev = np.linalg.eigvalsh(C)[::-1]
cum = np.cumsum(ev) / nw
k95 = int(np.searchsorted(cum, 0.95)) + 1
print("\neffective dimension of the %d sub-bands: %d components carry 95%% of "
      "the variance\n  (condition number %.3g)" % (nw, k95, ev[0] / max(ev[-1], 1e-300)))
