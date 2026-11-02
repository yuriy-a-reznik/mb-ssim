#!/usr/bin/env python3
"""Summarise the band-count sweep: held-out five-fold PLCC/SRCC/RMSE per N.

Reads runs/fit_n{3,4,5,6}.log (and, for reference, the Wilson six-channel fit
runs/fit_win.log) and prints one row per bank, plus the LaTeX table body.
"""
import re, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def folds(path):
    d, on = {}, False
    for L in open(path):
        if L.startswith("5-fold"): on = True; continue
        if "mean held-out" in L or L.startswith("final"): break
        if not on: continue
        m = re.match(r"\s+(\S+)\s+n=\s*(\d+)\s+PLCC\s+([\d.]+)\s+SRCC\s+([\d.]+)"
                     r"(?:\s+KRCC\s+([\d.]+)\s+RMSE\s+([\d.]+))?", L)
        if m:
            d.setdefault(m.group(1), []).append(
                tuple(float(m.group(i)) if m.group(i) else float('nan')
                      for i in (3, 4, 5, 6)))
    return d

def mean(d):
    return {k: tuple(np.nanmean([v[i] for v in vs]) for i in range(4))
            for k, vs in d.items()}

def header(path):
    """channels actually used at each condition, and the shared weight count"""
    nb, nw = {}, None
    for L in open(path):
        m = re.match(r"\s+(\S+)\s+ppd\s+[\d.]+\s+(\d+) bands\s+(\d+) exponents", L)
        if m: nb[m.group(1)] = int(m.group(2))
        m2 = re.match(r"master bank : (\d+) bands, (\d+) shared", L)
        if m2: nw = (int(m2.group(1)), int(m2.group(2)))
    return nb, nw

RUNS = [("3", f"{ROOT}/runs/fit_n3.log"),
        ("4", f"{ROOT}/runs/fit_n4.log"),
        ("5", f"{ROOT}/runs/fit_n5.log"),
        ("6", f"{ROOT}/runs/fit_n6.log")]
REF  = ("6 (Wilson)", f"{ROOT}/runs/fit_win.log")

rows = []
for name, p in RUNS + [REF]:
    if not os.path.exists(p): continue
    f = folds(p)
    if "pooled" not in f: continue
    m = mean(f)
    nb, nw = header(p)
    per = [v[0] for v in f["pooled"]]
    rows.append(dict(name=name, m=m, nb=nb, nw=nw, per=np.array(per)))

print("%-11s %-7s %-9s   %-7s %-7s %-7s   %-7s %-7s" %
      ("bank", "weights", "ch 50/100", "PLCC50", "PLCC100", "PLCCpool",
       "SRCCpool", "RMSEpool"))
for r in rows:
    nb = r["nb"]
    print("N = %-7s %-7s %-9s   %.4f  %.4f  %.4f    %.4f  %.4f" %
          (r["name"], r["nw"][1] if r["nw"] else "?",
           "%d/%d" % (nb.get("50cm", 0), nb.get("100cm", 0)),
           r["m"]["50cm"][0], r["m"]["100cm"][0], r["m"]["pooled"][0],
           r["m"]["pooled"][1], r["m"]["pooled"][3]))

if len(rows) > 1:
    from scipy import stats
    base = [r for r in rows if r["name"] == "6 (Wilson)"]
    if base:
        b = base[0]
        print("\npaired against the Wilson six-channel bank (pooled PLCC, 5 folds):")
        for r in rows:
            if r is b: continue
            d = r["per"] - b["per"]
            t, p = stats.ttest_rel(r["per"], b["per"])
            print("   N = %-11s diff %+.4f +- %.4f   p = %.3f"
                  % (r["name"], d.mean(), d.std(ddof=1), p))

print("\n%% LaTeX table body")
for r in rows:
    nb = r["nb"]
    print(r"%s & %s & %d\,/\,%d & %.4f & %.4f & %.4f \\" %
          (r["name"], r["nw"][1] if r["nw"] else "?",
           nb.get("50cm", 0), nb.get("100cm", 0),
           r["m"]["50cm"][0], r["m"]["100cm"][0], r["m"]["pooled"][0]))

# ---- emit the LaTeX table ------------------------------------------------
import io
OUT = os.path.join(ROOT, "runs", "table_sweep.tex")
L = []
A = L.append
A(r"\begin{table}[tp]")
A(r"\caption{Band-count sweep on CID:IQ. Held-out five-fold cross-validation,")
A(r"folds split by reference image, one set of sub-band weights fitted jointly")
A(r"across both viewing distances, all weights free; PLCC, higher is better.")
A(r"The four banks of Figure~\ref{fig:sweep} span the same range with $N$")
A(r"channels each; the last row is the Wilson bank of Section~\ref{sec:bank},")
A(r"which differs only in where its intermediate peaks sit. `channels' counts")
A(r"those below the display Nyquist rate at each distance.}")
A(r"\label{tab:sweep}")
A(r"\small")
A(r"\centering")
A(r"\begin{tabular}{@{}lccccc@{}}")
A(r"\toprule")
A(r"bank & weights & channels, $50/100$\,cm & 50\,cm & 100\,cm & pooled \\")
A(r"\midrule")
for r in rows:
    nb = r["nb"]
    nm = "$N=%s$" % r["name"] if r["name"].isdigit() else \
         r"$N=6$, Wilson peaks"
    A(r"%s & %s & %d\,/\,%d & %.4f & %.4f & %.4f \\" %
      (nm, r["nw"][1] if r["nw"] else "?",
       nb.get("50cm", 0), nb.get("100cm", 0),
       r["m"]["50cm"][0], r["m"]["100cm"][0], r["m"]["pooled"][0]))
A(r"\bottomrule")
A(r"\end{tabular}")
A(r"\end{table}")
open(OUT, "w").write("\n".join(L) + "\n")
print("\nwrote", OUT)
