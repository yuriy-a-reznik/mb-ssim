#!/usr/bin/env python3
"""Summarise the LIVE Release 2 ppd sweep into LaTeX table rows."""
import re, sys, glob, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PPDS = [15, 20, 25, 30, 40, 55]

rows = []
for P in PPDS:
    log = f"{ROOT}/runs/live_f{P}.log"
    if not os.path.exists(log):
        rows.append((P, None, None, None)); continue
    nb, plcc, srcc = None, None, None
    for L in open(log):
        m = re.search(r"(\d+) bands", L)
        if m and nb is None: nb = int(m.group(1))
        m = re.match(r"\s+\S+\s+n=\s*\d+\s+PLCC\s+([\d.]+)\s+SRCC\s+([\d.]+)", L)
        if m: plcc, srcc = float(m.group(1)), float(m.group(2))   # last = final fit
    rows.append((P, nb, plcc, srcc))

def row(vals, fmt="%s"):
    return " & ".join(fmt % v if v is not None else "--" for v in vals)

print("\\renewcommand{\\LIVEppdrow}{%s}"  % row([r[0] for r in rows], "%d"))
print("\\renewcommand{\\LIVEnbrow}{%s}"   % row([r[1] for r in rows], "%d"))
print("\\renewcommand{\\LIVEplccrow}{%s}" % row([r[2] for r in rows], "%.4f"))
print("\\renewcommand{\\LIVEsrccrow}{%s}" % row([r[3] for r in rows], "%.4f"))
best = max((r for r in rows if r[2] is not None), key=lambda r: r[2], default=None)
if best:
    print("%% best ppd = %d (PLCC %.4f)" % (best[0], best[2]))
