#!/usr/bin/env python3
"""Extended ablation sweep for MB-SSIM: perturb the filterbank, not just its size.

The paper's band-count sweep (runs/sweep_bands.sh) varies only N.  This driver
varies, around the Wilson six-channel reference, the things a reviewer will
ask about next:

  shift     the whole bank slid by -1/3, -1/6, +1/6, +1/3 octave
  jitter    each peak moved independently by U(-1/3, +1/3) octave  (J draws)
  bwscale   all bandwidths scaled by 0.6, 0.8, 1.25, 1.5
  bwconst   every channel given the same 1.5-octave bandwidth
  spacing   peaks linear in frequency between 0.8 and 16 c/deg
  shape     log-Gabor without taper, taper onset 0.8, difference-of-Gaussians
  window    statistics window sigma 1.0 (7x7) and 2.5 (17x17) instead of 1.5
  random    peak jitter and bandwidth scale drawn jointly  (R draws)

Every configuration runs exactly the paper's protocol: two log-xi tensors
(50 and 100 cm, stride 16), one joint five-fold fit with the same search
budget as runs/sweep_bands.sh, and the held-out numbers are read from the
fit log.  Nothing else differs between configurations.

Usage
  scripts/ablation_sweep.py list                      # show the configurations
  scripts/ablation_sweep.py check                     # realise every bank with nbsbank (no images needed)
  scripts/ablation_sweep.py run  [--shard=I/N] [--dry-run]
  scripts/ablation_sweep.py summarise                 # CSV + LaTeX table, paired tests vs the reference

Options: --jitter=J (default 20) --random=R (default 20) --seed=S (default 1)
         --only=glob[,glob...]  run only matching configurations, in that order
         --budget=fast|full (full = the paper's search budget, default;
                 fast = 1 restart, 2 basin rounds, ~4x cheaper; logs are ablf_*)
         --jobs=2 (build the two tensors of a configuration concurrently)
         --families=shift,bwscale,... (run only these perturbation families)
         --root=DIR (repo root; default: parent of this script)
         --data50=PATH --data100=PATH (default data/cidiq_{50,100}cm.csv)

`run` is resumable: a configuration whose fit log already contains the
"mean held-out" block is skipped, and the two tensors are deleted after the
fit so that disk use stays at two tensors per running shard.  Cost per
configuration is roughly 35 core-minutes (README timings), so the default
55-configuration sweep is about 32 core-hours; --shard splits it.
"""
import os, re, sys, math, random, subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WIL_PEAKS = [0.8, 1.7, 2.8, 4.0, 8.0, 16.0]
WIL_BW    = [2.5, 2.2, 1.9, 1.7, 1.5, 1.25]
PPD = {"50cm": 31.612, "100cm": 63.223}

FIT_FULL = ("--combine=mean --objective=plcc --folds=5 --restarts=2 --basin=6 "
            "--maxit=3000 --eval-stride=8")          # the paper's budget (runs/sweep_bands.sh)
FIT_FAST = ("--combine=mean --objective=plcc --folds=5 --restarts=1 --basin=2 "
            "--maxit=3000 --eval-stride=8")          # ~4x cheaper; use for every config incl. the reference
FIT = FIT_FULL
TAG = "abl"
CACHE = "--transform=dft --weight-sym=off --stride=16 --progress --maxit=1 --restarts=1 --folds=0"


def interp_bw(peak):
    """Wilson's bandwidth interpolated (in log frequency) at an arbitrary peak."""
    lp = [math.log2(p) for p in WIL_PEAKS]
    x = math.log2(peak)
    if x <= lp[0]: return WIL_BW[0]
    if x >= lp[-1]: return WIL_BW[-1]
    for k in range(len(lp) - 1):
        if lp[k] <= x <= lp[k + 1]:
            t = (x - lp[k]) / (lp[k + 1] - lp[k])
            return WIL_BW[k] + t * (WIL_BW[k + 1] - WIL_BW[k])


def bank_args(peaks, bws, bank="loggabor", taper="1.0"):
    s = ("--bank=%s --bands=%d --peaks=%s --first=bandpass --last=bandpass "
         "--bw=%s --band-norm=peak" % (bank, len(peaks),
                                       ",".join("%.4f" % p for p in peaks),
                                       ",".join("%.4f" % b for b in bws)))
    if bank == "loggabor":
        s += " --taper=%s" % taper
    return s


def configs(J=20, R=20, seed=1):
    rng = random.Random(seed)
    C = []
    def add(name, family, peaks, bws, bank="loggabor", taper="1.0", extra=""):
        C.append(dict(name=name, family=family, peaks=list(peaks), bws=list(bws),
                      args=bank_args(peaks, bws, bank, taper) + (" " + extra if extra else "")))
    add("ref_wilson", "reference", WIL_PEAKS, WIL_BW)
    for s in (-1/3, -1/6, 1/6, 1/3):
        add("shift_%+.3foct" % s, "shift", [p * 2 ** s for p in WIL_PEAKS], WIL_BW)
    for f in (0.6, 0.8, 1.25, 1.5):
        add("bwscale_%.2f" % f, "bwscale", WIL_PEAKS, [b * f for b in WIL_BW])
    add("bwconst_1.5", "bwconst", WIL_PEAKS, [1.5] * 6)
    lin = [0.8 + k * (16 - 0.8) / 5 for k in range(6)]
    add("spacing_linear", "spacing", lin, [interp_bw(p) for p in lin])
    add("shape_notaper", "shape", WIL_PEAKS, WIL_BW, taper="off")
    add("shape_taper0.8", "shape", WIL_PEAKS, WIL_BW, taper="0.8")
    add("shape_dog", "shape", WIL_PEAKS, WIL_BW, bank="dog")
    add("window_s1.0", "window", WIL_PEAKS, WIL_BW, extra="--nbs-win=7 --win-sigma=1.0")
    add("window_s2.5", "window", WIL_PEAKS, WIL_BW, extra="--nbs-win=17 --win-sigma=2.5")
    for k in range(J):
        pk = sorted(p * 2 ** rng.uniform(-1/3, 1/3) for p in WIL_PEAKS)
        add("jitter_%02d" % k, "jitter", pk, [interp_bw(p) for p in pk])
    for k in range(R):
        pk = sorted(p * 2 ** rng.uniform(-1/3, 1/3) for p in WIL_PEAKS)
        f = rng.uniform(0.7, 1.4)
        add("random_%02d" % k, "random", pk, [interp_bw(p) * f for p in pk])
    return C


def fit_log(c):   return os.path.join(ROOT, "runs", "%s_%s.log" % (TAG, c["name"]))
def cache(c, cond): return os.path.join(ROOT, "runs", "abl_%s_%s.bin" % (c["name"], cond))
def wfile(c):     return os.path.join(ROOT, "runs", "%s_w_%s.txt" % (TAG, c["name"]))


def done(c):
    p = fit_log(c)
    return os.path.exists(p) and "mean held-out" in open(p).read()


def commands(c, data):
    cmds = []
    for cond in ("50cm", "100cm"):
        if os.path.exists(cache(c, cond)):
            cmds.append("true  # tensor %s already built" % cache(c, cond)); continue
        cmds.append("%s/bin/nbstune %s %s --ppd=%s --cache-save=%s --out=/dev/null %s > %s/runs/abl_%s_%s.log 2>&1"
                    % (ROOT, c["args"], CACHE, PPD[cond], cache(c, cond), data[cond],
                       ROOT, c["name"], "c" + cond))
    cmds.append("%s/bin/nbsjoint --label=50cm --cond=%s --label=100cm --cond=%s %s --out=%s > %s 2>&1"
                % (ROOT, cache(c, "50cm"), cache(c, "100cm"), FIT, wfile(c), fit_log(c)))
    cmds.append("rm -f %s %s" % (cache(c, "50cm"), cache(c, "100cm")))
    return cmds


def parse_folds(path):
    per, on = [], False
    for L in open(path):
        if L.startswith("5-fold"): on = True; continue
        if "mean held-out" in L: break
        if not on: continue
        m = re.match(r"\s+pooled\s+n=\s*\d+\s+PLCC\s+([\d.]+)", L)
        if m: per.append(float(m.group(1)))
    mean = {}
    on = False
    for L in open(path):
        if "mean held-out" in L: on = True; continue
        if on:
            m = re.match(r"\s+(\S+)\s+PLCC\s+([\d.]+)\s+SRCC\s+([\d.]+)(?:\s+KRCC\s+[\d.]+)?\s+RMSE\s+([\d.]+)", L)
            if m: mean[m.group(1)] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
            elif L.startswith("final"): break
    return per, mean


def paired_t(a, b):
    d = [x - y for x, y in zip(a, b)]
    n = len(d)
    if n < 2: return float("nan"), float("nan")
    m = sum(d) / n
    s = math.sqrt(sum((x - m) ** 2 for x in d) / (n - 1))
    if s == 0: return m, float("nan")
    t = m / (s / math.sqrt(n))
    try:
        from scipy.stats import t as tdist
        p = 2 * (1 - tdist.cdf(abs(t), n - 1))
    except Exception:
        p = float("nan")
    return m, p


def main():
    args = sys.argv[1:]
    if not args or args[0] not in ("list", "check", "run", "summarise"):
        print(__doc__); sys.exit(1)
    mode = args[0]
    opt = dict(jitter=20, random=20, seed=1, shard="1/1",
               data50=os.path.join(ROOT, "data", "cidiq_50cm.csv"),
               data100=os.path.join(ROOT, "data", "cidiq_100cm.csv"))
    dry = False
    for a in args[1:]:
        if a == "--dry-run": dry = True; continue
        k, v = a.lstrip("-").split("=", 1)
        opt[k] = int(v) if k in ("jitter", "random", "seed") else v
    global FIT, TAG
    if opt.get("budget", "full") == "fast":
        FIT, TAG = FIT_FAST, "ablf"
    C = configs(opt["jitter"], opt["random"], opt["seed"])
    data = {"50cm": opt["data50"], "100cm": opt["data100"]}

    if mode == "list":
        print("%-18s %-9s  %-46s %s" % ("name", "family", "peaks [c/deg]", "bandwidths [oct]"))
        for c in C:
            print("%-18s %-9s  %-46s %s" % (c["name"], c["family"],
                  " ".join("%.2f" % p for p in c["peaks"]), " ".join("%.2f" % b for b in c["bws"])))
        print("\n%d configurations" % len(C))
        return

    if mode == "check":
        bad = 0
        for c in C:
            for cond, ppd in PPD.items():
                cmd = "%s/bin/nbsbank %s --ppd=%s --points=50" % (ROOT, c["args"].replace("--nbs-win=7 ", "").replace("--nbs-win=17 ", ""), ppd)
                r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                if r.returncode != 0:
                    bad += 1; print("FAIL", c["name"], cond, r.stderr.strip()[:200])
        print("%d configurations x 2 conditions realised, %d failures" % (len(C), bad))
        return

    if mode == "run":
        i, n = (int(x) for x in opt["shard"].split("/"))
        fams = opt.get("families")
        if "only" in opt:                      # --only=glob,glob : restrict and order by pattern
            import fnmatch
            sel = []
            for pat in opt["only"].split(","):
                sel += [c for c in C if fnmatch.fnmatch(c["name"], pat) and c not in sel]
            C = sel
        mine = [c for k, c in enumerate(C) if k % n == i - 1
                and (not fams or c["family"] in fams.split(","))]
        jobs = int(opt.get("jobs", 1))
        os.makedirs(os.path.join(ROOT, "runs"), exist_ok=True)
        for c in mine:
            if done(c):
                print("skip  ", c["name"]); continue
            print("run   ", c["name"], flush=True)
            cmds = commands(c, data)
            if dry:
                for cmd in cmds: print("   ", cmd)
                continue
            ok = True
            if jobs >= 2:   # the two tensors are independent: build them side by side
                procs = [subprocess.Popen(cmd, shell=True) for cmd in cmds[:2]]
                ok = all(p.wait() == 0 for p in procs)
                rest = cmds[2:]
            else:
                rest = cmds
            for cmd in rest:
                if not ok: break
                if subprocess.run(cmd, shell=True).returncode != 0:
                    ok = False
            if not ok: print("FAILED:", c["name"])
        return

    if mode == "summarise":
        rows = []
        for c in C:
            if not done(c): continue
            per, mean = parse_folds(fit_log(c))
            rows.append((c, per, mean))
        if not rows:
            print("no finished configurations"); return
        ref = next((r for r in rows if r[0]["name"] == "ref_wilson"), None)
        if ref is None:
            ref = next((r for r in rows if "fit_win" in r[0]["name"]), rows[0])
        out = os.path.join(ROOT, "runs", "%s_sweep.csv" % TAG)
        with open(out, "w") as f:
            f.write("name,family,peaks,bandwidths,plcc50,plcc100,plcc_pooled,srcc_pooled,rmse_pooled,delta_vs_ref,p_paired\n")
            print("%-18s %-9s %8s %8s %8s %8s %8s %9s %6s" % ("name", "family", "PLCC50", "PLCC100", "pooled", "SRCC", "RMSE", "d vs ref", "p"))
            for c, per, mean in rows:
                d, p = paired_t(per, ref[1])
                g = lambda k, i: mean.get(k, (float('nan'),) * 3)[i]
                f.write("%s,%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.3f\n" % (
                    c["name"], c["family"], " ".join("%.3f" % x for x in c["peaks"]),
                    " ".join("%.3f" % x for x in c["bws"]), g("50cm", 0), g("100cm", 0),
                    g("pooled", 0), g("pooled", 1), g("pooled", 2), d, p))
                print("%-18s %-9s %8.4f %8.4f %8.4f %8.4f %8.4f %+9.4f %6.2f" % (
                    c["name"], c["family"], g("50cm", 0), g("100cm", 0), g("pooled", 0),
                    g("pooled", 1), g("pooled", 2), d, p))
        # family summary for the paper: best / worst / mean pooled PLCC per family
        print("\nper family (held-out pooled PLCC):")
        fam = {}
        for c, per, mean in rows:
            fam.setdefault(c["family"], []).append(mean.get("pooled", (float('nan'),))[0])
        tex = [r"\begin{tabular}{@{}lrccc@{}}", r"\toprule",
               r"perturbation & configs & min & mean & max \\", r"\midrule"]
        for k, v in fam.items():
            v = [x for x in v if x == x]
            if not v: continue
            print("  %-10s n=%2d  min %.4f  mean %.4f  max %.4f" % (k, len(v), min(v), sum(v) / len(v), max(v)))
            tex.append(r"%s & %d & %.4f & %.4f & %.4f \\" % (k, len(v), min(v), sum(v) / len(v), max(v)))
        tex += [r"\bottomrule", r"\end{tabular}"]
        open(os.path.join(ROOT, "runs", "table_%s_sweep.tex" % TAG), "w").write("\n".join(tex) + "\n")
        print("\nwrote", out, "and runs/table_%s_sweep.tex" % TAG)


if __name__ == "__main__":
    main()
