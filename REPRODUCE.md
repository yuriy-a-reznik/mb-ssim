# Reproducing the EI 2027 paper

Every number and figure in *"MB-SSIM: A Multi-Band Variant of SSIM,
Tunable to Viewing Conditions"* comes from the commands below. Timings are for one core of a 2.8 GHz Xeon; the whole
sequence takes about three hours on two cores.

## 0. Build and self-test

```sh
make
make check          # 333 assertions, must report 0 failures
```

The assertions that matter for the paper's claims:

| assertion | what it protects |
|---|---|
| `analytic band responses` | the realised bank *is* the Wilson bank: peaks within 0.01 c/deg, unit gain, bandwidths within 0.01 octave; the taper moves neither |
| `dft core against the spatial core` | the two cores agree to 5e-3 over the same pixels, and the DFT answer is bit-identical when `--radius-mult` is halved while the spatial answer is not |
| `dft log-xi cache` | the tensor the fitters optimise against reproduces the direct metric to 1e-12 |
| `weight file across band counts` | a six-channel exponent vector read by a five-channel bank keeps each exponent's `(i,j)` identity |
| `free (asymmetric) exponents` | a mirrored 36-entry matrix reproduces the 21-entry symmetric metric exactly, and transposing it does change the score |
| `EI-2023 two-band identity` | the two-band reading of SSIM, worst \|Δ\| = 0.015 |

## 1. The database

CID:IQ is available from the Norwegian Colour and Visual Computing
Laboratory. Unpack it and build the two index files — one per viewing
distance. The `ppd` values are not ours: they come from the database's own
`RunMetrics.m`, which fixes the display at 92 dpi and computes
`ppd = dpi × distance_inches × tan(1°)`.

```sh
python3 scripts/make_index.py cidiq /path/to/CID-IQ --distance=50  -o data/cidiq_50cm.csv
python3 scripts/make_index.py cidiq /path/to/CID-IQ --distance=100 -o data/cidiq_100cm.csv
```

690 pairs each, 23 references, 800×800 pixels; ppd 31.612 and 63.223;
image subtends 23.8° and 12.5°.

## 2. Baselines and the untuned metric  (~30 min, two cores)

```sh
for D in 50cm:31.612 100cm:63.223; do
  N=${D%%:*}; P=${D##*:}
  bin/nbsbench --preset=wilson --transform=dft --taper=1.0 --ppd=$P \
               --progress --out=runs/cidiq_$N.csv data/cidiq_$N.csv &
done; wait
```

Reported in the paper as the SSIM and MS-SSIM columns of Table 2, and as the
untuned MB-SSIM starting point.

## 3. The log-xi tensors  (~25 min, two cores)

The metric is `mean_p exp(Σ_k w_k L_k(p))` with `L_k = log ξ_k` independent
of the exponents, so the tensor is computed once and every candidate vector
afterwards costs only dot products.

```sh
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --stride=16 --folds=0 --maxit=1 --restarts=1"
bin/nbstune $C --ppd=31.612 --cache-save=runs/cache_50cm.bin  --out=/dev/null data/cidiq_50cm.csv &
bin/nbstune $C --ppd=63.223 --cache-save=runs/cache_100cm.bin --out=/dev/null data/cidiq_100cm.csv &
wait
```

## 4. The joint fit  (~60 min)

One exponent vector for both viewing distances, five-fold cross-validation
with folds split by reference image:

```sh
bin/nbsjoint --label=50cm  --cond=runs/cache_50cm.bin \
             --label=100cm --cond=runs/cache_100cm.bin \
             --combine=mean --objective=plcc \
             --folds=5 --restarts=2 --basin=6 --maxit=3000 --eval-stride=8 \
             --out=results/w_cidiq_joint.txt
```

The per-fold held-out numbers are Table 2. Swap `--combine=pooled` for the
robustness check quoted in the text (pooled PLCC 0.8343 against 0.8334), and
drop `--weight-sym=off` from step 3 for the 21-parameter symmetric row of
the ablation table (that needs its own tensor, since the packing differs).

## 5. Final scores

```sh
bin/nbsjoint --label=50cm  --cond=runs/cache_50cm.bin \
             --label=100cm --cond=runs/cache_100cm.bin \
             --combine=mean --folds=0 \
             --load-weights=results/w_cidiq_joint.txt \
             --dump-scores=runs/scores_final.csv
```

The per-image scores of `runs/scores_final.csv` are what every figure of
the paper is drawn from; the figure scripts themselves are not part of
this package.

## 6. Design measurements quoted in the text

```sh
# realised bank, both viewing distances
bin/nbsbank --preset=wilson --taper=1.0 --ppd=63.223
bin/nbsbank --preset=wilson --taper=1.0 --ppd=31.612    # reports the dropped channel

# what truncation does to the lowest channel (Figure 2b)
bin/nbsbank --preset=wilson --taper=1.0 --ppd=63.223 --radius-mult=4 --points=4000
bin/nbsbank --preset=wilson --taper=1.0 --ppd=63.223 --radius-mult=2 --points=4000

# runtime, one 800x800 pair
time bin/nbsmetric --preset=wilson --ppd=63.223 --transform=spatial --quiet REF DIST
time bin/nbsmetric --preset=wilson --ppd=63.223 --transform=dft --taper=1.0 --quiet REF DIST
```

## 7. Ablations and the sub-band redundancy analysis

The symmetric row needs a second pair of tensors, built without
`--weight-sym=off`, then the same joint fit. The uniform row, the CSF rows
and the one-parameter row are produced by the scripts of section 7c below,
under the same fold protocol as the fitted rows (the `baseline` block that
`nbsjoint` prints is the uniform set on all data, 0.6589 pooled, not the
fold mean).

The band-count rows -- three channels in place of six, obtained by merging
Wilson's channels in pairs (peak at the geometric mean of the pair,
bandwidth at their sum) -- need their own pair of tensors and one more joint
fit, all of it in `runs/b3_run.sh` (and `runs/b3_sym.sh` for the symmetric
row, which is the same with `--weight-sym=on`):

```sh
B="--bank=loggabor --bands=3 --peaks=1.1662,3.3466,11.3137 \
   --first=bandpass --last=bandpass --bw=4.7,3.6,2.75 --band-norm=peak"
C="$B --transform=dft --taper=1.0 --weight-sym=off --stride=16 \
   --maxit=1 --restarts=1 --folds=0"
bin/nbstune $C --ppd=31.612 --cache-save=runs/b3_50cm.bin  \
            --out=/dev/null data/cidiq_50cm.csv
bin/nbstune $C --ppd=63.223 --cache-save=runs/b3_100cm.bin \
            --out=/dev/null data/cidiq_100cm.csv
bin/nbsjoint --label=50cm --cond=runs/b3_50cm.bin \
             --label=100cm --cond=runs/b3_100cm.bin \
             --combine=mean --objective=plcc --folds=5 --restarts=2 \
             --basin=6 --maxit=3000 --eval-stride=8 --out=runs/w_b3_mean.txt
```

The band-count sweep of Table 5 is `runs/sweep_bands.sh`, which repeats the
same construction for N = 3, 4, 5 and 6 with the first and last peaks pinned
at 0.8 and 16 c/deg, the intermediate peaks spaced log-uniformly between
them, and each bandwidth interpolated from Wilson's values and scaled by
6/N. `scripts/sweep_summary.py` reads the four logs, runs the paired tests
and writes `runs/table_sweep.tex`.

The search budget is deliberately the same as the six-channel fit of
step 4, so that only the bank differs. Tensor construction takes about
20 minutes; the fit about 40.


The redundancy analysis -- the mirror-pair correlations of Figure 7(b), the
condition number and the effective dimension -- comes from

```sh
python3 scripts/subband_corr.py runs/cache_100cm.bin
```

which reads the tensor directly and correlates the per-image mean log-xi of
every sub-band against every other.

## 8. LIVE Release 2  (~1 h, two cores)

LIVE publishes no viewing geometry, so ppd is swept. The primary
experiment transfers the CID:IQ exponents verbatim -- nothing is fitted on
LIVE at all:

```sh
bin/nbsbench --metrics=ssim,msssim --progress --out=runs/live_base.csv \
             data/live_r2_realigned.csv
for P in 15 20 25 30 40 55 80 110; do
  bin/nbsbench --preset=wilson --transform=dft --taper=1.0 --weight-sym=off \
               --weight-file=results/cidiq/w_cidiq_joint.txt --metrics=nbssim \
               --ppd=$P --progress --out=runs/liveU_$P.csv \
               data/live_r2_realigned.csv
done
python3 scripts/live_untuned_summary.py
```

The refitted comparison additionally needs a tensor at the selected ppd and
the usual five-fold joint fit with a single condition.

## 9. Full-resolution confirmation  (~30 min, two cores)

The tensors are built on every sixteenth pixel. To confirm that this does
not bias the reported correlations, rescore every pixel of all 1380 images
with the fitted exponents:

```sh
for D in 50cm:31.612 100cm:63.223; do
  N=${D%%:*}; P=${D##*:}
  bin/nbsbench --preset=wilson --transform=dft --taper=1.0 --weight-sym=off \
               --weight-file=results/w_cidiq_joint.txt --metrics=nbssim \
               --ppd=$P --progress --out=runs/final_$N.csv data/cidiq_$N.csv &
done; wait
```

Pooled PLCC moves by 0.0012.

## 7b. Extended ablation: peak placement, bandwidth, shape, window (Table 7)

`scripts/ablation_sweep.py` perturbs the six-channel Wilson bank along the
axes a reader will ask about next -- a global shift of all peaks (±1/6,
±1/3 octave), independent per-channel jitter, bandwidth scaling (0.6x to
1.5x) and a constant-bandwidth bank, linear instead of log spacing, the
filter shape (log-Gabor with and without taper, difference-of-Gaussians),
the statistics window (sigma 1.0 and 2.5), and joint random draws -- under
exactly the protocol and search budget of `runs/sweep_bands.sh`.

```sh
python3 scripts/ablation_sweep.py list            # 56 configurations by default
python3 scripts/ablation_sweep.py check           # realise every bank with nbsbank; needs no images
python3 scripts/ablation_sweep.py run --shard=1/8 # one of eight shards; resumable
python3 scripts/ablation_sweep.py summarise       # runs/ablation_sweep.csv + table_ablation_sweep.tex
```

Each configuration costs two tensors (~20 min) and one joint fit (~15 min
at this budget), about 35 core-minutes; the default sweep is ~32 core-hours.
`summarise` reports held-out pooled PLCC per configuration with a paired
fold-wise test against the Wilson reference, and a min/mean/max per
perturbation family in the paper's table format.

The 28 configurations reported in Table 7 of the paper were run with
`--budget=fast`; their fit logs are `runs/ablf_*.log`, the fitted weights
`runs/ablf_w_*.txt`, and the summary `runs/ablf_sweep.csv` (the script
also writes a LaTeX version next to it). At the paper's full budget the reference
configuration reproduces Table 2 exactly (`runs/abl_ref_wilson.log`).

## 7c. Contrast-sensitivity weights and the one-parameter form (Table 6)  (~15 min, two cores)

Table 6's parameter-free and one-parameter rows need the two six-channel
tensors once more, this time kept on disk (the joint fit of step 4 builds
them in memory):

```sh
mkdir -p runs/csf
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --stride=16 \
   --maxit=1 --restarts=1 --folds=0 --out=/dev/null"
bin/nbstune $C --ppd=31.612 --cache-save=runs/csf/cache_50cm.bin  data/cidiq_50cm.csv
bin/nbstune $C --ppd=63.223 --cache-save=runs/csf/cache_100cm.bin data/cidiq_100cm.csv
python3 scripts/csf_weights.py            # uniform, Barten CSF (S, sqrt S, S^2; L = 20/40/80 cd/m2;
                                          # per-condition field), square-root compression
python3 scripts/csf_lowest.py             # CSF + lowest sub-band, beta chosen per training fold
python3 scripts/fold_protocol_baselines.py  # SSIM / MS-SSIM under both protocols
python3 scripts/csf_refit.py 1.0 0.75 0.5 # all 36 weights refitted under exp(-sum w (-log xi)^gamma)
```

`csf_weights.py` reads the `NBSCACHE3` tensors directly and reproduces the
`baseline` block of `nbsjoint` (0.7009 / 0.6328 / 0.6589 on all data) and
the all-data score of the released weights (0.8196) to four digits. Every
MB-SSIM number in the paper is a mean over the five content folds (group
index mod 5) with the logistic fitted on the held-out fold; the scripts
print both that and the all-data protocol. Barten's CSF is eq. 3.26 of his
1999 book with its published constants, evaluated at 40 cd/m2 and a 20
degree field (the CID:IQ display is 80 cd/m2 at the white; the 20-80 cd/m2
range moves the pooled result by 0.014 and the actual 12.6/24.9 degree
fields by 0.001). The outputs that went into the paper are in
`results/cidiq/csf/`: `csf_weights.txt`, `csf_lowest.txt`,
`fold_protocol_baselines.txt` and the refit logs `refit_gamma*.log`
(the gamma = 1 refit lands at 0.8315 pooled, within 0.003 of the
`nbsjoint` search; gamma = 3/4 at 0.8196; gamma = 1/2, the square-root
integral's compression, at 0.7969; gamma = 2 about 0.03 below gamma = 1
on the three folds run).

## 7d. The normaliser of the ratio (supplementary to the hidden-gamma paper)  (~2 h, two cores)

`--xi-norm=P` replaces the sum of the two band energies in the ratio by
twice their power mean of order P (`1` is SSIM's arithmetic mean and the
default, `geo` or `0` the geometric mean, i.e. a correlation coefficient,
`max` the larger energy). `scripts/norm/run_norm.sh` builds the two CID:IQ
tensors and runs the fast-budget joint fit for P = geo, max, 0.5 and 2;
`scripts/norm/normeval.py runs/norm/p_<P>` then scores the uniform, CSF and
CSF + lowest-sub-band weight sets on those tensors. Fit logs and weights
are in `results/cidiq/norm/`. The P = 1 reference numbers are those of
sections 7b (fast-budget fit, `runs/ablf_ref_wilson.log`) and 7c.

Results (held-out pooled PLCC, 36 weights refitted): geometric mean 0.773,
P = 0.5 0.822, P = 1 0.831, P = 2 0.830, max 0.761 -- the arithmetic mean
sits at the optimum.

## 7e. The stabilising constants (Table 8)  (~4 h, two cores)

`--kc=a,b,...` sets the per-channel constants k (C = (k L)^2, indexed by
the finer channel of each sub-band); the paper uses SSIM's (0.01, 0.03)
everywhere except in Table 8. `scripts/constants/run_constants.sh` and
`run_constants_more.sh` build the CID:IQ tensors at (K1, K2) scaled by
1/10, 1/5, 1/3, 1/2 and 3, and at CSF-shaped constants, score the
parameter-free weight sets (`scripts/norm/normeval.py`) and run the
fast-budget joint fit at each setting; `scripts/constants/cv_select.py`
chooses the constant level and the lowest-sub-band weight on the training
folds and scores the held-out folds (it needs the section 7c tensors in
`runs/csf/` for the x1 level). Outputs for every setting are in
`results/cidiq/constants/` (`eval.txt`, `fit.log`, `w.txt`) with the
cross-validated selection in `cv_select.txt`.

## 10. Cross-verification on the Brightcove-Kingston multi-screen dataset (Table 4)

The BC-KU Multi-Screen dataset (Barman, Reznik & Martini, QoMEX 2023;
<https://github.com/NabajeetBarman/Multiscreen-Dataset>) rates 27 HEVC
renditions of three 1080p sources on a 55" UHD TV at 3H, a 12.4" tablet at
18 in and a 6.1" phone at 12-14 in. On the 1080p grid the effective angular
sampling densities are 56.55 (TV, video x2), 57.4 (tablet, x1.458) and
97.9 c/deg (phone, 1:1); `scripts/multiscreen/msds_pipeline.py` has the
derivation in its docstring.

```sh
S=scripts/multiscreen; W=/path/to/work        # W holds Results.xlsx
python3 $S/msds_pipeline.py frames /path/source_sequences /path/test_sequences $W/frames --step=20
python3 $S/msds_build_index.py $W            # msds_{mobile,tablet,tv}.csv, mos.json
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --stride=16 --maxit=1 --restarts=1 --folds=0"
for D in tv:56.55 tablet:57.4 mobile:97.9; do N=${D%%:*}; P=${D##*:}
  bin/nbstune $C --ppd=$P --cache-save=$W/cache_$N.bin --out=/dev/null $W/msds_$N.csv; done
bin/nbsbench --metrics=ssim,msssim --out=$W/bench_ssim.csv $W/msds_tv.csv
python3 $S/msds_eval.py $W                   # clip-level table; LOCO refits via bin/nbscachesel
```

`msds_eval.py` evaluates SSIM, MS-SSIM, uniform-weight MB-SSIM and MB-SSIM
with the CID:IQ weights transferred unchanged, all at the clip level after
averaging over frames; this is Table 4. Passing `--refit` additionally runs
leave-one-content-out refits (joint across devices, and per device, using
`bin/nbscachesel` to slice the cached tensors by content); these are not
reported in the paper -- three contents do not identify 36 weights, and
the refits stay below the transferred weights -- but the code is kept so
the comparison can be repeated on larger sets.
`results/multiscreen/` holds the clip scores and the SSIM/MS-SSIM
frame scores.
