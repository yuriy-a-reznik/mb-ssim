# MB-SSIM: a multi-band variant of SSIM, tunable to viewing conditions

[![licence: MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](LICENSE)

A C89 implementation of MB-SSIM, a full-reference image quality metric that
keeps SSIM's per-pixel product of similarity ratios but replaces its implicit
two-band decomposition with a psychophysically placed filterbank specified in
**cycles per degree** — so the viewing geometry is an input to the metric
rather than a correction applied afterwards.

> Y. Reznik, *"MB-SSIM: A Multi-Band Variant of SSIM, Tunable to
> Viewing Conditions"*, IS&T Electronic Imaging 2027, Image Quality and
> System Performance. **(this repository:
> <https://github.com/yuriy-a-reznik/mb-ssim>)**

The tools call the metric `nbssim` (N-band SSIM) throughout, which is what
it was called while the design was being worked out; `MB-SSIM` is the name
it goes by in the paper. They are the same metric. Likewise, what this
document calls the per-sub-band *exponents* `w_ij` the paper calls the
**sub-band weights** — again the same quantities, under the name that reads
better in prose.

Building on:

> Y. Reznik, *"Another look at SSIM image quality metric"*, IS&T Electronic
> Imaging 2023, Image Quality and System Performance XX, pp. 305-1–305-7.
>
> Y. Reznik, *"SSIM Revisited: From Structural Similarity to Error
> Visibility"*, ACM MHV'24, Denver, CO, pp. 129–130.

together with faithful reference implementations of SSIM and MS-SSIM, a
testbench that follows the evaluation methodology of the original SSIM and
MS-SSIM papers, and a tuner for the metric's free parameters.

Everything is strict C89 with no dependencies beyond the standard library —
including the PNG decoder, which carries its own DEFLATE implementation.
The Python scripts are only used to prepare dataset index files and to
cross-check the C code; nothing in the metric path needs them.

---

## Quick start

```sh
make && make check                 # strict C89, no dependencies; 333 assertions
bin/nbsmetric --preset=wilson --transform=dft --taper=1.0 --ppd=63.223 ref.png dist.png
```

The only viewing-dependent input is `--ppd`, the angular sampling density in
pixels per degree of visual angle, `ppd = 1 / (2 atan(1 / (2 d rho)))` for
viewing distance `d` and display pixel density `rho` (eq. 7 of the paper).
Typical values: a 1080p television at three picture heights is 57, a 4K
television at the same distance 113, a phone at 13 inches about 104.

Two further options matter for the paper's last two ablations:
`--kc=a,b,...` sets the per-channel stabilising constants (the paper's
Table 8 scales SSIM's `0.01,0.03`; `--kc=0.0033,0.01` is the level that
scored best on CID:IQ), and `--xi-norm=P` changes the normaliser of the
similarity ratio from the arithmetic mean of the two band energies (`1`,
SSIM's choice and the default) to the power mean of order `P`
(`geo` for the geometric mean, i.e. a correlation coefficient).

## Where each number in the paper comes from

| paper | where |
|---|---|
| Table 1 (viewing setups) | `docs/DATASETS.md` (the geometry of each database) |
| Table 2 (CID:IQ), Figures 5, 6 | `REPRODUCE.md` §2–5; `results/cidiq/` |
| Table 3 (LIVE Release 2) | `REPRODUCE.md` §8; `results/live_r2/` |
| Table 4 (multi-screen cross-verification) | `REPRODUCE.md` §10; `results/multiscreen/` |
| Table 5 (channel-count sweep), Figure 7 | `REPRODUCE.md` §7; `results/sweep/` |
| Table 6 (weight parameterisations, CSF rows) | `REPRODUCE.md` §7, §7c; `results/cidiq/csf/` |
| Table 7 (bank perturbations) | `REPRODUCE.md` §7b; `runs/ablf_*` |
| Table 8 (stabilising constants) | `REPRODUCE.md` §7e; `results/cidiq/constants/` |
| Figure 8 (fitted weights) | `REPRODUCE.md` §5, §7 |

## What is in the repository

| path | contents |
|---|---|
| `src/`, `include/` | the metric, the SSIM/MS-SSIM references, DFT and spatial cores, log-xi cache, logistic fit and correlation measures |
| `tools/` | `nbsmetric` (one pair), `nbsbench` (a database), `nbstune` and `nbsjoint` (fit the sub-band weights, single and multi-condition), `nbsbank` (dump the realised filterbank) |
| `test/` | `nbstest`, run by `make check` |
| `scripts/` | dataset index builders, summaries, the EI-2023 reference check, `ablation_sweep.py`, the contrast-sensitivity-weights experiment (`csf_weights.py`, `csf_lowest.py`, `csf_refit.py`), `norm/`, `constants/` and `multiscreen/` |
| `runs/` | the exact command lines and logs behind every number in the paper |
| `results/` | fitted weights and scores for CID:IQ (including the CSF-weight rows of Table 7), LIVE Release 2, the channel-count sweep and the multi-screen cross-verification |
| `docs/` | how to obtain the subjective databases and the viewing-geometry reasoning behind each `--ppd` |
| `REPRODUCE.md` | step-by-step reproduction of every table and figure, with timings |

The raw subjective databases (CID:IQ, LIVE) are not redistributed; see
`docs/DATASETS.md`.

## Reproducing and extending the paper

`REPRODUCE.md` walks through every experiment. The two headline commands
are the joint fit across both CID:IQ viewing distances (Table 3 of the
paper) and the perturbation study of the filterbank (Table 7):

```sh
# one set of weights for both viewing distances, held-out five-fold CV
bin/nbsjoint --label=50cm --cond=runs/cache_50cm.bin --label=100cm --cond=runs/cache_100cm.bin \
             --combine=mean --objective=plcc --folds=5 --restarts=2 --basin=6 --maxit=3000
# 56 perturbations of the bank: peak placement, bandwidth, shape, window; sharded and resumable
python3 scripts/ablation_sweep.py list
python3 scripts/ablation_sweep.py run --budget=fast --jobs=2 --shard=1/4
python3 scripts/ablation_sweep.py summarise --budget=fast
```

Every design parameter of the bank -- `--bands`, `--peaks`, `--bw`,
`--bank={loggabor,dog,cumgauss}`, `--taper`, `--win`, `--win-sigma`,
`--kc` -- is a command-line argument, so speed/accuracy trade-offs can be
explored rather than fixed at design time.

---

## What the metric is

The 2023 paper observes that SSIM is, to a good approximation, the product of
one and the same distance operator applied to a low-pass and a high-pass
version of the two images:

```
SSIM(x, y)  ≈  ξ(x_L, y_L) · ξ(x_H, y_H)

              2·E[a·b] + C
ξ(a, b)  =  ─────────────────────
            E[a²] + E[b²] + C
```

with `E[·]` a Gaussian-weighted local expectation (11×11, σ = 1.5) and the
band split produced by a Gaussian of σ = 3 px. The MHV'24 abstract proposes
generalising this to N bands arranged like the spatial channels of the
contrast sensitivity function, combined by a weighted product.

This implementation realises that as a **separable two-dimensional bank**:
N one-dimensional bands per axis give N×N two-dimensional subbands, and

```
MB-SSIM(x, y)  =  mean_p  ∏_{i,j} ξ_ij(p) ^ w_ij
```

The subband scores are combined *per spatial location* and only then averaged
over the image, as the abstract specifies. The exponents are constrained to be
symmetric under transposition of the two axes, `w_ij == w_ji`, which reflects
the fact that there is no reason to treat horizontal and vertical detail
differently — for the default N = 6 that leaves **21 free exponents**, not 36.

### The filterbank

Band frequencies are specified in **cycles per degree of visual angle**, which
means the metric needs to know the viewing geometry: `--ppd` (pixels per
degree). A Gaussian of pixel standard deviation σ has corner frequency
`f = ppd / (2πσ)` — the convention used in the 2023 paper, where σ = 3 px maps
to ≈ 4.24 cpd at 80 cpd sampling.

Three bank constructions are available:

| `--bank=` | construction | notes |
|---|---|---|
| `dog` *(default)* | each band-pass is a difference of two Gaussians whose corners straddle the requested peak | peaks land **exactly** where asked; `--rho` sets the bandwidth. For N = 2 it degenerates to exactly the two-band model of the 2023 paper. |
| `cumgauss` | nested Gaussian low-passes, bands taken as successive differences | sums to unity exactly, but band-pass peaks cannot be spaced closer than about one octave |
| `loggabor` | Gaussian magnitude on a log-frequency axis, by frequency sampling | the literal reading of slide 17 of the MHV'24 deck |

**Why the default is `dog`, not `cumgauss`.** The nested-Gaussian construction
is the most elegant option — it is a partition of unity, so the bands sum back
to the original image exactly, and it is the direct generalisation of the
two-band model. But it has a hard structural limit: the difference of two
nested Gaussians with corners `f₁ < f₂` peaks at

```
f*² = 4·ln(f₂/f₁) / (1/f₁² − 1/f₂²)      →      f* > √2 · f₁
```

so each band-pass peak must exceed √2 times the *previous cut-off*, and the
cut-off itself sits well above the previous peak. The requested channel set
(0.8, 1.7, 2.8, 4, 8 cpd — spacings of 1.09, 0.72, 0.52 and 1.0 octaves) is
**not realisable** that way: after placing the 1.7 cpd peak the next reachable
peak is already above 3.4 cpd. Independent difference-of-Gaussians bands have
no such constraint, at the cost of the partition-of-unity property, which the
weighted-product combination does not actually require. `cumgauss` remains
available and is unit-tested for the partition property.

Band-pass corners are **calibrated against the discrete response**, not the
continuous formula: for σ below about one pixel, sampling aliases the Gaussian
enough to shift a DoG peak by 15 % or more. `nbsmetric --bands-report` prints
the realised peak of every band so this is visible rather than assumed.

### The stabilising constants

`C_ij = (kc[max(i,j)] · L)²`, indexed by the *coarser* of the two 1-D bands, so
the DC-carrying band keeps a `C1`-like constant (default k = 0.01) and the
detail bands a `C2`-like one (default k = 0.03). Tunable via `--kc`.

---

## Build

```sh
make            # bin/nbsmetric, nbsbench, nbstune, nbsjoint, nbsbank, nbstest
make check      # runs the unit tests
```

`-std=c89 -Wall -Wextra -pedantic`, warning-free. Tested with gcc.

---

## Tools

### `nbsmetric` — one image pair

```sh
bin/nbsmetric --preset=wilson --transform=dft --taper=1.0 --ppd=63.223 \
              --bands-report ref.bmp distorted.bmp
```

```
image   : 800 x 800
SSIM    : 0.98583497
MS-SSIM : 0.97647277
MB-SSIM : 0.92740567
```

`bin/nbsbank --preset=wilson --ppd=63.223 --taper=1.0` reports the realised
bank; at ppd 31.612 it reports five channels and names the one it dropped.

Accepts PNG (own decoder, 8/16-bit, colour types 0/2/3/4/6, non-interlaced),
BMP (8/24/32-bit uncompressed) and PNM. Luma is extracted with BT.601
coefficients by default (matching MATLAB's `rgb2gray`, which is what the SSIM
reference code uses); `--luma=709` and `--luma=mean` are available.

### `nbsbench` — a whole database

```sh
bin/nbsbench --out=scores.csv --progress data/live.csv
```

Reports PLCC, SRCC, KRCC, RMSE and MAE. PLCC/RMSE/MAE are computed **after**
the VQEG five-parameter logistic fit

```
f(x) = β₁·(0.5 − 1/(1+exp(β₂(x−β₃)))) + β₄·x + β₅
```

exactly as in the SSIM and MS-SSIM papers, so a metric is not penalised for a
monotone nonlinearity. SRCC/KRCC signs are flipped to be positive-when-agreeing
regardless of whether the database uses MOS or DMOS.

`--shard=I/N` splits the work across processes; `--report=scores.csv` merges
and re-reports without recomputation.

### `nbstune` — fit the exponents

```sh
bin/nbstune --stride=12 --folds=5 --out=weights.txt --progress data/live.csv
```

The metric is `mean_p exp(Σ_k w_k·L_k(p))` where `L_k = log ξ_k` depends only
on the images, so for a fixed filterbank the `L_k` maps are computed **once**
and cached on a subsampled pixel grid. Every subsequent exponent evaluation is
then a handful of dot products, which is what makes Nelder-Mead over 21
parameters tractable. Changing a band frequency invalidates the cache.

Cross-validation folds are split **by reference image**, so no source content
is shared between train and test. The held-out numbers are the ones to quote;
the all-data fit is reported alongside precisely so the gap is visible.

### `nbsbank` — dump the realised responses

```sh
bin/nbsbank --ppd=30 --points=400 > bank.csv
```

Emits the exact DTFT of every kernel the metric uses, plus their sum. Useful
for reproducing the figure on slide 17 and for confirming a bank is doing what
you think.

### `nbsjoint` — fit one exponent vector across several viewing conditions

```sh
bin/nbstune --preset=wilson --transform=dft --taper=1.0 --stride=16 \
            --ppd=31.612 --cache-save=cache_50cm.bin  data/cidiq_50cm.csv
bin/nbstune --preset=wilson --transform=dft --taper=1.0 --stride=16 \
            --ppd=63.223 --cache-save=cache_100cm.bin data/cidiq_100cm.csv

bin/nbsjoint --label=50cm  --cond=cache_50cm.bin \
             --label=100cm --cond=cache_100cm.bin \
             --combine=pooled --folds=5 --basin=4 --out=w.txt
```

Fitting each viewing distance separately cannot answer the question the metric
is built to answer: any per-condition improvement could be the filterbank
tracking the geometry, or the fit absorbing the difference. `nbsjoint` forces
**one** exponent vector to serve all conditions at once.

Conditions may have different band counts — at 50 cm the 16 c/deg channel is
above the display Nyquist rate and is dropped, leaving five channels and 15
exponents against six and 21 at 100 cm. They are not independent parameter
sets: the surviving channels are the same channels at the same cycles per
degree, so each condition indexes into a master matrix by `(i,j)` and a
dropped channel simply contributes no factor to the product.

`--combine=mean` averages the per-condition accuracies (each with its own
logistic fit); `--combine=pooled` runs one logistic fit over all conditions,
which additionally demands that the same MB-SSIM value mean the same quality
at every distance; `--combine=min` optimises the worst condition.

---

## The Wilson filterbank preset

`--preset=wilson` configures six band-pass channels peaking at
0.8, 1.7, 2.8, 4, 8 and 16 c/deg with half-amplitude bandwidths of
2.5, 2.2, 1.9, 1.7, 1.5 and 1.25 octaves, after Wilson, McFarlane & Phillips
(1983) and Wilson & Gelb (1984). Channels are log-Gabor: a difference of
Gaussians cannot go below about 1.77 octaves, so it cannot reach the narrow
high-frequency mechanisms at all.

`--taper=F` gives each channel compact support in frequency by rolling it off
to zero at `F × bandwidth` octaves from the centre. Without it a 2.5-octave
log-Gabor at 0.8 c/deg is still at 2% of peak by 6 c/deg and cannot be
critically sampled. The taper leaves the realised peaks and half-amplitude
bandwidths unchanged to three decimal places.

Channels at or above the display Nyquist rate `ppd/2` are dropped, and the
tools say so. This is the model asserting that a distortion the display
cannot present cannot be seen.

---

## Two cores: `--transform=dft` and `--transform=spatial`

The spatial core convolves with truncated, Hamming-windowed FIR
approximations to the channel transfer functions. It is the obvious
implementation and it has three problems, all of which change the metric
rather than merely approximating it:

0. **The patch window inside ξ is not optional.** `--nbs-win=1` evaluates ξ
   on the sub-band samples themselves, with no local expectation. It is a
   reasonable thing to try — the filterbank has already band-limited each
   channel — and it costs 0.040 PLCC pooled on CID:IQ once the exponents
   are refitted. Band-limiting a channel and estimating a local second
   moment of it are different operations. The default (11×11 Gaussian,
   σ = 1.5, as in SSIM) is the right one.

1. **Truncation changes the channel.** A 2.5-octave channel at 0.8 c/deg with
   ppd 63 needs 203 taps. Shortening it lets DC through — 21% at ±4σ, 77% at
   ±2σ — so the lowest "band-pass" channel silently becomes part low-pass, by
   an amount set by an arbitrary support parameter.
2. **The margin discards much of the image.** The 203-tap kernel plus the
   statistics window forces a 106-pixel margin at every edge, so an 800×800
   image is scored on 54% of its area against SSIM's 98%.
3. **The statistics window is the wrong size for the low channels.** An 11×11
   window at ppd 63 spans 0.17°; one cycle of the 0.8 c/deg channel is 1.25°.
   The window sees a seventh of a cycle of the band it is measuring.

`--transform=dft` samples the transfer functions directly on the DFT grid.
Peaks, gains and bandwidths come out exact, no calibration is needed, the
answer does not depend on `--radius-mult` at all (there is a unit test for
that), and there is no margin. `--band-decimate=auto` then inverse-transforms
each sub-band at its own critical rate, so the 11×11 window covers a fixed
number of *cycles* of each channel rather than a fixed number of pixels.

It is also 5.6× faster: 1.33 s against 7.39 s for an 800×800 pair at ppd 63.2
on one core of a 2.8 GHz Xeon, while scoring 1.9× more pixels. Both images
are carried in one complex plane (reference in the real part, distorted in the
imaginary part), which halves every transform stage, and the decimated
sub-bands cut the N² statistics passes by the square of the mean decimation
factor.

---

## Datasets

See `docs/DATASETS.md` for download instructions, score-file semantics and the
viewing-geometry reasoning behind each default `--ppd`.

```sh
python3 scripts/make_index.py live    /path/to/databaserelease2 -o data/live.csv
python3 scripts/make_index.py tid2013 /path/to/tid2013          -o data/tid2013.csv
python3 scripts/make_index.py csiq    /path/to/CSIQ             -o data/csiq.csv
```

The index format is a comma-separated list of
`reference,distorted,score,reference_id` with a small `#`-prefixed header
carrying the score orientation and the viewing geometry. It is documented in
`include/nbs_data.h` and is trivial to generate for any other database.

---

## Verification

`make check` runs the unit tests: kernel normalisation and symmetry, DC
preservation under convolution, exact band-peak placement, the
partition-of-unity property of the `cumgauss` bank, ppd-invariance of the
cut-off solution, `metric(x,x) == 1` for all three metrics, symmetry and
monotonicity of MB-SSIM, PNG/BMP/PNM decoding, correlation measures against
closed-form cases, logistic-fit recovery of a monotone nonlinearity, and
consistency between the direct metric path and the cached log-ξ path used by
the tuner.

Two checks are worth calling out:

**SSIM and MS-SSIM against an independent reference.** `scripts/ref_check.py`
is a numpy transcription of Wang's `ssim_index.m` and `ms_ssim.m`, written
against the papers rather than against this code. Agreement over noise, blur
and three image sizes, including the automatic-downsampling path:

```
max |C − numpy|  SSIM: 4.5e-13   MS-SSIM: 4.5e-13
```

**The 2023 paper's central claim.** The test suite computes the isotropic
two-band model ξ(x_L,y_L)·ξ(x_H,y_H) with σ = 3 directly and compares it to
reference SSIM across noise and blur distortions: worst |Δ| = 0.015,
rms = 0.008 — the same order as Tables 1–4 of the EI 2023 paper.

---

## Status and caveats

- The exponents shipped as defaults are **uniform** (`w_ij = 2/N²`, total
  exponent 2, matching SSIM's two unit-exponent factors). They are a starting
  point, not a result. Run `nbstune` or `nbsjoint` against a real subjective
  database before quoting anything. The vector fitted on CID:IQ for the 2027
  paper is in `results/w_cidiq_joint.txt`.
- The fitted exponents are **database- and viewing-condition-specific**. They
  do not transfer between databases in our experience, and the failure is
  directional rather than a matter of scale.
- The bank is radial only. Wilson's mechanisms are orientation-selective; the
  separable N×N construction here gives a crude four-way orientation sampling,
  not oriented channels.
- The DFT core evaluates the whole frame by default (`--nbs-margin=auto` means
  zero margin there, since there is no filter support to fit inside the
  image). The spatial core still crops `longest filter radius + window radius`
  from each edge. Pass an explicit `--nbs-margin=N` to compare the two cores
  over identical pixels — there is a unit test that does exactly this.
- ξ is clamped to `[--xi-floor, 1]` before exponentiation. Band-pass ξ can go
  negative when E[ab] < 0, which has no real-valued fractional power.
- An all-band-pass bank carries no mean luminance, so a pure luminance shift
  is invisible to it. The fitted weight on the 0.8×0.8 c/deg sub-band appears
  to be doing the work SSIM's luminance term does; a bank with an explicit
  low-pass channel has not been evaluated.

## Citing

If you use this code, please cite the paper and the two it builds on:

```bibtex
@inproceedings{reznik2027mbssim,
  author    = {Yuriy Reznik},
  title     = {{MB-SSIM}: A Multi-Band Variant of {SSIM}, Tunable to Viewing Conditions},
  booktitle = {Proc. IS\&T Int'l. Symp. on Electronic Imaging: Image Quality and System Performance},
  year      = {2027}
}
```

See `CITATION.cff` for a machine-readable version.

## Licence

MIT -- see `LICENSE`. The reference SSIM/MS-SSIM algorithms are due to
Z. Wang et al.; please cite the original papers. Subjective databases
carry their own licence terms -- see `docs/DATASETS.md`.
