/* nbssim_dft.c -- DFT-domain evaluation of the N-band SSIM.
 *
 * Strict C89.
 *
 * WHY THIS EXISTS
 *
 * The spatial core convolves with truncated, Hamming-windowed FIR
 * approximations to the band transfer functions.  Three things go wrong
 * with that as the bank gets wide in octaves and low in frequency:
 *
 *   1. Truncation moves the realised peak and widens the realised
 *      bandwidth, so the bank has to be numerically re-calibrated and the
 *      answer still depends on an arbitrary support parameter.  Worse, a
 *      2.5-octave channel centred at 0.8 cpd truncated to +-2 sigma passes
 *      78%% of DC -- it stops being band-pass at all.
 *
 *   2. The longest kernel dictates a margin, and at ppd 63 that margin is
 *      106 px at every edge, so an 800x800 image is scored on 54%% of its
 *      area against SSIM's 98%%.
 *
 *   3. Every sub-band is carried at full resolution, so the 11x11
 *      statistics window spans 0.17 deg for all of them.  One cycle of the
 *      0.8 cpd channel is 1.25 deg.  The window cannot see a seventh of a
 *      cycle of the band whose structural similarity it reports.
 *
 * Sampling the transfer function on the DFT grid fixes (1) exactly and
 * (2) by construction.  Inverse-transforming each sub-band at its own
 * critical rate fixes (3): the window then covers a fixed number of cycles
 * of each channel, which is what MS-SSIM achieves by decimating between
 * scales.  It also cuts the n^2 statistics work by the square of the mean
 * decimation factor.
 *
 * STRUCTURE
 *
 * The mask H_i(nu_x) H_j(nu_y) is separable, so the transform is done as
 * rows-then-columns and the row stage is shared across j, mirroring the
 * spatial core's row-filtered cache:
 *
 *     mirror-extend to Wp x Hp (powers of two)
 *     FFT rows                                        -> X[v][u]
 *     for each row band i:
 *         multiply by H_i, fold u to Wd_i, inverse FFT rows   -> P_i
 *         FFT columns of P_i                                  -> C_i
 *         for each column band j:
 *             multiply by H_j, fold v to Hd_j, inverse FFT    -> S_ij
 *             windowed statistics and log xi on the Wd_i x Hd_j grid
 *             bilinear-resample onto the output grid and accumulate
 *
 * Folding a masked spectrum and inverse-transforming at the reduced length
 * is the exact DFT decimation identity, so the only approximation in the
 * decimated path is whatever the band's transfer function leaves above its
 * own critical rate -- which is why the log-Gabor sections are given
 * compact support by a raised-cosine taper (see nbs_bank_cfg.taper_mult).
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nbs.h"
#include "nbs_fft.h"
#include "nbs_metrics.h"

/* ------------------------------------------------------------------ */
/* Helpers                                                             */
/* ------------------------------------------------------------------ */

/* Reflect an out-of-range coordinate back inside [0,n), mirroring without
 * repeating the edge sample -- the same rule as NBS_BORDER_SYMMETRIC. */
static int refl(int s, int n)
{
    int p;
    if (n <= 1) return 0;
    p = 2 * n - 2;
    s = s % p;
    if (s < 0) s += p;
    return (s < n) ? s : p - s;
}

/* Cache-blocked transpose.  src is h rows of w; dst is w rows of h. */
static void transpose(const double *src, double *dst, int w, int h)
{
    const int T = 32;
    int i0, j0, i, j;
    for (i0 = 0; i0 < h; i0 += T) {
        int i1 = (i0 + T < h) ? i0 + T : h;
        for (j0 = 0; j0 < w; j0 += T) {
            int j1 = (j0 + T < w) ? j0 + T : w;
            for (i = i0; i < i1; ++i)
                for (j = j0; j < j1; ++j)
                    dst[(size_t)j * h + i] = src[(size_t)i * w + j];
        }
    }
}

/* Windowed second moments of a sub-band pair.
 *
 * The generic separable convolution allocates its own scratch on every
 * call, which for n^2 sub-bands means a hundred multi-megabyte first-touch
 * allocations per image pair.  This version takes its scratch from the
 * caller, folds the symmetric window so each tap is used twice, and forms
 * xx, yy and xy in one traversal so the sub-band samples are read once.
 * Border handling is symmetric reflection, matching the spatial core. */
static void win_stats(const double *sx, const double *sy, int w, int h,
                      const nbs_kernel *win,
                      double *exx, double *eyy, double *exy,
                      double *t1, double *t2, double *t3)
{
    int r = win->r, x, y, t;
    const double *kp = win->tap + r;

    for (y = 0; y < h; ++y) {
        const double *px = sx + (size_t)y * w, *py = sy + (size_t)y * w;
        double *a = t1 + (size_t)y * w, *b = t2 + (size_t)y * w, *c = t3 + (size_t)y * w;
        for (x = 0; x < w; ++x) {
            double s1 = kp[0] * px[x] * px[x];
            double s2 = kp[0] * py[x] * py[x];
            double s3 = kp[0] * px[x] * py[x];
            if (x >= r && x + r < w) {
                for (t = 1; t <= r; ++t) {
                    double kk = kp[t];
                    double xl = px[x - t], xr = px[x + t];
                    double yl = py[x - t], yr = py[x + t];
                    s1 += kk * (xl * xl + xr * xr);
                    s2 += kk * (yl * yl + yr * yr);
                    s3 += kk * (xl * yl + xr * yr);
                }
            } else {
                for (t = 1; t <= r; ++t) {
                    double kk = kp[t];
                    int il = refl(x - t, w), ir = refl(x + t, w);
                    double xl = px[il], xr = px[ir];
                    double yl = py[il], yr = py[ir];
                    s1 += kk * (xl * xl + xr * xr);
                    s2 += kk * (yl * yl + yr * yr);
                    s3 += kk * (xl * yl + xr * yr);
                }
            }
            a[x] = s1; b[x] = s2; c[x] = s3;
        }
    }

    for (y = 0; y < h; ++y) {
        double *oa = exx + (size_t)y * w, *ob = eyy + (size_t)y * w, *oc = exy + (size_t)y * w;
        const double *a0 = t1 + (size_t)y * w, *b0 = t2 + (size_t)y * w, *c0 = t3 + (size_t)y * w;
        for (x = 0; x < w; ++x) { oa[x] = kp[0] * a0[x]; ob[x] = kp[0] * b0[x]; oc[x] = kp[0] * c0[x]; }
        for (t = 1; t <= r; ++t) {
            double kk = kp[t];
            int yl = (y - t >= 0) ? y - t : refl(y - t, h);
            int yr = (y + t < h)  ? y + t : refl(y + t, h);
            const double *al = t1 + (size_t)yl * w, *ar = t1 + (size_t)yr * w;
            const double *bl = t2 + (size_t)yl * w, *br = t2 + (size_t)yr * w;
            const double *cl = t3 + (size_t)yl * w, *cr = t3 + (size_t)yr * w;
            for (x = 0; x < w; ++x) {
                oa[x] += kk * (al[x] + ar[x]);
                ob[x] += kk * (bl[x] + br[x]);
                oc[x] += kk * (cl[x] + cr[x]);
            }
        }
    }
}

/* Transform length for an image dimension: a power of two with room for a
 * guard band of about an eighth of the image on each side, so that the
 * periodic extension the DFT assumes is a mirror rather than a jump. */
static int pad_size(int n)
{
    int want = n + n / 4;
    int p = nbs_fft_size(want);
    if (p < n) return 0;
    return p;
}

/* Per-band decimated length: the smallest power of two that can carry the
 * band's support without aliasing, clamped to [min_len, full]. */
static int decim_len(const nbs_nbssim_cfg *cfg, int b, int full, int min_len)
{
    double eps, edge;
    int need, d;

    if (!cfg->decimate) return full;
    eps = (cfg->decim_eps > 0.0) ? cfg->decim_eps : 0.02;
    edge = nbs_band_edge(&cfg->bank, b, eps);      /* cycles/pixel */
    need = (int)ceil(2.0 * edge * (double)full);
    d = nbs_fft_size(need > 1 ? need : 1);
    if (d < min_len) d = min_len;
    if (d > full) d = full;
    return d;
}

int nbs_nbssim_dft_grid(const nbs_nbssim_cfg *cfg, int w, int h,
                        int *out_wp, int *out_hp, int *out_wd, int *out_hd)
{
    int Wp = pad_size(w), Hp = pad_size(h), b, n = cfg->bank.nbands;
    int minw, minh;

    if (Wp == 0 || Hp == 0) { nbs_set_error("nbssim/dft: image too large"); return 1; }
    minw = nbs_fft_size(cfg->win > 8 ? cfg->win * 2 : 16);
    minh = minw;
    if (out_wp) *out_wp = Wp;
    if (out_hp) *out_hp = Hp;
    for (b = 0; b < n; ++b) {
        if (out_wd) out_wd[b] = decim_len(cfg, b, Wp, minw);
        if (out_hd) out_hd[b] = decim_len(cfg, b, Hp, minh);
    }
    return 0;
}

/* Bilinear resampling tables from the output grid onto a decimated grid.
 * Output pixel o maps to padded coordinate o+pad, hence to decimated
 * coordinate (o+pad)/d.  Indices wrap, because the decimated sub-band is
 * periodic. */
static void resample_table(int nout, int pad, int full, int len,
                           int *i0, int *i1, double *fr)
{
    double d = (double)full / (double)len;
    int o;
    for (o = 0; o < nout; ++o) {
        double g = ((double)(o + pad)) / d;
        int    a = (int)floor(g);
        double t = g - (double)a;
        int    aa = a % len, bb;
        if (aa < 0) aa += len;
        bb = (aa + 1) % len;
        i0[o] = aa; i1[o] = bb; fr[o] = t;
    }
}

/* ------------------------------------------------------------------ */
/* Core                                                                */
/* ------------------------------------------------------------------ */

int nbs_nbssim_dft(const nbs_plane *xin, const nbs_plane *yin,
                   const nbs_nbssim_cfg *cfg, int mode, int stride,
                   double *score, double *cache, int *out_npix)
{
    const nbs_bank_cfg *bank = &cfg->bank;
    nbs_plane xs, ys;
    const nbs_plane *px, *py;
    nbs_fft_plan pw, ph, prow[NBS_MAX_BANDS], pcol[NBS_MAX_BANDS];
    /* One complex plane carries both images: the reference in the real
     * part, the distorted in the imaginary part.  Every stage from here to
     * the sub-bands is linear with real coefficients, so the two never mix
     * and the whole transform pipeline costs half what two real transforms
     * would.  Only the pointwise statistics need them apart. */
    double *ZR = NULL, *ZI = NULL;
    double *CZR = NULL, *CZI = NULL;
    double *SXb = NULL, *SYb = NULL, *E1 = NULL, *E2 = NULL, *E3 = NULL;
    double *TMP = NULL, *T2 = NULL, *T3 = NULL, *LOGB = NULL, *acc = NULL;
    double *A[NBS_MAX_BANDS], *B[NBS_MAX_BANDS];
    double *tr = NULL, *ti = NULL;
    int *ix0[NBS_MAX_BANDS], *ix1[NBS_MAX_BANDS];
    int *iy0[NBS_MAX_BANDS], *iy1[NBS_MAX_BANDS];
    double *fx[NBS_MAX_BANDS], *fy[NBS_MAX_BANDS];
    int wd[NBS_MAX_BANDS], hd[NBS_MAX_BANDS];
    nbs_kernel win;
    int n, i, j, k, u, v, b, rc = 1;
    int W, H, Wp = 0, Hp = 0, padl = 0, padt = 0, m0, f;
    int nw, nsamp = 0, ox, oy;
    size_t big;

    xs.d = NULL; ys.d = NULL; win.tap = NULL;
    pw.n = 0; pw.cs = pw.sn = NULL; pw.rev = NULL;
    ph.n = 0; ph.cs = ph.sn = NULL; ph.rev = NULL;
    for (k = 0; k < NBS_MAX_BANDS; ++k) {
        A[k] = B[k] = NULL; fx[k] = fy[k] = NULL;
        ix0[k] = ix1[k] = iy0[k] = iy1[k] = NULL;
        prow[k].n = 0; prow[k].cs = prow[k].sn = NULL; prow[k].rev = NULL;
        pcol[k].n = 0; pcol[k].cs = pcol[k].sn = NULL; pcol[k].rev = NULL;
    }

    n  = bank->nbands;
    nw = nbs_w_count_ex(n, cfg->w_sym);

    if (xin->w != yin->w || xin->h != yin->h) {
        nbs_set_error("nbssim/dft: size mismatch"); return 1;
    }
    if (n < 2 || n > NBS_MAX_BANDS) {
        nbs_set_error("nbssim/dft: bad band count"); return 1;
    }
    if (stride < 1) stride = 1;

    f = cfg->downsample;
    if (f < 0) {
        int mn = (xin->w < xin->h) ? xin->w : xin->h;
        f = (int)floor((double)mn / 256.0 + 0.5);
        if (f < 1) f = 1;
    }
    if (f > 1) {
        if (nbs_downsample(xin, &xs, f, NBS_BORDER_SYMMETRIC) != 0) goto done;
        if (nbs_downsample(yin, &ys, f, NBS_BORDER_SYMMETRIC) != 0) goto done;
        px = &xs; py = &ys;
    } else {
        px = xin; py = yin;
    }
    W = px->w; H = px->h;

    /* The DFT core has no filter support to fit inside the image, so the
     * automatic margin is zero: unlike the spatial core it scores the whole
     * frame.  An explicit margin is still honoured, which is what lets the
     * two paths be compared over the same pixels. */
    m0 = (cfg->margin > 0) ? cfg->margin : 0;
    if (W - 2 * m0 <= 0 || H - 2 * m0 <= 0) {
        nbs_set_error("nbssim/dft: margin larger than the image");
        goto done;
    }

    if (nbs_nbssim_dft_grid(cfg, W, H, &Wp, &Hp, wd, hd) != 0) goto done;
    padl = (Wp - W) / 2;
    padt = (Hp - H) / 2;

    if (nbs_fft_plan_init(&pw, Wp) != 0) { nbs_set_error("nbssim/dft: plan"); goto done; }
    if (nbs_fft_plan_init(&ph, Hp) != 0) { nbs_set_error("nbssim/dft: plan"); goto done; }
    for (b = 0; b < n; ++b) {
        if (nbs_fft_plan_init(&prow[b], wd[b]) != 0) { nbs_set_error("nbssim/dft: plan"); goto done; }
        if (nbs_fft_plan_init(&pcol[b], hd[b]) != 0) { nbs_set_error("nbssim/dft: plan"); goto done; }
    }

    big = (size_t)Wp * (size_t)Hp;
    ZR = (double *)nbs_malloc(big * sizeof(double));
    ZI = (double *)nbs_malloc(big * sizeof(double));
    CZR = (double *)nbs_malloc(big * sizeof(double));
    CZI = (double *)nbs_malloc(big * sizeof(double));
    SXb = (double *)nbs_malloc(big * sizeof(double));
    SYb = (double *)nbs_malloc(big * sizeof(double));
    E1  = (double *)nbs_malloc(big * sizeof(double));
    E2  = (double *)nbs_malloc(big * sizeof(double));
    E3  = (double *)nbs_malloc(big * sizeof(double));
    TMP = (double *)nbs_malloc(big * sizeof(double));
    T2  = (double *)nbs_malloc(big * sizeof(double));
    T3  = (double *)nbs_malloc(big * sizeof(double));
    LOGB = (double *)nbs_malloc(big * sizeof(double));
    tr = (double *)nbs_malloc((size_t)(Wp > Hp ? Wp : Hp) * sizeof(double));
    ti = (double *)nbs_malloc((size_t)(Wp > Hp ? Wp : Hp) * sizeof(double));
    if (ZR == NULL || ZI == NULL || CZR == NULL || CZI == NULL ||
        SXb == NULL || SYb == NULL || E1 == NULL || E2 == NULL ||
        E3 == NULL || TMP == NULL || T2 == NULL || T3 == NULL || LOGB == NULL || tr == NULL || ti == NULL)
        goto done;

    /* Sampled transfer functions, one array per band per axis, indexed by
     * DFT bin.  Normalising by the analytic peak makes "unit peak gain"
     * exact rather than a numerical fit. */
    for (b = 0; b < n; ++b) {
        double gp = bank->norm_peak ? nbs_band_peak(bank, b, NULL) : 1.0;
        if (!(gp > 1e-12)) gp = 1.0;
        A[b] = (double *)nbs_malloc((size_t)Wp * sizeof(double));
        B[b] = (double *)nbs_malloc((size_t)Hp * sizeof(double));
        if (A[b] == NULL || B[b] == NULL) goto done;
        for (u = 0; u < Wp; ++u) {
            int uu = (u <= Wp / 2) ? u : Wp - u;
            A[b][u] = nbs_band_H(bank, b, (double)uu / (double)Wp) / gp;
        }
        for (v = 0; v < Hp; ++v) {
            int vv = (v <= Hp / 2) ? v : Hp - v;
            B[b][v] = nbs_band_H(bank, b, (double)vv / (double)Hp) / gp;
        }
    }

    /* Resampling tables from the output grid onto each band's grid. */
    for (b = 0; b < n; ++b) {
        ix0[b] = (int *)nbs_malloc((size_t)W * sizeof(int));
        ix1[b] = (int *)nbs_malloc((size_t)W * sizeof(int));
        fx[b]  = (double *)nbs_malloc((size_t)W * sizeof(double));
        iy0[b] = (int *)nbs_malloc((size_t)H * sizeof(int));
        iy1[b] = (int *)nbs_malloc((size_t)H * sizeof(int));
        fy[b]  = (double *)nbs_malloc((size_t)H * sizeof(double));
        if (ix0[b] == NULL || ix1[b] == NULL || fx[b] == NULL ||
            iy0[b] == NULL || iy1[b] == NULL || fy[b] == NULL) goto done;
        resample_table(W, padl, Wp, wd[b], ix0[b], ix1[b], fx[b]);
        resample_table(H, padt, Hp, hd[b], iy0[b], iy1[b], fy[b]);
    }

    if (nbs_ssim_window(&win, cfg->win, cfg->win_sigma) != 0) goto done;

    /* Output accumulator / sample count. */
    if (mode == 0) {
        acc = (double *)nbs_calloc(big, sizeof(double));
        if (acc == NULL) goto done;
    } else {
        nsamp = 0;
        for (oy = m0; oy < H - m0; oy += stride)
            for (ox = m0; ox < W - m0; ox += stride)
                nsamp++;
        if (out_npix) *out_npix = nsamp;
        if (cache == NULL) { rc = 0; goto done; }
        memset(cache, 0, (size_t)nsamp * (size_t)nw * sizeof(double));
    }

    /* ---- mirror-extend and transform rows ---- */
    for (v = 0; v < Hp; ++v) {
        int sy = refl(v - padt, H);
        double *zr = ZR + (size_t)v * Wp, *zi = ZI + (size_t)v * Wp;
        for (u = 0; u < Wp; ++u) {
            int sx = refl(u - padl, W);
            zr[u] = NBS_AT(px, sx, sy);
            zi[u] = NBS_AT(py, sx, sy);
        }
        nbs_fft(&pw, zr, zi, 0);
    }

    /* ---- band loop ---- */
    for (i = 0; i < n; ++i) {
        int Wd = wd[i];
        int du = Wp / Wd;
        double su = 1.0 / (double)du;

        /* Rows: mask, fold to Wd, inverse transform, then transform the
         * columns of the row-decimated plane. */
        for (v = 0; v < Hp; ++v) {
            const double *zr = ZR + (size_t)v * Wp, *zi = ZI + (size_t)v * Wp;
            double *cr = CZR + (size_t)v * Wd, *ci = CZI + (size_t)v * Wd;
            int a;
            for (k = 0; k < Wd; ++k) { cr[k] = 0.0; ci[k] = 0.0; }
            for (a = 0; a < du; ++a) {
                const double *g = A[i] + (size_t)a * Wd;
                const double *sr = zr + (size_t)a * Wd, *si = zi + (size_t)a * Wd;
                for (k = 0; k < Wd; ++k) { cr[k] += g[k] * sr[k]; ci[k] += g[k] * si[k]; }
            }
            for (k = 0; k < Wd; ++k) { cr[k] *= su; ci[k] *= su; }
            nbs_fft(&prow[i], cr, ci, 1);
        }
        /* CXR/CYR now hold the row-decimated planes, Hp rows of Wd real
         * samples.  Everything from here on runs down columns, so transpose
         * once and keep the transposed layout: a strided FFT over a 1024-row
         * column walks 8 KB per element and spends its whole time in cache
         * misses, which measured slower than the spatial core it replaced. */
        transpose(CZR, SXb, Wd, Hp);   /* SXb, SYb double as scratch here */
        transpose(CZI, SYb, Wd, Hp);
        memcpy(CZR, SXb, (size_t)Wd * (size_t)Hp * sizeof(double));
        memcpy(CZI, SYb, (size_t)Wd * (size_t)Hp * sizeof(double));
        for (k = 0; k < Wd; ++k)
            nbs_fft(&ph, CZR + (size_t)k * Hp, CZI + (size_t)k * Hp, 0);

        for (j = 0; j < n; ++j) {
            int Hd = hd[j];
            int dv = Hp / Hd;
            double sv = 1.0 / (double)dv;
            int wi = nbs_w_index_ex(i, j, n, cfg->w_sym);
            double kk = cfg->kc[(i > j) ? i : j];
            double C;
            size_t sz = (size_t)Wd * (size_t)Hd;

            C = kk * cfg->L; C *= C;

            /* Columns (the fast axis in the transposed layout): mask, fold
             * to Hd, inverse transform.  Sub-bands come out as Wd rows of
             * Hd samples, i.e. transposed. */
            for (k = 0; k < Wd; ++k) {
                const double *czr = CZR + (size_t)k * Hp, *czi = CZI + (size_t)k * Hp;
                double *sx = SXb + (size_t)k * Hd, *sy = SYb + (size_t)k * Hd;
                int q, a;
                for (q = 0; q < Hd; ++q) { tr[q] = 0.0; ti[q] = 0.0; }
                for (a = 0; a < dv; ++a) {
                    const double *br = czr + (size_t)a * Hd, *bi = czi + (size_t)a * Hd;
                    const double *g = B[j] + (size_t)a * Hd;
                    for (q = 0; q < Hd; ++q) { tr[q] += g[q] * br[q]; ti[q] += g[q] * bi[q]; }
                }
                for (q = 0; q < Hd; ++q) { tr[q] *= sv; ti[q] *= sv; }
                nbs_fft(&pcol[j], tr, ti, 1);
                /* real part -> reference sub-band, imaginary -> distorted */
                for (q = 0; q < Hd; ++q) { sx[q] = tr[q]; sy[q] = ti[q]; }
            }

            /* Windowed second moments on the sub-band's own grid.  The
             * window is separable and isotropic, so it does not care that
             * the plane is transposed. */
            win_stats(SXb, SYb, Hd, Wd, &win, E1, E2, E3, TMP, T2, T3);

            for (k = 0; (size_t)k < sz; ++k) {
                double xi_ = nbs_xi_ratio(E3[k], E1[k], E2[k], C, cfg->norm_p);
                if (xi_ < cfg->xi_floor) xi_ = cfg->xi_floor;
                if (xi_ > 1.0) xi_ = 1.0;
                TMP[k] = log(xi_);
            }
            transpose(TMP, LOGB, Hd, Wd);   /* back to Hd rows of Wd */

            /* Resample onto the output grid and accumulate. */
            {
                const int *px0 = ix0[i], *px1 = ix1[i];
                const int *py0 = iy0[j], *py1 = iy1[j];
                const double *pfx = fx[i], *pfy = fy[j];
                double wgt = (mode == 0) ? cfg->w[wi] : 1.0;
                int p = 0;
                int y_lo = m0, y_hi = H - m0, x_lo = m0, x_hi = W - m0;
                int st = (mode == 0) ? 1 : stride;

                for (oy = y_lo; oy < y_hi; oy += st) {
                    const double *r0 = LOGB + (size_t)py0[oy] * Wd;
                    const double *r1 = LOGB + (size_t)py1[oy] * Wd;
                    double ty = pfy[oy], sy1 = 1.0 - ty;
                    for (ox = x_lo; ox < x_hi; ox += st) {
                        int a0 = px0[ox], a1 = px1[ox];
                        double tx = pfx[ox], sx1 = 1.0 - tx;
                        double val = sy1 * (sx1 * r0[a0] + tx * r0[a1])
                                   + ty  * (sx1 * r1[a0] + tx * r1[a1]);
                        if (mode == 0) acc[(size_t)oy * W + ox] += wgt * val;
                        else { cache[(size_t)p * (size_t)nw + (size_t)wi] += val; p++; }
                    }
                }
            }
        }
    }

    if (mode == 0) {
        double s = 0.0;
        long cnt = 0;
        for (oy = m0; oy < H - m0; ++oy)
            for (ox = m0; ox < W - m0; ++ox) { s += exp(acc[(size_t)oy * W + ox]); cnt++; }
        *score = (cnt > 0) ? s / (double)cnt : 0.0;
    }
    rc = 0;

done:
    for (k = 0; k < NBS_MAX_BANDS; ++k) {
        free(A[k]); free(B[k]);
        free(ix0[k]); free(ix1[k]); free(fx[k]);
        free(iy0[k]); free(iy1[k]); free(fy[k]);
        nbs_fft_plan_free(&prow[k]);
        nbs_fft_plan_free(&pcol[k]);
    }
    nbs_fft_plan_free(&pw); nbs_fft_plan_free(&ph);
    free(ZR); free(ZI); free(CZR); free(CZI);
    free(SXb); free(SYb); free(E1); free(E2); free(E3);
    free(TMP); free(T2); free(T3); free(LOGB); free(acc); free(tr); free(ti);
    nbs_kernel_free(&win);
    nbs_plane_free(&xs); nbs_plane_free(&ys);
    return rc;
}
