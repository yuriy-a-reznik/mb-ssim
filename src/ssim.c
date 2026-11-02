/* ssim.c -- reference implementation of the SSIM index.
 *
 * Faithful port of Z. Wang's ssim_index.m / ssim.m:
 *   Z. Wang, A. C. Bovik, H. R. Sheikh and E. P. Simoncelli,
 *   "Image quality assessment: from error visibility to structural
 *    similarity", IEEE Trans. Image Processing, 13(4):600-612, 2004.
 *
 * Windowing is an 11x11 circular-symmetric Gaussian with sigma = 1.5,
 * normalised to unit sum, and the correlation is evaluated over the
 * "valid" region only, exactly as filter2(...,'valid') does.
 *
 * Strict C89.
 */

#include "nbs.h"
#include "nbs_metrics.h"

#include <stdlib.h>
#include <math.h>

void nbs_ssim_cfg_default(nbs_ssim_cfg *c)
{
    c->win        = 11;
    c->win_sigma  = 1.5;
    c->K1         = 0.01;
    c->K2         = 0.03;
    c->L          = 255.0;
    c->downsample = -1;   /* automatic, as in ssim.m */
}

/* Separable "valid" convolution with a symmetric kernel.
 * dst is (h-2r) x (w-2r) and must be pre-allocated by the caller. */
int nbs_conv_valid_sep(const nbs_plane *src, const nbs_kernel *k, nbs_plane *dst)
{
    int r = k->r;
    int w = src->w, h = src->h;
    int vw = w - 2 * r, vh = h - 2 * r;
    int x, y, t;
    double *tmp;

    if (vw <= 0 || vh <= 0) { nbs_set_error("valid conv: image smaller than window"); return 1; }
    if (dst->w != vw || dst->h != vh) { nbs_set_error("valid conv: dst size mismatch"); return 1; }

    /* Horizontal pass into a (vw x h) scratch. */
    tmp = (double *)nbs_malloc((size_t)vw * (size_t)h * sizeof(double));
    if (tmp == NULL) return 1;

    for (y = 0; y < h; ++y) {
        const double *srow = src->d + (size_t)y * (size_t)w;
        double *drow = tmp + (size_t)y * (size_t)vw;
        for (x = 0; x < vw; ++x) {
            double acc = 0.0;
            for (t = 0; t <= 2 * r; ++t) acc += k->tap[t] * srow[x + t];
            drow[x] = acc;
        }
    }
    /* Vertical pass. */
    for (y = 0; y < vh; ++y) {
        double *drow = dst->d + (size_t)y * (size_t)vw;
        for (x = 0; x < vw; ++x) {
            double acc = 0.0;
            for (t = 0; t <= 2 * r; ++t)
                acc += k->tap[t] * tmp[(size_t)(y + t) * (size_t)vw + (size_t)x];
            drow[x] = acc;
        }
    }
    free(tmp);
    return 0;
}

/* Build the 11x11-equivalent 1-D Gaussian window (unit sum). */
int nbs_ssim_window(nbs_kernel *k, int win, double sigma)
{
    int i, r, n;
    double s, t, inv2s2;

    if (win < 1 || (win % 2) == 0) { nbs_set_error("ssim: window size must be odd"); return 1; }
    r = (win - 1) / 2;
    n = win;
    if (nbs_kernel_alloc(k, r) != 0) return 1;

    /* MATLAB fspecial('gaussian', N, sigma) samples at -r..r and normalises
     * the *2-D* kernel.  For a separable Gaussian this is identical to
     * normalising the 1-D kernel and taking the outer product. */
    inv2s2 = 1.0 / (2.0 * sigma * sigma);
    s = 0.0;
    for (i = -r; i <= r; ++i) {
        t = exp(-((double)i * (double)i) * inv2s2);
        k->tap[i + r] = t;
        s += t;
    }
    for (i = 0; i < n; ++i) k->tap[i] /= s;
    return 0;
}

int nbs_ssim_ex(const nbs_plane *x, const nbs_plane *y,
                const nbs_ssim_cfg *cfg,
                double *out_mssim, double *out_mcs)
{
    nbs_kernel win;
    nbs_plane xs, ys;          /* possibly downsampled inputs */
    nbs_plane xx, yy, xy;      /* products */
    nbs_plane mx, my, sxx, syy, sxy;
    const nbs_plane *px, *py;
    double C1, C2, ssum, cssum;
    int r, vw, vh, i, n, f, rc = 1;

    xs.d = ys.d = NULL;
    xx.d = yy.d = xy.d = NULL;
    mx.d = my.d = sxx.d = syy.d = sxy.d = NULL;
    win.tap = NULL;

    if (x->w != y->w || x->h != y->h) { nbs_set_error("ssim: size mismatch"); return 1; }

    /* Automatic downsampling, per ssim.m. */
    f = cfg->downsample;
    if (f < 0) {
        int mn = (x->w < x->h) ? x->w : x->h;
        f = (int)floor((double)mn / 256.0 + 0.5);
        if (f < 1) f = 1;
    }
    if (f > 1) {
        if (nbs_downsample(x, &xs, f, NBS_BORDER_SYMMETRIC) != 0) goto done;
        if (nbs_downsample(y, &ys, f, NBS_BORDER_SYMMETRIC) != 0) goto done;
        px = &xs; py = &ys;
    } else {
        px = x; py = y;
    }

    if (nbs_ssim_window(&win, cfg->win, cfg->win_sigma) != 0) goto done;
    r  = win.r;
    vw = px->w - 2 * r;
    vh = px->h - 2 * r;
    if (vw <= 0 || vh <= 0) { nbs_set_error("ssim: image smaller than window"); goto done; }

    if (nbs_plane_alloc(&xx, px->w, px->h) != 0) goto done;
    if (nbs_plane_alloc(&yy, px->w, px->h) != 0) goto done;
    if (nbs_plane_alloc(&xy, px->w, px->h) != 0) goto done;
    n = px->w * px->h;
    for (i = 0; i < n; ++i) {
        xx.d[i] = px->d[i] * px->d[i];
        yy.d[i] = py->d[i] * py->d[i];
        xy.d[i] = px->d[i] * py->d[i];
    }

    if (nbs_plane_alloc(&mx,  vw, vh) != 0) goto done;
    if (nbs_plane_alloc(&my,  vw, vh) != 0) goto done;
    if (nbs_plane_alloc(&sxx, vw, vh) != 0) goto done;
    if (nbs_plane_alloc(&syy, vw, vh) != 0) goto done;
    if (nbs_plane_alloc(&sxy, vw, vh) != 0) goto done;

    if (nbs_conv_valid_sep(px,  &win, &mx)  != 0) goto done;
    if (nbs_conv_valid_sep(py,  &win, &my)  != 0) goto done;
    if (nbs_conv_valid_sep(&xx, &win, &sxx) != 0) goto done;
    if (nbs_conv_valid_sep(&yy, &win, &syy) != 0) goto done;
    if (nbs_conv_valid_sep(&xy, &win, &sxy) != 0) goto done;

    C1 = cfg->K1 * cfg->L; C1 *= C1;
    C2 = cfg->K2 * cfg->L; C2 *= C2;

    ssum = 0.0; cssum = 0.0;
    n = vw * vh;
    for (i = 0; i < n; ++i) {
        double mx2  = mx.d[i] * mx.d[i];
        double my2  = my.d[i] * my.d[i];
        double mxmy = mx.d[i] * my.d[i];
        double vx   = sxx.d[i] - mx2;
        double vy   = syy.d[i] - my2;
        double vxy  = sxy.d[i] - mxmy;
        double lum  = (2.0 * mxmy + C1) / (mx2 + my2 + C1);
        double cs   = (2.0 * vxy  + C2) / (vx  + vy  + C2);
        ssum  += lum * cs;
        cssum += cs;
    }

    if (out_mssim) *out_mssim = ssum  / (double)n;
    if (out_mcs)   *out_mcs   = cssum / (double)n;
    rc = 0;

done:
    nbs_kernel_free(&win);
    nbs_plane_free(&xs); nbs_plane_free(&ys);
    nbs_plane_free(&xx); nbs_plane_free(&yy); nbs_plane_free(&xy);
    nbs_plane_free(&mx); nbs_plane_free(&my);
    nbs_plane_free(&sxx); nbs_plane_free(&syy); nbs_plane_free(&sxy);
    return rc;
}

int nbs_ssim(const nbs_plane *x, const nbs_plane *y,
             const nbs_ssim_cfg *cfg, double *out_mssim)
{
    return nbs_ssim_ex(x, y, cfg, out_mssim, NULL);
}
