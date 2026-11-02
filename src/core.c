/* core.c -- planes, kernels, separable convolution. Strict C89. */

#include "nbs.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ------------------------------------------------------------------ */

static char nbs_err[512] = "";

const char *nbs_last_error(void)
{
    return nbs_err;
}

void nbs_set_error(const char *msg)
{
    if (msg == NULL) {
        nbs_err[0] = '\0';
    } else {
        size_t n = strlen(msg);
        if (n >= sizeof(nbs_err)) n = sizeof(nbs_err) - 1;
        memcpy(nbs_err, msg, n);
        nbs_err[n] = '\0';
    }
}

void *nbs_malloc(size_t n)
{
    void *p;
    if (n == 0) n = 1;
    p = malloc(n);
    if (p == NULL) nbs_set_error("out of memory");
    return p;
}

void *nbs_calloc(size_t n, size_t sz)
{
    void *p;
    if (n == 0 || sz == 0) { n = 1; sz = 1; }
    p = calloc(n, sz);
    if (p == NULL) nbs_set_error("out of memory");
    return p;
}

/* ------------------------------------------------------------------ */
/* Planes                                                              */
/* ------------------------------------------------------------------ */

int nbs_plane_alloc(nbs_plane *p, int w, int h)
{
    if (p == NULL) return 1;
    p->w = 0; p->h = 0; p->d = NULL;
    if (w <= 0 || h <= 0) { nbs_set_error("plane: bad dimensions"); return 1; }
    p->d = (double *)nbs_malloc((size_t)w * (size_t)h * sizeof(double));
    if (p->d == NULL) return 1;
    p->w = w; p->h = h;
    return 0;
}

void nbs_plane_free(nbs_plane *p)
{
    if (p == NULL) return;
    free(p->d);
    p->d = NULL;
    p->w = 0;
    p->h = 0;
}

int nbs_plane_copy(nbs_plane *dst, const nbs_plane *src)
{
    if (nbs_plane_alloc(dst, src->w, src->h) != 0) return 1;
    memcpy(dst->d, src->d, (size_t)src->w * (size_t)src->h * sizeof(double));
    return 0;
}

int nbs_plane_mul(const nbs_plane *a, const nbs_plane *b, nbs_plane *out)
{
    size_t i, n;
    if (a->w != b->w || a->h != b->h) { nbs_set_error("mul: size mismatch"); return 1; }
    if (out->d == NULL && nbs_plane_alloc(out, a->w, a->h) != 0) return 1;
    n = (size_t)a->w * (size_t)a->h;
    for (i = 0; i < n; ++i) out->d[i] = a->d[i] * b->d[i];
    return 0;
}

int nbs_plane_sub(const nbs_plane *a, const nbs_plane *b, nbs_plane *out)
{
    size_t i, n;
    if (a->w != b->w || a->h != b->h) { nbs_set_error("sub: size mismatch"); return 1; }
    if (out->d == NULL && nbs_plane_alloc(out, a->w, a->h) != 0) return 1;
    n = (size_t)a->w * (size_t)a->h;
    for (i = 0; i < n; ++i) out->d[i] = a->d[i] - b->d[i];
    return 0;
}

double nbs_plane_mean(const nbs_plane *p)
{
    size_t i, n = (size_t)p->w * (size_t)p->h;
    double s = 0.0;
    for (i = 0; i < n; ++i) s += p->d[i];
    return (n > 0) ? s / (double)n : 0.0;
}

/* ------------------------------------------------------------------ */
/* Kernels                                                             */
/* ------------------------------------------------------------------ */

int nbs_kernel_alloc(nbs_kernel *k, int r)
{
    if (k == NULL) return 1;
    k->r = 0; k->tap = NULL;
    if (r < 0) { nbs_set_error("kernel: negative radius"); return 1; }
    k->tap = (double *)nbs_calloc((size_t)(2 * r + 1), sizeof(double));
    if (k->tap == NULL) return 1;
    k->r = r;
    return 0;
}

void nbs_kernel_free(nbs_kernel *k)
{
    if (k == NULL) return;
    free(k->tap);
    k->tap = NULL;
    k->r = 0;
}

int nbs_kernel_gaussian(nbs_kernel *k, double sigma, double radius_mult)
{
    int i, r;
    double s, t, inv2s2;

    if (sigma <= 0.0) {
        /* Degenerate: identity kernel. */
        if (nbs_kernel_alloc(k, 0) != 0) return 1;
        k->tap[0] = 1.0;
        return 0;
    }
    if (radius_mult <= 0.0) radius_mult = 4.0;
    r = (int)ceil(radius_mult * sigma);
    if (r < 1) r = 1;
    if (nbs_kernel_alloc(k, r) != 0) return 1;

    inv2s2 = 1.0 / (2.0 * sigma * sigma);
    s = 0.0;
    for (i = -r; i <= r; ++i) {
        t = exp(-((double)i * (double)i) * inv2s2);
        k->tap[i + r] = t;
        s += t;
    }
    for (i = 0; i <= 2 * r; ++i) k->tap[i] /= s;
    return 0;
}

int nbs_kernel_box(nbs_kernel *k, int n)
{
    int i, r;
    if (n < 1 || (n % 2) == 0) { nbs_set_error("box kernel: n must be odd >= 1"); return 1; }
    r = (n - 1) / 2;
    if (nbs_kernel_alloc(k, r) != 0) return 1;
    for (i = 0; i < n; ++i) k->tap[i] = 1.0 / (double)n;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Convolution                                                         */
/* ------------------------------------------------------------------ */

/* Map an out-of-range coordinate into [0,n). */
static int reflect(int i, int n, nbs_border b)
{
    if (n == 1) return 0;
    switch (b) {
    case NBS_BORDER_REPLICATE:
        if (i < 0) return 0;
        if (i >= n) return n - 1;
        return i;
    case NBS_BORDER_SYMMETRIC:
    default:
        /* Mirror with edge repeated: -1 -> 0, -2 -> 1, n -> n-1, ...
         * This matches MATLAB imfilter(...,'symmetric'). */
        while (i < 0 || i >= n) {
            if (i < 0) i = -i - 1;
            if (i >= n) i = 2 * n - i - 1;
        }
        return i;
    }
}

int nbs_conv_sep(const nbs_plane *src, nbs_plane *dst,
                 const nbs_kernel *kx, const nbs_kernel *ky,
                 nbs_border border)
{
    int x, y, t, w, h, rx, ry;
    double *tmp;
    double acc;

    if (src == NULL || dst == NULL) return 1;
    w = src->w; h = src->h;
    if (dst->d == NULL && nbs_plane_alloc(dst, w, h) != 0) return 1;
    if (dst->w != w || dst->h != h) { nbs_set_error("conv: dst size mismatch"); return 1; }
    if (dst->d == src->d) { nbs_set_error("conv: src and dst alias"); return 1; }

    tmp = (double *)nbs_malloc((size_t)w * (size_t)h * sizeof(double));
    if (tmp == NULL) return 1;

    /* Horizontal pass. */
    if (kx != NULL && !(kx->r == 0 && kx->tap[0] == 1.0)) {
        rx = kx->r;
        for (y = 0; y < h; ++y) {
            const double *srow = src->d + (size_t)y * (size_t)w;
            double *drow = tmp + (size_t)y * (size_t)w;
            for (x = 0; x < w; ++x) {
                acc = 0.0;
                if (x >= rx && x + rx < w) {
                    for (t = -rx; t <= rx; ++t)
                        acc += kx->tap[t + rx] * srow[x + t];
                } else if (border == NBS_BORDER_ZERO) {
                    for (t = -rx; t <= rx; ++t) {
                        int xi = x + t;
                        if (xi >= 0 && xi < w) acc += kx->tap[t + rx] * srow[xi];
                    }
                } else {
                    for (t = -rx; t <= rx; ++t)
                        acc += kx->tap[t + rx] * srow[reflect(x + t, w, border)];
                }
                drow[x] = acc;
            }
        }
    } else {
        memcpy(tmp, src->d, (size_t)w * (size_t)h * sizeof(double));
    }

    /* Vertical pass. */
    if (ky != NULL && !(ky->r == 0 && ky->tap[0] == 1.0)) {
        ry = ky->r;
        for (y = 0; y < h; ++y) {
            double *drow = dst->d + (size_t)y * (size_t)w;
            for (x = 0; x < w; ++x) {
                acc = 0.0;
                if (y >= ry && y + ry < h) {
                    for (t = -ry; t <= ry; ++t)
                        acc += ky->tap[t + ry] * tmp[(size_t)(y + t) * (size_t)w + (size_t)x];
                } else if (border == NBS_BORDER_ZERO) {
                    for (t = -ry; t <= ry; ++t) {
                        int yi = y + t;
                        if (yi >= 0 && yi < h)
                            acc += ky->tap[t + ry] * tmp[(size_t)yi * (size_t)w + (size_t)x];
                    }
                } else {
                    for (t = -ry; t <= ry; ++t)
                        acc += ky->tap[t + ry]
                             * tmp[(size_t)reflect(y + t, h, border) * (size_t)w + (size_t)x];
                }
                drow[x] = acc;
            }
        }
    } else {
        memcpy(dst->d, tmp, (size_t)w * (size_t)h * sizeof(double));
    }

    free(tmp);
    return 0;
}

int nbs_gauss_blur(const nbs_plane *src, nbs_plane *dst, double sigma,
                   nbs_border border)
{
    nbs_kernel k;
    int rc;
    if (nbs_kernel_gaussian(&k, sigma, 4.0) != 0) return 1;
    rc = nbs_conv_sep(src, dst, &k, &k, border);
    nbs_kernel_free(&k);
    return rc;
}

int nbs_downsample(const nbs_plane *src, nbs_plane *dst, int factor,
                   nbs_border border)
{
    nbs_plane lp;
    nbs_kernel k;
    int nw, nh, x, y, i;
    int r;

    if (factor <= 1) return nbs_plane_copy(dst, src);

    /* Wang's ms_ssim.m / ssim.m use ones(f)/f^2, i.e. an even-length box
     * filter, then take every f-th sample starting at index 0.  An even
     * box has no integer centre; we emulate it exactly below rather than
     * going through nbs_conv_sep. */
    (void)k;
    lp.d = NULL;
    if (nbs_plane_alloc(&lp, src->w, src->h) != 0) return 1;

    r = factor;
    for (y = 0; y < src->h; ++y) {
        for (x = 0; x < src->w; ++x) {
            double acc = 0.0;
            int dy, dx;
            /* MATLAB imfilter with an even f x f kernel places the origin at
             * floor(f/2) (0-based), so taps cover offsets
             * [-floor((f-1)/2) .. floor(f/2)]. */
            int lo = -((r - 1) / 2);
            int hi = r / 2;
            for (dy = lo; dy <= hi; ++dy) {
                int yy = reflect(y + dy, src->h, border);
                for (dx = lo; dx <= hi; ++dx) {
                    int xx = reflect(x + dx, src->w, border);
                    acc += src->d[(size_t)yy * (size_t)src->w + (size_t)xx];
                }
            }
            lp.d[(size_t)y * (size_t)src->w + (size_t)x] = acc / ((double)r * (double)r);
        }
    }

    nw = (src->w + factor - 1) / factor;
    nh = (src->h + factor - 1) / factor;
    if (nbs_plane_alloc(dst, nw, nh) != 0) { nbs_plane_free(&lp); return 1; }
    for (y = 0; y < nh; ++y) {
        for (x = 0; x < nw; ++x) {
            i = 0; (void)i;
            dst->d[(size_t)y * (size_t)nw + (size_t)x] =
                lp.d[(size_t)(y * factor) * (size_t)src->w + (size_t)(x * factor)];
        }
    }
    nbs_plane_free(&lp);
    return 0;
}
