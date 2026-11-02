/* fft.c -- iterative radix-2 complex FFT, strict C89.
 *
 * Nothing clever here: decimation in time, bit-reversal permutation up
 * front, precomputed twiddles.  The metric spends most of its time in
 * the windowed statistics rather than in the transform, so a textbook
 * kernel is fast enough and much easier to trust.
 */

#include <stdlib.h>
#include <math.h>
#include "nbs_fft.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int nbs_fft_size(int n)
{
    int s = 1;
    if (n < 1) return 0;
    while (s < n) {
        if (s > (1 << 29)) return 0;
        s <<= 1;
    }
    return s;
}

int nbs_fft_plan_init(nbs_fft_plan *p, int n)
{
    int k, half, j, i, bit;

    p->n = 0; p->cs = NULL; p->sn = NULL; p->rev = NULL;
    if (n < 1 || (n & (n - 1)) != 0) return 1;   /* powers of two only */

    half = (n > 1) ? n / 2 : 1;
    p->cs  = (double *)malloc((size_t)half * sizeof(double));
    p->sn  = (double *)malloc((size_t)half * sizeof(double));
    p->rev = (int *)malloc((size_t)n * sizeof(int));
    if (p->cs == NULL || p->sn == NULL || p->rev == NULL) {
        nbs_fft_plan_free(p);
        return 1;
    }
    for (k = 0; k < half; ++k) {
        double a = -2.0 * M_PI * (double)k / (double)n;
        p->cs[k] = cos(a);
        p->sn[k] = sin(a);
    }
    /* Bit-reversal permutation, built incrementally. */
    j = 0;
    for (i = 0; i < n; ++i) {
        p->rev[i] = j;
        bit = n >> 1;
        while (bit > 0 && (j & bit) != 0) { j ^= bit; bit >>= 1; }
        j |= bit;
    }
    p->n = n;
    return 0;
}

void nbs_fft_plan_free(nbs_fft_plan *p)
{
    if (p == NULL) return;
    free(p->cs); free(p->sn); free(p->rev);
    p->cs = NULL; p->sn = NULL; p->rev = NULL; p->n = 0;
}

void nbs_fft_strided(const nbs_fft_plan *p, double *re, double *im,
                     int stride, int inverse)
{
    int n = p->n, i, j, len, half, step, k, a, b;
    double wr, wi, tr, ti;

    if (n <= 1) return;

    for (i = 0; i < n; ++i) {
        j = p->rev[i];
        if (j > i) {
            double t;
            t = re[(size_t)i * stride]; re[(size_t)i * stride] = re[(size_t)j * stride]; re[(size_t)j * stride] = t;
            t = im[(size_t)i * stride]; im[(size_t)i * stride] = im[(size_t)j * stride]; im[(size_t)j * stride] = t;
        }
    }

    for (len = 2; len <= n; len <<= 1) {
        half = len >> 1;
        step = n / len;
        for (i = 0; i < n; i += len) {
            for (k = 0; k < half; ++k) {
                wr = p->cs[k * step];
                wi = p->sn[k * step];
                if (inverse) wi = -wi;
                a = i + k;
                b = a + half;
                tr = re[(size_t)b * stride] * wr - im[(size_t)b * stride] * wi;
                ti = re[(size_t)b * stride] * wi + im[(size_t)b * stride] * wr;
                re[(size_t)b * stride] = re[(size_t)a * stride] - tr;
                im[(size_t)b * stride] = im[(size_t)a * stride] - ti;
                re[(size_t)a * stride] += tr;
                im[(size_t)a * stride] += ti;
            }
        }
    }

    if (inverse) {
        double s = 1.0 / (double)n;
        for (i = 0; i < n; ++i) {
            re[(size_t)i * stride] *= s;
            im[(size_t)i * stride] *= s;
        }
    }
}

void nbs_fft(const nbs_fft_plan *p, double *re, double *im, int inverse)
{
    nbs_fft_strided(p, re, im, 1, inverse);
}
