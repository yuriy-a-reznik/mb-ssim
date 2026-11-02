/* nbs_fft.h -- minimal self-contained FFT for the DFT-domain filterbank.
 *
 * Strict C89.  No external dependencies beyond libm.
 *
 * The transform sizes used by the metric are always powers of two: the
 * caller mirror-pads the image up to nbs_fft_size() before transforming,
 * which both removes the wrap-around of circular convolution and lets a
 * plain iterative radix-2 kernel do all the work.
 *
 * Complex data is stored as two separate double arrays (re[], im[])
 * rather than interleaved.  That keeps the row/column passes of the 2-D
 * transform simple and lets the real-input path leave im[] alone.
 */

#ifndef NBS_FFT_H
#define NBS_FFT_H

/* Smallest power of two >= n (n >= 1).  Returns 0 on overflow. */
int nbs_fft_size(int n);

/* Twiddle-factor table for a given transform length.  Build once, reuse
 * for every row and column pass of the same length. */
typedef struct {
    int     n;        /* transform length, a power of two */
    double *cs;       /* cos(-2*pi*k/n), k < n/2 */
    double *sn;       /* sin(-2*pi*k/n), k < n/2 */
    int    *rev;      /* bit-reversal permutation */
} nbs_fft_plan;

int  nbs_fft_plan_init(nbs_fft_plan *p, int n);
void nbs_fft_plan_free(nbs_fft_plan *p);

/* In-place complex FFT of length p->n.  inverse == 0 gives the forward
 * transform with no scaling; inverse != 0 gives the inverse transform
 * including the 1/n factor, so forward followed by inverse is the
 * identity. */
void nbs_fft(const nbs_fft_plan *p, double *re, double *im, int inverse);

/* Strided variants, for transforming the columns of a row-major plane
 * without copying.  Element k lives at re[k*stride]. */
void nbs_fft_strided(const nbs_fft_plan *p, double *re, double *im,
                     int stride, int inverse);

#endif /* NBS_FFT_H */
