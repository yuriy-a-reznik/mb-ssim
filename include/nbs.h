/* nbs.h -- common types and utilities for the N-band SSIM project.
 *
 * Strict C89.  No dependencies beyond the C89 standard library.
 *
 * Copyright (c) 2026.  Provided as research code.
 */

#ifndef NBS_H
#define NBS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ */
/* Planes                                                              */
/* ------------------------------------------------------------------ */

/* A single-channel image of doubles, row-major, stride == w. */
typedef struct {
    int     w;
    int     h;
    double *d;
} nbs_plane;

int  nbs_plane_alloc(nbs_plane *p, int w, int h);
void nbs_plane_free(nbs_plane *p);
int  nbs_plane_copy(nbs_plane *dst, const nbs_plane *src);

#define NBS_AT(p, x, y) ((p)->d[(size_t)(y) * (size_t)(p)->w + (size_t)(x)])

/* ------------------------------------------------------------------ */
/* Image loading                                                       */
/* ------------------------------------------------------------------ */

/* Luma coefficient sets. */
typedef enum {
    NBS_LUMA_BT601 = 0,   /* 0.299 / 0.587 / 0.114  (MATLAB rgb2gray)   */
    NBS_LUMA_BT709 = 1,   /* 0.2126 / 0.7152 / 0.0722                   */
    NBS_LUMA_MEAN  = 2    /* (R+G+B)/3                                  */
} nbs_luma_mode;

/* An 8-bit/channel image as loaded from disk.  `nch` is 1 or 3. */
typedef struct {
    int            w;
    int            h;
    int            nch;
    unsigned char *px;    /* w*h*nch, row-major, interleaved */
} nbs_image8;

void nbs_image8_free(nbs_image8 *im);

/* Load a PNM (P2/P5/P6), BMP (24/32-bit uncompressed, or 8-bit paletted)
 * or PNG (bit depth 8, colour types 0/2/3/4/6, non-interlaced) file.
 * Format is detected from the file's magic bytes, not its extension.
 * Returns 0 on success, non-zero on failure (message in nbs_last_error()). */
int nbs_image8_load(const char *path, nbs_image8 *out);

/* Extract a luma plane from an 8-bit image.  Values are in [0,255]. */
int nbs_image8_to_luma(const nbs_image8 *im, nbs_luma_mode mode, nbs_plane *out);

/* Convenience: load + luma in one call. */
int nbs_load_luma(const char *path, nbs_luma_mode mode, nbs_plane *out);

const char *nbs_last_error(void);
void        nbs_set_error(const char *fmt_free_text);

/* ------------------------------------------------------------------ */
/* Filtering                                                           */
/* ------------------------------------------------------------------ */

/* Border handling for separable convolution. */
typedef enum {
    NBS_BORDER_SYMMETRIC = 0, /* mirror without repeating the edge sample  */
    NBS_BORDER_REPLICATE = 1,
    NBS_BORDER_ZERO      = 2
} nbs_border;

/* A symmetric FIR kernel of odd length (2*r+1), tap[0] is the centre. */
typedef struct {
    int     r;      /* radius */
    double *tap;    /* length 2*r+1, index 0 .. 2*r, centre at r */
} nbs_kernel;

int  nbs_kernel_alloc(nbs_kernel *k, int r);
void nbs_kernel_free(nbs_kernel *k);

/* Sampled, unit-sum Gaussian with standard deviation `sigma` (in pixels).
 * The radius is chosen as ceil(radius_mult*sigma), clamped to >= 1. */
int nbs_kernel_gaussian(nbs_kernel *k, double sigma, double radius_mult);

/* Unit-sum box kernel of length n (n odd). */
int nbs_kernel_box(nbs_kernel *k, int n);

/* Separable convolution, "same" size.  src and dst must not alias. */
int nbs_conv_sep(const nbs_plane *src, nbs_plane *dst,
                 const nbs_kernel *kx, const nbs_kernel *ky,
                 nbs_border border);

/* In-place-safe wrapper that allocates its own scratch. */
int nbs_gauss_blur(const nbs_plane *src, nbs_plane *dst, double sigma,
                   nbs_border border);

/* Downsample by an integer factor after box-averaging over factor x factor
 * blocks (the low-pass used by Wang's ms_ssim.m and ssim.m). */
int nbs_downsample(const nbs_plane *src, nbs_plane *dst, int factor,
                   nbs_border border);

/* ------------------------------------------------------------------ */
/* Pointwise helpers                                                   */
/* ------------------------------------------------------------------ */

int nbs_plane_mul(const nbs_plane *a, const nbs_plane *b, nbs_plane *out);
int nbs_plane_sub(const nbs_plane *a, const nbs_plane *b, nbs_plane *out);
double nbs_plane_mean(const nbs_plane *p);

/* ------------------------------------------------------------------ */
/* Misc                                                                */
/* ------------------------------------------------------------------ */

void *nbs_malloc(size_t n);
void *nbs_calloc(size_t n, size_t sz);

#ifdef __cplusplus
}
#endif

#endif /* NBS_H */
