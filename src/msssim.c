/* msssim.c -- reference implementation of the MS-SSIM index.
 *
 * Faithful port of Z. Wang's ms_ssim.m:
 *   Z. Wang, E. P. Simoncelli and A. C. Bovik, "Multiscale structural
 *   similarity for image quality assessment", 37th Asilomar Conf. on
 *   Signals, Systems & Computers, 2003, vol. 2, pp. 1398-1402.
 *
 *   MS-SSIM = [ prod_{j=1..M-1} cs_j^beta_j ] * (l_M * cs_M)^beta_M
 *
 * with beta = [0.0448 0.2856 0.3001 0.2363 0.1333] and a 2x2 box
 * low-pass followed by decimation by two between scales.
 *
 * Strict C89.
 */

#include "nbs.h"
#include "nbs_metrics.h"

#include <stdlib.h>
#include <math.h>

void nbs_msssim_cfg_default(nbs_msssim_cfg *c)
{
    int i;
    nbs_ssim_cfg_default(&c->base);
    c->base.downsample = 0;   /* ms_ssim.m does not pre-downsample */
    c->levels = 5;
    c->weight[0] = 0.0448;
    c->weight[1] = 0.2856;
    c->weight[2] = 0.3001;
    c->weight[3] = 0.2363;
    c->weight[4] = 0.1333;
    for (i = 5; i < NBS_MSSSIM_MAX_LEVELS; ++i) c->weight[i] = 0.0;
}

int nbs_msssim(const nbs_plane *x, const nbs_plane *y,
               const nbs_msssim_cfg *cfg, double *out)
{
    nbs_plane a, b, na, nb;
    double mssim[NBS_MSSSIM_MAX_LEVELS];
    double mcs[NBS_MSSSIM_MAX_LEVELS];
    double acc;
    int l, rc = 1, minsz;
    nbs_ssim_cfg sc;

    a.d = b.d = na.d = nb.d = NULL;

    if (x->w != y->w || x->h != y->h) { nbs_set_error("ms-ssim: size mismatch"); return 1; }
    if (cfg->levels < 1 || cfg->levels > NBS_MSSSIM_MAX_LEVELS) {
        nbs_set_error("ms-ssim: bad level count"); return 1;
    }

    /* ms_ssim.m requires min(M,N) > (win-1)*2^(levels-1) so that the
     * coarsest scale is still at least as large as the window. */
    minsz = (x->w < x->h) ? x->w : x->h;
    if (minsz < (cfg->base.win - 1) * (1 << (cfg->levels - 1)) + 1) {
        nbs_set_error("ms-ssim: image too small for the requested number of scales");
        return 1;
    }

    sc = cfg->base;
    sc.downsample = 0;

    if (nbs_plane_copy(&a, x) != 0) goto done;
    if (nbs_plane_copy(&b, y) != 0) goto done;

    for (l = 0; l < cfg->levels; ++l) {
        if (nbs_ssim_ex(&a, &b, &sc, &mssim[l], &mcs[l]) != 0) goto done;
        if (l + 1 < cfg->levels) {
            na.d = NULL; nb.d = NULL;
            if (nbs_downsample(&a, &na, 2, NBS_BORDER_SYMMETRIC) != 0) goto done;
            if (nbs_downsample(&b, &nb, 2, NBS_BORDER_SYMMETRIC) != 0) goto done;
            nbs_plane_free(&a); nbs_plane_free(&b);
            a = na; b = nb;
            na.d = NULL; nb.d = NULL;
        }
    }

    acc = 0.0;
    for (l = 0; l < cfg->levels - 1; ++l) {
        double v = mcs[l];
        if (v < 1e-12) v = 1e-12;
        acc += cfg->weight[l] * log(v);
    }
    {
        double v = mssim[cfg->levels - 1];
        if (v < 1e-12) v = 1e-12;
        acc += cfg->weight[cfg->levels - 1] * log(v);
    }
    *out = exp(acc);
    rc = 0;

done:
    nbs_plane_free(&a); nbs_plane_free(&b);
    nbs_plane_free(&na); nbs_plane_free(&nb);
    return rc;
}
