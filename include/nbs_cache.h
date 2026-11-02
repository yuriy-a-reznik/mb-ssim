/* nbs_cache.h -- on-disk cache of the log-xi tensor.
 *
 * Strict C89.
 *
 * Building the tensor costs one full pass over a database; fitting a weight
 * vector against it costs a few dot products.  Persisting it is what makes
 * trying a second objective, a different regularisation, or -- the reason
 * this moved out of nbstune.c into its own module -- a *joint* fit across
 * two viewing distances, essentially free.
 *
 * The metric for any exponent vector w is
 *     mean_p exp( sum_k w[k] * L[p][k] )
 * so the tensor is stored as floats: it is summed in double and the extra
 * precision buys nothing against a 24-bit mantissa on a log-ratio bounded
 * by log(xi_floor).
 */

#ifndef NBS_CACHE_H
#define NBS_CACHE_H

#include "nbs.h"
#include "nbs_metrics.h"

typedef struct {
    float  *L;        /* concatenated: entry e occupies np[e]*nw floats */
    long   *off;
    int    *np;
    double *subj;     /* subjective score                              */
    int    *group;    /* reference-image id, for content-disjoint folds */
    int     n;        /* number of entries with data                   */
    int     nw;       /* nbs_w_count(nbands)                           */

    /* Header: what the tensor depends on.  A joint fit needs these to
     * check that two caches describe the same bank at different ppd. */
    int     nbands;
    int     stride;
    int     w_sym;    /* 1 = w_ij tied to w_ji, 0 = all n^2 free       */
    double  ppd;
    double  f[NBS_MAX_BANDS];
    double  rho[NBS_MAX_BANDS];
    double  kc[NBS_MAX_BANDS];
} nbs_logxi_cache;

int  nbs_cache_save(const nbs_logxi_cache *c, const char *path);
int  nbs_cache_load(nbs_logxi_cache *c, const char *path);
void nbs_cache_free(nbs_logxi_cache *c);

/* Fill the header fields from a metric configuration. */
void nbs_cache_stamp(nbs_logxi_cache *c, const nbs_nbssim_cfg *cfg, int stride);

/* Refuse a cache that was not built with this configuration.  Returns 0 if
 * it matches, non-zero with a message in nbs_last_error() if not. */
int  nbs_cache_check(const nbs_logxi_cache *c, const nbs_nbssim_cfg *cfg,
                     int stride);

#endif /* NBS_CACHE_H */
