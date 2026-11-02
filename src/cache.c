/* cache.c -- read/write the log-xi tensor.  Strict C89. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nbs_cache.h"

/* The layout is unchanged from the version nbstune wrote inline; what is
 * new is that the header is read back out rather than merely checked, so a
 * second tool can discover a cache's ppd and band set instead of having to
 * be told what they are. */
#define NBSC_MAGIC "NBSCACHE3"

void nbs_cache_stamp(nbs_logxi_cache *c, const nbs_nbssim_cfg *cfg, int stride)
{
    int i;
    c->nbands = cfg->bank.nbands;
    c->ppd    = cfg->bank.ppd;
    c->stride = stride;
    c->w_sym  = cfg->w_sym;
    c->nw     = nbs_w_count_ex(cfg->bank.nbands, cfg->w_sym);
    for (i = 0; i < NBS_MAX_BANDS; ++i) {
        c->f[i]   = cfg->bank.f[i];
        c->rho[i] = cfg->bank.rho[i];
        c->kc[i]  = cfg->kc[i];
    }
}

int nbs_cache_check(const nbs_logxi_cache *c, const nbs_nbssim_cfg *cfg,
                    int stride)
{
    int i;
    if (c->nbands != cfg->bank.nbands ||
        c->w_sym != cfg->w_sym ||
        c->nw != nbs_w_count_ex(cfg->bank.nbands, cfg->w_sym) ||
        c->stride != stride ||
        fabs(c->ppd - cfg->bank.ppd) > 1e-9) {
        nbs_set_error("cache was built with a different band count, exponent "
                      "symmetry, stride or ppd -- rebuild it");
        return 1;
    }
    for (i = 0; i < c->nbands; ++i) {
        if (fabs(c->f[i] - cfg->bank.f[i]) > 1e-9 ||
            fabs(c->rho[i] - cfg->bank.rho[i]) > 1e-9 ||
            fabs(c->kc[i] - cfg->kc[i]) > 1e-12) {
            nbs_set_error("cache was built with different band frequencies, "
                          "rho or kc -- rebuild it");
            return 1;
        }
    }
    return 0;
}

int nbs_cache_save(const nbs_logxi_cache *c, const char *path)
{
    FILE *f = fopen(path, "wb");
    int i, ok = 1;
    long total = 0;

    if (f == NULL) { nbs_set_error("cannot write cache file"); return 1; }
    for (i = 0; i < c->n; ++i) total += (long)c->np[i] * (long)c->nw;

    fwrite(NBSC_MAGIC, 1, 9, f);
    fwrite(&c->n, sizeof(int), 1, f);
    fwrite(&c->nw, sizeof(int), 1, f);
    fwrite(&c->stride, sizeof(int), 1, f);
    fwrite(&total, sizeof(long), 1, f);
    fwrite(&c->nbands, sizeof(int), 1, f);
    fwrite(&c->w_sym, sizeof(int), 1, f);
    fwrite(&c->ppd, sizeof(double), 1, f);
    fwrite(c->f, sizeof(double), NBS_MAX_BANDS, f);
    fwrite(c->rho, sizeof(double), NBS_MAX_BANDS, f);
    fwrite(c->kc, sizeof(double), NBS_MAX_BANDS, f);
    ok &= (fwrite(c->off, sizeof(long), (size_t)c->n, f) == (size_t)c->n);
    ok &= (fwrite(c->np, sizeof(int), (size_t)c->n, f) == (size_t)c->n);
    ok &= (fwrite(c->subj, sizeof(double), (size_t)c->n, f) == (size_t)c->n);
    ok &= (fwrite(c->group, sizeof(int), (size_t)c->n, f) == (size_t)c->n);
    ok &= (fwrite(c->L, sizeof(float), (size_t)total, f) == (size_t)total);
    fclose(f);
    if (!ok) { nbs_set_error("short write on cache file"); return 1; }
    return 0;
}

int nbs_cache_load(nbs_logxi_cache *c, const char *path)
{
    FILE *f = fopen(path, "rb");
    char magic[9];
    int ok = 1;
    long total;

    memset(c, 0, sizeof(*c));
    if (f == NULL) { nbs_set_error("cannot open cache file"); return 1; }
    if (fread(magic, 1, 9, f) != 9 || memcmp(magic, NBSC_MAGIC, 9) != 0) {
        fclose(f); nbs_set_error("cache file: bad magic"); return 1;
    }
    ok &= (fread(&c->n, sizeof(int), 1, f) == 1);
    ok &= (fread(&c->nw, sizeof(int), 1, f) == 1);
    ok &= (fread(&c->stride, sizeof(int), 1, f) == 1);
    ok &= (fread(&total, sizeof(long), 1, f) == 1);
    ok &= (fread(&c->nbands, sizeof(int), 1, f) == 1);
    ok &= (fread(&c->w_sym, sizeof(int), 1, f) == 1);
    ok &= (fread(&c->ppd, sizeof(double), 1, f) == 1);
    ok &= (fread(c->f, sizeof(double), NBS_MAX_BANDS, f) == NBS_MAX_BANDS);
    ok &= (fread(c->rho, sizeof(double), NBS_MAX_BANDS, f) == NBS_MAX_BANDS);
    ok &= (fread(c->kc, sizeof(double), NBS_MAX_BANDS, f) == NBS_MAX_BANDS);
    if (!ok || total <= 0 || c->n <= 0) {
        fclose(f); nbs_set_error("cache file: bad header"); return 1;
    }

    c->off   = (long *)  nbs_malloc((size_t)c->n * sizeof(long));
    c->np    = (int *)   nbs_malloc((size_t)c->n * sizeof(int));
    c->subj  = (double *)nbs_malloc((size_t)c->n * sizeof(double));
    c->group = (int *)   nbs_malloc((size_t)c->n * sizeof(int));
    c->L     = (float *) nbs_malloc((size_t)total * sizeof(float));
    if (!c->off || !c->np || !c->subj || !c->group || !c->L) {
        fclose(f); nbs_cache_free(c); return 1;
    }
    ok &= (fread(c->off, sizeof(long), (size_t)c->n, f) == (size_t)c->n);
    ok &= (fread(c->np, sizeof(int), (size_t)c->n, f) == (size_t)c->n);
    ok &= (fread(c->subj, sizeof(double), (size_t)c->n, f) == (size_t)c->n);
    ok &= (fread(c->group, sizeof(int), (size_t)c->n, f) == (size_t)c->n);
    ok &= (fread(c->L, sizeof(float), (size_t)total, f) == (size_t)total);
    fclose(f);
    if (!ok) { nbs_set_error("cache file is truncated"); nbs_cache_free(c); return 1; }
    return 0;
}

void nbs_cache_free(nbs_logxi_cache *c)
{
    if (c == NULL) return;
    free(c->L); free(c->off); free(c->np); free(c->subj); free(c->group);
    memset(c, 0, sizeof(*c));
}
