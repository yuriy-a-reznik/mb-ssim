/* nbs_cli.h -- shared command line option handling.  Strict C89. */

#ifndef NBS_CLI_H
#define NBS_CLI_H

#include "nbs.h"
#include "nbs_metrics.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    nbs_luma_mode  luma;
    nbs_ssim_cfg   ssim;
    nbs_msssim_cfg msssim;
    nbs_nbssim_cfg nbs;
    int            ppd_set;      /* 1 if --ppd was given explicitly */
} nbs_opts;

void nbs_opts_default(nbs_opts *o);

/* Returns 1 if `arg` was consumed, 0 if unrecognised, -1 on error. */
int  nbs_opts_parse_one(nbs_opts *o, const char *arg);

/* Re-derive dependent state (band cut-offs) after ppd or peaks changed. */
int  nbs_opts_finalise(nbs_opts *o);

void nbs_opts_print_help(void);

/* Load a weight vector from a text file (whitespace separated). */
int  nbs_load_weights(const char *path, double *w, int nbands);
int  nbs_load_weights_ex(const char *path, double *w, int nbands, int sym);
int  nbs_save_weights(const char *path, const nbs_nbssim_cfg *cfg);

#ifdef __cplusplus
}
#endif

#endif /* NBS_CLI_H */
