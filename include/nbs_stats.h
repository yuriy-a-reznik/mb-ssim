/* nbs_stats.h -- correlation measures, VQEG logistic fitting, Nelder-Mead.
 * Strict C89.
 */

#ifndef NBS_STATS_H
#define NBS_STATS_H

#ifdef __cplusplus
extern "C" {
#endif

double nbs_pearson (const double *a, const double *b, int n);
double nbs_spearman(const double *a, const double *b, int n);
double nbs_kendall (const double *a, const double *b, int n);   /* tau-b */

/* VQEG five-parameter logistic:
 *   f(x) = p[0] * (0.5 - 1/(1+exp(p[1]*(x-p[2])))) + p[3]*x + p[4]
 * Fits p to map objective scores `x` onto subjective scores `y`.
 * Returns 0 on success; `p` receives the fitted parameters. */
int nbs_logistic_fit(const double *x, const double *y, int n, double *p);

/* Warm-startable variant: `p` is used as the starting point when
 * `warm` is non-zero, and `passes`/`maxit` bound the search.  The default fit
 * above is nbs_logistic_fit_ex(x, y, n, p, 0, 4, 4000). */
int nbs_logistic_fit_ex(const double *x, const double *y, int n, double *p,
                        int warm, int passes, int maxit);

double nbs_logistic_eval(const double *p, double x);

/* Full report: PLCC/RMSE/MAE are computed after the logistic mapping. */
typedef struct {
    double plcc;
    double srcc;
    double krcc;
    double rmse;
    double mae;
    int    n;
    int    fit_ok;
    double p[5];
} nbs_fitreport;

int nbs_evaluate(const double *obj, const double *subj, int n, nbs_fitreport *rep);

/* Generic Nelder-Mead minimiser.
 *   nd    -- dimensionality
 *   x     -- in: start point, out: best point found
 *   step  -- initial simplex edge lengths (one per dimension)
 *   fn    -- objective
 *   ctx   -- opaque user pointer
 *   maxit -- iteration cap
 *   tol   -- relative convergence tolerance on the function spread
 * Returns the best function value. */
double nbs_nelder_mead(int nd, double *x, const double *step,
                       double (*fn)(const double *, void *), void *ctx,
                       int maxit, double tol);

#ifdef __cplusplus
}
#endif

#endif /* NBS_STATS_H */
