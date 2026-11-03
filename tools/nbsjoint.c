/* nbsjoint -- fit ONE exponent vector across several viewing conditions.
 *
 * Strict C89.
 *
 * WHY
 *
 * CID:IQ shows the same 690 distorted images at 50 cm and at 100 cm and
 * collects a separate MOS for each.  The mean opinion score rises by 0.61
 * when the observer steps back, and 392 of the 690 pairs move by more than
 * half a point.  SSIM and MS-SSIM cannot represent any of that: they take
 * no viewing geometry, so they emit one number for two different answers.
 *
 * NB-SSIM does take viewing geometry -- it is the only thing that changes
 * between the two conditions here, through the cycles-per-degree to
 * cycles-per-pixel mapping of the filterbank.  If a *single* exponent
 * vector, fitted once, predicts both conditions, then the bank is carrying
 * the viewing-distance dependence.  If instead each condition needs its own
 * exponents, the fit is absorbing what the model should have explained.
 * Tuning each distance separately cannot tell those two apart.  This tool
 * exists so the question can be asked properly.
 *
 * HOW
 *
 * At 100 cm (ppd 63.2) the bank has six channels; at 50 cm (ppd 31.6) the
 * 16 c/deg channel is above the display Nyquist rate and is dropped, so the
 * bank has five.  The conditions therefore have different numbers of free
 * exponents -- 21 and 15.  They are not independent parameter sets: the
 * five surviving channels are the same five, at the same cycles per degree,
 * so w_ij for i,j <= 4 is shared and only the 16 c/deg row and column are
 * exclusive to the near condition.  The master vector is the largest
 * condition's; each condition indexes into it through a map built from
 * nbs_w_index at its own band count.  A dropped channel simply contributes
 * no factor to the product, which is the physically honest reading: the
 * display cannot show it, so it cannot carry a distortion.
 *
 * COMBINING
 *
 *   mean    average of the per-condition accuracies.  Each condition gets
 *           its own logistic fit, so this asks only that the metric rank
 *           well within each condition.
 *   pooled  one logistic fit over both conditions at once.  This is the
 *           strict test: the same NB-SSIM value must mean the same quality
 *           whether it came from 50 cm or 100 cm.  SSIM and MS-SSIM are
 *           structurally incapable of scoring well here.
 *   min     the worse of the two, for a fit that cannot trade one viewing
 *           distance away against the other.
 */

#include "nbs.h"
#include "nbs_metrics.h"
#include "nbs_stats.h"
#include "nbs_cache.h"
#include "nbs_cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_COND 8

typedef struct {
    nbs_logxi_cache c;       /* full tensor, used for every reported number */
    nbs_logxi_cache s;       /* physically compacted subsample, for the search */
    int    have_s;
    int    map[NBS_MAX_W];   /* condition weight index -> master index */
    char   label[64];
} condition;

/* Sample subsampling used *only* while searching.  Scoring one entry costs
 * np*nw multiply-adds and the simplex wants tens of thousands of
 * evaluations; at stride 16 over 800x800 that is 2500 samples per image and
 * a full CV sweep runs for the better part of a day.  Every reported number
 * is recomputed over all samples. */
static int    g_eval_stride = 1;

static double g_wmax = 20.0;
static int    g_allow_neg = 0;
static double g_l2 = 0.0;

/* The one place the parameter vector becomes usable exponents.  Objective,
 * reports and the saved file all pass through here, so what is written out
 * is exactly what was scored. */
static void clamp_w(const double *v, double *w, int nw)
{
    int k;
    for (k = 0; k < nw; ++k) {
        double t = v[k];
        if (!g_allow_neg && t < 0.0) t = 0.0;
        if (t >  g_wmax) t =  g_wmax;
        if (t < -g_wmax) t = -g_wmax;
        w[k] = t;
    }
}

typedef struct {
    condition *cond;
    int        ncond;
    int        nw;          /* master exponent count        */
    int        objective;   /* 0 plcc, 1 srcc, 2 mix        */
    int        combine;     /* 0 mean, 1 pooled, 2 min      */
    const int *use;         /* per-group mask, NULL for all */
    double    *scr;         /* scratch, total entries       */
    double    *sub;
    double    *pred;
    double     pfit[MAX_COND + 1][5];
    int        have_fit[MAX_COND + 1];
} objctx;

/* exp() dominates the search: one objective evaluation exponentiates once
 * per cached pixel per image, and the simplex wants tens of thousands of
 * evaluations.  Range reduction to 2^n * exp(r) with |r| <= ln2/2 and a
 * degree-7 Taylor is good to about 1e-13 here -- far tighter than the
 * float storage of the tensor itself -- and several times quicker than the
 * library call.  The exponent is a weighted sum of log-xi values, each of
 * which is negative and bounded below by log(xi_floor), so the scale factor
 * is always in range for the small table. */
static double fexp_pow2[64];
static int    fexp_ready = 0;

static void fexp_init(void)
{
    int i;
    for (i = 0; i < 64; ++i) fexp_pow2[i] = ldexp(1.0, i - 62);
    fexp_ready = 1;
}

static double fexp(double x)
{
    static const double LOG2E = 1.4426950408889634074;
    static const double LN2   = 0.69314718055994530942;
    double r, p;
    int n;

    if (x >= 0.0) return (x > 0.0) ? exp(x) : 1.0;   /* xi <= 1, so x <= 0 */
    n = (int)(x * LOG2E - 0.5);                       /* x < 0 */
    if (n <= -62) return 0.0;
    r = x - (double)n * LN2;
    /* |r| <= ln2 here; Taylor to degree 7 by Horner. */
    p = 1.0 + r * (1.0 + r * (0.5 + r * (1.0 / 6.0 + r * (1.0 / 24.0
        + r * (1.0 / 120.0 + r * (1.0 / 720.0 + r * (1.0 / 5040.0
        + r * (1.0 / 40320.0))))))));
    return p * fexp_pow2[n + 62];
}

/* Take every g_eval_stride-th cached pixel and store the result densely.
 * Striding through the full tensor leaves the search touching one row in
 * eight across 250 MB; compacting first drops the working set to something
 * that stays in cache and makes the inner loop sequential. */
static int cache_compact(const nbs_logxi_cache *src, nbs_logxi_cache *dst, int stride)
{
    long total = 0;
    int e, p, k;

    *dst = *src;
    dst->L = NULL; dst->off = NULL; dst->np = NULL;
    dst->subj = NULL; dst->group = NULL;
    dst->off   = (long *)  nbs_malloc((size_t)src->n * sizeof(long));
    dst->np    = (int *)   nbs_malloc((size_t)src->n * sizeof(int));
    dst->subj  = (double *)nbs_malloc((size_t)src->n * sizeof(double));
    dst->group = (int *)   nbs_malloc((size_t)src->n * sizeof(int));
    if (!dst->off || !dst->np || !dst->subj || !dst->group) return 1;
    for (e = 0; e < src->n; ++e) {
        int keep = (src->np[e] + stride - 1) / stride;
        dst->off[e] = total;
        dst->np[e]  = keep;
        dst->subj[e]  = src->subj[e];
        dst->group[e] = src->group[e];
        total += (long)keep * (long)src->nw;
    }
    dst->L = (float *)nbs_malloc((size_t)total * sizeof(float));
    if (dst->L == NULL) return 1;
    for (e = 0; e < src->n; ++e) {
        const float *si = src->L + src->off[e];
        float *di = dst->L + dst->off[e];
        int q = 0;
        for (p = 0; p < src->np[e]; p += stride) {
            for (k = 0; k < src->nw; ++k)
                di[(size_t)q * src->nw + k] = si[(size_t)p * src->nw + k];
            q++;
        }
    }
    return 0;
}

/* Score one condition's entries for a master weight vector. */
static int score_cond(const condition *cd, const double *wmaster,
                      const int *use, double *scr, double *sub)
{
    const nbs_logxi_cache *c =
        (g_eval_stride > 1 && cd->have_s) ? &cd->s : &cd->c;
    double wv[NBS_MAX_W];
    int    wi[NBS_MAX_W];
    int e, p, k, m = 0, nz = 0;

    if (!fexp_ready) fexp_init();

    /* A converged fit zeroes most of the 21 exponents, and a zero exponent
     * contributes nothing to the product.  Walking only the live ones cuts
     * the inner loop by two thirds late in the search, where most of the
     * evaluations happen. */
    for (k = 0; k < c->nw; ++k) {
        double t = wmaster[cd->map[k]];
        if (t != 0.0) { wv[nz] = t; wi[nz] = k; nz++; }
    }

    for (e = 0; e < c->n; ++e) {
        double acc = 0.0;
        const float *L;
        int used = 0;
        if (use != NULL && !use[c->group[e]]) continue;
        L = c->L + c->off[e];
        for (p = 0; p < c->np[e]; ++p) {
            const float *row = L + (size_t)p * (size_t)c->nw;
            double s = 0.0;
            for (k = 0; k < nz; ++k) s += wv[k] * (double)row[wi[k]];
            acc += fexp(s);
            used++;
        }
        scr[m] = acc / (double)(used ? used : 1);
        sub[m] = c->subj[e];
        m++;
    }
    return m;
}

/* Accuracy of one score/subject set: PLCC after a logistic fit, SRCC, or
 * the mean of the two.  `slot` selects a warm-start slot for the fit. */
static double accuracy(objctx *o, double *scr, double *sub, int m, int slot)
{
    double plcc = 0.0, srcc = 0.0;

    if (m < 6) return 0.0;
    if (o->objective != 0) srcc = fabs(nbs_spearman(scr, sub, m));
    if (o->objective == 1) return srcc;

    if (!o->have_fit[slot]) {
        nbs_logistic_fit_ex(scr, sub, m, o->pfit[slot], 0, 4, 4000);
        o->have_fit[slot] = 1;
    } else {
        nbs_logistic_fit_ex(scr, sub, m, o->pfit[slot], 1, 2, 400);
    }
    {
        int i;
        for (i = 0; i < m; ++i) o->pred[i] = nbs_logistic_eval(o->pfit[slot], scr[i]);
        plcc = fabs(nbs_pearson(o->pred, sub, m));
    }
    if (o->objective == 2) return 0.5 * plcc + 0.5 * srcc;
    return plcc;
}

static double combined(objctx *o, const double *wmaster, double *per_cond)
{
    int i, m, tot = 0;
    double sum = 0.0, worst = 1e30;

    for (i = 0; i < o->ncond; ++i) {
        m = score_cond(&o->cond[i], wmaster, o->use, o->scr + tot, o->sub + tot);
        if (m < 6) return 0.0;
        if (o->combine != 1) {
            double a = accuracy(o, o->scr + tot, o->sub + tot, m, i);
            if (per_cond) per_cond[i] = a;
            sum += a;
            if (a < worst) worst = a;
        }
        tot += m;
    }
    if (o->combine == 1) return accuracy(o, o->scr, o->sub, tot, MAX_COND);
    if (o->combine == 2) return worst;
    return sum / (double)o->ncond;
}

static double objective_fn(const double *v, void *ctx)
{
    objctx *o = (objctx *)ctx;
    double w[NBS_MAX_W];
    double val, pen = 0.0;
    int k, nw = o->nw;      /* always the MASTER width.  Deriving it from
                             * the conditions present breaks the
                             * single-condition comparison below, where the
                             * near-distance bank has five channels but its
                             * exponents still live in the six-channel
                             * master matrix -- clamping only the first
                             * fifteen entries would leave the rest of the
                             * vector unwritten. */

    clamp_w(v, w, nw);
    for (k = 0; k < nw; ++k) pen += w[k] * w[k];
    pen *= g_l2;

    val = combined(o, w, NULL);
    if (!(val == val)) return 1e9;
    return -val + pen;
}

static unsigned long g_rng = 987654321UL;
static double rnd01(void)
{
    g_rng = g_rng * 1103515245UL + 12345UL;
    return (double)((g_rng >> 16) & 0x7fffUL) / 32767.0;
}

/* Nelder-Mead parks on a sparse local optimum in 21 dimensions; basin
 * hopping around it -- perturb, re-run the simplex, keep only improvements
 * -- escapes those parks without needing derivatives of an objective whose
 * inner logistic fit has none. */
static double optimise(double *w, int nw, objctx *oc, int restarts, int maxit,
                       int basin, double spread)
{
    double step[NBS_MAX_W], best[NBS_MAX_W], cand[NBS_MAX_W];
    double fbest, f;
    int rs, b, k;

    memcpy(best, w, (size_t)nw * sizeof(double));
    for (rs = 0; rs < restarts; ++rs) {
        for (k = 0; k < nw; ++k)
            step[k] = (0.4 * (fabs(best[k]) + 0.05)) * pow(0.4, (double)rs);
        nbs_nelder_mead(nw, best, step, objective_fn, oc, maxit, 1e-10);
    }
    fbest = objective_fn(best, oc);

    for (b = 0; b < basin; ++b) {
        memcpy(cand, best, (size_t)nw * sizeof(double));
        for (k = 0; k < nw; ++k) {
            double jitter = spread * (2.0 * rnd01() - 1.0);
            cand[k] = cand[k] * (1.0 + jitter) + spread * 0.5 * rnd01();
        }
        for (rs = 0; rs < restarts; ++rs) {
            for (k = 0; k < nw; ++k)
                step[k] = (0.4 * (fabs(cand[k]) + 0.05)) * pow(0.4, (double)rs);
            nbs_nelder_mead(nw, cand, step, objective_fn, oc, maxit, 1e-10);
        }
        f = objective_fn(cand, oc);
        if (f < fbest) { fbest = f; memcpy(best, cand, (size_t)nw * sizeof(double)); }
    }
    clamp_w(best, w, nw);
    return fbest;
}

static void report(objctx *o, const double *v, const int *use, const char *tag)
{
    double w[NBS_MAX_W];
    int i, m, tot = 0, nw = o->nw, saved_stride = g_eval_stride;
    nbs_fitreport r;
    const int *saved = o->use;

    g_eval_stride = 1;

    clamp_w(v, w, nw);
    o->use = use;

    printf("  %s\n", tag);
    for (i = 0; i < o->ncond; ++i) {
        m = score_cond(&o->cond[i], w, use, o->scr + tot, o->sub + tot);
        if (m >= 6) {
            nbs_evaluate(o->scr + tot, o->sub + tot, m, &r);
            printf("    %-14s n=%4d   PLCC %7.4f   SRCC %7.4f   KRCC %7.4f   RMSE %7.4f\n",
                   o->cond[i].label, m, fabs(r.plcc), fabs(r.srcc), fabs(r.krcc), r.rmse);
        }
        tot += m;
    }
    if (tot >= 6) {
        nbs_evaluate(o->scr, o->sub, tot, &r);
        printf("    %-14s n=%4d   PLCC %7.4f   SRCC %7.4f   KRCC %7.4f   RMSE %7.4f\n",
               "pooled", tot, fabs(r.plcc), fabs(r.srcc), fabs(r.krcc), r.rmse);
    }
    o->use = saved;
    g_eval_stride = saved_stride;
}

static void usage(void)
{
    printf("usage: nbsjoint [options] --cond=CACHE [--cond=CACHE ...]\n\n");
    printf("Fits one exponent vector across several viewing conditions, each\n");
    printf("given as a log-xi cache built by nbstune --cache-save at its own\n");
    printf("ppd.  Conditions may have different band counts; channels above a\n");
    printf("condition's display Nyquist rate simply do not participate.\n\n");
    printf("  --cond=PATH        add a condition (repeatable, up to %d)\n", MAX_COND);
    printf("  --label=NAME       name for the next --cond (default: its ppd)\n");
    printf("  --combine=mean|pooled|min   how the per-condition accuracies\n");
    printf("                     become one number (default mean)\n");
    printf("  --objective=plcc|srcc|mix   accuracy measure (default plcc)\n");
    printf("  --folds=K          content-disjoint cross-validation folds,\n");
    printf("                     split by reference image (default 5)\n");
    printf("  --restarts=N       Nelder-Mead restarts (default 4)\n");
    printf("  --basin=N          basin-hopping rounds (default 0)\n");
    printf("  --spread=F         basin perturbation scale (default 0.6)\n");
    printf("  --maxit=N          iterations per restart (default 6000)\n");
    printf("  --eval-stride=N    subsample cached pixels during the search\n");
    printf("                     only; reports always use every sample\n");
    printf("  --allow-negative   permit negative exponents\n");
    printf("  --wmax=F           hard bound on any exponent (default 20)\n");
    printf("  --l2=F             ridge penalty on the exponent vector\n");
    printf("  --out=PATH         write the fitted exponents here\n");
    printf("  --per-cond         also fit each condition on its own, for\n");
    printf("                     comparison against the shared fit\n");
    printf("  --dump-scores=PATH write per-pair scores for the final vector,\n");
    printf("                     as condition,index,group,subjective,score\n\n");
    nbs_opts_print_help();
}

int main(int argc, char **argv)
{
    nbs_opts o;
    condition cond[MAX_COND];
    objctx oc;
    const char *outp = NULL, *pending_label = NULL;
    const char *dumpp = NULL, *loadw = NULL;
    int ncond = 0, i, k, rc = 1, folds = 5, objective = 0, combine = 0;
    int master_sym = 1;
    int restarts = 4, maxit = 6000, basin = 0, allow_neg = 0, per_cond = 0;
    double spread = 0.6;
    double w[NBS_MAX_W];
    int master_nb = 0, master_nw = 0, ngroup = 0, ntot = 0;

    memset(cond, 0, sizeof(cond));
    nbs_opts_default(&o);

    for (i = 1; i < argc; ++i) {
        int r;
        if (strcmp(argv[i], "--help") == 0) { usage(); return 0; }
        if (strncmp(argv[i], "--cond=", 7) == 0) {
            if (ncond >= MAX_COND) { fprintf(stderr, "too many conditions\n"); return 2; }
            if (nbs_cache_load(&cond[ncond].c, argv[i] + 7) != 0) {
                fprintf(stderr, "error: %s: %s\n", argv[i] + 7, nbs_last_error());
                return 2;
            }
            if (pending_label != NULL) {
                strncpy(cond[ncond].label, pending_label, sizeof(cond[0].label) - 1);
                pending_label = NULL;
            } else {
                sprintf(cond[ncond].label, "ppd %.4g", cond[ncond].c.ppd);
            }
            ncond++;
            continue;
        }
        if (strncmp(argv[i], "--label=", 8) == 0) { pending_label = argv[i] + 8; continue; }
        if (strncmp(argv[i], "--out=", 6) == 0) { outp = argv[i] + 6; continue; }
        if (strncmp(argv[i], "--dump-scores=", 14) == 0) { dumpp = argv[i] + 14; continue; }
        if (strncmp(argv[i], "--load-weights=", 15) == 0) { loadw = argv[i] + 15; continue; }
        if (strncmp(argv[i], "--folds=", 8) == 0) { folds = atoi(argv[i] + 8); continue; }
        if (strncmp(argv[i], "--restarts=", 11) == 0) { restarts = atoi(argv[i] + 11); continue; }
        if (strncmp(argv[i], "--basin=", 8) == 0) { basin = atoi(argv[i] + 8); continue; }
        if (strncmp(argv[i], "--spread=", 9) == 0) { spread = atof(argv[i] + 9); continue; }
        if (strncmp(argv[i], "--maxit=", 8) == 0) { maxit = atoi(argv[i] + 8); continue; }
        if (strncmp(argv[i], "--eval-stride=", 14) == 0) {
            g_eval_stride = atoi(argv[i] + 14);
            if (g_eval_stride < 1) g_eval_stride = 1;
            continue;
        }
        if (strcmp(argv[i], "--allow-negative") == 0) { allow_neg = 1; continue; }
        if (strcmp(argv[i], "--per-cond") == 0) { per_cond = 1; continue; }
        if (strncmp(argv[i], "--wmax=", 7) == 0) { g_wmax = atof(argv[i] + 7); continue; }
        if (strncmp(argv[i], "--l2=", 5) == 0) { g_l2 = atof(argv[i] + 5); continue; }
        if (strncmp(argv[i], "--objective=", 12) == 0) {
            const char *v = argv[i] + 12;
            if (strcmp(v, "plcc") == 0) objective = 0;
            else if (strcmp(v, "srcc") == 0) objective = 1;
            else if (strcmp(v, "mix") == 0) objective = 2;
            else { fprintf(stderr, "bad --objective\n"); return 2; }
            continue;
        }
        if (strncmp(argv[i], "--combine=", 10) == 0) {
            const char *v = argv[i] + 10;
            if (strcmp(v, "mean") == 0) combine = 0;
            else if (strcmp(v, "pooled") == 0) combine = 1;
            else if (strcmp(v, "min") == 0) combine = 2;
            else { fprintf(stderr, "bad --combine\n"); return 2; }
            continue;
        }
        r = nbs_opts_parse_one(&o, argv[i]);
        if (r < 0) { fprintf(stderr, "error: %s\n", nbs_last_error()); return 2; }
        if (r == 1) continue;
        fprintf(stderr, "unknown argument: %s\n", argv[i]);
        return 2;
    }
    if (ncond < 1) { usage(); return 2; }
    if (!(g_wmax > 0.0)) g_wmax = 20.0;
    g_allow_neg = allow_neg;

    /* The master band set is the richest condition's.  Every other
     * condition must be a prefix of it at the same cycles per degree --
     * that is what makes a shared exponent meaningful. */
    for (i = 0; i < ncond; ++i)
        if (cond[i].c.nbands > master_nb) master_nb = cond[i].c.nbands;
    master_sym = cond[0].c.w_sym;
    for (i = 1; i < ncond; ++i) {
        if (cond[i].c.w_sym != master_sym) {
            fprintf(stderr, "conditions disagree on exponent symmetry\n");
            goto done;
        }
    }
    master_nw = nbs_w_count_ex(master_nb, master_sym);

    for (i = 0; i < ncond; ++i) {
        int a, b, t = 0, nb = cond[i].c.nbands;
        for (k = 0; k < nb; ++k) {
            int ref = 0;
            for (a = 0; a < ncond; ++a)
                if (cond[a].c.nbands == master_nb) { ref = a; break; }
            if (fabs(cond[i].c.f[k] - cond[ref].c.f[k]) > 1e-9) {
                fprintf(stderr,
                        "condition %d is not a prefix of the master bank: "
                        "band %d is %.6g cpd, master has %.6g cpd\n",
                        i, k, cond[i].c.f[k], cond[ref].c.f[k]);
                goto done;
            }
        }
        for (a = 0; a < nb; ++a)
            for (b = (master_sym ? a : 0); b < nb; ++b)
                cond[i].map[t++] = nbs_w_index_ex(a, b, master_nb, master_sym);
        if (t != cond[i].c.nw) {
            fprintf(stderr, "condition %d: weight count mismatch\n", i);
            goto done;
        }
        for (k = 0; k < cond[i].c.n; ++k) {
            if (cond[i].c.group[k] + 1 > ngroup) ngroup = cond[i].c.group[k] + 1;
        }
        ntot += cond[i].c.n;
    }

    printf("conditions:\n");
    for (i = 0; i < ncond; ++i)
        printf("  %-14s ppd %8.4g   %d bands   %2d exponents   %4d pairs   stride %d\n",
               cond[i].label, cond[i].c.ppd, cond[i].c.nbands,
               cond[i].c.nw, cond[i].c.n, cond[i].c.stride);
    printf("master bank : %d bands, %d shared exponents (%s)\n",
           master_nb, master_nw,
           master_sym ? "w_ij == w_ji" : "all free");
    printf("bands       : [");
    for (i = 0; i < master_nb; ++i) {
        int ref = 0, a;
        for (a = 0; a < ncond; ++a) if (cond[a].c.nbands == master_nb) { ref = a; break; }
        printf("%s%.4g", (i ? " " : ""), cond[ref].c.f[i]);
    }
    printf("] cpd\n");
    printf("objective   : %s of per-condition %s\n",
           combine == 0 ? "mean" : combine == 1 ? "pooled fit" : "worst",
           objective == 0 ? "PLCC" : objective == 1 ? "SRCC" : "PLCC/SRCC mix");
    printf("\n");
    fflush(stdout);

    oc.cond = cond; oc.ncond = ncond; oc.nw = master_nw;
    oc.objective = objective; oc.combine = combine;
    oc.use = NULL;
    oc.scr = (double *)nbs_malloc((size_t)ntot * sizeof(double));
    oc.sub = (double *)nbs_malloc((size_t)ntot * sizeof(double));
    oc.pred = (double *)nbs_malloc((size_t)ntot * sizeof(double));
    if (!oc.scr || !oc.sub || !oc.pred) goto done;

    if (g_eval_stride > 1) {
        for (i = 0; i < ncond; ++i) {
            if (cache_compact(&cond[i].c, &cond[i].s, g_eval_stride) != 0) {
                fprintf(stderr, "out of memory compacting the search cache\n");
                goto done;
            }
            cond[i].have_s = 1;
        }
        printf("search tensor subsampled by %d\n\n", g_eval_stride);
    }
    for (i = 0; i <= MAX_COND; ++i) oc.have_fit[i] = 0;

    /* Uniform start, total exponent 2 as in the original SSIM. */
    for (k = 0; k < master_nw; ++k) w[k] = 2.0 / (double)(master_nb * master_nb);
    o.nbs.w_sym = master_sym;

    printf("baseline (uniform exponents):\n");
    report(&oc, w, NULL, "all data");
    printf("\n");
    fflush(stdout);

    if (folds > 1) {
        double sp[MAX_COND], ss[MAX_COND], sk[MAX_COND], se[MAX_COND];
        double pp = 0.0, ps = 0.0, pk = 0.0, pe = 0.0;
        int f, nf = 0, search_stride;
        int *tr = (int *)nbs_malloc((size_t)ngroup * sizeof(int));
        int *te = (int *)nbs_malloc((size_t)ngroup * sizeof(int));
        if (!tr || !te) { free(tr); free(te); goto done; }
        for (i = 0; i < ncond; ++i) { sp[i] = 0.0; ss[i] = 0.0; sk[i] = 0.0; se[i] = 0.0; }
        printf("%d-fold cross-validation, folds split by reference image,\n"
               "one shared exponent vector fitted per fold:\n", folds);
        for (f = 0; f < folds; ++f) {
            double wf[NBS_MAX_W];
            int g, m, tot = 0;
            nbs_fitreport r;
            for (g = 0; g < ngroup; ++g) { te[g] = ((g % folds) == f); tr[g] = !te[g]; }
            memcpy(wf, w, (size_t)master_nw * sizeof(double));
            oc.use = tr;
            for (i = 0; i <= MAX_COND; ++i) oc.have_fit[i] = 0;
            g_rng = 987654321UL + (unsigned long)f * 7919UL;
            optimise(wf, master_nw, &oc, restarts, maxit, basin, spread);

            oc.use = te;
            search_stride = g_eval_stride;
            g_eval_stride = 1;          /* held-out numbers use every sample */
            printf("  fold %d:\n", f);
            for (i = 0; i < ncond; ++i) {
                m = score_cond(&cond[i], wf, te, oc.scr + tot, oc.sub + tot);
                if (m >= 6) {
                    nbs_evaluate(oc.scr + tot, oc.sub + tot, m, &r);
                    printf("    %-14s n=%4d   PLCC %7.4f   SRCC %7.4f   "
                           "KRCC %7.4f   RMSE %7.4f\n",
                           cond[i].label, m, fabs(r.plcc), fabs(r.srcc),
                           fabs(r.krcc), r.rmse);
                    sp[i] += fabs(r.plcc); ss[i] += fabs(r.srcc);
                    sk[i] += fabs(r.krcc); se[i] += r.rmse;
                }
                tot += m;
            }
            if (tot >= 6) {
                nbs_evaluate(oc.scr, oc.sub, tot, &r);
                printf("    %-14s n=%4d   PLCC %7.4f   SRCC %7.4f   "
                       "KRCC %7.4f   RMSE %7.4f\n",
                       "pooled", tot, fabs(r.plcc), fabs(r.srcc),
                       fabs(r.krcc), r.rmse);
                pp += fabs(r.plcc); ps += fabs(r.srcc);
                pk += fabs(r.krcc); pe += r.rmse;
            }
            g_eval_stride = search_stride;
            nf++;
            fflush(stdout);
        }
        if (nf > 0) {
            printf("  mean held-out:\n");
            for (i = 0; i < ncond; ++i)
                printf("    %-14s          PLCC %7.4f   SRCC %7.4f   "
                       "KRCC %7.4f   RMSE %7.4f\n",
                       cond[i].label, sp[i] / nf, ss[i] / nf, sk[i] / nf, se[i] / nf);
            printf("    %-14s          PLCC %7.4f   SRCC %7.4f   "
                   "KRCC %7.4f   RMSE %7.4f\n",
                   "pooled", pp / nf, ps / nf, pk / nf, pe / nf);
        }
        free(tr); free(te);
        oc.use = NULL;
        printf("\n");
    }

    if (loadw != NULL) {
        double tmp[NBS_MAX_W];
        int a, b, t = 0;
        if (nbs_load_weights_ex(loadw, tmp, master_nb, master_sym) != 0) {
            fprintf(stderr, "error: %s\n", nbs_last_error());
            goto done;
        }
        for (a = 0; a < master_nb; ++a)
            for (b = (master_sym ? a : 0); b < master_nb; ++b)
                w[t++] = tmp[nbs_w_index_ex(a, b, master_nb, master_sym)];
        printf("exponents loaded from %s, no search:\n", loadw);
    } else {
        printf("final shared fit on all data:\n");
        oc.use = NULL;
        for (i = 0; i <= MAX_COND; ++i) oc.have_fit[i] = 0;
        g_rng = 987654321UL;
        optimise(w, master_nw, &oc, restarts, maxit, basin, spread);
    }
    report(&oc, w, NULL, "all data (shared fit)");
    printf("\n");

    if (dumpp != NULL) {
        FILE *fp = fopen(dumpp, "w");
        if (fp == NULL) {
            fprintf(stderr, "cannot write %s\n", dumpp);
        } else {
            int saved_stride = g_eval_stride, e;
            g_eval_stride = 1;
            fprintf(fp, "condition,index,group,subjective,score\n");
            for (i = 0; i < ncond; ++i) {
                int m = score_cond(&cond[i], w, NULL, oc.scr, oc.sub);
                (void)m;
                for (e = 0; e < cond[i].c.n; ++e)
                    fprintf(fp, "%s,%d,%d,%.6f,%.10f\n", cond[i].label, e,
                            cond[i].c.group[e], cond[i].c.subj[e], oc.scr[e]);
            }
            g_eval_stride = saved_stride;
            fclose(fp);
            printf("per-pair scores written to %s\n", dumpp);
        }
    }

    /* How much is given up by insisting on one vector?  Fit each condition
     * alone and show the gap. */
    if (per_cond) {
        printf("for comparison, each condition fitted on its own:\n");
        for (i = 0; i < ncond; ++i) {
            objctx s;
            double wi[NBS_MAX_W];
            int m;
            nbs_fitreport r;
            s = oc; s.cond = &cond[i]; s.ncond = 1; s.combine = 0; s.use = NULL;
            for (k = 0; k <= MAX_COND; ++k) s.have_fit[k] = 0;
            for (k = 0; k < master_nw; ++k) wi[k] = 2.0 / (double)(master_nb * master_nb);
            g_rng = 987654321UL;
            optimise(wi, master_nw, &s, restarts, maxit, basin, spread);
            m = score_cond(&cond[i], wi, NULL, oc.scr, oc.sub);
            if (m >= 6) {
                nbs_evaluate(oc.scr, oc.sub, m, &r);
                printf("    %-14s n=%4d   PLCC %7.4f   SRCC %7.4f  (own exponents)\n",
                       cond[i].label, m, fabs(r.plcc), fabs(r.srcc));
            }
        }
        printf("\n");
    }

    printf("fitted exponents (%s):\n",
           master_sym ? "upper triangle, w_ij == w_ji"
                      : "row i = horizontal channel, column j = vertical");
    {
        int a, b, t = 0;
        for (a = 0; a < master_nb; ++a) {
            printf("  ");
            if (master_sym) for (b = 0; b < a; ++b) printf("%9s", "");
            for (b = (master_sym ? a : 0); b < master_nb; ++b)
                printf("%9.5f", w[t++]);
            printf("\n");
        }
    }
    if (!master_sym) {
        /* How well does the fit respect the isotropy the symmetric variant
         * assumes?  Report it rather than assume it. */
        double num = 0.0, den = 0.0;
        int a, b;
        for (a = 0; a < master_nb; ++a)
            for (b = a + 1; b < master_nb; ++b) {
                double p = w[nbs_w_index_ex(a, b, master_nb, 0)];
                double q = w[nbs_w_index_ex(b, a, master_nb, 0)];
                num += fabs(p - q);
                den += 0.5 * (fabs(p) + fabs(q));
            }
        printf("\noff-diagonal asymmetry: mean |w_ij - w_ji| / mean |w| = %.3f\n",
               (den > 0.0) ? num / den : 0.0);
        {
            double dsum = 0.0, osum = 0.0;
            int nd = 0, no = 0;
            for (a = 0; a < master_nb; ++a)
                for (b = 0; b < master_nb; ++b) {
                    double v = w[nbs_w_index_ex(a, b, master_nb, 0)];
                    if (a == b) { dsum += v; nd++; } else { osum += v; no++; }
                }
            printf("mean exponent on the oblique diagonal (i==j): %.5f over %d\n",
                   dsum / (nd ? nd : 1), nd);
            printf("mean exponent off  the diagonal  (i!=j)     : %.5f over %d\n",
                   osum / (no ? no : 1), no);
        }
    }

    if (outp != NULL) {
        int a, b, t = 0;
        o.nbs.bank.nbands = master_nb;
        o.nbs.w_sym = master_sym;
        for (a = 0; a < master_nb; ++a)
            for (b = (master_sym ? a : 0); b < master_nb; ++b)
                o.nbs.w[nbs_w_index_ex(a, b, master_nb, master_sym)] = w[t++];
        if (nbs_save_weights(outp, &o.nbs) != 0)
            fprintf(stderr, "error: %s\n", nbs_last_error());
        else
            printf("\nwritten to %s\n", outp);
    }
    rc = 0;

done:
    free(oc.scr); free(oc.sub); free(oc.pred);
    for (i = 0; i < ncond; ++i) {
        nbs_cache_free(&cond[i].c);
        if (cond[i].have_s) nbs_cache_free(&cond[i].s);
    }
    return rc;
}
