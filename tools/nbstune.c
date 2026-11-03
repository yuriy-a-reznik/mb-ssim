/* nbstune -- fit the N-band SSIM exponents (and optionally the filterbank
 * cut-offs) to maximise agreement with subjective scores.
 *
 * Strategy
 * --------
 * The metric is
 *
 *     S(x,y) = mean_p exp( sum_k w_k * L_k(p) ),   L_k(p) = log xi_k(p)
 *
 * summed over the N(N+1)/2 symmetric band pairs k.  For a fixed filterbank
 * the maps L_k depend only on the images, so they are computed once and
 * cached on a subsampled pixel grid.  Every subsequent weight evaluation is
 * then a handful of dot products, which makes Nelder-Mead over 21 exponents
 * cheap.  Changing a cut-off invalidates the cache, so cut-offs are handled
 * by a coarse outer pattern search with far fewer evaluations.
 *
 * Overfitting is measured, not assumed: folds are split by reference image,
 * so no source content is shared between train and test.
 *
 * Strict C89.
 */

#include "nbs.h"
#include "nbs_metrics.h"
#include "nbs_data.h"
#include "nbs_stats.h"
#include "nbs_cli.h"
#include "nbs_cache.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static double g_wmax = 20.0;   /* hard bound on any single exponent */
static int    g_allow_neg = 0;
static double g_l2 = 0.0;      /* ridge penalty on the exponent vector */

/* The one place the parameter vector is turned into usable exponents.  The
 * objective, the reports and the saved file all go through this, so what is
 * written out is exactly what was scored. */
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
    const nbs_logxi_cache *c;
    const int   *use;      /* 1 if entry participates, NULL for all */
    int          objective;/* 0 = plcc, 1 = srcc, 2 = mix           */
    int          allow_neg;
    double      *scr;      /* scratch, size c->n */
    double      *sub;      /* scratch, size c->n */
    double       pfit[5];  /* warm start for the logistic fit */
    int          have_fit;
} objctx;

/* score all selected entries for a weight vector */
static int score_all(const nbs_logxi_cache *c, const double *w, const int *use,
                     double *scr, double *sub, const double *subj)
{
    int e, p, k, m = 0;
    for (e = 0; e < c->n; ++e) {
        double acc = 0.0;
        const float *L;
        if (use != NULL && !use[e]) continue;
        L = c->L + c->off[e];
        for (p = 0; p < c->np[e]; ++p) {
            double s = 0.0;
            const float *row = L + (size_t)p * (size_t)c->nw;
            for (k = 0; k < c->nw; ++k) s += w[k] * (double)row[k];
            acc += exp(s);
        }
        scr[m] = acc / (double)c->np[e];
        sub[m] = subj[e];
        m++;
    }
    return m;
}

static double objective_fn(const double *v, void *ctx)
{
    objctx *o = (objctx *)ctx;
    double w[NBS_MAX_W];
    int k, m;
    double val, pen = 0.0;

    clamp_w(v, w, o->c->nw);
    for (k = 0; k < o->c->nw; ++k) pen += w[k] * w[k];
    pen *= g_l2;
    m = score_all(o->c, w, o->use, o->scr, o->sub, o->c->subj);
    if (m < 6) return 1e9;

    if (o->objective == 1) {
        val = fabs(nbs_spearman(o->scr, o->sub, m));
    } else {
        /* PLCC after a warm-started logistic fit.  The search only needs a
         * good relative ordering, so a short refit from the previous
         * solution replaces the full multi-pass fit; the reported numbers
         * always use the full fit. */
        double plcc;
        int i;
        if (!o->have_fit) {
            nbs_logistic_fit_ex(o->scr, o->sub, m, o->pfit, 0, 4, 4000);
            o->have_fit = 1;
        } else {
            nbs_logistic_fit_ex(o->scr, o->sub, m, o->pfit, 1, 2, 400);
        }
        {
            double *pred = (double *)nbs_malloc((size_t)m * sizeof(double));
            if (pred == NULL) return 1e9;
            for (i = 0; i < m; ++i) pred[i] = nbs_logistic_eval(o->pfit, o->scr[i]);
            plcc = fabs(nbs_pearson(pred, o->sub, m));
            free(pred);
        }
        if (o->objective == 2)
            val = 0.5 * plcc + 0.5 * fabs(nbs_spearman(o->scr, o->sub, m));
        else
            val = plcc;
    }
    if (!(val == val)) return 1e9;
    return -val + pen;
}

static void eval_report(const nbs_logxi_cache *c, const double *v, const int *use,
                        double *scr, double *sub, const char *label)
{
    nbs_fitreport r;
    double w[NBS_MAX_W];
    int m;
    clamp_w(v, w, c->nw);
    m = score_all(c, w, use, scr, sub, c->subj);
    if (m < 6) { printf("  %-18s (too few samples)\n", label); return; }
    nbs_evaluate(scr, sub, m, &r);
    printf("  %-18s n=%4d   PLCC %7.4f   SRCC %7.4f   KRCC %7.4f   RMSE %8.4f\n",
           label, m, fabs(r.plcc), fabs(r.srcc), fabs(r.krcc), r.rmse);
}

/* Nelder-Mead degrades badly above roughly ten dimensions: with 21 exponents
 * it reliably parks on a sparse local optimum.  Basin hopping wraps it --
 * perturb the incumbent, re-run the simplex, keep the result only if it
 * improved -- which is enough to escape those parks without needing
 * derivatives of an objective whose inner logistic fit has none. */
static unsigned long g_rng = 987654321UL;

static double rnd01(void)
{
    g_rng = g_rng * 1103515245UL + 12345UL;
    return (double)((g_rng >> 16) & 0x7fffUL) / 32767.0;
}

static double optimise(double *w, int nw, objctx *oc, int restarts, int maxit,
                       int basin, double spread)
{
    double step[NBS_MAX_W];
    double best[NBS_MAX_W], cand[NBS_MAX_W];
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
            /* Perturb multiplicatively where a weight is alive, additively
             * where it has collapsed to zero -- otherwise a zeroed exponent
             * can never come back. */
            double jitter = spread * (2.0 * rnd01() - 1.0);
            cand[k] = cand[k] * (1.0 + jitter) + spread * 0.5 * rnd01();
        }
        for (rs = 0; rs < restarts; ++rs) {
            for (k = 0; k < nw; ++k)
                step[k] = (0.4 * (fabs(cand[k]) + 0.05)) * pow(0.4, (double)rs);
            nbs_nelder_mead(nw, cand, step, objective_fn, oc, maxit, 1e-10);
        }
        f = objective_fn(cand, oc);
        if (f < fbest) {
            fbest = f;
            memcpy(best, cand, (size_t)nw * sizeof(double));
        }
    }
    clamp_w(best, w, nw);
    return fbest;
}

static void usage(void)
{
    printf("usage: nbstune [options] INDEX.csv\n\n");
    printf("  --out=PATH          write the fitted exponents here\n");
    printf("  --stride=N          pixel subsampling for the cache (default 16)\n");
    printf("  --folds=K           content-disjoint cross-validation folds\n");
    printf("                      (default 5; 0 or 1 disables CV)\n");
    printf("  --objective=plcc|srcc|mix   default plcc\n");
    printf("  --restarts=N        Nelder-Mead restarts (default 4)\n");
    printf("  --basin=N           basin-hopping rounds around the simplex\n");
    printf("                      result (default 0); escapes local optima\n");
    printf("  --spread=F          basin perturbation scale (default 0.6)\n");
    printf("  --maxit=N           iterations per restart (default 6000)\n");
    printf("  --allow-negative    permit negative exponents\n");
    printf("  --wmax=F            hard bound on any exponent (default 20)\n");
    printf("  --l2=F              ridge penalty on the exponent vector\n");
    printf("  --cache-save=PATH   write the log-xi cache for reuse\n");
    printf("  --cache-load=PATH   reuse a cache instead of rescanning images\n");

    printf("  --subset=N          use only every N-th dataset entry\n");
    printf("  --progress          progress to stderr\n\n");
    nbs_opts_print_help();
}

int main(int argc, char **argv)
{
    nbs_opts o;
    nbs_dataset ds;
    nbs_logxi_cache c;
    objctx oc;
    const char *index = NULL, *outp = NULL, *csave = NULL, *cload = NULL;
    int i, k, rc = 1, stride = 16, folds = 5, objective = 0;
    int restarts = 4, maxit = 6000, allow_neg = 0, subset = 1, progress = 0;
    int basin = 0;
    double spread = 0.6;
    double *scr = NULL, *sub = NULL;
    int *use = NULL;
    double w[NBS_MAX_W];
    double step[NBS_MAX_W];
    long total = 0;
    int nw;

    memset(&ds, 0, sizeof(ds));
    memset(&c, 0, sizeof(c));
    nbs_opts_default(&o);

    for (i = 1; i < argc; ++i) {
        int r;
        if (strcmp(argv[i], "--help") == 0) { usage(); return 0; }
        if (strncmp(argv[i], "--out=", 6) == 0) { outp = argv[i] + 6; continue; }
        if (strncmp(argv[i], "--stride=", 9) == 0) { stride = atoi(argv[i] + 9); continue; }
        if (strncmp(argv[i], "--folds=", 8) == 0) { folds = atoi(argv[i] + 8); continue; }
        if (strncmp(argv[i], "--restarts=", 11) == 0) { restarts = atoi(argv[i] + 11); continue; }
        if (strncmp(argv[i], "--basin=", 8) == 0) { basin = atoi(argv[i] + 8); continue; }
        if (strncmp(argv[i], "--spread=", 9) == 0) { spread = atof(argv[i] + 9); continue; }
        if (strncmp(argv[i], "--maxit=", 8) == 0) { maxit = atoi(argv[i] + 8); continue; }
        if (strncmp(argv[i], "--subset=", 9) == 0) { subset = atoi(argv[i] + 9); continue; }
        if (strcmp(argv[i], "--allow-negative") == 0) { allow_neg = 1; continue; }
        if (strncmp(argv[i], "--wmax=", 7) == 0) { g_wmax = atof(argv[i] + 7); continue; }
        if (strncmp(argv[i], "--l2=", 5) == 0) { g_l2 = atof(argv[i] + 5); continue; }
        if (strncmp(argv[i], "--cache-save=", 13) == 0) { csave = argv[i] + 13; continue; }
        if (strncmp(argv[i], "--cache-load=", 13) == 0) { cload = argv[i] + 13; continue; }
        if (strcmp(argv[i], "--progress") == 0) { progress = 1; continue; }
        if (strncmp(argv[i], "--objective=", 12) == 0) {
            const char *v = argv[i] + 12;
            if (strcmp(v, "plcc") == 0) objective = 0;
            else if (strcmp(v, "srcc") == 0) objective = 1;
            else if (strcmp(v, "mix") == 0) objective = 2;
            else { fprintf(stderr, "bad --objective\n"); return 2; }
            continue;
        }
        r = nbs_opts_parse_one(&o, argv[i]);
        if (r < 0) { fprintf(stderr, "error: %s\n", nbs_last_error()); return 2; }
        if (r == 1) continue;
        if (argv[i][0] == '-') { fprintf(stderr, "unknown option: %s\n", argv[i]); return 2; }
        index = argv[i];
    }
    if (index == NULL) { usage(); return 2; }
    if (stride < 1) stride = 1;
    if (subset < 1) subset = 1;
    if (!(g_wmax > 0.0)) g_wmax = 20.0;
    g_allow_neg = allow_neg;

    if (nbs_dataset_load(index, &ds) != 0) {
        fprintf(stderr, "error: %s\n", nbs_last_error()); return 2;
    }
    if (!o.ppd_set) o.nbs.bank.ppd = ds.ppd;
    if (nbs_opts_finalise(&o) != 0) { fprintf(stderr, "error: %s\n", nbs_last_error()); return 2; }

    nw = nbs_w_count_ex(o.nbs.bank.nbands, o.nbs.w_sym);
    c.nw = nw;

    printf("dataset : %s  (%d pairs, %d groups)\n", ds.name, ds.n, ds.ngroups);
    printf("ppd     : %.6g   bands: %d   free exponents: %d (%s)\n",
           o.nbs.bank.ppd, o.nbs.bank.nbands, nw,
           o.nbs.w_sym ? "w_ij == w_ji" : "all free");
    printf("bands   : [");
    for (i = 0; i < o.nbs.bank.nbands; ++i)
        printf("%s%.4g%s", (i ? " " : ""), o.nbs.bank.f[i],
               o.nbs.bank.role[i] == NBS_BAND_LP ? "L" :
               o.nbs.bank.role[i] == NBS_BAND_HP ? "H" : "");
    printf("] cpd\n");
    if (cload != NULL) printf("loading cached log-xi from %s ...\n", cload);
    else               printf("caching log-xi with stride %d ...\n", stride);
    fflush(stdout);

    c.off   = (long *)  nbs_malloc((size_t)ds.n * sizeof(long));
    c.np    = (int *)   nbs_malloc((size_t)ds.n * sizeof(int));
    c.subj  = (double *)nbs_malloc((size_t)ds.n * sizeof(double));
    c.group = (int *)   nbs_malloc((size_t)ds.n * sizeof(int));
    if (!c.off || !c.np || !c.subj || !c.group) goto done;

    /* Pass 1: build the cache (or load one built earlier). */
    c.n = 0;
    if (cload != NULL) {
        free(c.off); free(c.np); free(c.subj); free(c.group);
        c.off = NULL; c.np = NULL; c.subj = NULL; c.group = NULL;
        if (nbs_cache_load(&c, cload) != 0) goto done;
        if (nbs_cache_check(&c, &o.nbs, stride) != 0) goto done;
        printf("loaded cache: %d pairs\n", c.n);
        total = 0;
        for (i = 0; i < c.n; ++i) total += (long)c.np[i] * (long)c.nw;
    } else {
        size_t cap = 0;
        for (i = 0; i < ds.n; ++i) {
            nbs_plane x, y;
            char pr[2048], pd[2048];
            int npix = 0;
            double *tmp;

            if ((i % subset) != 0) continue;
            x.d = NULL; y.d = NULL;
            nbs_path_join(pr, (int)sizeof(pr), ds.root, ds.e[i].ref);
            nbs_path_join(pd, (int)sizeof(pd), ds.root, ds.e[i].dist);
            if (nbs_load_luma(pr, o.luma, &x) != 0) {
                fprintf(stderr, "skip %s: %s\n", pr, nbs_last_error()); continue;
            }
            if (nbs_load_luma(pd, o.luma, &y) != 0) {
                fprintf(stderr, "skip %s: %s\n", pd, nbs_last_error());
                nbs_plane_free(&x); continue;
            }
            if (x.w != y.w || x.h != y.h) {
                fprintf(stderr, "skip %s: size mismatch\n", pd);
                nbs_plane_free(&x); nbs_plane_free(&y); continue;
            }
            if (nbs_nbssim_logxi(&x, &y, &o.nbs, stride, NULL, &npix) != 0 || npix <= 0) {
                fprintf(stderr, "skip %s: %s\n", pd, nbs_last_error());
                nbs_plane_free(&x); nbs_plane_free(&y); continue;
            }
            tmp = (double *)nbs_malloc((size_t)npix * (size_t)nw * sizeof(double));
            if (tmp == NULL) { nbs_plane_free(&x); nbs_plane_free(&y); goto done; }
            if (nbs_nbssim_logxi(&x, &y, &o.nbs, stride, tmp, &npix) != 0) {
                free(tmp); nbs_plane_free(&x); nbs_plane_free(&y); continue;
            }
            nbs_plane_free(&x); nbs_plane_free(&y);

            if ((size_t)total + (size_t)npix * (size_t)nw > cap) {
                size_t nc = cap ? cap * 2 : (size_t)1 << 20;
                float *p;
                while (nc < (size_t)total + (size_t)npix * (size_t)nw) nc *= 2;
                p = (float *)realloc(c.L, nc * sizeof(float));
                if (p == NULL) { free(tmp); nbs_set_error("out of memory"); goto done; }
                c.L = p; cap = nc;
            }
            for (k = 0; k < npix * nw; ++k) c.L[total + k] = (float)tmp[k];
            free(tmp);

            c.off[c.n]   = total;
            c.np[c.n]    = npix;
            c.subj[c.n]  = ds.e[i].score;
            c.group[c.n] = ds.e[i].group;
            total += (long)npix * (long)nw;
            c.n++;

            if (progress && (c.n % 5) == 0) {
                fprintf(stderr, "\rcached %d (%.1f MB)", c.n,
                        (double)total * sizeof(float) / 1048576.0);
                fflush(stderr);
            }
        }
    }
    if (progress && cload == NULL) fprintf(stderr, "\n");
    if (c.n < 10) { fprintf(stderr, "too few usable pairs (%d)\n", c.n); goto done; }
    if (cload == NULL)
        printf("cached %d pairs, %.1f MB\n", c.n, (double)total * sizeof(float) / 1048576.0);
    if (csave != NULL) {
        nbs_cache_stamp(&c, &o.nbs, stride);
        if (nbs_cache_save(&c, csave) != 0)
            fprintf(stderr, "warning: %s\n", nbs_last_error());
        else
            printf("cache written to %s\n", csave);
    }
    printf("\n");
    fflush(stdout);

    scr = (double *)nbs_malloc((size_t)c.n * sizeof(double));
    sub = (double *)nbs_malloc((size_t)c.n * sizeof(double));
    use = (int *)   nbs_malloc((size_t)c.n * sizeof(int));
    if (!scr || !sub || !use) goto done;

    oc.c = &c; oc.objective = objective; oc.allow_neg = allow_neg;
    oc.scr = scr; oc.sub = sub; oc.have_fit = 0;

    /* Baseline: the starting (uniform) exponents. */
    {
        int a, b, t = 0;
        for (a = 0; a < o.nbs.bank.nbands; ++a)
            for (b = (o.nbs.w_sym ? a : 0); b < o.nbs.bank.nbands; ++b)
                w[t++] = o.nbs.w[nbs_w_index_ex(a, b, o.nbs.bank.nbands, o.nbs.w_sym)];
    }
    for (i = 0; i < c.n; ++i) use[i] = 1;
    printf("baseline (starting exponents):\n");
    eval_report(&c, w, NULL, scr, sub, "all data");
    fflush(stdout);

    /* Cross-validation. */
    if (folds > 1) {
        double sum_plcc = 0.0, sum_srcc = 0.0;
        int f, nfold = 0;
        printf("\n%d-fold cross-validation (folds split by reference image):\n", folds);
        for (f = 0; f < folds; ++f) {
            double wf[NBS_MAX_W];
            int *tr = (int *)nbs_malloc((size_t)c.n * sizeof(int));
            int *te = (int *)nbs_malloc((size_t)c.n * sizeof(int));
            nbs_fitreport r;
            int m, rs;
            if (!tr || !te) { free(tr); free(te); goto done; }
            for (i = 0; i < c.n; ++i) {
                te[i] = ((c.group[i] % folds) == f);
                tr[i] = !te[i];
            }
            memcpy(wf, w, (size_t)nw * sizeof(double));
            oc.use = tr;
            oc.have_fit = 0;
            g_rng = 987654321UL + (unsigned long)f * 7919UL;
            optimise(wf, nw, &oc, restarts, maxit, basin, spread);
            (void)rs; (void)step;
            m = score_all(&c, wf, te, scr, sub, c.subj);
            if (m >= 6) {
                nbs_evaluate(scr, sub, m, &r);
                printf("  fold %d: test n=%4d   PLCC %7.4f   SRCC %7.4f\n",
                       f, m, fabs(r.plcc), fabs(r.srcc));
                sum_plcc += fabs(r.plcc); sum_srcc += fabs(r.srcc);
                nfold++;
            }
            free(tr); free(te);
            fflush(stdout);
        }
        if (nfold > 0)
            printf("  mean held-out:   PLCC %7.4f   SRCC %7.4f\n",
                   sum_plcc / nfold, sum_srcc / nfold);
    }

    /* Final fit on everything. */
    printf("\nfinal fit on all data:\n");
    oc.use = NULL;
    oc.have_fit = 0;
    g_rng = 987654321UL;
    optimise(w, nw, &oc, restarts, maxit, basin, spread);
    eval_report(&c, w, NULL, scr, sub, "all data (fitted)");

    printf("\nfitted exponents (%s):\n",
           o.nbs.w_sym ? "upper triangle, w_ij == w_ji"
                       : "row i = horizontal channel, column j = vertical");
    {
        int a, b, t = 0;
        for (a = 0; a < o.nbs.bank.nbands; ++a) {
            printf("  ");
            if (o.nbs.w_sym) for (b = 0; b < a; ++b) printf("%9s", "");
            for (b = (o.nbs.w_sym ? a : 0); b < o.nbs.bank.nbands; ++b)
                printf("%9.5f", w[t++]);
            printf("\n");
        }
    }

    if (outp != NULL) {
        int a, b, t = 0;
        for (a = 0; a < o.nbs.bank.nbands; ++a)
            for (b = (o.nbs.w_sym ? a : 0); b < o.nbs.bank.nbands; ++b)
                o.nbs.w[nbs_w_index_ex(a, b, o.nbs.bank.nbands, o.nbs.w_sym)] = w[t++];
        if (nbs_save_weights(outp, &o.nbs) != 0) {
            fprintf(stderr, "error: %s\n", nbs_last_error());
        } else {
            printf("\nwritten to %s\n", outp);
        }
    }
    rc = 0;

done:
    free(scr); free(sub); free(use);
    nbs_cache_free(&c);
    nbs_dataset_free(&ds);
    if (rc != 0 && nbs_last_error()[0] != '\0')
        fprintf(stderr, "error: %s\n", nbs_last_error());
    return rc;
}
