/* stats.c -- correlations, VQEG logistic fit and Nelder-Mead.  Strict C89. */

#include "nbs.h"
#include "nbs_stats.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ------------------------------------------------------------------ */

double nbs_pearson(const double *a, const double *b, int n)
{
    double ma = 0.0, mb = 0.0, sa = 0.0, sb = 0.0, sab = 0.0;
    int i;
    if (n < 2) return 0.0;
    for (i = 0; i < n; ++i) { ma += a[i]; mb += b[i]; }
    ma /= (double)n; mb /= (double)n;
    for (i = 0; i < n; ++i) {
        double da = a[i] - ma, db = b[i] - mb;
        sa += da * da; sb += db * db; sab += da * db;
    }
    if (sa <= 0.0 || sb <= 0.0) return 0.0;
    return sab / sqrt(sa * sb);
}

/* ---- ranking with average ties ----------------------------------- */

typedef struct { double v; int i; } rkitem;

static int rkcmp(const void *pa, const void *pb)
{
    const rkitem *a = (const rkitem *)pa;
    const rkitem *b = (const rkitem *)pb;
    if (a->v < b->v) return -1;
    if (a->v > b->v) return  1;
    return 0;
}

static int rank_avg(const double *x, int n, double *rank)
{
    rkitem *t = (rkitem *)nbs_malloc((size_t)n * sizeof(rkitem));
    int i, j;
    if (t == NULL) return 1;
    for (i = 0; i < n; ++i) { t[i].v = x[i]; t[i].i = i; }
    qsort(t, (size_t)n, sizeof(rkitem), rkcmp);
    i = 0;
    while (i < n) {
        double r;
        j = i;
        while (j + 1 < n && t[j + 1].v == t[i].v) j++;
        r = 0.5 * ((double)(i + 1) + (double)(j + 1));
        for (; i <= j; ++i) rank[t[i].i] = r;
    }
    free(t);
    return 0;
}

double nbs_spearman(const double *a, const double *b, int n)
{
    double *ra, *rb, r;
    if (n < 2) return 0.0;
    ra = (double *)nbs_malloc((size_t)n * sizeof(double));
    rb = (double *)nbs_malloc((size_t)n * sizeof(double));
    if (ra == NULL || rb == NULL) { free(ra); free(rb); return 0.0; }
    if (rank_avg(a, n, ra) != 0 || rank_avg(b, n, rb) != 0) { free(ra); free(rb); return 0.0; }
    r = nbs_pearson(ra, rb, n);
    free(ra); free(rb);
    return r;
}

double nbs_kendall(const double *a, const double *b, int n)
{
    long conc = 0, disc = 0, ta = 0, tb = 0;
    int i, j;
    double n0, den;
    if (n < 2) return 0.0;
    for (i = 0; i < n; ++i) {
        for (j = i + 1; j < n; ++j) {
            double da = a[i] - a[j];
            double db = b[i] - b[j];
            if (da == 0.0 && db == 0.0) { ta++; tb++; }
            else if (da == 0.0) ta++;
            else if (db == 0.0) tb++;
            else if ((da > 0.0) == (db > 0.0)) conc++;
            else disc++;
        }
    }
    n0 = 0.5 * (double)n * (double)(n - 1);
    den = sqrt((n0 - (double)ta) * (n0 - (double)tb));
    if (den <= 0.0) return 0.0;
    return (double)(conc - disc) / den;
}

/* ------------------------------------------------------------------ */
/* Nelder-Mead                                                         */
/* ------------------------------------------------------------------ */

double nbs_nelder_mead(int nd, double *x, const double *step,
                       double (*fn)(const double *, void *), void *ctx,
                       int maxit, double tol)
{
    double *sx, *fv, *xc, *xr, *xe;
    int i, j, it, lo, hi, nhi;
    double best;

    if (nd < 1) return fn(x, ctx);

    sx = (double *)nbs_malloc((size_t)(nd + 1) * (size_t)nd * sizeof(double));
    fv = (double *)nbs_malloc((size_t)(nd + 1) * sizeof(double));
    xc = (double *)nbs_malloc((size_t)nd * sizeof(double));
    xr = (double *)nbs_malloc((size_t)nd * sizeof(double));
    xe = (double *)nbs_malloc((size_t)nd * sizeof(double));
    if (!sx || !fv || !xc || !xr || !xe) {
        free(sx); free(fv); free(xc); free(xr); free(xe);
        return fn(x, ctx);
    }

    for (i = 0; i <= nd; ++i) {
        for (j = 0; j < nd; ++j) sx[i * nd + j] = x[j];
        if (i > 0) {
            double s = step[i - 1];
            if (s == 0.0) s = (x[i - 1] != 0.0) ? 0.05 * x[i - 1] : 0.05;
            sx[i * nd + (i - 1)] += s;
        }
        fv[i] = fn(&sx[i * nd], ctx);
    }

    for (it = 0; it < maxit; ++it) {
        double frange, fr, fe, fc2;

        lo = 0; hi = 0; nhi = 0;
        for (i = 1; i <= nd; ++i) {
            if (fv[i] < fv[lo]) lo = i;
            if (fv[i] > fv[hi]) hi = i;
        }
        nhi = (hi == 0) ? 1 : 0;
        for (i = 0; i <= nd; ++i)
            if (i != hi && fv[i] > fv[nhi]) nhi = i;

        frange = fabs(fv[hi]) + fabs(fv[lo]);
        if (frange < 1e-300) frange = 1e-300;
        if (2.0 * fabs(fv[hi] - fv[lo]) / frange < tol) break;

        for (j = 0; j < nd; ++j) {
            double s = 0.0;
            for (i = 0; i <= nd; ++i) if (i != hi) s += sx[i * nd + j];
            xc[j] = s / (double)nd;
        }

        for (j = 0; j < nd; ++j) xr[j] = xc[j] + 1.0 * (xc[j] - sx[hi * nd + j]);
        fr = fn(xr, ctx);

        if (fr < fv[lo]) {
            for (j = 0; j < nd; ++j) xe[j] = xc[j] + 2.0 * (xc[j] - sx[hi * nd + j]);
            fe = fn(xe, ctx);
            if (fe < fr) { memcpy(&sx[hi * nd], xe, (size_t)nd * sizeof(double)); fv[hi] = fe; }
            else         { memcpy(&sx[hi * nd], xr, (size_t)nd * sizeof(double)); fv[hi] = fr; }
        } else if (fr < fv[nhi]) {
            memcpy(&sx[hi * nd], xr, (size_t)nd * sizeof(double)); fv[hi] = fr;
        } else {
            for (j = 0; j < nd; ++j) xe[j] = xc[j] + 0.5 * (sx[hi * nd + j] - xc[j]);
            fc2 = fn(xe, ctx);
            if (fc2 < fv[hi]) {
                memcpy(&sx[hi * nd], xe, (size_t)nd * sizeof(double)); fv[hi] = fc2;
            } else {
                for (i = 0; i <= nd; ++i) {
                    if (i == lo) continue;
                    for (j = 0; j < nd; ++j)
                        sx[i * nd + j] = 0.5 * (sx[i * nd + j] + sx[lo * nd + j]);
                    fv[i] = fn(&sx[i * nd], ctx);
                }
            }
        }
    }

    lo = 0;
    for (i = 1; i <= nd; ++i) if (fv[i] < fv[lo]) lo = i;
    memcpy(x, &sx[lo * nd], (size_t)nd * sizeof(double));
    best = fv[lo];

    free(sx); free(fv); free(xc); free(xr); free(xe);
    return best;
}

/* ------------------------------------------------------------------ */
/* VQEG logistic                                                       */
/* ------------------------------------------------------------------ */

double nbs_logistic_eval(const double *p, double x)
{
    double t = p[1] * (x - p[2]);
    double e;
    if (t > 700.0) t = 700.0;
    if (t < -700.0) t = -700.0;
    e = 1.0 / (1.0 + exp(-t));
    return p[0] * (e - 0.5) + p[3] * x + p[4];
}

typedef struct {
    const double *x;
    const double *y;
    int n;
} fitctx;

static double fit_sse(const double *p, void *ctx)
{
    const fitctx *c = (const fitctx *)ctx;
    double s = 0.0;
    int i;
    for (i = 0; i < c->n; ++i) {
        double d = nbs_logistic_eval(p, c->x[i]) - c->y[i];
        s += d * d;
    }
    if (!(s == s)) return 1e300;
    return s;
}

int nbs_logistic_fit(const double *x, const double *y, int n, double *p)
{
    return nbs_logistic_fit_ex(x, y, n, p, 0, 4, 4000);
}

int nbs_logistic_fit_ex(const double *x, const double *y, int n, double *p,
                        int warm, int passes, int maxit)
{
    fitctx c;
    double step[5];
    double mx = 0.0, my = 0.0, sx = 0.0, sy = 0.0, r;
    int i, pass;

    if (n < 6) return 1;
    for (i = 0; i < n; ++i) { mx += x[i]; my += y[i]; }
    mx /= (double)n; my /= (double)n;
    for (i = 0; i < n; ++i) {
        sx += (x[i] - mx) * (x[i] - mx);
        sy += (y[i] - my) * (y[i] - my);
    }
    sx = sqrt(sx / (double)n);
    sy = sqrt(sy / (double)n);
    if (sx <= 0.0) return 1;

    r = nbs_pearson(x, y, n);

    if (!warm) {
        p[0] = (r >= 0.0 ? 1.0 : -1.0) * 4.0 * sy;
        p[1] = 2.0 / sx;
        p[2] = mx;
        p[3] = 0.0;
        p[4] = my;
    }

    c.x = x; c.y = y; c.n = n;
    for (pass = 0; pass < passes; ++pass) {
        step[0] = 0.5 * fabs(p[0]) + 0.1 * sy;
        step[1] = 0.5 * fabs(p[1]) + 0.1 / (sx > 0.0 ? sx : 1.0);
        step[2] = 0.5 * sx;
        step[3] = 0.5 * (sy / (sx > 0.0 ? sx : 1.0)) + 1e-3;
        step[4] = 0.5 * sy + 1e-3;
        for (i = 0; i < 5; ++i) step[i] *= pow(0.35, (double)pass + (warm ? 2.0 : 0.0));
        nbs_nelder_mead(5, p, step, fit_sse, &c, maxit, 1e-12);
    }
    return 0;
}

int nbs_evaluate(const double *obj, const double *subj, int n, nbs_fitreport *rep)
{
    double *pred;
    int i;
    double sse = 0.0, sae = 0.0;

    memset(rep, 0, sizeof(*rep));
    rep->n = n;
    if (n < 2) return 1;

    rep->srcc = nbs_spearman(obj, subj, n);
    rep->krcc = nbs_kendall(obj, subj, n);

    pred = (double *)nbs_malloc((size_t)n * sizeof(double));
    if (pred == NULL) return 1;

    if (nbs_logistic_fit(obj, subj, n, rep->p) == 0) {
        rep->fit_ok = 1;
        for (i = 0; i < n; ++i) pred[i] = nbs_logistic_eval(rep->p, obj[i]);
    } else {
        rep->fit_ok = 0;
        for (i = 0; i < n; ++i) pred[i] = obj[i];
    }

    rep->plcc = nbs_pearson(pred, subj, n);
    for (i = 0; i < n; ++i) {
        double d = pred[i] - subj[i];
        sse += d * d;
        sae += fabs(d);
    }
    rep->rmse = sqrt(sse / (double)n);
    rep->mae  = sae / (double)n;

    free(pred);
    return 0;
}
