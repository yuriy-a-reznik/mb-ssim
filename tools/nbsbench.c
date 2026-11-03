/* nbsbench -- run SSIM, MS-SSIM and NB-SSIM over a dataset index and report
 * correlation with the subjective scores, using the VQEG methodology
 * (five-parameter logistic fit, then PLCC / RMSE / MAE, plus SRCC and KRCC).
 *
 * Strict C89.
 */

#include "nbs.h"
#include "nbs_metrics.h"
#include "nbs_data.h"
#include "nbs_stats.h"
#include "nbs_cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NMET 3
static const char *met_name[NMET] = { "SSIM", "MS-SSIM", "NB-SSIM" };

static int g_sign = 1;   /* +1 if a larger subjective score is better */
static int g_want[NMET] = { 1, 1, 1 };

/* --resume support.  Long runs get interrupted; rows are flushed per pair, so
 * an interrupted output file is a valid partial result and there is no reason
 * to recompute it. */
typedef struct { char **key; int n, cap; } donelist;

static void done_add(donelist *d, const char *ref, const char *dist)
{
    size_t n = strlen(ref) + strlen(dist) + 2;
    char *k;
    if (d->n == d->cap) {
        int nc = d->cap ? d->cap * 2 : 256;
        char **p = (char **)realloc(d->key, (size_t)nc * sizeof(char *));
        if (p == NULL) return;
        d->key = p; d->cap = nc;
    }
    k = (char *)nbs_malloc(n);
    if (k == NULL) return;
    strcpy(k, ref); strcat(k, "\001"); strcat(k, dist);
    d->key[d->n++] = k;
}

static int done_has(const donelist *d, const char *ref, const char *dist)
{
    int i;
    size_t lr = strlen(ref);
    for (i = 0; i < d->n; ++i) {
        const char *k = d->key[i];
        if (strncmp(k, ref, lr) == 0 && k[lr] == '\001' &&
            strcmp(k + lr + 1, dist) == 0) return 1;
    }
    return 0;
}

static void done_load(donelist *d, const char *path)
{
    FILE *f = fopen(path, "r");
    char line[4096];
    if (f == NULL) return;
    while (fgets(line, (int)sizeof(line), f) != NULL) {
        char *c1, *c2;
        if (line[0] == '#' || strncmp(line, "ref,", 4) == 0) continue;
        c1 = strchr(line, ',');
        if (c1 == NULL) continue;
        *c1 = '\0';
        c2 = strchr(c1 + 1, ',');
        if (c2 == NULL) continue;
        *c2 = '\0';
        done_add(d, line, c1 + 1);
    }
    fclose(f);
}

static void usage(void)
{
    printf("usage: nbsbench [options] INDEX.csv\n\n");
    printf("  --out=PATH        write per-image scores to PATH (CSV)\n");
    printf("  --shard=I/N       process only every N-th entry starting at I\n");
    printf("  --metrics=LIST    comma list of ssim,msssim,nbssim (default all).\n");
    printf("                    Skipped metrics are written as -1.  SSIM and\n");
    printf("                    MS-SSIM do not depend on ppd, so a ppd sweep\n");
    printf("                    only needs --metrics=nbssim.\n");
    printf("  --report=PATH     skip computation; read scores from PATH and\n");
    printf("                    report statistics only (use to merge shards)\n");
    printf("  --resume          append to --out, skipping pairs already in it\n");
    printf("  --by-distortion   also report per-reference-group breakdown\n");
    printf("  --progress        print a progress line to stderr\n\n");
    nbs_opts_print_help();
}

static int read_scores(const char *path, double **obj, double **subj,
                       int *n, int nmet)
{
    FILE *f = fopen(path, "r");
    char line[4096];
    int cap = 1024, k = 0, m;
    double *o[NMET], *s;

    if (f == NULL) { nbs_set_error("cannot open score file"); return 1; }
    for (m = 0; m < nmet; ++m) o[m] = (double *)nbs_malloc((size_t)cap * sizeof(double));
    s = (double *)nbs_malloc((size_t)cap * sizeof(double));

    while (fgets(line, (int)sizeof(line), f) != NULL) {
        char *p = line, *c;
        int col = 0;
        double v[NMET + 1];
        if (line[0] == '#') {
            if (strncmp(line, "# lower_better:", 15) == 0)
                g_sign = atoi(line + 15) ? -1 : 1;
            continue;
        }
        if (strncmp(line, "ref,", 4) == 0) continue;
        /* fields: ref,dist,refid,subj,ssim,msssim,nbssim */
        for (col = 0; col < 3; ++col) {
            c = strchr(p, ','); if (c == NULL) break; p = c + 1;
        }
        if (col < 3) continue;
        for (col = 0; col <= nmet; ++col) {
            v[col] = atof(p);
            c = strchr(p, ',');
            if (c == NULL) { if (col < nmet) col = -1; break; }
            p = c + 1;
        }
        if (col < 0) continue;
        if (k == cap) {
            cap *= 2;
            for (m = 0; m < nmet; ++m) o[m] = (double *)realloc(o[m], (size_t)cap * sizeof(double));
            s = (double *)realloc(s, (size_t)cap * sizeof(double));
        }
        s[k] = v[0];
        for (m = 0; m < nmet; ++m) o[m][k] = v[m + 1];
        k++;
    }
    fclose(f);
    for (m = 0; m < nmet; ++m) obj[m] = o[m];
    *subj = s;
    *n = k;
    return 0;
}

static void report(const char *label, const double *obj, const double *subj, int n)
{
    nbs_fitreport r;
    if (nbs_evaluate(obj, subj, n, &r) != 0) {
        printf("  %-10s  (evaluation failed)\n", label);
        return;
    }
    printf("  %-10s  PLCC %7.4f   SRCC %7.4f   KRCC %7.4f   RMSE %8.4f   MAE %8.4f\n",
           label, fabs(r.plcc), g_sign * r.srcc, g_sign * r.krcc, r.rmse, r.mae);
}

int main(int argc, char **argv)
{
    nbs_opts o;
    nbs_dataset ds;
    const char *index = NULL, *outp = NULL, *reportp = NULL;
    int i, m, shard_i = 0, shard_n = 1, progress = 0, rc = 1, resume = 0;
    donelist done;
    double *obj[NMET], *subj = NULL;
    int nvalid = 0;
    FILE *fout = NULL;

    for (m = 0; m < NMET; ++m) obj[m] = NULL;
    memset(&ds, 0, sizeof(ds));
    memset(&done, 0, sizeof(done));
    nbs_opts_default(&o);

    for (i = 1; i < argc; ++i) {
        int r;
        if (strcmp(argv[i], "--help") == 0) { usage(); return 0; }
        if (strncmp(argv[i], "--out=", 6) == 0) { outp = argv[i] + 6; continue; }
        if (strncmp(argv[i], "--report=", 9) == 0) { reportp = argv[i] + 9; continue; }
        if (strncmp(argv[i], "--metrics=", 10) == 0) {
            const char *v = argv[i] + 10;
            for (m = 0; m < NMET; ++m) g_want[m] = 0;
            if (strstr(v, "ssim")   != NULL) g_want[0] = 1;
            if (strstr(v, "msssim") != NULL) g_want[1] = 1;
            if (strstr(v, "nbssim") != NULL) g_want[2] = 1;
            /* "ssim" is a substring of the other two names; require an exact
             * token match for the plain SSIM slot. */
            {
                const char *p2 = v;
                int found = 0;
                while (*p2) {
                    if (strncmp(p2, "ssim", 4) == 0 &&
                        (p2 == v || p2[-1] == ',') &&
                        (p2[4] == '\0' || p2[4] == ',')) { found = 1; break; }
                    p2++;
                }
                g_want[0] = found;
            }
            if (!g_want[0] && !g_want[1] && !g_want[2]) {
                fprintf(stderr, "--metrics: nothing selected\n"); return 2;
            }
            continue;
        }
        if (strcmp(argv[i], "--progress") == 0) { progress = 1; continue; }
        if (strcmp(argv[i], "--resume") == 0) { resume = 1; continue; }
        if (strncmp(argv[i], "--shard=", 8) == 0) {
            const char *s = argv[i] + 8;
            const char *sl = strchr(s, '/');
            if (sl == NULL) { fprintf(stderr, "--shard needs I/N\n"); return 2; }
            shard_i = atoi(s);
            shard_n = atoi(sl + 1);
            if (shard_n < 1 || shard_i < 0 || shard_i >= shard_n) {
                fprintf(stderr, "--shard: bad range\n"); return 2;
            }
            continue;
        }
        r = nbs_opts_parse_one(&o, argv[i]);
        if (r < 0) { fprintf(stderr, "error: %s\n", nbs_last_error()); return 2; }
        if (r == 1) continue;
        if (argv[i][0] == '-') { fprintf(stderr, "unknown option: %s\n", argv[i]); return 2; }
        index = argv[i];
    }

    if (reportp != NULL) {
        int n = 0;
        if (read_scores(reportp, obj, &subj, &n, NMET) != 0) {
            fprintf(stderr, "error: %s\n", nbs_last_error()); return 2;
        }
        printf("scores file: %s   (%d entries)\n", reportp, n);
        for (m = 0; m < NMET; ++m) report(met_name[m], obj[m], subj, n);
        for (m = 0; m < NMET; ++m) free(obj[m]);
        free(subj);
        return 0;
    }

    if (index == NULL) { usage(); return 2; }
    if (nbs_dataset_load(index, &ds) != 0) {
        fprintf(stderr, "error: %s\n", nbs_last_error()); return 2;
    }
    if (!o.ppd_set) o.nbs.bank.ppd = ds.ppd;
    if (nbs_opts_finalise(&o) != 0) { fprintf(stderr, "error: %s\n", nbs_last_error()); return 2; }

    for (m = 0; m < NMET; ++m) {
        obj[m] = (double *)nbs_malloc((size_t)ds.n * sizeof(double));
        if (obj[m] == NULL) goto done;
    }
    subj = (double *)nbs_malloc((size_t)ds.n * sizeof(double));
    if (subj == NULL) goto done;

    if (outp != NULL) {
        if (resume) done_load(&done, outp);
        fout = fopen(outp, resume && done.n > 0 ? "a" : "w");
        if (fout == NULL) { fprintf(stderr, "cannot write %s\n", outp); goto done; }
        if (!(resume && done.n > 0)) {
            fprintf(fout, "# dataset: %s\n# ppd: %.6g\n# lower_better: %d\n",
                    ds.name, o.nbs.bank.ppd, ds.lower_better);
            fprintf(fout, "ref,dist,refid,subj,ssim,msssim,nbssim\n");
        } else {
            printf("resuming: %d pairs already present in %s\n", done.n, outp);
        }
    }

    printf("dataset : %s  (%d pairs, %d reference groups)\n", ds.name, ds.n, ds.ngroups);
    g_sign = ds.lower_better ? -1 : 1;
    printf("score   : %s (%s is better)\n", ds.score_name, ds.lower_better ? "lower" : "higher");
    printf("note    : SRCC/KRCC signs are flipped so that positive means agreement\n");
    printf("ppd     : %.6g\n", o.nbs.bank.ppd);
    printf("bands   : %d  [", o.nbs.bank.nbands);
    for (i = 0; i < o.nbs.bank.nbands; ++i)
        printf("%s%.4g%s", (i ? " " : ""), o.nbs.bank.f[i],
               o.nbs.bank.role[i] == NBS_BAND_LP ? "L" :
               o.nbs.bank.role[i] == NBS_BAND_HP ? "H" : "");
    printf("] cpd\n\n");
    fflush(stdout);

    for (i = 0; i < ds.n; ++i) {
        nbs_plane x, y;
        char pr[2048], pd[2048];
        double v[NMET];
        int bad = 0;

        if (shard_n > 1 && (i % shard_n) != shard_i) continue;
        if (resume && done.n > 0 && done_has(&done, ds.e[i].ref, ds.e[i].dist))
            continue;

        x.d = NULL; y.d = NULL;
        nbs_path_join(pr, (int)sizeof(pr), ds.root, ds.e[i].ref);
        nbs_path_join(pd, (int)sizeof(pd), ds.root, ds.e[i].dist);

        if (nbs_load_luma(pr, o.luma, &x) != 0) {
            fprintf(stderr, "skip %s: %s\n", pr, nbs_last_error()); bad = 1;
        } else if (nbs_load_luma(pd, o.luma, &y) != 0) {
            fprintf(stderr, "skip %s: %s\n", pd, nbs_last_error()); bad = 1;
        } else if (x.w != y.w || x.h != y.h) {
            fprintf(stderr, "skip %s: dimension mismatch\n", pd); bad = 1;
        }
        v[0] = v[1] = v[2] = -1.0;
        if (!bad) {
            if (g_want[0] && nbs_ssim(&x, &y, &o.ssim, &v[0]) != 0) bad = 1;
            if (!bad && g_want[1] && nbs_msssim(&x, &y, &o.msssim, &v[1]) != 0) v[1] = 0.0;
            if (!bad && g_want[2] && nbs_nbssim(&x, &y, &o.nbs, &v[2]) != 0) bad = 1;
        }
        nbs_plane_free(&x); nbs_plane_free(&y);
        if (bad) continue;

        for (m = 0; m < NMET; ++m) obj[m][nvalid] = v[m];
        subj[nvalid] = ds.e[i].score;
        nvalid++;

        if (fout != NULL) {
            fprintf(fout, "%s,%s,%s,%.8g,%.8f,%.8f,%.8f\n",
                    ds.e[i].ref, ds.e[i].dist, ds.e[i].refid, ds.e[i].score,
                    v[0], v[1], v[2]);
            fflush(fout);
        }
        if (progress && (nvalid % 10) == 0) {
            fprintf(stderr, "\r%d/%d", i + 1, ds.n);
            fflush(stderr);
        }
    }
    if (progress) fprintf(stderr, "\n");

    printf("evaluated %d pairs\n\n", nvalid);
    if (resume && done.n > 0 && outp != NULL) {
        /* The arrays hold only the pairs computed in THIS run, so summarising
         * them would silently describe a subset.  Re-read the complete file
         * instead. */
        double *fobj[NMET], *fsubj = NULL;
        int fn = 0;
        if (fout != NULL) { fclose(fout); fout = NULL; }
        if (read_scores(outp, fobj, &fsubj, &fn, NMET) == 0) {
            printf("summary over the complete file (%d pairs, %d new):\n", fn, nvalid);
            for (m = 0; m < NMET; ++m) if (g_want[m]) report(met_name[m], fobj[m], fsubj, fn);
            for (m = 0; m < NMET; ++m) free(fobj[m]);
            free(fsubj);
            rc = 0;
            goto done;
        }
    }
    if (shard_n > 1) {
        printf("(shard %d/%d -- run with --report to combine)\n", shard_i, shard_n);
    }
    for (m = 0; m < NMET; ++m) if (g_want[m]) report(met_name[m], obj[m], subj, nvalid);
    rc = 0;

done:
    if (fout != NULL) fclose(fout);
    for (m = 0; m < NMET; ++m) free(obj[m]);
    free(subj);
    nbs_dataset_free(&ds);
    return rc;
}
