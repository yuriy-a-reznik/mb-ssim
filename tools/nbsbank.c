/* nbsbank -- print the realised filterbank frequency responses as CSV.
 *
 * Columns: f_cpd, nu_cyc_per_px, H_0 .. H_{N-1}, sum
 * The responses are the exact DTFTs of the FIR kernels the metric actually
 * uses, not idealised shapes, so this doubles as a sanity check on the
 * kernel construction.
 *
 * Strict C89.
 */

#include "nbs.h"
#include "nbs_metrics.h"
#include "nbs_cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(int argc, char **argv)
{
    nbs_opts o;
    nbs_kernel ks[NBS_MAX_BANDS];
    int i, b, n, npts = 400;
    double fmin = 0.05, fmax;

    nbs_opts_default(&o);
    for (i = 1; i < argc; ++i) {
        int r;
        if (strcmp(argv[i], "--help") == 0) {
            printf("usage: nbsbank [options] > bank.csv\n\n");
            printf("  --points=N     number of log-spaced frequency samples\n\n");
            nbs_opts_print_help();
            return 0;
        }
        if (strncmp(argv[i], "--points=", 9) == 0) { npts = atoi(argv[i] + 9); continue; }
        r = nbs_opts_parse_one(&o, argv[i]);
        if (r < 0) { fprintf(stderr, "error: %s\n", nbs_last_error()); return 2; }
        if (r == 0) { fprintf(stderr, "unknown option: %s\n", argv[i]); return 2; }
    }
    if (nbs_opts_finalise(&o) != 0) { fprintf(stderr, "error: %s\n", nbs_last_error()); return 2; }
    if (npts < 2) npts = 2;

    n = o.nbs.bank.nbands;
    fmax = 0.5 * o.nbs.bank.ppd;

    if (nbs_bank_kernels(&o.nbs.bank, ks) != 0) {
        fprintf(stderr, "error: %s\n", nbs_last_error());
        return 2;
    }

    fprintf(stderr, "# ppd = %.6g, display Nyquist = %.6g cpd\n",
            o.nbs.bank.ppd, fmax);
    for (b = 0; b < n; ++b) {
        double best = -1e30, bf = 0.0, nu;
        const char *role = (o.nbs.bank.role[b] == NBS_BAND_LP) ? "LP" :
                           (o.nbs.bank.role[b] == NBS_BAND_HP) ? "HP" : "BP";
        for (nu = 1e-5; nu < 0.5; nu += 1e-5) {
            double r = nbs_kernel_response(&ks[b], nu);
            if (r > best) { best = r; bf = nu * o.nbs.bank.ppd; }
        }
        fprintf(stderr, "# band %d  %s  nominal %8.4f cpd  realised peak %8.4f cpd"
                        "  gain %.4f  taps %d\n",
                b, role, o.nbs.bank.f[b], bf, best, 2 * ks[b].r + 1);
    }

    printf("f_cpd,nu");
    for (b = 0; b < n; ++b) printf(",H%d", b);
    printf(",sum\n");

    for (i = 0; i < npts; ++i) {
        double t = (double)i / (double)(npts - 1);
        double f = fmin * pow(fmax / fmin, t);
        double nu = f / o.nbs.bank.ppd;
        double s = 0.0;
        printf("%.6g,%.6g", f, nu);
        for (b = 0; b < n; ++b) {
            double h = nbs_kernel_response(&ks[b], nu);
            s += h;
            printf(",%.6g", h);
        }
        printf(",%.6g\n", s);
    }

    for (b = 0; b < n; ++b) nbs_kernel_free(&ks[b]);
    return 0;
}
