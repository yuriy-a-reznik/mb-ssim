/* nbsmetric -- compute SSIM, MS-SSIM and N-band SSIM for one image pair.
 * Strict C89.
 */

#include "nbs.h"
#include "nbs_metrics.h"
#include "nbs_cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(void)
{
    printf("usage: nbsmetric [options] REFERENCE DISTORTED\n\n");
    printf("Computes SSIM, MS-SSIM and the proposed N-band SSIM between two\n");
    printf("8-bit/channel images (PNG, BMP or PNM).  Both are converted to\n");
    printf("luma first.\n\n");
    printf("  --quiet          print only the three numbers\n");
    printf("  --bands-report   list the filterbank cut-offs and peaks\n\n");
    nbs_opts_print_help();
}

int main(int argc, char **argv)
{
    nbs_opts o;
    nbs_plane x, y;
    const char *pref = NULL, *pdis = NULL;
    double s_ssim = 0.0, s_ms = 0.0, s_nb = 0.0;
    int i, quiet = 0, report = 0, rc = 1;

    x.d = NULL; y.d = NULL;
    nbs_opts_default(&o);

    for (i = 1; i < argc; ++i) {
        int r;
        if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) { usage(); return 0; }
        if (strcmp(argv[i], "--quiet") == 0) { quiet = 1; continue; }
        if (strcmp(argv[i], "--bands-report") == 0) { report = 1; continue; }
        r = nbs_opts_parse_one(&o, argv[i]);
        if (r < 0) { fprintf(stderr, "error: %s\n", nbs_last_error()); return 2; }
        if (r == 1) continue;
        if (argv[i][0] == '-') { fprintf(stderr, "unknown option: %s\n", argv[i]); return 2; }
        if (pref == NULL) pref = argv[i];
        else if (pdis == NULL) pdis = argv[i];
        else { fprintf(stderr, "too many positional arguments\n"); return 2; }
    }
    if (pref == NULL || pdis == NULL) { usage(); return 2; }
    if (nbs_opts_finalise(&o) != 0) { fprintf(stderr, "error: %s\n", nbs_last_error()); return 2; }

    if (nbs_load_luma(pref, o.luma, &x) != 0) {
        fprintf(stderr, "error loading '%s': %s\n", pref, nbs_last_error());
        goto done;
    }
    if (nbs_load_luma(pdis, o.luma, &y) != 0) {
        fprintf(stderr, "error loading '%s': %s\n", pdis, nbs_last_error());
        goto done;
    }
    if (x.w != y.w || x.h != y.h) {
        fprintf(stderr, "error: image dimensions differ (%dx%d vs %dx%d)\n",
                x.w, x.h, y.w, y.h);
        goto done;
    }

    if (report) {
        nbs_kernel ks[NBS_MAX_BANDS];
        int b, n = o.nbs.bank.nbands;
        printf("filterbank: %d bands, ppd = %.4g (Nyquist %.4g cpd)\n",
               n, o.nbs.bank.ppd, 0.5 * o.nbs.bank.ppd);
        if (nbs_bank_kernels(&o.nbs.bank, ks) == 0) {
            for (b = 0; b < n; ++b) {
                double best = -1e30, bf = 0.0, nu;
                const char *role = (o.nbs.bank.role[b] == NBS_BAND_LP) ? "LP" :
                                   (o.nbs.bank.role[b] == NBS_BAND_HP) ? "HP" : "BP";
                for (nu = 1e-5; nu < 0.5; nu += 1e-5) {
                    double r2 = nbs_kernel_response(&ks[b], nu);
                    if (r2 > best) { best = r2; bf = nu * o.nbs.bank.ppd; }
                }
                printf("  band %d %s  nominal %8.4f cpd  realised %8.4f cpd"
                       "  gain %.4f  taps %d\n",
                       b, role, o.nbs.bank.f[b], bf, best, 2 * ks[b].r + 1);
                nbs_kernel_free(&ks[b]);
            }
        }
        printf("\n");
    }

    if (nbs_ssim(&x, &y, &o.ssim, &s_ssim) != 0) {
        fprintf(stderr, "SSIM failed: %s\n", nbs_last_error()); goto done;
    }
    if (nbs_msssim(&x, &y, &o.msssim, &s_ms) != 0) {
        fprintf(stderr, "MS-SSIM failed: %s\n", nbs_last_error());
        s_ms = -1.0;
    }
    if (nbs_nbssim(&x, &y, &o.nbs, &s_nb) != 0) {
        fprintf(stderr, "NB-SSIM failed: %s\n", nbs_last_error()); goto done;
    }

    if (quiet) {
        printf("%.12f %.12f %.12f\n", s_ssim, s_ms, s_nb);
    } else {
        printf("image   : %d x %d\n", x.w, x.h);
        printf("SSIM    : %.8f\n", s_ssim);
        if (s_ms >= 0.0) printf("MS-SSIM : %.8f\n", s_ms);
        else             printf("MS-SSIM : n/a (image too small for %d scales)\n", o.msssim.levels);
        printf("NB-SSIM : %.8f\n", s_nb);
    }
    rc = 0;

done:
    nbs_plane_free(&x);
    nbs_plane_free(&y);
    return rc;
}
