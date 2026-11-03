/* nbstest -- unit tests.  Strict C89. */

#include "nbs.h"
#include "nbs_metrics.h"
#include "nbs_stats.h"
#include "nbs_cli.h"
#include "nbs_fft.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int failures = 0;
static int checks   = 0;

static void ck(int cond, const char *what)
{
    checks++;
    if (!cond) { failures++; printf("  FAIL: %s\n", what); }
}

static void ck_near(double a, double b, double tol, const char *what)
{
    checks++;
    if (!(fabs(a - b) <= tol)) {
        failures++;
        printf("  FAIL: %s (%.10g vs %.10g, tol %g)\n", what, a, b, tol);
    }
}

/* Deterministic pseudo-random generator (no rand() dependence). */
static unsigned long rng_state = 12345UL;
static double rnd(void)
{
    rng_state = rng_state * 1103515245UL + 12345UL;
    return (double)((rng_state >> 16) & 0x7fffUL) / 32767.0;
}

static void fill_texture(nbs_plane *p, double seed)
{
    int x, y;
    for (y = 0; y < p->h; ++y)
        for (x = 0; x < p->w; ++x)
            NBS_AT(p, x, y) = 128.0
                + 60.0 * sin((x + seed) * 0.11) * cos((y - seed) * 0.07)
                + 25.0 * sin((x + 2 * y) * 0.31)
                + 10.0 * sin(x * 0.9 + y * 1.3);
}

/* ---------------------------------------------------------------- */

static void test_kernels(void)
{
    nbs_kernel k;
    int i;
    double s = 0.0;

    printf("kernels\n");
    ck(nbs_kernel_gaussian(&k, 1.5, 4.0) == 0, "gaussian alloc");
    for (i = 0; i <= 2 * k.r; ++i) s += k.tap[i];
    ck_near(s, 1.0, 1e-12, "gaussian unit sum");
    ck(k.tap[k.r] == k.tap[k.r], "gaussian centre finite");
    for (i = 0; i < k.r; ++i)
        ck(fabs(k.tap[i] - k.tap[2 * k.r - i]) < 1e-15, "gaussian symmetry");
    nbs_kernel_free(&k);
}

static void test_conv_dc(void)
{
    nbs_plane a, b;
    int i, n;

    printf("convolution\n");
    nbs_plane_alloc(&a, 40, 30);
    b.d = NULL;
    n = a.w * a.h;
    for (i = 0; i < n; ++i) a.d[i] = 7.25;
    ck(nbs_gauss_blur(&a, &b, 3.0, NBS_BORDER_SYMMETRIC) == 0, "blur ok");
    for (i = 0; i < n; ++i)
        if (fabs(b.d[i] - 7.25) > 1e-9) { ck(0, "blur preserves DC"); break; }
    checks++;
    nbs_plane_free(&a); nbs_plane_free(&b);
}

static void test_bank_dog(void)
{
    nbs_bank_cfg bank;
    nbs_kernel ks[NBS_MAX_BANDS];
    static const double freq[5] = { 0.8, 1.7, 2.8, 4.0, 8.0 };
    int b;

    printf("filterbank (DoG)\n");
    memset(&bank, 0, sizeof(bank));
    bank.kind = NBS_BANK_DOG;
    bank.nbands = 6;
    bank.ppd = 30.0;
    bank.radius_mult = 4.0;
    bank.border = NBS_BORDER_SYMMETRIC;
    for (b = 0; b < NBS_MAX_BANDS; ++b) { bank.rho[b] = 3.0; bank.bw_oct[b] = 0.0; }
    ck(nbs_bank_set(&bank, freq, 5, 0, 0) == 0, "dog bank accepts the requested peaks");
    ck(bank.role[0] == NBS_BAND_LP, "band 0 is low-pass");
    ck(bank.role[5] == NBS_BAND_HP, "band 5 is high-pass");
    ck_near(bank.f[5], 8.0 * sqrt(2.0), 1e-12, "derived high-pass corner");

    ck(nbs_bank_kernels(&bank, ks) == 0, "dog kernels build");

    /* Band-pass peaks must land where requested. */
    for (b = 1; b <= 4; ++b) {
        double best = -1e30, bf = 0.0, nu;
        for (nu = 1e-5; nu < 0.5; nu += 1e-5) {
            double r = nbs_kernel_response(&ks[b], nu);
            if (r > best) { best = r; bf = nu * bank.ppd; }
        }
        ck_near(bf, freq[b], 0.03 * freq[b], "band-pass peak at the requested cpd");
        ck(best > 0.2, "band-pass gain is not degenerate");
    }
    /* Low-pass has unit DC gain, high-pass has zero DC gain. */
    ck_near(nbs_kernel_response(&ks[0], 0.0), 1.0, 1e-9, "low-pass DC gain 1");
    ck_near(nbs_kernel_response(&ks[5], 0.0), 0.0, 1e-9, "high-pass DC gain 0");
    for (b = 1; b <= 4; ++b)
        ck_near(nbs_kernel_response(&ks[b], 0.0), 0.0, 1e-9, "band-pass DC gain 0");

    for (b = 0; b < 6; ++b) nbs_kernel_free(&ks[b]);
}

static void test_bank_wilson(void)
{
    nbs_bank_cfg bank;
    nbs_kernel ks[NBS_MAX_BANDS];
    static const double wf[6] = { 0.8, 1.7, 2.8, 4.0, 8.0, 16.0 };
    static const double wb[6] = { 2.5, 2.2, 1.9, 1.7, 1.5, 1.25 };
    int b, i;

    printf("Wilson band-pass bank\n");

    /* At CID:IQ's 100 cm distance every channel is representable. */
    memset(&bank, 0, sizeof(bank));
    bank.kind = NBS_BANK_LOGGABOR; bank.nbands = 6; bank.ppd = 63.223;
    bank.radius_mult = 4.0; bank.border = NBS_BORDER_SYMMETRIC;
    bank.norm_peak = 1;
    for (i = 0; i < 6; ++i) bank.bw_oct[i] = wb[i];
    ck(nbs_bank_set(&bank, wf, 6, 1, 1) == 0, "wilson bank builds at ppd 63.2");
    ck(bank.nbands == 6, "all six channels survive at ppd 63.2");
    ck(bank.nyquist_dropped == 0, "nothing dropped at ppd 63.2");

    ck(nbs_bank_kernels(&bank, ks) == 0, "wilson kernels build");
    for (b = 0; b < bank.nbands; ++b) {
        double best = -1e30, bf = 0.0, nu, half, lo = 0.0, hi = 0.0, prev = 0.0, pn = 0.0;
        for (nu = 1e-5; nu < 0.5; nu += 1e-5) {
            double r = nbs_kernel_response(&ks[b], nu);
            if (r > best) { best = r; bf = nu * bank.ppd; }
        }
        ck_near(bf, wf[b], 0.06 * wf[b], "wilson peak at the requested cpd");
        ck_near(best, 1.0, 0.02, "peak-normalised band has unit gain");
        /* measured half-height bandwidth */
        half = 0.5 * best;
        for (nu = 1e-5; nu < 0.5; nu += 1e-5) {
            double r = nbs_kernel_response(&ks[b], nu);
            if (lo == 0.0 && r >= half) lo = nu;
            if (lo != 0.0 && hi == 0.0 && nu > bf / bank.ppd && r <= half) hi = nu;
            prev = r; pn = nu;
        }
        (void)prev; (void)pn;
        if (lo > 0.0 && hi > 0.0)
            ck_near(log(hi / lo) / log(2.0), wb[b], 0.25, "bandwidth within 0.25 octave of target");
    }
    for (b = 0; b < bank.nbands; ++b) nbs_kernel_free(&ks[b]);

    /* At 50 cm the 16 cpd channel is above the display Nyquist and must be
     * dropped, not silently mis-realised. */
    memset(&bank, 0, sizeof(bank));
    bank.kind = NBS_BANK_LOGGABOR; bank.nbands = 6; bank.ppd = 31.612;
    bank.radius_mult = 4.0; bank.border = NBS_BORDER_SYMMETRIC;
    bank.norm_peak = 1;
    for (i = 0; i < 6; ++i) bank.bw_oct[i] = wb[i];
    ck(nbs_bank_set(&bank, wf, 6, 1, 1) == 0, "wilson bank builds at ppd 31.6");
    ck(bank.nbands == 5, "the 16 cpd channel is dropped at ppd 31.6");
    ck(bank.nyquist_dropped == 1, "the drop is recorded");
    ck_near(bank.dropped_f[0], 16.0, 1e-9, "the dropped channel is reported as 16 cpd");
}

static void test_bank_cumgauss_partition(void)
{
    nbs_bank_cfg bank;
    nbs_kernel ks[NBS_MAX_BANDS];
    /* Peaks must be spaced by more than sqrt(2) x the previous cut-off, so
     * use a wider set than the CSF-motivated default. */
    static const double peaks[4] = { 0.8, 2.0, 5.0, 12.0 };
    int b, i;
    double worst = 0.0;

    printf("filterbank (cumulative Gaussian)\n");
    memset(&bank, 0, sizeof(bank));
    bank.kind = NBS_BANK_CUMGAUSS;
    bank.nbands = 5;
    bank.ppd = 40.0;
    bank.radius_mult = 5.0;
    bank.border = NBS_BORDER_SYMMETRIC;
    ck(nbs_bank_set_from_peaks(&bank, peaks) == 0, "cumgauss peak placement solves");
    for (b = 0; b + 1 < bank.nbands - 1; ++b)
        ck(bank.fc[b] < bank.fc[b + 1], "cut-offs ascending");

    ck(nbs_bank_kernels(&bank, ks) == 0, "cumgauss kernels build");
    for (i = 0; i <= 500; ++i) {
        double nu = 0.5 * (double)i / 500.0;
        double s = 0.0;
        for (b = 0; b < bank.nbands; ++b) s += nbs_kernel_response(&ks[b], nu);
        if (fabs(s - 1.0) > worst) worst = fabs(s - 1.0);
    }
    ck_near(worst, 0.0, 1e-12, "cumgauss bank is a partition of unity");

    /* The convention H(f_corner) = exp(-1/2) for the low-pass stage. */
    ck_near(nbs_kernel_response(&ks[0], bank.fc[0] / bank.ppd), exp(-0.5), 2e-3,
            "low-pass corner convention");

    for (b = 0; b < bank.nbands; ++b) nbs_kernel_free(&ks[b]);
}

static void test_bank_ppd_invariance(void)
{
    nbs_bank_cfg a, b;
    static const double peaks[4] = { 0.8, 2.0, 5.0, 12.0 };
    int i;

    memset(&a, 0, sizeof(a)); memset(&b, 0, sizeof(b));
    a.kind = b.kind = NBS_BANK_CUMGAUSS;
    a.nbands = b.nbands = 5;
    a.ppd = 30.0; b.ppd = 64.0;
    nbs_bank_set_from_peaks(&a, peaks);
    nbs_bank_set_from_peaks(&b, peaks);
    for (i = 0; i < 4; ++i)
        ck_near(a.fc[i], b.fc[i], 1e-9, "cut-offs in cpd independent of ppd");
}

static void test_ssim_identity(void)
{
    nbs_plane x, y;
    nbs_ssim_cfg cfg;
    nbs_msssim_cfg mcfg;
    nbs_nbssim_cfg ncfg;
    double v;

    printf("metric identities\n");
    nbs_plane_alloc(&x, 256, 256);
    fill_texture(&x, 0.0);
    nbs_plane_copy(&y, &x);

    nbs_ssim_cfg_default(&cfg);
    cfg.downsample = 0;
    ck(nbs_ssim(&x, &y, &cfg, &v) == 0, "ssim runs");
    ck_near(v, 1.0, 1e-12, "SSIM(x,x) == 1");

    nbs_msssim_cfg_default(&mcfg);
    ck(nbs_msssim(&x, &y, &mcfg, &v) == 0, "ms-ssim runs");
    ck_near(v, 1.0, 1e-12, "MS-SSIM(x,x) == 1");

    nbs_nbssim_cfg_default(&ncfg);
    ck(nbs_nbssim(&x, &y, &ncfg, &v) == 0, "nb-ssim runs");
    ck_near(v, 1.0, 1e-9, "NB-SSIM(x,x) == 1");

    nbs_plane_free(&x); nbs_plane_free(&y);
}

static void test_symmetry_and_monotonicity(void)
{
    nbs_plane x, y, z;
    nbs_nbssim_cfg ncfg;
    nbs_ssim_cfg scfg;
    double a, b, s1, s2;
    int i, n;

    nbs_plane_alloc(&x, 256, 256);
    fill_texture(&x, 0.0);
    nbs_plane_copy(&y, &x);
    nbs_plane_copy(&z, &x);
    n = x.w * x.h;
    for (i = 0; i < n; ++i) y.d[i] += 6.0 * (rnd() - 0.5);
    for (i = 0; i < n; ++i) z.d[i] += 30.0 * (rnd() - 0.5);

    nbs_nbssim_cfg_default(&ncfg);
    nbs_ssim_cfg_default(&scfg);
    scfg.downsample = 0;

    ck(nbs_nbssim(&x, &y, &ncfg, &a) == 0, "nbssim x,y");
    ck(nbs_nbssim(&y, &x, &ncfg, &b) == 0, "nbssim y,x");
    ck_near(a, b, 1e-12, "NB-SSIM is symmetric");

    ck(nbs_nbssim(&x, &z, &ncfg, &b) == 0, "nbssim x,z");
    ck(b < a, "NB-SSIM decreases with more noise");

    nbs_ssim(&x, &y, &scfg, &s1);
    nbs_ssim(&x, &z, &scfg, &s2);
    ck(s2 < s1, "SSIM decreases with more noise");
    ck(a < 1.0 && a > 0.0, "NB-SSIM in (0,1) for noisy input");

    nbs_plane_free(&x); nbs_plane_free(&y); nbs_plane_free(&z);
}

/* With two bands and unit exponents on all four 2-D sub-bands the metric
 * should track SSIM closely: this is the N = 2 case of the 2023 paper,
 * separably extended to two dimensions. */
static void test_two_band_tracks_ssim(void)
{
    nbs_plane x, y;
    nbs_nbssim_cfg ncfg;
    nbs_ssim_cfg scfg;
    double sa[6], sb[6];
    int t, i, n, ok = 1;

    printf("two-band model vs SSIM\n");
    nbs_plane_alloc(&x, 256, 256);
    fill_texture(&x, 0.0);
    nbs_plane_alloc(&y, 256, 256);
    n = x.w * x.h;

    nbs_nbssim_cfg_default(&ncfg);
    ncfg.bank.nbands = 2;
    {
        double f2[2];
        f2[0] = ncfg.bank.ppd / (2.0 * 3.14159265358979323846 * 3.0);
        f2[1] = f2[0];
        ncfg.bank.role[0] = NBS_BAND_LP;
        ncfg.bank.role[1] = NBS_BAND_HP;
        ncfg.bank.f[0] = f2[0];
        ncfg.bank.f[1] = f2[0];
    }
    for (i = 0; i < NBS_MAX_W; ++i) ncfg.w[i] = 0.0;
    ncfg.w[nbs_w_index(0, 0, 2)] = 1.0;
    ncfg.w[nbs_w_index(0, 1, 2)] = 0.5;
    ncfg.w[nbs_w_index(1, 1, 2)] = 1.0;

    nbs_ssim_cfg_default(&scfg);
    scfg.downsample = 0;

    for (t = 0; t < 6; ++t) {
        double amp = 2.0 * (double)(t + 1) * (double)(t + 1);
        for (i = 0; i < n; ++i) y.d[i] = x.d[i] + amp * (rnd() - 0.5);
        nbs_ssim(&x, &y, &scfg, &sa[t]);
        nbs_nbssim(&x, &y, &ncfg, &sb[t]);
    }
    for (t = 1; t < 6; ++t) {
        if (!(sa[t] < sa[t - 1]) || !(sb[t] < sb[t - 1])) ok = 0;
    }
    ck(ok, "both metrics decrease monotonically with distortion");
    ck(nbs_spearman(sa, sb, 6) > 0.99, "two-band NB-SSIM rank-tracks SSIM");

    nbs_plane_free(&x); nbs_plane_free(&y);
}

/* Reproduce the central claim of Reznik, EI 2023: patch-level SSIM is well
 * approximated by xi(x_L,y_L) . xi(x_H,y_H), with x_L a Gaussian low-pass of
 * sigma = 3 px and x_H = x - x_L.  This is the isotropic two-band model, so
 * it is computed here directly rather than through the separable bank. */
static double two_band_model(const nbs_plane *x, const nbs_plane *y)
{
    nbs_plane xl, yl, xh, yh, t, exx, eyy, exy;
    nbs_plane fxx, fyy, fxy, gxx, gyy, gxy;
    nbs_kernel win;
    double C1, C2, s = 0.0;
    int i, n, vn;

    xl.d = yl.d = xh.d = yh.d = t.d = NULL;
    exx.d = eyy.d = exy.d = NULL;
    fxx.d = fyy.d = fxy.d = gxx.d = gyy.d = gxy.d = NULL;

    nbs_gauss_blur(x, &xl, 3.0, NBS_BORDER_SYMMETRIC);
    nbs_gauss_blur(y, &yl, 3.0, NBS_BORDER_SYMMETRIC);
    nbs_plane_alloc(&xh, x->w, x->h);
    nbs_plane_alloc(&yh, x->w, x->h);
    n = x->w * x->h;
    for (i = 0; i < n; ++i) {
        xh.d[i] = x->d[i] - xl.d[i];
        yh.d[i] = y->d[i] - yl.d[i];
    }

    nbs_ssim_window(&win, 11, 1.5);
    vn = (x->w - 2 * win.r) * (x->h - 2 * win.r);

    nbs_plane_alloc(&t, x->w, x->h);
    nbs_plane_alloc(&exx, x->w - 2 * win.r, x->h - 2 * win.r);
    nbs_plane_alloc(&eyy, x->w - 2 * win.r, x->h - 2 * win.r);
    nbs_plane_alloc(&exy, x->w - 2 * win.r, x->h - 2 * win.r);
    nbs_plane_alloc(&fxx, x->w - 2 * win.r, x->h - 2 * win.r);
    nbs_plane_alloc(&fyy, x->w - 2 * win.r, x->h - 2 * win.r);
    nbs_plane_alloc(&fxy, x->w - 2 * win.r, x->h - 2 * win.r);

    for (i = 0; i < n; ++i) t.d[i] = xl.d[i] * xl.d[i];
    nbs_conv_valid_sep(&t, &win, &exx);
    for (i = 0; i < n; ++i) t.d[i] = yl.d[i] * yl.d[i];
    nbs_conv_valid_sep(&t, &win, &eyy);
    for (i = 0; i < n; ++i) t.d[i] = xl.d[i] * yl.d[i];
    nbs_conv_valid_sep(&t, &win, &exy);

    for (i = 0; i < n; ++i) t.d[i] = xh.d[i] * xh.d[i];
    nbs_conv_valid_sep(&t, &win, &fxx);
    for (i = 0; i < n; ++i) t.d[i] = yh.d[i] * yh.d[i];
    nbs_conv_valid_sep(&t, &win, &fyy);
    for (i = 0; i < n; ++i) t.d[i] = xh.d[i] * yh.d[i];
    nbs_conv_valid_sep(&t, &win, &fxy);

    C1 = 0.01 * 255.0; C1 *= C1;
    C2 = 0.03 * 255.0; C2 *= C2;
    for (i = 0; i < vn; ++i) {
        double xiL = (2.0 * exy.d[i] + C1) / (exx.d[i] + eyy.d[i] + C1);
        double xiH = (2.0 * fxy.d[i] + C2) / (fxx.d[i] + fyy.d[i] + C2);
        s += xiL * xiH;
    }
    s /= (double)vn;

    nbs_kernel_free(&win);
    nbs_plane_free(&xl); nbs_plane_free(&yl);
    nbs_plane_free(&xh); nbs_plane_free(&yh);
    nbs_plane_free(&t);
    nbs_plane_free(&exx); nbs_plane_free(&eyy); nbs_plane_free(&exy);
    nbs_plane_free(&fxx); nbs_plane_free(&fyy); nbs_plane_free(&fxy);
    nbs_plane_free(&gxx); nbs_plane_free(&gyy); nbs_plane_free(&gxy);
    return s;
}

static void test_ei2023_two_band(void)
{
    nbs_plane x, y;
    nbs_ssim_cfg scfg;
    double worst = 0.0, rms = 0.0;
    int t, i, n;

    printf("EI-2023 two-band identity\n");
    nbs_plane_alloc(&x, 256, 256);
    fill_texture(&x, 0.0);
    nbs_plane_alloc(&y, 256, 256);
    n = x.w * x.h;
    nbs_ssim_cfg_default(&scfg);
    scfg.downsample = 0;

    for (t = 0; t < 7; ++t) {
        double ssim_ref, model, d;
        if (t < 4) {
            double amp = 4.0 * (double)(t + 1) * (double)(t + 1);
            for (i = 0; i < n; ++i) y.d[i] = x.d[i] + amp * (rnd() - 0.5);
        } else {
            nbs_plane b;
            b.d = NULL;
            nbs_gauss_blur(&x, &b, 0.5 * (double)(t - 3), NBS_BORDER_SYMMETRIC);
            for (i = 0; i < n; ++i) y.d[i] = b.d[i];
            nbs_plane_free(&b);
        }
        nbs_ssim(&x, &y, &scfg, &ssim_ref);
        model = two_band_model(&x, &y);
        d = fabs(ssim_ref - model);
        rms += d * d;
        if (d > worst) worst = d;
    }
    rms = sqrt(rms / 7.0);
    printf("  worst |SSIM - two-band| = %.4f, rms = %.4f\n", worst, rms);
    ck(worst < 0.05, "two-band model matches SSIM to better than 0.05");
    ck(rms < 0.03, "two-band model rms agreement with SSIM better than 0.03");

    nbs_plane_free(&x); nbs_plane_free(&y);
}

static void test_stats(void)
{
    double a[10], b[10], c[10];
    int i;

    printf("statistics\n");
    for (i = 0; i < 10; ++i) { a[i] = (double)i; b[i] = 2.0 * i + 1.0; c[i] = -3.0 * i; }
    ck_near(nbs_pearson(a, b, 10), 1.0, 1e-12, "pearson +1");
    ck_near(nbs_pearson(a, c, 10), -1.0, 1e-12, "pearson -1");
    ck_near(nbs_spearman(a, b, 10), 1.0, 1e-12, "spearman +1");
    ck_near(nbs_kendall(a, b, 10), 1.0, 1e-12, "kendall +1");
    ck_near(nbs_kendall(a, c, 10), -1.0, 1e-12, "kendall -1");

    /* Ties. */
    {
        double t1[6], t2[6];
        t1[0]=1; t1[1]=1; t1[2]=2; t1[3]=3; t1[4]=3; t1[5]=4;
        t2[0]=1; t2[1]=1; t2[2]=2; t2[3]=3; t2[4]=3; t2[5]=4;
        ck_near(nbs_spearman(t1, t2, 6), 1.0, 1e-12, "spearman with ties");
        ck_near(nbs_kendall(t1, t2, 6), 1.0, 1e-12, "kendall tau-b with ties");
    }

    /* Logistic fit must undo a monotone nonlinearity. */
    {
        double x[40], y[40];
        nbs_fitreport r;
        for (i = 0; i < 40; ++i) {
            x[i] = 0.02 + 0.024 * i;
            y[i] = 100.0 / (1.0 + exp(8.0 * (x[i] - 0.5))) + 0.5 * (rnd() - 0.5);
        }
        nbs_evaluate(x, y, 40, &r);
        ck(fabs(r.srcc) > 0.999, "srcc on monotone data");
        ck(fabs(r.plcc) > 0.995, "logistic fit recovers linearity");
    }
}

static void test_png_roundtrip(void)
{
    /* A 2x2 truecolour PNG, stored (uncompressed deflate blocks) so we can
     * embed it literally.  Pixels: red, green, blue, white. */
    static const unsigned char png[] = {
        0x89,0x50,0x4e,0x47,0x0d,0x0a,0x1a,0x0a,0x00,0x00,0x00,0x0d,
        0x49,0x48,0x44,0x52,0x00,0x00,0x00,0x02,0x00,0x00,0x00,0x02,
        0x08,0x02,0x00,0x00,0x00,0xfd,0xd4,0x9a,0x73,0x00,0x00,0x00,
        0x12,0x49,0x44,0x41,0x54,0x78,0xda,0x63,0xf8,0xcf,0xc0,0xc0,
        0x00,0xc2,0x0c,0xff,0x81,0x00,0x00,0x1f,0xee,0x05,0xfb,0xf1,
        0xab,0xba,0x77,0x00,0x00,0x00,0x00,0x49,0x45,0x4e,0x44,0xae,
        0x42,0x60,0x82
    };
    nbs_image8 im;
    FILE *f;
    const char *path = "/tmp/nbstest_tmp.png";

    printf("png decoder\n");
    f = fopen(path, "wb");
    if (f == NULL) { ck(0, "cannot create temp png"); return; }
    fwrite(png, 1, sizeof(png), f);
    fclose(f);

    if (nbs_image8_load(path, &im) != 0) {
        printf("  (png load failed: %s)\n", nbs_last_error());
        ck(0, "png loads");
        remove(path);
        return;
    }
    ck(im.w == 2 && im.h == 2 && im.nch == 3, "png dimensions");
    ck(im.px[0] == 255 && im.px[1] == 0 && im.px[2] == 0, "png pixel 0 red");
    ck(im.px[3] == 0 && im.px[4] == 255 && im.px[5] == 0, "png pixel 1 green");
    ck(im.px[6] == 0 && im.px[7] == 0 && im.px[8] == 255, "png pixel 2 blue");
    ck(im.px[9] == 255 && im.px[10] == 255 && im.px[11] == 255, "png pixel 3 white");
    nbs_image8_free(&im);
    remove(path);
}

static void test_pnm_bmp(void)
{
    const char *pp = "/tmp/nbstest_tmp.ppm";
    const char *bp = "/tmp/nbstest_tmp.bmp";
    FILE *f;
    nbs_image8 im;
    unsigned char bmp[54 + 3 * 4];
    int i;

    printf("pnm / bmp decoders\n");
    f = fopen(pp, "wb");
    if (f == NULL) { ck(0, "cannot create temp ppm"); return; }
    fprintf(f, "P6\n2 2\n255\n");
    fputc(255, f); fputc(0, f); fputc(0, f);
    fputc(0, f); fputc(255, f); fputc(0, f);
    fputc(0, f); fputc(0, f); fputc(255, f);
    fputc(255, f); fputc(255, f); fputc(255, f);
    fclose(f);
    ck(nbs_image8_load(pp, &im) == 0, "ppm loads");
    ck(im.w == 2 && im.h == 2 && im.nch == 3, "ppm dimensions");
    ck(im.px[0] == 255 && im.px[1] == 0, "ppm first pixel");
    nbs_image8_free(&im);
    remove(pp);

    /* 2x1 24-bit BMP, bottom-up, 8-byte padded rows. */
    memset(bmp, 0, sizeof(bmp));
    bmp[0] = 'B'; bmp[1] = 'M';
    bmp[2] = (unsigned char)(54 + 8);
    bmp[10] = 54;
    bmp[14] = 40;
    bmp[18] = 2;                 /* width  = 2 */
    bmp[22] = 1;                 /* height = 1 */
    bmp[26] = 1;                 /* planes */
    bmp[28] = 24;                /* bpp */
    /* BGR: pixel0 = red, pixel1 = blue */
    bmp[54 + 0] = 0;   bmp[54 + 1] = 0;   bmp[54 + 2] = 255;
    bmp[54 + 3] = 255; bmp[54 + 4] = 0;   bmp[54 + 5] = 0;
    f = fopen(bp, "wb");
    if (f == NULL) { ck(0, "cannot create temp bmp"); return; }
    fwrite(bmp, 1, 54, f);
    for (i = 0; i < 8; ++i) fputc(bmp[54 + i], f);
    fclose(f);
    ck(nbs_image8_load(bp, &im) == 0, "bmp loads");
    if (im.px != NULL) {
        ck(im.w == 2 && im.h == 1 && im.nch == 3, "bmp dimensions");
        ck(im.px[0] == 255 && im.px[1] == 0 && im.px[2] == 0, "bmp pixel 0 red");
        ck(im.px[3] == 0 && im.px[4] == 0 && im.px[5] == 255, "bmp pixel 1 blue");
    }
    nbs_image8_free(&im);
    remove(bp);
}

static void test_logxi_consistency(void)
{
    nbs_plane x, y;
    nbs_nbssim_cfg cfg;
    double direct, viacache = 0.0;
    double *buf;
    int npix = 0, p, k, n, i, nw;

    printf("log-xi cache consistency\n");
    nbs_plane_alloc(&x, 200, 180);
    fill_texture(&x, 3.0);
    nbs_plane_copy(&y, &x);
    n = x.w * x.h;
    for (i = 0; i < n; ++i) y.d[i] += 12.0 * (rnd() - 0.5);

    nbs_nbssim_cfg_default(&cfg);
    cfg.bank.ppd = 20.0;
    nw = nbs_w_count(cfg.bank.nbands);

    ck(nbs_nbssim(&x, &y, &cfg, &direct) == 0, "direct score");
    ck(nbs_nbssim_logxi(&x, &y, &cfg, 1, NULL, &npix) == 0, "logxi size query");
    buf = (double *)malloc((size_t)npix * (size_t)nw * sizeof(double));
    ck(buf != NULL, "logxi buffer");
    if (buf == NULL) { nbs_plane_free(&x); nbs_plane_free(&y); return; }
    ck(nbs_nbssim_logxi(&x, &y, &cfg, 1, buf, &npix) == 0, "logxi fill");

    for (p = 0; p < npix; ++p) {
        double s = 0.0;
        for (k = 0; k < nw; ++k) s += cfg.w[k] * buf[(size_t)p * (size_t)nw + (size_t)k];
        viacache += exp(s);
    }
    viacache /= (double)npix;
    ck_near(direct, viacache, 1e-12, "cached path reproduces the direct score");

    free(buf);
    nbs_plane_free(&x); nbs_plane_free(&y);
}

/* nbs_save_weights writes a '#' header; nbs_load_weights must skip it.
 * These two lived in the same file and disagreed on the format. */
static void test_weight_file_roundtrip(void)
{
    nbs_nbssim_cfg a, b;
    const char *path = "/tmp/nbstest_w.txt";
    int i, j, n, ok = 1;

    printf("weight file round-trip\n");
    nbs_nbssim_cfg_default(&a);
    nbs_nbssim_cfg_default(&b);
    n = a.bank.nbands;
    for (i = 0; i < n; ++i)
        for (j = i; j < n; ++j)
            a.w[nbs_w_index(i, j, n)] = 0.1 + 0.37 * (double)(i * n + j);

    ck(nbs_save_weights(path, &a) == 0, "weights save");
    ck(nbs_load_weights(path, b.w, n) == 0, "weights load");
    for (i = 0; i < n; ++i)
        for (j = i; j < n; ++j) {
            int k = nbs_w_index(i, j, n);
            if (fabs(a.w[k] - b.w[k]) > 1e-9) ok = 0;
        }
    ck(ok, "saved weights survive a load unchanged");
    remove(path);
}


/* ---------------------------------------------------------------- */
/* DFT-domain core                                                    */
/* ---------------------------------------------------------------- */

static void test_fft(void)
{
    int sizes[6];
    int s, n, k, t;
    printf("fft\n");

    ck(nbs_fft_size(1) == 1, "next pow2 of 1");
    ck(nbs_fft_size(3) == 4, "next pow2 of 3");
    ck(nbs_fft_size(800) == 1024, "next pow2 of 800");
    ck(nbs_fft_size(1024) == 1024, "next pow2 of 1024");

    sizes[0] = 2; sizes[1] = 4; sizes[2] = 8;
    sizes[3] = 16; sizes[4] = 64; sizes[5] = 256;
    for (s = 0; s < 6; ++s) {
        double *re, *im, *r0, *i0;
        double emax = 0.0, rmax = 0.0;
        nbs_fft_plan p;
        n = sizes[s];
        re = (double *)malloc((size_t)n * sizeof(double));
        im = (double *)malloc((size_t)n * sizeof(double));
        r0 = (double *)malloc((size_t)n * sizeof(double));
        i0 = (double *)malloc((size_t)n * sizeof(double));
        ck(nbs_fft_plan_init(&p, n) == 0, "fft plan");
        for (k = 0; k < n; ++k) {
            re[k] = r0[k] = sin(0.7 * k) + 0.3 * cos(2.1 * k);
            im[k] = i0[k] = cos(0.4 * k);
        }
        nbs_fft(&p, re, im, 0);
        for (k = 0; k < n; ++k) {
            double ar = 0.0, ai = 0.0;
            for (t = 0; t < n; ++t) {
                double a = -2.0 * 3.14159265358979323846 * k * t / n;
                ar += r0[t] * cos(a) - i0[t] * sin(a);
                ai += r0[t] * sin(a) + i0[t] * cos(a);
            }
            if (fabs(ar - re[k]) > emax) emax = fabs(ar - re[k]);
            if (fabs(ai - im[k]) > emax) emax = fabs(ai - im[k]);
        }
        nbs_fft(&p, re, im, 1);
        for (k = 0; k < n; ++k) {
            if (fabs(re[k] - r0[k]) > rmax) rmax = fabs(re[k] - r0[k]);
            if (fabs(im[k] - i0[k]) > rmax) rmax = fabs(im[k] - i0[k]);
        }
        ck(emax < 1e-9, "fft matches the naive DFT");
        ck(rmax < 1e-12, "fft round-trip is the identity");
        nbs_fft_plan_free(&p);
        free(re); free(im); free(r0); free(i0);
    }
}

/* The point of the DFT core: the bank it realises is the bank that was
 * asked for.  The spatial core needs a numerical calibration to get close;
 * here the peaks, gains and bandwidths come out of the definition. */
static void test_band_analytic(void)
{
    nbs_opts o;
    static const double tgt_f[6]  = { 0.8, 1.7, 2.8, 4.0, 8.0, 16.0 };
    static const double tgt_bw[6] = { 2.5, 2.2, 1.9, 1.7, 1.5, 1.25 };
    int i;

    printf("analytic band responses\n");
    nbs_opts_default(&o);
    ck(nbs_opts_parse_one(&o, "--preset=wilson") == 1, "wilson preset");
    ck(nbs_opts_parse_one(&o, "--ppd=63.223") == 1, "ppd");
    ck(nbs_opts_finalise(&o) == 0, "finalise");
    ck(o.nbs.bank.nbands == 6, "six channels below Nyquist at ppd 63");

    for (i = 0; i < o.nbs.bank.nbands; ++i) {
        double at = 0.0, pk = nbs_band_peak(&o.nbs.bank, i, &at);
        ck_near(at * o.nbs.bank.ppd, tgt_f[i], 0.01, "analytic peak lands on target");
        ck_near(pk, 1.0, 1e-6, "unit peak gain after normalisation");
        ck_near(nbs_band_bw_oct(&o.nbs.bank, i), tgt_bw[i], 0.01,
                "analytic bandwidth matches Wilson");
    }

    /* The taper gives the channel compact support without moving it. */
    ck(nbs_opts_parse_one(&o, "--taper=1.0") == 1, "taper");
    ck(nbs_opts_finalise(&o) == 0, "finalise with taper");
    for (i = 0; i < o.nbs.bank.nbands; ++i) {
        double at = 0.0;
        nbs_band_peak(&o.nbs.bank, i, &at);
        ck_near(at * o.nbs.bank.ppd, tgt_f[i], 0.01, "taper leaves the peak alone");
        ck_near(nbs_band_bw_oct(&o.nbs.bank, i), tgt_bw[i], 0.02,
                "taper leaves the half-amplitude bandwidth alone");
        ck(nbs_band_edge(&o.nbs.bank, i, 1e-9)
           <= nbs_band_edge(&o.nbs.bank, i, 0.02) * 1.35 + 1e-9,
           "tapered band has compact support");
    }
}

/* The two cores are different implementations of one definition.  With
 * decimation off and the same evaluation region they must agree to within
 * the spatial core's truncation error, and the DFT core must be
 * insensitive to the support parameter that the spatial one depends on. */
static void test_dft_matches_spatial(void)
{
    nbs_plane x, y;
    nbs_nbssim_cfg cfg;
    double sp = 0.0, df = 0.0, dfdec = 0.0, sp2 = 0.0, df2 = 0.0;
    int i, n;

    printf("dft core against the spatial core\n");
    nbs_plane_alloc(&x, 192, 160);
    fill_texture(&x, 5.0);
    nbs_plane_copy(&y, &x);
    n = x.w * x.h;
    for (i = 0; i < n; ++i) y.d[i] += 9.0 * (rnd() - 0.5);

    nbs_nbssim_cfg_default(&cfg);
    cfg.bank.ppd = 20.0;
    cfg.margin = 40;

    cfg.transform = NBS_XFORM_SPATIAL;
    ck(nbs_nbssim(&x, &y, &cfg, &sp) == 0, "spatial score");
    cfg.transform = NBS_XFORM_DFT;
    cfg.decimate = 0;
    ck(nbs_nbssim(&x, &y, &cfg, &df) == 0, "dft score");
    ck_near(sp, df, 5e-3, "dft agrees with spatial over the same pixels");

    /* Decimation changes the pooling scale by design, so it is allowed to
     * move the score -- but not far, and not erratically. */
    cfg.decimate = 1;
    ck(nbs_nbssim(&x, &y, &cfg, &dfdec) == 0, "dft decimated score");
    ck_near(df, dfdec, 2e-2, "decimated score stays close to full rate");

    /* Halving the spatial support moves the spatial answer; it must not
     * move the DFT answer at all, because there is no support to halve. */
    cfg.transform = NBS_XFORM_SPATIAL;
    cfg.bank.radius_mult = 2.0;
    nbs_bank_set(&cfg.bank, NULL, 0, 0, 0);
    ck(nbs_nbssim(&x, &y, &cfg, &sp2) == 0, "spatial score, short support");
    cfg.transform = NBS_XFORM_DFT;
    cfg.decimate = 0;
    ck(nbs_nbssim(&x, &y, &cfg, &df2) == 0, "dft score, short support");
    ck_near(df, df2, 1e-12, "dft is independent of --radius-mult");
    ck(fabs(sp - sp2) > fabs(df - df2), "spatial is not");

    nbs_plane_free(&x); nbs_plane_free(&y);
}

/* Identical inputs must score exactly 1 whatever the exponents. */
static void test_dft_identity(void)
{
    nbs_plane x;
    nbs_nbssim_cfg cfg;
    double s = 0.0;

    printf("dft identity\n");
    nbs_plane_alloc(&x, 130, 97);
    fill_texture(&x, 1.0);
    nbs_nbssim_cfg_default(&cfg);
    cfg.bank.ppd = 24.0;
    cfg.transform = NBS_XFORM_DFT;
    ck(nbs_nbssim(&x, &x, &cfg, &s) == 0, "dft identity score");
    ck_near(s, 1.0, 1e-9, "identical images score 1");
    nbs_plane_free(&x);
}

/* The log-xi cache must reproduce the DFT score too, since the tuner
 * evaluates every candidate exponent vector through it. */
static void test_dft_logxi(void)
{
    nbs_plane x, y;
    nbs_nbssim_cfg cfg;
    double direct = 0.0, viacache = 0.0;
    double *buf;
    int npix = 0, p, k, n, i, nw;

    printf("dft log-xi cache\n");
    nbs_plane_alloc(&x, 150, 120);
    fill_texture(&x, 7.0);
    nbs_plane_copy(&y, &x);
    n = x.w * x.h;
    for (i = 0; i < n; ++i) y.d[i] += 15.0 * (rnd() - 0.5);

    nbs_nbssim_cfg_default(&cfg);
    cfg.bank.ppd = 22.0;
    cfg.transform = NBS_XFORM_DFT;
    nw = nbs_w_count(cfg.bank.nbands);

    ck(nbs_nbssim(&x, &y, &cfg, &direct) == 0, "dft direct score");
    ck(nbs_nbssim_logxi(&x, &y, &cfg, 1, NULL, &npix) == 0, "dft logxi size");
    buf = (double *)malloc((size_t)npix * (size_t)nw * sizeof(double));
    ck(buf != NULL, "dft logxi buffer");
    if (buf == NULL) { nbs_plane_free(&x); nbs_plane_free(&y); return; }
    ck(nbs_nbssim_logxi(&x, &y, &cfg, 1, buf, &npix) == 0, "dft logxi fill");
    for (p = 0; p < npix; ++p) {
        double s = 0.0;
        for (k = 0; k < nw; ++k) s += cfg.w[k] * buf[(size_t)p * (size_t)nw + (size_t)k];
        viacache += exp(s);
    }
    viacache /= (double)npix;
    ck_near(direct, viacache, 1e-12, "dft cached path reproduces the direct score");
    free(buf);
    nbs_plane_free(&x); nbs_plane_free(&y);
}


/* A vector fitted with six channels must be readable by the five-channel
 * bank that a shorter viewing distance leaves behind, and the surviving
 * exponents must be the same (i,j) entries -- not the first fifteen numbers
 * in the file, which is a different matrix entirely. */
static void test_weight_file_band_mapping(void)
{
    nbs_nbssim_cfg six, five;
    const char *path = "/tmp/nbstest_w6.txt";
    int i, j;

    printf("weight file across band counts\n");
    nbs_nbssim_cfg_default(&six);
    six.bank.nbands = 6;
    for (i = 0; i < 6; ++i)
        for (j = i; j < 6; ++j)
            six.w[nbs_w_index(i, j, 6)] = 100.0 * (i + 1) + (j + 1);
    ck(nbs_save_weights(path, &six) == 0, "save six-band exponents");

    nbs_nbssim_cfg_default(&five);
    five.bank.nbands = 5;
    for (i = 0; i < NBS_MAX_W; ++i) five.w[i] = -1.0;
    ck(nbs_load_weights(path, five.w, 5) == 0, "load into a five-band bank");
    for (i = 0; i < 5; ++i)
        for (j = i; j < 5; ++j)
            ck_near(five.w[nbs_w_index(i, j, 5)], 100.0 * (i + 1) + (j + 1), 1e-12,
                    "surviving exponent keeps its (i,j) identity");

    /* Going the other way is not defined and must be refused. */
    ck(nbs_save_weights(path, &five) == 0, "save five-band exponents");
    ck(nbs_load_weights(path, six.w, 6) != 0,
       "a five-band file cannot fill a six-band bank");
}


/* With all N^2 exponents free, setting w_ij = w_ji by hand must reproduce
 * the symmetric metric exactly -- otherwise the two packings disagree
 * about which sub-band an exponent belongs to.  And the metric must NOT be
 * invariant to transposing the exponent matrix, because sub-band (i,j) has
 * horizontal frequency f_i and vertical f_j: if transposing changed
 * nothing, the extra parameters would be measuring nothing. */
static void test_free_exponents(void)
{
    nbs_plane x, y;
    nbs_nbssim_cfg sym, freec;
    double a = 0.0, b = 0.0, c = 0.0;
    int i, j, n, np;

    printf("free (asymmetric) exponents\n");
    ck(nbs_w_count_ex(6, 1) == 21, "symmetric packing has 21 exponents");
    ck(nbs_w_count_ex(6, 0) == 36, "free packing has 36 exponents");
    for (i = 0; i < 6; ++i)
        for (j = 0; j < 6; ++j) {
            ck(nbs_w_index_ex(i, j, 6, 0) == i * 6 + j, "free index is row-major");
            ck(nbs_w_index_ex(i, j, 6, 1) == nbs_w_index_ex(j, i, 6, 1),
               "symmetric index ignores the order");
        }

    nbs_plane_alloc(&x, 160, 128);
    fill_texture(&x, 11.0);
    nbs_plane_copy(&y, &x);
    np = x.w * x.h;
    for (i = 0; i < np; ++i) y.d[i] += 11.0 * (rnd() - 0.5);

    nbs_nbssim_cfg_default(&sym);
    sym.bank.ppd = 20.0;
    sym.transform = NBS_XFORM_DFT;
    n = sym.bank.nbands;
    for (i = 0; i < n; ++i)
        for (j = i; j < n; ++j)
            sym.w[nbs_w_index_ex(i, j, n, 1)] = 0.02 + 0.011 * (double)(i + 2 * j);

    freec = sym;
    freec.w_sym = 0;
    for (i = 0; i < NBS_MAX_W; ++i) freec.w[i] = 0.0;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            freec.w[nbs_w_index_ex(i, j, n, 0)] = sym.w[nbs_w_index_ex(i, j, n, 1)];

    ck(nbs_nbssim(&x, &y, &sym, &a) == 0, "symmetric score");
    ck(nbs_nbssim(&x, &y, &freec, &b) == 0, "free score, mirrored exponents");
    ck_near(a, b, 1e-12, "a mirrored free matrix reproduces the symmetric metric");

    /* Now make it genuinely asymmetric and check it matters. */
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            freec.w[nbs_w_index_ex(i, j, n, 0)] =
                sym.w[nbs_w_index_ex(i, j, n, 1)] * (i < j ? 1.7 : 0.3);
    ck(nbs_nbssim(&x, &y, &freec, &c) == 0, "free score, asymmetric exponents");
    ck(fabs(c - a) > 1e-6,
       "transposing weight between (i,j) and (j,i) changes the score");

    nbs_plane_free(&x); nbs_plane_free(&y);
}

/* The free packing must survive a file round-trip, and a symmetric file
 * must be readable by a free bank (by mirroring) but not the reverse. */
static void test_free_weight_file(void)
{
    nbs_nbssim_cfg a, b;
    const char *pf = "/tmp/nbstest_free.txt";
    const char *ps = "/tmp/nbstest_sym.txt";
    int i, j, n = 6, bad;

    printf("free exponent file round-trip\n");
    nbs_nbssim_cfg_default(&a);
    a.bank.nbands = n; a.w_sym = 0;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            a.w[nbs_w_index_ex(i, j, n, 0)] = 10.0 * (i + 1) + (j + 1);
    ck(nbs_save_weights(pf, &a) == 0, "save free exponents");

    nbs_nbssim_cfg_default(&b);
    b.bank.nbands = n; b.w_sym = 0;
    for (i = 0; i < NBS_MAX_W; ++i) b.w[i] = -1.0;
    ck(nbs_load_weights_ex(pf, b.w, n, 0) == 0, "load free exponents");
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            ck_near(b.w[nbs_w_index_ex(i, j, n, 0)], 10.0 * (i + 1) + (j + 1),
                    1e-12, "free exponent survives the round trip");

    /* Symmetric file into a free bank: must mirror. */
    nbs_nbssim_cfg_default(&a);
    a.bank.nbands = n; a.w_sym = 1;
    for (i = 0; i < n; ++i)
        for (j = i; j < n; ++j)
            a.w[nbs_w_index_ex(i, j, n, 1)] = 100.0 * (i + 1) + (j + 1);
    ck(nbs_save_weights(ps, &a) == 0, "save symmetric exponents");
    for (i = 0; i < NBS_MAX_W; ++i) b.w[i] = -1.0;
    ck(nbs_load_weights_ex(ps, b.w, n, 0) == 0, "symmetric file into a free bank");
    bad = 0;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j) {
            int lo = (i < j) ? i : j, hi = (i < j) ? j : i;
            if (fabs(b.w[nbs_w_index_ex(i, j, n, 0)]
                     - (100.0 * (lo + 1) + (hi + 1))) > 1e-12) bad++;
        }
    ck(bad == 0, "symmetric file mirrors into both triangles");

    /* Free file into a symmetric bank: must be refused. */
    ck(nbs_load_weights_ex(pf, b.w, n, 1) != 0,
       "an asymmetric file cannot fill a symmetric bank");
}

int main(void)
{
    printf("nbssim unit tests\n\n");
    test_kernels();
    test_conv_dc();
    test_bank_dog();
    test_bank_wilson();
    test_bank_cumgauss_partition();
    test_bank_ppd_invariance();
    test_ssim_identity();
    test_symmetry_and_monotonicity();
    test_two_band_tracks_ssim();
    test_ei2023_two_band();
    test_stats();
    test_png_roundtrip();
    test_pnm_bmp();
    test_logxi_consistency();
    test_weight_file_roundtrip();
    test_weight_file_band_mapping();
    test_free_exponents();
    test_free_weight_file();
    test_fft();
    test_band_analytic();
    test_dft_matches_spatial();
    test_dft_identity();
    test_dft_logxi();

    printf("\n%d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
