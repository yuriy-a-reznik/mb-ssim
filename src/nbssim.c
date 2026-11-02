/* nbssim.c -- the proposed N-band filterbank SSIM-like metric.
 *
 * Background
 * ----------
 * Reznik (Electronic Imaging 2023, "Another look at SSIM image quality
 * metric") showed that SSIM is well approximated by
 *
 *     SSIM(x,y) ~= xi(x_L, y_L) . xi(x_H, y_H),
 *     xi(a,b)    = (2 E[ab] + C) / (E[a^2] + E[b^2] + C)
 *
 * i.e. one and the same distance operator applied to a low-pass and a
 * high-pass version of the two images, with E[.] taken over an 11x11
 * Gaussian window.  Reznik (MHV'24, "SSIM Revisited: From Structural
 * Similarity to Error Visibility") proposed the N-band generalisation
 *
 *     NB-SSIM(x,y) = prod_i xi^{w_i}(x_i, y_i)
 *
 * with a sub-band structure resembling the spatial channels of the contrast
 * sensitivity function, and CSF-motivated exponents w_i.
 *
 * This file implements that metric with a separable two-dimensional bank:
 * N one-dimensional bands per axis give N x N two-dimensional sub-bands,
 * and the exponents are constrained to be symmetric under transposition of
 * the axes, w_ij == w_ji, leaving N(N+1)/2 free parameters.  The final
 * score is
 *
 *     NB-SSIM = mean_p  prod_{i,j} xi_ij(p)^{w_ij}
 *
 * i.e. the sub-band scores are combined per spatial location and only then
 * averaged over the image, as the MHV'24 abstract specifies.
 *
 * Band construction
 * -----------------
 * All banks are built from sampled Gaussians.  With the convention of the
 * 2023 paper a Gaussian of pixel standard deviation sigma has corner
 * frequency f = ppd / (2 pi sigma), i.e. sigma = ppd / (2 pi f), and
 * H(f_corner) = exp(-1/2).
 *
 *   NBS_BAND_LP:  h = G(sigma(f))
 *   NBS_BAND_HP:  h = delta - G(sigma(f))
 *   NBS_BAND_BP:  h = G(sigma(a)) - G(sigma(b)),  a < b,
 *                 chosen so the response peaks exactly at f:
 *
 *                     f^2 = 4 ln(rho) / (1/a^2 - 1/b^2),  b = rho a
 *                 =>  a   = f sqrt( (rho^2 - 1) / (4 rho^2 ln rho) )
 *
 * so the peak location is exact and rho controls the bandwidth.  For N = 2
 * the bank is exactly {G(sigma), delta - G(sigma)}: the two-band model of
 * the 2023 paper.
 *
 * Strict C89.
 */

#include "nbs.h"
#include "nbs_metrics.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/* from ssim.c */
int nbs_ssim_window(nbs_kernel *k, int win, double sigma);

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Weight indexing                                                     */
/*                                                                      */
/* Sub-band (i,j) is filtered with channel i along the rows and channel */
/* j along the columns, so i indexes horizontal frequency and j         */
/* vertical.  Reflecting across the main diagonal exchanges the two     */
/* cardinal axes; the diagonal entries w_ii are the sub-bands whose 2-D */
/* centre frequency lies at 45 degrees, i.e. the oblique ones.          */
/* ------------------------------------------------------------------ */

double nbs_xi_ratio(double exy, double exx, double eyy, double C, double p)
{
    double m;
    if (p == 1.0) return (2.0 * exy + C) / (exx + eyy + C);
    if (exx < 0.0) exx = 0.0;
    if (eyy < 0.0) eyy = 0.0;
    if (p == 0.0) m = sqrt(exx * eyy);
    else if (p > 1e6) m = (exx > eyy) ? exx : eyy;
    else if (p < -1e6) m = (exx < eyy) ? exx : eyy;
    else m = pow(0.5 * (pow(exx, p) + pow(eyy, p)), 1.0 / p);
    return (2.0 * exy + C) / (2.0 * m + C);
}

int nbs_w_count_ex(int nbands, int sym)
{
    return sym ? nbands * (nbands + 1) / 2 : nbands * nbands;
}

/* Symmetric: compact upper-triangular row-major packing,
 * (0,0),(0,1),...,(0,n-1),(1,1),...,(1,n-1),(2,2),...
 * Free: plain row-major, i*n + j. */
int nbs_w_index_ex(int i, int j, int n, int sym)
{
    int a, b;
    if (!sym) return i * n + j;
    a = (i < j) ? i : j;
    b = (i < j) ? j : i;
    return a * n - (a * (a - 1)) / 2 + (b - a);
}

int nbs_w_count(int nbands)          { return nbs_w_count_ex(nbands, 1); }
int nbs_w_index(int i, int j, int n) { return nbs_w_index_ex(i, j, n, 1); }

/* ------------------------------------------------------------------ */
/* Bank configuration                                                  */
/* ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ */
/* Discrete Gaussian response and band-pass calibration                */
/* ------------------------------------------------------------------ */

static double sigma_of(double ppd, double f)
{
    return ppd / (2.0 * M_PI * f);
}

/* DTFT of a unit-sum sampled Gaussian, via Poisson summation.  Sampling
 * aliases the continuous Gaussian, and for sigma below about one pixel the
 * difference is large enough to shift a difference-of-Gaussians peak
 * noticeably -- which is exactly why the band-pass corners are calibrated
 * against this function rather than the continuous formula. */
static double gauss_dtft(double nu, double sigma)
{
    double a = 2.0 * M_PI * M_PI * sigma * sigma;
    double num = 0.0, den = 0.0;
    int k;
    if (sigma <= 0.0) return 1.0;
    for (k = -3; k <= 3; ++k) {
        double t = nu + (double)k;
        num += exp(-a * t * t);
        den += exp(-a * (double)k * (double)k);
    }
    return num / den;
}

/* Peak location (in cycles/pixel) of the discrete DoG with corners a < b. */
static double dog_peak_discrete(double ppd, double a, double rho)
{
    double sa = sigma_of(ppd, a);
    double sb = sigma_of(ppd, rho * a);
    double best = -1e30, bnu = 0.0;
    int i;
    for (i = 1; i <= 2000; ++i) {
        double nu = 0.5 * (double)i / 2000.0;
        double h  = gauss_dtft(nu, sb) - gauss_dtft(nu, sa);
        if (h > best) { best = h; bnu = nu; }
    }
    /* Parabolic refinement around the grid maximum. */
    {
        double d = 0.5 / 2000.0;
        double n0 = bnu - d, n2 = bnu + d;
        double y0, y1, y2, den;
        if (n0 <= 0.0) n0 = 1e-6;
        if (n2 >= 0.5) n2 = 0.5;
        y0 = gauss_dtft(n0, sb) - gauss_dtft(n0, sa);
        y1 = best;
        y2 = gauss_dtft(n2, sb) - gauss_dtft(n2, sa);
        den = y0 - 2.0 * y1 + y2;
        if (fabs(den) > 1e-300) bnu = bnu - 0.5 * d * (y2 - y0) / den;
    }
    if (bnu < 1e-6) bnu = 1e-6;
    if (bnu > 0.5) bnu = 0.5;
    return bnu;
}

/* Analytic (continuous-domain) starting point:
 * a = f sqrt((rho^2-1) / (4 rho^2 ln rho)). */
static double dog_lower_corner(double f, double rho)
{
    double r2 = rho * rho;
    return f * sqrt((r2 - 1.0) / (4.0 * r2 * log(rho)));
}

/* Choose the lower corner so the *discrete* response peaks exactly at f. */
static double calibrate_dog(double ppd, double f, double rho)
{
    double target = f / ppd;
    double lo, hi, mid;
    int it;

    if (!(rho > 1.0)) rho = 3.0;
    mid = dog_lower_corner(f, rho);
    if (target >= 0.5) return mid;

    lo = mid * 0.2;
    hi = mid * 5.0;
    for (it = 0; it < 60; ++it) {
        double p;
        mid = 0.5 * (lo + hi);
        p = dog_peak_discrete(ppd, mid, rho);
        if (p < target) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

/* Half-height full bandwidth, in octaves, of the discrete DoG with lower
 * corner a and corner ratio rho.  Returns 0 if the band is not resolvable
 * (peak too close to DC or to Nyquist for both crossings to exist). */
static double dog_bandwidth_oct_n(double ppd, double a, double rho, int N)
{
    double sa = sigma_of(ppd, a), sb = sigma_of(ppd, rho * a);
    double best = -1e30, bnu = 0.0, half, lo = 0.0, hi = 0.0;
    int i;
    double prev = 0.0, prevnu = 0.0;

    for (i = 1; i <= N; ++i) {
        double nu = 0.5 * (double)i / (double)N;
        double h = gauss_dtft(nu, sb) - gauss_dtft(nu, sa);
        if (h > best) { best = h; bnu = nu; }
    }
    if (best <= 0.0) return 0.0;
    half = 0.5 * best;

    prev = 0.0; prevnu = 0.0;
    for (i = 1; i <= N; ++i) {
        double nu = 0.5 * (double)i / (double)N;
        double h = gauss_dtft(nu, sb) - gauss_dtft(nu, sa);
        if (lo == 0.0 && nu <= bnu && h >= half) {
            lo = (prev < half && nu > prevnu)
               ? prevnu + (nu - prevnu) * (half - prev) / (h - prev) : nu;
        }
        if (lo != 0.0 && hi == 0.0 && nu > bnu && h <= half) {
            hi = (prev > half) ? prevnu + (nu - prevnu) * (prev - half) / (prev - h) : nu;
        }
        prev = h; prevnu = nu;
    }
    if (lo <= 0.0 || hi <= 0.0 || hi <= lo) return 0.0;
    return log(hi / lo) / log(2.0);
}

static double dog_bandwidth_oct(double ppd, double a, double rho)
{
    return dog_bandwidth_oct_n(ppd, a, rho, 4000);
}

/* Solve for the corner ratio rho giving a target bandwidth at a target peak.
 * Bandwidth grows monotonically with rho, so bisection is safe.  The lower
 * corner is re-calibrated for the peak at every trial rho. */
static double solve_rho_for_bw(double ppd, double f, double target_bw)
{
    double lo = 1.02, hi = 40.0, mid = 3.0;
    int it;
    if (!(target_bw > 0.0)) return 3.0;
    /* Bisection only needs the bandwidth to a fraction of a percent, and
     * this loop is the whole cost of configuring a bank -- a 4000-point
     * sweep inside it cost five seconds of start-up per process.  Sweep
     * coarsely here; the reported bandwidth is measured at full resolution
     * afterwards. */
    for (it = 0; it < 50; ++it) {
        double a, bw;
        mid = 0.5 * (lo + hi);
        /* Analytic lower corner here, not the discrete calibration: the
         * bandwidth-versus-rho curve is what is being inverted, and the
         * discrete correction shifts it by far less than the bisection
         * tolerance.  Nesting a 2000-point peak search inside this loop
         * cost several seconds every time a bank was configured. */
        a = dog_lower_corner(f, mid);
        bw = dog_bandwidth_oct_n(ppd, a, mid, 512);
        if (bw <= 0.0) { hi = mid; continue; }   /* unresolvable -> narrow it */
        if (bw < target_bw) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

int nbs_bank_set(nbs_bank_cfg *bank, const double *freq, int nfreq,
                 int first_is_bandpass, int last_is_bandpass)
{
    int n = bank->nbands, i;

    if (n < 2 || n > NBS_MAX_BANDS) { nbs_set_error("bank: bad band count"); return 1; }
    if (nfreq != n && nfreq != n - 1) {
        nbs_set_error("bank: expected nbands or nbands-1 frequencies");
        return 1;
    }
    for (i = 0; i < nfreq; ++i) {
        if (!(freq[i] > 0.0)) { nbs_set_error("bank: frequencies must be positive"); return 1; }
        if (i > 0 && !(freq[i] > freq[i - 1])) {
            nbs_set_error("bank: frequencies must be strictly ascending");
            return 1;
        }
        bank->f[i] = freq[i];
    }
    if (nfreq == n - 1) {
        /* Derive the high-pass corner: place it just above the last peak. */
        bank->f[n - 1] = freq[n - 2] * sqrt(2.0);
    }

    for (i = 0; i < n; ++i) {
        bank->role[i] = NBS_BAND_BP;
        if (bank->rho[i] <= 1.0) bank->rho[i] = 3.0;
    }
    if (!first_is_bandpass) bank->role[0]     = NBS_BAND_LP;
    if (!last_is_bandpass)  bank->role[n - 1] = NBS_BAND_HP;

    if (!(bank->ppd > 0.0)) { nbs_set_error("bank: ppd must be set before the bands"); return 1; }

    /* Drop any band the display cannot present.  Keeping it would mean
     * "measuring" a channel whose centre frequency is beyond what the pixel
     * grid can carry, which is not a modelling choice but an error. */
    bank->nyquist_dropped = 0;
    {
        int keep = 0;
        for (i = 0; i < n; ++i) {
            if (bank->f[i] >= 0.5 * bank->ppd) {
                bank->dropped_f[bank->nyquist_dropped++] = bank->f[i];
                continue;
            }
            bank->f[keep]      = bank->f[i];
            bank->role[keep]   = bank->role[i];
            bank->rho[keep]    = bank->rho[i];
            bank->bw_oct[keep] = bank->bw_oct[i];
            keep++;
        }
        if (keep < 2) {
            nbs_set_error("bank: fewer than two bands remain below the display "
                          "Nyquist rate -- raise ppd or lower the band frequencies");
            return 1;
        }
        if (keep != n) {
            /* The last surviving band takes over the high-pass role if the
             * original high-pass was dropped. */
            if (!last_is_bandpass) bank->role[keep - 1] = NBS_BAND_HP;
            n = keep;
            bank->nbands = n;
        }
    }

    for (i = 0; i < n; ++i) {
        bank->a_lo[i] = 0.0;
        if (bank->kind != NBS_BANK_DOG || bank->role[i] != NBS_BAND_BP)
            continue;
        /* A difference of two Gaussians cannot be made narrower than about
         * 1.77 octaves -- as the corner ratio approaches 1 the shape tends to
         * the Laplacian of a Gaussian, whose bandwidth is fixed while its gain
         * collapses.  Narrower channels (Wilson's high-frequency mechanisms
         * are 1.25-1.5 octaves) need a log-Gabor bank instead. */
        if (bank->bw_oct[i] > 0.0)
            bank->rho[i] = solve_rho_for_bw(bank->ppd, bank->f[i], bank->bw_oct[i]);
        bank->a_lo[i] = calibrate_dog(bank->ppd, bank->f[i], bank->rho[i]);
    }
    return 0;
}

/* --- CUMGAUSS peak placement -------------------------------------- */

static double dog_peak_from_corners(double flo, double fhi)
{
    double num, den;
    if (!(fhi > flo) || flo <= 0.0) return 0.0;
    num = 4.0 * log(fhi / flo);
    den = 1.0 / (flo * flo) - 1.0 / (fhi * fhi);
    if (den <= 0.0) return 0.0;
    return sqrt(num / den);
}

static int dog_solve(double flo, double target, double *fhi)
{
    double lo, hi, mid, v;
    int it;
    if (flo <= 0.0) return 1;
    if (target <= sqrt(2.0) * flo * (1.0 + 1e-9)) return 1;
    lo = flo * (1.0 + 1e-12);
    hi = flo * 2.0;
    for (it = 0; it < 200; ++it) {
        v = dog_peak_from_corners(flo, hi);
        if (v >= target) break;
        hi *= 2.0;
        if (hi > 1e9) return 1;
    }
    for (it = 0; it < 200; ++it) {
        mid = 0.5 * (lo + hi);
        v = dog_peak_from_corners(flo, mid);
        if (v < target) lo = mid; else hi = mid;
    }
    *fhi = 0.5 * (lo + hi);
    return 0;
}

int nbs_bank_set_from_peaks(nbs_bank_cfg *bank, const double *peaks)
{
    int k, n = bank->nbands;
    if (n < 2 || n > NBS_MAX_BANDS) { nbs_set_error("bank: bad band count"); return 1; }
    if (peaks[0] <= 0.0) { nbs_set_error("bank: low-pass corner must be positive"); return 1; }
    bank->fc[0] = peaks[0];
    for (k = 1; k <= n - 2; ++k) {
        double f;
        if (dog_solve(bank->fc[k - 1], peaks[k], &f) != 0) {
            nbs_set_error("cumgauss bank: requested peak unreachable "
                          "(must exceed sqrt(2) x the previous cut-off); "
                          "use the default DoG bank for closely spaced peaks");
            return 1;
        }
        bank->fc[k] = f;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Kernel construction                                                 */
/* ------------------------------------------------------------------ */

/* Add kernel `src` (radius rs) scaled by s into dst (radius rd >= rs). */
static void kadd(nbs_kernel *dst, const nbs_kernel *src, double s)
{
    int i;
    for (i = -src->r; i <= src->r; ++i)
        dst->tap[i + dst->r] += s * src->tap[i + src->r];
}

static int build_gauss_band(nbs_kernel *out, double ppd, double f, double rm,
                            nbs_band_role role, double rho, double a_lo)
{
    nbs_kernel ga, gb;
    double sa, sb, a, b;
    int r;

    ga.tap = NULL; gb.tap = NULL;

    if (role == NBS_BAND_LP) {
        return nbs_kernel_gaussian(out, sigma_of(ppd, f), rm);
    }
    if (role == NBS_BAND_HP) {
        if (nbs_kernel_gaussian(&ga, sigma_of(ppd, f), rm) != 0) return 1;
        if (nbs_kernel_alloc(out, ga.r) != 0) { nbs_kernel_free(&ga); return 1; }
        kadd(out, &ga, -1.0);
        out->tap[out->r] += 1.0;
        nbs_kernel_free(&ga);
        return 0;
    }

    /* Band-pass: difference of two Gaussians with corners a < b = rho*a. */
    if (!(rho > 1.0)) rho = 3.0;
    a = (a_lo > 0.0) ? a_lo : dog_lower_corner(f, rho);
    b = rho * a;
    sa = sigma_of(ppd, a);
    sb = sigma_of(ppd, b);
    if (nbs_kernel_gaussian(&ga, sa, rm) != 0) return 1;
    if (nbs_kernel_gaussian(&gb, sb, rm) != 0) { nbs_kernel_free(&ga); return 1; }
    r = (ga.r > gb.r) ? ga.r : gb.r;
    if (nbs_kernel_alloc(out, r) != 0) { nbs_kernel_free(&ga); nbs_kernel_free(&gb); return 1; }
    kadd(out, &gb, 1.0);    /* narrower sigma -> wider pass band */
    kadd(out, &ga, -1.0);
    nbs_kernel_free(&ga);
    nbs_kernel_free(&gb);
    return 0;
}

static double gauss_lp(double f, double f0)
{
    if (f0 <= 0.0) return 0.0;
    return exp(-(f * f) / (2.0 * f0 * f0));
}

/* Raised-cosine taper in log frequency: 1 out to 0.6*B octaves from the
 * centre, smoothly to 0 at B octaves.  B <= 0 disables it. */
static double log_taper(double t, double B)
{
    double t0, u;
    if (B <= 0.0) return 1.0;
    t = (t < 0.0) ? -t : t;
    t0 = 0.6 * B;
    if (t <= t0) return 1.0;
    if (t >= B)  return 0.0;
    u = (t - t0) / (B - t0);
    return 0.5 * (1.0 + cos(M_PI * u));
}

/* Log-Gabor magnitude, Gaussian on a log2 frequency axis, optionally
 * given compact support by the raised-cosine taper. */
static double loggabor_resp(double f, double f0, double bw_oct, double taper_mult)
{
    double s, t;
    if (f <= 0.0 || f0 <= 0.0) return 0.0;
    if (bw_oct <= 0.0) bw_oct = 1.0;
    s = bw_oct / (2.0 * sqrt(2.0 * log(2.0)));
    t = log(f / f0) / log(2.0);
    return exp(-(t * t) / (2.0 * s * s))
         * log_taper(t, (taper_mult > 0.0) ? taper_mult * bw_oct : 0.0);
}

/* The band's transfer function with an explicit centre/corner, so the
 * calibration search can move it without mutating the bank. */
static double band_H_at(const nbs_bank_cfg *bank, int b, double nu, double f0)
{
    double f = nu * bank->ppd;

    if (bank->kind == NBS_BANK_CUMGAUSS) {
        int n = bank->nbands;
        if (b == 0)     return gauss_lp(f, bank->fc[0]);
        if (b == n - 1) return 1.0 - gauss_lp(f, bank->fc[n - 2]);
        return gauss_lp(f, bank->fc[b]) - gauss_lp(f, bank->fc[b - 1]);
    }
    if (bank->role[b] == NBS_BAND_LP) return gauss_lp(f, f0);
    if (bank->role[b] == NBS_BAND_HP) return 1.0 - gauss_lp(f, f0);
    if (bank->kind == NBS_BANK_LOGGABOR)
        return loggabor_resp(f, f0, bank->bw_oct[b], bank->taper_mult);
    {
        double rho = (bank->rho[b] > 1.0) ? bank->rho[b] : 3.0;
        double a = (bank->a_lo[b] > 0.0) ? bank->a_lo[b]
                                         : dog_lower_corner(f0, rho);
        return gauss_lp(f, rho * a) - gauss_lp(f, a);
    }
}

static int build_loggabor_band_m(nbs_kernel *out, const nbs_bank_cfg *bank,
                                 int b, int r, int m, double f_use)
{
    int i, t;

    if (nbs_kernel_alloc(out, r) != 0) return 1;
    for (i = -r; i <= r; ++i) {
        double acc = 0.0;
        for (t = 0; t < m; ++t) {
            double nu = 0.5 * ((double)t + 0.5) / (double)m;
            acc += band_H_at(bank, b, nu, f_use) * cos(2.0 * M_PI * nu * (double)i);
        }
        acc *= 2.0 * (0.5 / (double)m);
        acc *= 0.54 + 0.46 * cos(M_PI * (double)i / (double)(r + 1));
        out->tap[i + r] = acc;
    }
    return 0;
}

static int build_loggabor_band(nbs_kernel *out, const nbs_bank_cfg *bank, int b, int r)
{
    double f_use = (bank->f_nom[b] > 0.0) ? bank->f_nom[b] : bank->f[b];
    return build_loggabor_band_m(out, bank, b, r, 4096, f_use);
}

/* Truncating and windowing a log-Gabor attenuates its low-frequency skirt,
 * which drags the realised peak upward -- by ~10% for a 2.5-octave channel at
 * 0.8 cpd when ppd is 63.  Solve for the nominal centre whose *realised* peak
 * lands where the model asks.  A coarse quadrature is enough here; the shipped
 * kernel is rebuilt at full resolution afterwards. */
static double calibrate_loggabor(const nbs_bank_cfg *bank, int b, int r)
{
    double lo, hi, mid = bank->f[b], target = bank->f[b];
    int it;
    lo = target * 0.35;
    hi = target * 1.2;
    for (it = 0; it < 30; ++it) {
        nbs_kernel k;
        double best = -1e30, bf = 0.0, nu;
        mid = 0.5 * (lo + hi);
        if (build_loggabor_band_m(&k, bank, b, r, 1024, mid) != 0) return target;
        for (nu = 1e-4; nu < 0.5; nu += 5e-4) {
            double rr = nbs_kernel_response(&k, nu);
            if (rr > best) { best = rr; bf = nu * bank->ppd; }
        }
        nbs_kernel_free(&k);
        if (bf < target) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

/* ------------------------------------------------------------------ */
/* Analytic transfer functions                                         */
/*                                                                      */
/* These are the definition of the bank.  Both the low-pass Gaussian    */
/* written as exp(-f^2/(2 f0^2)) in cycles/degree and the sampled       */
/* Gaussian of standard deviation ppd/(2 pi f0) used by the spatial     */
/* kernels have the same continuous response, so the two paths agree by */
/* construction and differ only by truncation and aliasing.             */
/* ------------------------------------------------------------------ */

double nbs_band_H(const nbs_bank_cfg *bank, int b, double nu)
{
    if (b < 0 || b >= bank->nbands) return 0.0;
    return band_H_at(bank, b, nu, bank->f[b]);
}

double nbs_band_peak(const nbs_bank_cfg *bank, int b, double *at_nu)
{
    double best = -1e30, bnu = 0.0, nu;
    int i, N = 20000;
    for (i = 0; i <= N; ++i) {
        double h;
        nu = 0.5 * (double)i / (double)N;
        h = nbs_band_H(bank, b, nu);
        if (h > best) { best = h; bnu = nu; }
    }
    if (at_nu) *at_nu = bnu;
    return best;
}

double nbs_band_edge(const nbs_bank_cfg *bank, int b, double eps)
{
    double peak = nbs_band_peak(bank, b, NULL), thr;
    int i, N = 20000;
    if (peak <= 0.0) return 0.5;
    thr = eps * peak;
    for (i = N; i >= 0; --i) {
        double nu = 0.5 * (double)i / (double)N;
        double h = nbs_band_H(bank, b, nu);
        if (h < 0.0) h = -h;
        if (h > thr) return nu;
    }
    return 0.0;
}

double nbs_band_bw_oct(const nbs_bank_cfg *bank, int b)
{
    double bnu = 0.0, peak = nbs_band_peak(bank, b, &bnu), half;
    double lo = 0.0, hi = 0.0, prev, prevnu, nu, h;
    int i, N = 20000;
    if (peak <= 0.0 || bnu <= 0.0) return 0.0;
    half = 0.5 * peak;
    prev = nbs_band_H(bank, b, 0.0); prevnu = 0.0;
    for (i = 1; i <= N; ++i) {
        nu = 0.5 * (double)i / (double)N;
        h  = nbs_band_H(bank, b, nu);
        if (lo == 0.0 && nu <= bnu && h >= half)
            lo = (prev < half) ? prevnu + (nu - prevnu) * (half - prev) / (h - prev) : nu;
        if (lo != 0.0 && hi == 0.0 && nu > bnu && h <= half)
            hi = (prev > half) ? prevnu + (nu - prevnu) * (prev - half) / (prev - h) : nu;
        prev = h; prevnu = nu;
    }
    if (lo <= 0.0 || hi <= lo) return 0.0;
    return log(hi / lo) / log(2.0);
}

int nbs_bank_kernels(const nbs_bank_cfg *bank, nbs_kernel *ks)
{
    int b, n = bank->nbands;
    double rm = (bank->radius_mult > 0.0) ? bank->radius_mult : 4.0;

    for (b = 0; b < n; ++b) { ks[b].tap = NULL; ks[b].r = 0; }

    if (bank->kind == NBS_BANK_CUMGAUSS) {
        nbs_kernel g[NBS_MAX_BANDS];
        int k;
        for (k = 0; k < NBS_MAX_BANDS; ++k) { g[k].tap = NULL; g[k].r = 0; }
        if (n < 2) { nbs_set_error("cumgauss bank needs at least 2 bands"); return 1; }
        for (k = 0; k < n - 1; ++k) {
            if (nbs_kernel_gaussian(&g[k], sigma_of(bank->ppd, bank->fc[k]), rm) != 0) {
                for (b = 0; b < k; ++b) nbs_kernel_free(&g[b]);
                return 1;
            }
        }
        for (b = 0; b < n; ++b) {
            int r;
            if (b == 0)            r = g[0].r;
            else if (b == n - 1)   r = g[n - 2].r;
            else                   r = (g[b].r > g[b - 1].r) ? g[b].r : g[b - 1].r;
            if (nbs_kernel_alloc(&ks[b], r) != 0) goto fail_cum;
            if (b == 0) {
                kadd(&ks[b], &g[0], 1.0);
            } else if (b == n - 1) {
                kadd(&ks[b], &g[n - 2], -1.0);
                ks[b].tap[ks[b].r] += 1.0;
            } else {
                kadd(&ks[b], &g[b], 1.0);
                kadd(&ks[b], &g[b - 1], -1.0);
            }
        }
        for (k = 0; k < NBS_MAX_BANDS; ++k) nbs_kernel_free(&g[k]);
        return 0;
    fail_cum:
        for (k = 0; k < NBS_MAX_BANDS; ++k) nbs_kernel_free(&g[k]);
        for (b = 0; b < n; ++b) nbs_kernel_free(&ks[b]);
        return 1;
    }

    if (bank->kind == NBS_BANK_LOGGABOR) {
        /* Per-band support.  The spatial envelope of a log-Gabor scales as
         * 1/(f0 * bw), so a channel at 16 cpd needs far fewer taps than one
         * at 0.8 cpd.  Sizing every band by the lowest one -- as an earlier
         * version did -- made the n^2 column passes cost 5x what they need
         * to.  The bandwidth factor is clamped at 1 so that the wide low
         * channels keep the support they had. */
        nbs_bank_cfg *mb = (nbs_bank_cfg *)bank;   /* fills the f_nom cache */
        for (b = 0; b < n; ++b) {
            double bwf = 1.0;
            int r;
            if (bank->role[b] == NBS_BAND_BP && bank->bw_oct[b] > 0.0) {
                bwf = 2.5 / bank->bw_oct[b];
                if (bwf < 1.0) bwf = 1.0;
            }
            r = (int)ceil(rm * sigma_of(bank->ppd, bank->f[b]) * 2.0 * bwf);
            if (r < 4) r = 4;
            if (r > 400) r = 400;
            if (mb->f_nom[b] <= 0.0 && bank->role[b] == NBS_BAND_BP)
                mb->f_nom[b] = calibrate_loggabor(bank, b, r);
            if (build_loggabor_band(&ks[b], bank, b, r) != 0) {
                for (--b; b >= 0; --b) nbs_kernel_free(&ks[b]);
                return 1;
            }
        }
        if (bank->norm_peak) {
            for (b = 0; b < n; ++b) {
                double best = 0.0, nu;
                int t;
                for (nu = 0.0; nu <= 0.5; nu += 1e-4) {
                    double rr = nbs_kernel_response(&ks[b], nu);
                    if (rr > best) best = rr;
                }
                if (best > 1e-9)
                    for (t = 0; t <= 2 * ks[b].r; ++t) ks[b].tap[t] /= best;
            }
        }
        return 0;
    }

    /* NBS_BANK_DOG */
    for (b = 0; b < n; ++b) {
        if (build_gauss_band(&ks[b], bank->ppd, bank->f[b], rm,
                             bank->role[b], bank->rho[b], bank->a_lo[b]) != 0) {
            for (--b; b >= 0; --b) nbs_kernel_free(&ks[b]);
            return 1;
        }
    }
    if (bank->norm_peak) {
        for (b = 0; b < n; ++b) {
            double best = 0.0, nu;
            int t;
            for (nu = 0.0; nu <= 0.5; nu += 1e-4) {
                double r = nbs_kernel_response(&ks[b], nu);
                if (r > best) best = r;
            }
            if (best > 1e-9) {
                for (t = 0; t <= 2 * ks[b].r; ++t) ks[b].tap[t] /= best;
            }
        }
    }
    return 0;
}

double nbs_kernel_response(const nbs_kernel *k, double nu)
{
    double acc = k->tap[k->r];
    int i;
    for (i = 1; i <= k->r; ++i)
        acc += (k->tap[k->r + i] + k->tap[k->r - i]) * cos(2.0 * M_PI * nu * (double)i);
    return acc;
}

/* ------------------------------------------------------------------ */
/* Configuration defaults                                              */
/* ------------------------------------------------------------------ */

void nbs_nbssim_cfg_default(nbs_nbssim_cfg *c)
{
    static const double freq[5] = { 0.8, 1.7, 2.8, 4.0, 8.0 };
    int i, nw;

    memset(c, 0, sizeof(*c));
    c->bank.kind        = NBS_BANK_DOG;
    c->bank.nbands      = 6;
    c->bank.ppd         = 30.0;
    c->bank.radius_mult = 4.0;
    c->bank.border      = NBS_BORDER_SYMMETRIC;
    c->bank.taper_mult  = 0.0;   /* no compact support unless asked */
    for (i = 0; i < NBS_MAX_BANDS; ++i) {
        c->bank.rho[i]    = 3.0;
        c->bank.bw_oct[i] = 1.0;
    }
    nbs_bank_set(&c->bank, freq, 5, 0, 0);

    c->win        = 11;
    c->win_sigma  = 1.5;
    c->L          = 255.0;
    c->xi_floor   = 1e-6;
    c->norm_p     = 1.0;
    c->downsample = 0;
    c->margin     = -1;        /* automatic */
    c->w_sym      = 1;
    c->transform  = NBS_XFORM_SPATIAL;
    c->decimate   = 1;
    c->decim_eps  = 0.02;

    c->kc[0] = 0.01;
    for (i = 1; i < NBS_MAX_BANDS; ++i) c->kc[i] = 0.03;

    /* Uniform start with total exponent 2, matching the two unit-exponent
     * factors of the original SSIM. */
    nw = c->bank.nbands * c->bank.nbands;
    for (i = 0; i < NBS_MAX_W; ++i) c->w[i] = 2.0 / (double)nw;
}

/* ------------------------------------------------------------------ */
/* Core evaluation                                                     */
/* ------------------------------------------------------------------ */

static int conv_axis(const nbs_plane *src, nbs_plane *dst,
                     const nbs_kernel *k, const nbs_kernel *id,
                     int axis, nbs_border border)
{
    if (axis == 0) return nbs_conv_sep(src, dst, k, id, border);
    return nbs_conv_sep(src, dst, id, k, border);
}

/* mode 0: mean weighted power product.
 * mode 1: fill the log-xi cache on a subsampled grid. */
static int nbssim_core(const nbs_plane *xin, const nbs_plane *yin,
                       const nbs_nbssim_cfg *cfg,
                       int mode, int stride,
                       double *score, double *cache, int *out_npix)
{
    nbs_kernel ks[NBS_MAX_BANDS];
    nbs_kernel win, id;
    nbs_plane rbx[NBS_MAX_BANDS], rby[NBS_MAX_BANDS];
    nbs_plane xs, ys, bxp, byp, exx, eyy, exy, tmp, acc;
    const nbs_plane *px, *py;
    int n, i, j, k, rc = 1, maxr = 0, f;
    int W, H, m0, sx, sy, nsamp = 0, nw, np;

    for (k = 0; k < NBS_MAX_BANDS; ++k) {
        ks[k].tap = NULL; ks[k].r = 0;
        rbx[k].d = NULL; rby[k].d = NULL;
    }
    win.tap = NULL; id.tap = NULL;
    xs.d = ys.d = bxp.d = byp.d = NULL;
    exx.d = eyy.d = exy.d = tmp.d = acc.d = NULL;

    n  = cfg->bank.nbands;
    nw = nbs_w_count_ex(n, cfg->w_sym);

    if (xin->w != yin->w || xin->h != yin->h) { nbs_set_error("nbssim: size mismatch"); return 1; }
    if (n < 2 || n > NBS_MAX_BANDS) { nbs_set_error("nbssim: bad band count"); return 1; }

    f = cfg->downsample;
    if (f < 0) {
        int mn = (xin->w < xin->h) ? xin->w : xin->h;
        f = (int)floor((double)mn / 256.0 + 0.5);
        if (f < 1) f = 1;
    }
    if (f > 1) {
        if (nbs_downsample(xin, &xs, f, NBS_BORDER_SYMMETRIC) != 0) goto done;
        if (nbs_downsample(yin, &ys, f, NBS_BORDER_SYMMETRIC) != 0) goto done;
        px = &xs; py = &ys;
    } else {
        px = xin; py = yin;
    }
    W = px->w; H = px->h;
    np = W * H;

    if (nbs_bank_kernels(&cfg->bank, ks) != 0) goto done;
    for (k = 0; k < n; ++k) if (ks[k].r > maxr) maxr = ks[k].r;
    if (nbs_ssim_window(&win, cfg->win, cfg->win_sigma) != 0) goto done;
    if (nbs_kernel_alloc(&id, 0) != 0) goto done;
    id.tap[0] = 1.0;

    /* Discard the margin contaminated by the border extension of the longest
     * analysis filter, plus the statistics window radius -- unless the caller
     * asked for a specific margin (0 = evaluate everywhere). */
    m0 = (cfg->margin >= 0) ? cfg->margin : (maxr + win.r);
    if (W - 2 * m0 <= 0 || H - 2 * m0 <= 0) {
        nbs_set_error("nbssim: image too small for this filterbank "
                      "(raise the low-pass corner, lower ppd, or lower "
                      "--radius-mult)");
        goto done;
    }

    /* Split along rows. */
    for (i = 0; i < n; ++i) {
        if (nbs_plane_alloc(&rbx[i], W, H) != 0) goto done;
        if (nbs_plane_alloc(&rby[i], W, H) != 0) goto done;
        if (conv_axis(px, &rbx[i], &ks[i], &id, 0, cfg->bank.border) != 0) goto done;
        if (conv_axis(py, &rby[i], &ks[i], &id, 0, cfg->bank.border) != 0) goto done;
    }

    if (nbs_plane_alloc(&bxp, W, H) != 0) goto done;
    if (nbs_plane_alloc(&byp, W, H) != 0) goto done;
    if (nbs_plane_alloc(&tmp, W, H) != 0) goto done;
    if (nbs_plane_alloc(&exx, W, H) != 0) goto done;
    if (nbs_plane_alloc(&eyy, W, H) != 0) goto done;
    if (nbs_plane_alloc(&exy, W, H) != 0) goto done;

    if (mode == 0) {
        if (nbs_plane_alloc(&acc, W, H) != 0) goto done;
        memset(acc.d, 0, (size_t)np * sizeof(double));
    } else {
        nsamp = 0;
        for (sy = m0; sy < H - m0; sy += stride)
            for (sx = m0; sx < W - m0; sx += stride)
                nsamp++;
        if (out_npix) *out_npix = nsamp;
        if (cache == NULL) { rc = 0; goto done; }
        memset(cache, 0, (size_t)nsamp * (size_t)nw * sizeof(double));
    }

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            double C, kk, wgt = 0.0;
            int wi = nbs_w_index_ex(i, j, n, cfg->w_sym);

            if (conv_axis(&rbx[i], &bxp, &ks[j], &id, 1, cfg->bank.border) != 0) goto done;
            if (conv_axis(&rby[i], &byp, &ks[j], &id, 1, cfg->bank.border) != 0) goto done;

            for (k = 0; k < np; ++k) tmp.d[k] = bxp.d[k] * bxp.d[k];
            if (nbs_conv_sep(&tmp, &exx, &win, &win, cfg->bank.border) != 0) goto done;
            for (k = 0; k < np; ++k) tmp.d[k] = byp.d[k] * byp.d[k];
            if (nbs_conv_sep(&tmp, &eyy, &win, &win, cfg->bank.border) != 0) goto done;
            for (k = 0; k < np; ++k) tmp.d[k] = bxp.d[k] * byp.d[k];
            if (nbs_conv_sep(&tmp, &exy, &win, &win, cfg->bank.border) != 0) goto done;

            kk = cfg->kc[(i > j) ? i : j];
            C  = kk * cfg->L; C *= C;
            if (mode == 0) wgt = cfg->w[wi];

            if (mode == 0) {
                for (sy = m0; sy < H - m0; ++sy) {
                    size_t base = (size_t)sy * (size_t)W;
                    for (sx = m0; sx < W - m0; ++sx) {
                        double xi = nbs_xi_ratio(exy.d[base + sx], exx.d[base + sx],
                                                 eyy.d[base + sx], C, cfg->norm_p);
                        if (xi < cfg->xi_floor) xi = cfg->xi_floor;
                        if (xi > 1.0) xi = 1.0;
                        acc.d[base + sx] += wgt * log(xi);
                    }
                }
            } else {
                int p = 0;
                for (sy = m0; sy < H - m0; sy += stride) {
                    size_t base = (size_t)sy * (size_t)W;
                    for (sx = m0; sx < W - m0; sx += stride) {
                        double xi = nbs_xi_ratio(exy.d[base + sx], exx.d[base + sx],
                                                 eyy.d[base + sx], C, cfg->norm_p);
                        if (xi < cfg->xi_floor) xi = cfg->xi_floor;
                        if (xi > 1.0) xi = 1.0;
                        cache[(size_t)p * (size_t)nw + (size_t)wi] += log(xi);
                        p++;
                    }
                }
            }
        }
    }

    if (mode == 0) {
        double s = 0.0;
        long cnt = 0;
        for (sy = m0; sy < H - m0; ++sy) {
            size_t base = (size_t)sy * (size_t)W;
            for (sx = m0; sx < W - m0; ++sx) { s += exp(acc.d[base + sx]); cnt++; }
        }
        *score = (cnt > 0) ? s / (double)cnt : 0.0;
    }
    rc = 0;

done:
    for (k = 0; k < NBS_MAX_BANDS; ++k) {
        nbs_kernel_free(&ks[k]);
        nbs_plane_free(&rbx[k]);
        nbs_plane_free(&rby[k]);
    }
    nbs_kernel_free(&win);
    nbs_kernel_free(&id);
    nbs_plane_free(&xs); nbs_plane_free(&ys);
    nbs_plane_free(&bxp); nbs_plane_free(&byp);
    nbs_plane_free(&exx); nbs_plane_free(&eyy); nbs_plane_free(&exy);
    nbs_plane_free(&tmp); nbs_plane_free(&acc);
    return rc;
}

int nbs_nbssim(const nbs_plane *x, const nbs_plane *y,
               const nbs_nbssim_cfg *cfg, double *out)
{
    if (cfg->transform == NBS_XFORM_DFT)
        return nbs_nbssim_dft(x, y, cfg, 0, 1, out, NULL, NULL);
    return nbssim_core(x, y, cfg, 0, 1, out, NULL, NULL);
}

int nbs_nbssim_logxi(const nbs_plane *x, const nbs_plane *y,
                     const nbs_nbssim_cfg *cfg, int stride,
                     double *out, int *out_npix)
{
    if (stride < 1) stride = 1;
    if (cfg->transform == NBS_XFORM_DFT)
        return nbs_nbssim_dft(x, y, cfg, 1, stride, NULL, out, out_npix);
    return nbssim_core(x, y, cfg, 1, stride, NULL, out, out_npix);
}
