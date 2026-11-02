/* nbs_metrics.h -- SSIM, MS-SSIM and the proposed N-band SSIM.
 *
 * Strict C89.
 */

#ifndef NBS_METRICS_H
#define NBS_METRICS_H

#include "nbs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================== */
/* Reference SSIM  (Wang, Bovik, Sheikh & Simoncelli, IEEE TIP 2004)   */
/* ================================================================== */

typedef struct {
    int    win;          /* window size, odd; 11 in the paper            */
    double win_sigma;    /* Gaussian window sigma; 1.5 in the paper      */
    double K1;           /* 0.01                                         */
    double K2;           /* 0.03                                         */
    double L;            /* dynamic range, 255 for 8-bit                 */
    int    downsample;   /* -1 = automatic (Wang's ssim.m rule),         */
                         /*  0 = none, n>0 = force factor n              */
} nbs_ssim_cfg;

void nbs_ssim_cfg_default(nbs_ssim_cfg *c);

int nbs_ssim(const nbs_plane *x, const nbs_plane *y,
             const nbs_ssim_cfg *cfg, double *out_mssim);

int nbs_ssim_ex(const nbs_plane *x, const nbs_plane *y,
                const nbs_ssim_cfg *cfg,
                double *out_mssim, double *out_mcs);

/* ================================================================== */
/* Reference MS-SSIM  (Wang, Simoncelli & Bovik, Asilomar 2003)        */
/* ================================================================== */

#define NBS_MSSSIM_MAX_LEVELS 8

typedef struct {
    int          levels;
    double       weight[NBS_MSSSIM_MAX_LEVELS];
    nbs_ssim_cfg base;
} nbs_msssim_cfg;

void nbs_msssim_cfg_default(nbs_msssim_cfg *c);

int nbs_msssim(const nbs_plane *x, const nbs_plane *y,
               const nbs_msssim_cfg *cfg, double *out);

/* ================================================================== */
/* Proposed N-band SSIM                                                */
/* ================================================================== */

#define NBS_MAX_BANDS 12
#define NBS_MAX_W (NBS_MAX_BANDS * NBS_MAX_BANDS)

/* Role of each 1-D band. */
typedef enum {
    NBS_BAND_LP = 0,   /* Gaussian low-pass, corner f            */
    NBS_BAND_BP = 1,   /* band-pass peaking at f                 */
    NBS_BAND_HP = 2    /* high-pass, corner f                    */
} nbs_band_role;

/* How band-pass sections are realised.
 *
 * NBS_BANK_DOG      each band-pass is a difference of two Gaussians whose
 *                   corner frequencies straddle the requested peak.  The
 *                   peak lands exactly where asked and the "sharpness"
 *                   rho (ratio of the two corners) sets the bandwidth.
 *                   For N = 2 this degenerates to exactly the two-band
 *                   low-pass / high-pass model of Reznik, EI 2023.
 *
 * NBS_BANK_CUMGAUSS nested Gaussian low-pass stages, bands taken as
 *                   successive differences.  Sums to unity exactly, but
 *                   band-pass peaks cannot be spaced closer than about
 *                   one octave, so it cannot realise arbitrary peak sets.
 *
 * NBS_BANK_LOGGABOR magnitude responses that are Gaussian on a log
 *                   frequency axis (the literal reading of slide 17),
 *                   realised by frequency sampling with a Hamming window.
 */
typedef enum {
    NBS_BANK_DOG      = 0,
    NBS_BANK_CUMGAUSS = 1,
    NBS_BANK_LOGGABOR = 2
} nbs_bank_kind;

typedef struct {
    nbs_bank_kind kind;
    int           nbands;
    double        ppd;                    /* pixels per degree            */

    nbs_band_role role[NBS_MAX_BANDS];
    double        f[NBS_MAX_BANDS];       /* corner or peak, cycles/degree */
    double        rho[NBS_MAX_BANDS];     /* DOG sharpness, > 1            */
    double        a_lo[NBS_MAX_BANDS];    /* calibrated lower DoG corner   */
    double        f_nom[NBS_MAX_BANDS];   /* calibrated log-Gabor centre   */
    /* Half-height bandwidth in octaves.  For LOGGABOR it is the Gaussian
     * width directly.  For DOG, a value > 0 overrides rho: the corner ratio
     * is solved numerically so the realised band has this bandwidth, which is
     * how psychophysical channel models specify their filters (Wilson,
     * McFarlane & Phillips 1983 give ~2-2.5 octaves at the low end narrowing
     * to ~1.25-1.5 at the high end).  0 means "use rho". */
    double        bw_oct[NBS_MAX_BANDS];

    /* Bands whose centre frequency reached or exceeded the display Nyquist
     * rate are dropped by nbs_bank_set, since the display cannot present
     * them; this records how many.  At CID:IQ's 50 cm distance (ppd 31.6,
     * Nyquist 15.8 cpd) a 16 cpd channel is one such band. */
    int           nyquist_dropped;
    double        dropped_f[NBS_MAX_BANDS];

    /* CUMGAUSS only: nbands-1 ascending cut-offs in cycles/degree. */
    double        fc[NBS_MAX_BANDS];

    /* Narrow band-pass sections have low peak gain (a 1-octave DoG peaks
     * near 0.2), so their sub-band energies sit far below the fixed
     * stabilising constant C and xi is driven towards 1 regardless of the
     * distortion -- the band stops carrying information.  Normalising each
     * band to unit peak gain removes that coupling.  Low-pass and high-pass
     * sections already have unit gain, so N = 2 is unaffected. */
    int           norm_peak;              /* 1 = scale each band to peak 1  */

    double        radius_mult;            /* Gaussian support in sigmas    */
    nbs_border    border;

    /* Frequency-domain taper, in units of bw_oct, applied to LOGGABOR
     * band-pass sections only.  A Gaussian on a log axis has heavy tails --
     * a 2.5-octave channel centred at 0.8 cpd is still at 2%% of peak by
     * 6 cpd -- so the channel cannot be critically sampled and the DFT core
     * cannot decimate it.  Multiplying by a raised cosine that reaches zero
     * at f0 * 2^(+-taper_mult*bw_oct) gives the channel compact support in
     * frequency, which is what Watson's cortex transform and Daly's VDP do
     * by construction.  0 disables the taper (infinite support).  */
    double        taper_mult;
} nbs_bank_cfg;

/* How the sub-band decomposition is carried out. */
typedef enum {
    NBS_XFORM_SPATIAL = 0,  /* truncated, windowed FIR kernels            */
    NBS_XFORM_DFT     = 1   /* exact transfer functions sampled on the
                             * DFT grid, sub-bands inverse-transformed at
                             * their own critical rate                    */
} nbs_transform;

typedef struct {
    nbs_bank_cfg bank;

    int    win;
    double win_sigma;

    /* C_ij = (kc[max(i,j)] * L)^2 -- indexed by the coarser 1-D band, so the
     * low-pass keeps a C1-like constant and the detail bands a C2-like one. */
    double kc[NBS_MAX_BANDS];
    double L;

    /* Exponents.  Packed for the actual band count; index with
     * nbs_w_index_ex(i, j, nbands, w_sym). */
    double w[NBS_MAX_W];

    /* 1 = tie w_ij to w_ji (isotropy across the cardinal axes), 0 = all
     * n^2 exponents free.  Leaving them free costs a larger search but
     * lets the data say whether the isotropy holds. */
    int    w_sym;

    double xi_floor;

    /* Normaliser of the similarity ratio: xi = (2 E[ab] + C) / (2 M_p + C)
     * with M_p the power mean of E[a^2] and E[b^2] of order p.  p = 1 is the
     * arithmetic mean, i.e. SSIM's ratio; p = 0 the geometric mean, which
     * makes xi a correlation coefficient blind to a gain on either band;
     * p -> inf the larger of the two energies.  Default 1. */
    double norm_p;
    int    downsample;

    /* Pixels discarded from each edge before averaging.  -1 selects the
     * automatic value (longest analysis filter radius + window radius), which
     * evaluates only where every filter fits entirely inside the image.  That
     * is the conservative choice, but it makes NB-SSIM average over far fewer
     * pixels than SSIM does -- 68%% of the image against SSIM's 92%% at
     * ppd 50 -- which is a confound when the two are compared.  Setting 0
     * evaluates everywhere and relies on the symmetric border extension, the
     * same compromise MS-SSIM makes in its inter-scale filtering. */
    int    margin;

    /* Sub-band decomposition method. */
    nbs_transform transform;

    /* DFT core only.  1 = inverse-transform each sub-band at its own
     * critical rate, so the 11x11 statistics window covers a fixed number
     * of cycles of the channel it is measuring rather than a fixed number
     * of pixels.  This matters: at ppd 63 an 11x11 window spans 0.17 deg,
     * while one cycle of the 0.8 cpd channel is 1.25 deg, so at full rate
     * the window sees a seventh of a cycle of the band whose similarity it
     * is supposed to be reporting.  MS-SSIM avoids this by decimating
     * between scales.  0 keeps every sub-band at full resolution, which
     * reproduces the spatial core's behaviour and is used to check the two
     * against each other. */
    int    decimate;

    /* Relative amplitude below which a band's transfer function is treated
     * as zero when choosing its critical rate. */
    double decim_eps;
} nbs_nbssim_cfg;

/* Packing of the exponent matrix.
 *
 * Sub-band (i,j) is filtered with channel i along the rows and channel j
 * along the columns, so i indexes HORIZONTAL frequency and j VERTICAL
 * frequency.  Reflecting across the main diagonal therefore exchanges the
 * two cardinal axes, and w_ij == w_ji is the statement that horizontal and
 * vertical detail count equally.  That is a strong assumption -- the
 * oblique effect says sensitivity falls away from the cardinal axes, and
 * the diagonal entries w_ii are exactly the obliquely-oriented sub-bands --
 * so it is a switch, not a law.
 *
 *   sym != 0   upper triangle only, n(n+1)/2 free exponents
 *   sym == 0   all n^2 exponents free, row-major
 */
int nbs_w_count_ex(int nbands, int sym);
int nbs_w_index_ex(int i, int j, int nbands, int sym);

/* Symmetric packing (sym = 1), kept for callers that assume it. */
int nbs_w_count(int nbands);
int nbs_w_index(int i, int j, int nbands);

void nbs_nbssim_cfg_default(nbs_nbssim_cfg *c);

/* Configure a DOG (or LOGGABOR) bank from a list of `nbands` frequencies in
 * cycles per degree.  freq[0] is the low-pass corner (unless
 * first_is_bandpass is set), freq[1..nbands-2] are band-pass peaks and
 * freq[nbands-1] is the high-pass corner.  If only nbands-1 values are
 * supplied the high-pass corner is taken as sqrt(2) times the last one.
 *
 * IMPORTANT: this also calibrates the band-pass corners against the *discrete*
 * filter response at bank->ppd, so it must be called after ppd, nbands and rho
 * are set, and re-called if any of them changes.  Mutating bank->ppd alone
 * leaves a stale calibration and mis-places the band-pass peaks.
 *
 * Returns 0 on success. */
int nbs_bank_set(nbs_bank_cfg *bank, const double *freq, int nfreq,
                 int first_is_bandpass, int last_is_bandpass);

/* CUMGAUSS helper: choose nested cut-offs so the interior DoG bands peak at
 * the requested frequencies.  peaks[0] is the low-pass corner. */
int nbs_bank_set_from_peaks(nbs_bank_cfg *bank, const double *peaks);

/* Build the N one-dimensional analysis kernels.  Caller frees each with
 * nbs_kernel_free.  ks must have room for bank->nbands entries. */
int nbs_bank_kernels(const nbs_bank_cfg *bank, nbs_kernel *ks);

/* Analytic 1-D transfer function of band b at normalised frequency nu
 * (cycles/pixel, 0 .. 0.5).  This is the *definition* of the band; the
 * spatial kernels are a truncated, windowed approximation to it, while the
 * DFT core samples it directly on the transform grid and so realises the
 * requested peak and bandwidth exactly. */
double nbs_band_H(const nbs_bank_cfg *bank, int b, double nu);

/* Peak value of |H_b| over 0 <= nu <= 0.5, and the frequency at which it
 * occurs (either may be NULL). */
double nbs_band_peak(const nbs_bank_cfg *bank, int b, double *at_nu);

/* Smallest nu such that |H_b(nu')| <= eps * peak for every nu' > nu.  Used
 * to pick each sub-band's critical rate. */
double nbs_band_edge(const nbs_bank_cfg *bank, int b, double eps);

/* Half-amplitude bandwidth of band b in octaves, measured on the analytic
 * response.  Returns 0 if the band has no two-sided half-amplitude pair
 * (low-pass and high-pass sections do not). */
double nbs_band_bw_oct(const nbs_bank_cfg *bank, int b);

/* |H(nu)| of an explicit kernel at normalised frequency nu (cycles/pixel).
 * Zero-phase symmetric kernels have a real response; the signed value is
 * returned. */
double nbs_kernel_response(const nbs_kernel *k, double nu);

/* The 11x11-equivalent unit-sum 1-D Gaussian window used by SSIM. */
int nbs_ssim_window(nbs_kernel *k, int win, double sigma);

/* Separable "valid" convolution; dst must be (h-2r) x (w-2r). */
int nbs_conv_valid_sep(const nbs_plane *src, const nbs_kernel *k, nbs_plane *dst);

/* Compute the metric. */
int nbs_nbssim(const nbs_plane *x, const nbs_plane *y,
               const nbs_nbssim_cfg *cfg, double *out);

/* Evaluate log(xi) per symmetric band pair on a subsampled grid.
 * `out` needs nbs_w_count(nbands) * (*out_npix) doubles, pixel-major.
 * Pass out == NULL to query the sample count only.
 *
 * With this cache, the metric for any exponent vector w is
 *     mean_p exp( sum_k w[k] * out[p*nw + k] ). */
double nbs_xi_ratio(double exy, double exx, double eyy, double C, double p);

int nbs_nbssim_logxi(const nbs_plane *x, const nbs_plane *y,
                     const nbs_nbssim_cfg *cfg, int stride,
                     double *out, int *out_npix);

/* DFT-domain evaluation, used by both of the above when cfg->transform is
 * NBS_XFORM_DFT.  mode 0 fills *score, mode 1 fills the log-xi cache. */
int nbs_nbssim_dft(const nbs_plane *x, const nbs_plane *y,
                   const nbs_nbssim_cfg *cfg, int mode, int stride,
                   double *score, double *cache, int *out_npix);

/* Report the sub-band grid the DFT core would use for an image of the
 * given size: transform size, and per-band decimated length.  Either
 * output array may be NULL.  Returns 0 on success. */
int nbs_nbssim_dft_grid(const nbs_nbssim_cfg *cfg, int w, int h,
                        int *out_wp, int *out_hp, int *out_wd, int *out_hd);

#ifdef __cplusplus
}
#endif

#endif /* NBS_METRICS_H */
