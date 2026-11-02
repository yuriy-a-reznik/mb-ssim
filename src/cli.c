/* cli.c -- shared option parsing.  Strict C89. */

#include "nbs.h"
#include "nbs_cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const double NBS_DEFAULT_PEAKS[5] = { 0.8, 1.7, 2.8, 4.0, 8.0 };
static double g_freq[NBS_MAX_BANDS];
static int    g_nfreq = 0;
static int    g_first_bp = 0;
static int    g_last_bp  = 0;
static int    g_have_cutoffs = 0;
static double g_wtmp[NBS_MAX_W];
static int    g_nwtmp = 0;
static char   g_wfile[1024];
static int    g_have_wfile = 0;

void nbs_opts_default(nbs_opts *o)
{
    memset(o, 0, sizeof(*o));
    o->luma = NBS_LUMA_BT601;
    nbs_ssim_cfg_default(&o->ssim);
    nbs_msssim_cfg_default(&o->msssim);
    nbs_nbssim_cfg_default(&o->nbs);
    o->ppd_set = 0;
    memcpy(g_freq, NBS_DEFAULT_PEAKS, sizeof(NBS_DEFAULT_PEAKS));
    g_nfreq = 5;
    g_first_bp = 0;
    g_last_bp = 0;
    g_have_cutoffs = 0;
    g_nwtmp = 0;
    g_have_wfile = 0;
}

static int parse_list(const char *s, double *out, int maxn)
{
    int n = 0;
    const char *p = s;
    while (*p != '\0' && n < maxn) {
        char *end;
        double v = strtod(p, &end);
        if (end == p) break;
        out[n++] = v;
        p = end;
        while (*p == ',' || *p == ' ' || *p == ';') p++;
    }
    return n;
}

static int starts(const char *a, const char *b, const char **rest)
{
    size_t n = strlen(b);
    if (strncmp(a, b, n) == 0) { *rest = a + n; return 1; }
    return 0;
}

int nbs_opts_parse_one(nbs_opts *o, const char *arg)
{
    const char *v;

    if (starts(arg, "--ppd=", &v)) {
        o->nbs.bank.ppd = atof(v);
        o->ppd_set = 1;
        if (o->nbs.bank.ppd <= 0.0) { nbs_set_error("--ppd must be positive"); return -1; }
        return 1;
    }
    if (starts(arg, "--bank=", &v)) {
        if (strcmp(v, "dog") == 0) o->nbs.bank.kind = NBS_BANK_DOG;
        else if (strcmp(v, "cumgauss") == 0) o->nbs.bank.kind = NBS_BANK_CUMGAUSS;
        else if (strcmp(v, "loggabor") == 0) o->nbs.bank.kind = NBS_BANK_LOGGABOR;
        else { nbs_set_error("--bank must be dog, cumgauss or loggabor"); return -1; }
        return 1;
    }
    if (starts(arg, "--peaks=", &v)) {
        g_nfreq = parse_list(v, g_freq, NBS_MAX_BANDS);
        if (g_nfreq < 1) { nbs_set_error("--peaks: could not parse list"); return -1; }
        g_have_cutoffs = 0;
        return 1;
    }
    if (starts(arg, "--bands=", &v)) {
        int n = atoi(v);
        if (n < 2 || n > NBS_MAX_BANDS) { nbs_set_error("--bands out of range"); return -1; }
        o->nbs.bank.nbands = n;
        return 1;
    }
    if (starts(arg, "--first=", &v)) {
        if (strcmp(v, "lowpass") == 0) g_first_bp = 0;
        else if (strcmp(v, "bandpass") == 0) g_first_bp = 1;
        else { nbs_set_error("--first must be lowpass or bandpass"); return -1; }
        return 1;
    }
    if (starts(arg, "--last=", &v)) {
        if (strcmp(v, "highpass") == 0) g_last_bp = 0;
        else if (strcmp(v, "bandpass") == 0) g_last_bp = 1;
        else { nbs_set_error("--last must be highpass or bandpass"); return -1; }
        return 1;
    }
    if (starts(arg, "--rho=", &v)) {
        double tmp[NBS_MAX_BANDS];
        int n = parse_list(v, tmp, NBS_MAX_BANDS), i;
        if (n < 1) { nbs_set_error("--rho: could not parse list"); return -1; }
        for (i = 0; i < NBS_MAX_BANDS; ++i)
            o->nbs.bank.rho[i] = (i < n) ? tmp[i] : tmp[n - 1];
        return 1;
    }
    if (starts(arg, "--band-norm=", &v)) {
        if (strcmp(v, "peak") == 0) o->nbs.bank.norm_peak = 1;
        else if (strcmp(v, "none") == 0) o->nbs.bank.norm_peak = 0;
        else { nbs_set_error("--band-norm must be peak or none"); return -1; }
        return 1;
    }
    if (starts(arg, "--bw-oct=", &v) || starts(arg, "--bw=", &v)) {
        double tmp[NBS_MAX_BANDS];
        int n = parse_list(v, tmp, NBS_MAX_BANDS), i;
        if (n < 1) { nbs_set_error("--bw: could not parse list"); return -1; }
        for (i = 0; i < NBS_MAX_BANDS; ++i)
            o->nbs.bank.bw_oct[i] = (i < n) ? tmp[i] : tmp[n - 1];
        return 1;
    }
    if (starts(arg, "--preset=", &v)) {
        /* Wilson, McFarlane & Phillips (1983) / Wilson & Gelb (1984): six
         * band-pass mechanisms about an octave apart, bandwidths narrowing
         * with frequency. */
        if (strcmp(v, "wilson") == 0) {
            static const double wf[6] = { 0.8, 1.7, 2.8, 4.0, 8.0, 16.0 };
            static const double wb[6] = { 2.5, 2.2, 1.9, 1.7, 1.5, 1.25 };
            int i;
            memcpy(g_freq, wf, sizeof(wf));
            g_nfreq = 6;
            g_first_bp = 1; g_last_bp = 1;
            g_have_cutoffs = 0;
            o->nbs.bank.nbands = 6;
            /* log-Gabor, not DoG: a difference of Gaussians bottoms out at
             * ~1.77 octaves and cannot reach Wilson's narrow high bands. */
            o->nbs.bank.kind = NBS_BANK_LOGGABOR;
            o->nbs.bank.norm_peak = 1;
            for (i = 0; i < 6; ++i) o->nbs.bank.bw_oct[i] = wb[i];
            for (i = 6; i < NBS_MAX_BANDS; ++i) o->nbs.bank.bw_oct[i] = 0.0;
        } else {
            nbs_set_error("--preset: only 'wilson' is defined");
            return -1;
        }
        return 1;
    }
    if (starts(arg, "--cutoffs=", &v)) {
        int n = parse_list(v, o->nbs.bank.fc, NBS_MAX_BANDS);
        if (n < 1) { nbs_set_error("--cutoffs: could not parse list"); return -1; }
        o->nbs.bank.nbands = n + 1;
        o->nbs.bank.kind = NBS_BANK_CUMGAUSS;
        g_have_cutoffs = 1;
        g_nfreq = 0;
        return 1;
    }
    if (starts(arg, "--weights=", &v)) {
        g_nwtmp = parse_list(v, g_wtmp, NBS_MAX_W);
        if (g_nwtmp < 1) { nbs_set_error("--weights: could not parse list"); return -1; }
        return 1;
    }
    if (starts(arg, "--weight-file=", &v)) {
        strncpy(g_wfile, v, sizeof(g_wfile) - 1);
        g_wfile[sizeof(g_wfile) - 1] = '\0';
        g_have_wfile = 1;
        return 1;
    }
    if (starts(arg, "--kc=", &v)) {
        double tmp[NBS_MAX_BANDS];
        int n = parse_list(v, tmp, NBS_MAX_BANDS), i;
        if (n < 1) { nbs_set_error("--kc: could not parse list"); return -1; }
        for (i = 0; i < NBS_MAX_BANDS; ++i)
            o->nbs.kc[i] = (i < n) ? tmp[i] : tmp[n - 1];
        return 1;
    }
    if (starts(arg, "--win=", &v)) {
        int n = atoi(v);
        if (n < 1 || (n % 2) == 0) { nbs_set_error("--win must be odd and >= 1"); return -1; }
        o->ssim.win = n; o->msssim.base.win = n; o->nbs.win = n;
        return 1;
    }
    if (starts(arg, "--nbs-win=", &v)) {
        int n = atoi(v);
        if (n < 1 || (n % 2) == 0) { nbs_set_error("--nbs-win must be odd and >= 1"); return -1; }
        o->nbs.win = n;
        return 1;
    }
    if (starts(arg, "--win-sigma=", &v)) {
        double s = atof(v);
        if (s <= 0.0) { nbs_set_error("--win-sigma must be positive"); return -1; }
        o->ssim.win_sigma = s; o->msssim.base.win_sigma = s; o->nbs.win_sigma = s;
        return 1;
    }
    if (starts(arg, "--radius-mult=", &v)) {
        o->nbs.bank.radius_mult = atof(v);
        return 1;
    }
    if (starts(arg, "--luma=", &v)) {
        if (strcmp(v, "601") == 0) o->luma = NBS_LUMA_BT601;
        else if (strcmp(v, "709") == 0) o->luma = NBS_LUMA_BT709;
        else if (strcmp(v, "mean") == 0) o->luma = NBS_LUMA_MEAN;
        else { nbs_set_error("--luma must be 601, 709 or mean"); return -1; }
        return 1;
    }
    if (starts(arg, "--ssim-downsample=", &v)) {
        if (strcmp(v, "auto") == 0) o->ssim.downsample = -1;
        else if (strcmp(v, "off") == 0) o->ssim.downsample = 0;
        else o->ssim.downsample = atoi(v);
        return 1;
    }
    if (starts(arg, "--nbs-downsample=", &v)) {
        if (strcmp(v, "auto") == 0) o->nbs.downsample = -1;
        else if (strcmp(v, "off") == 0) o->nbs.downsample = 0;
        else o->nbs.downsample = atoi(v);
        return 1;
    }
    if (starts(arg, "--levels=", &v)) {
        int n = atoi(v);
        if (n < 1 || n > NBS_MSSSIM_MAX_LEVELS) { nbs_set_error("--levels out of range"); return -1; }
        o->msssim.levels = n;
        return 1;
    }
    if (starts(arg, "--nbs-margin=", &v)) {
        if (strcmp(v, "auto") == 0) o->nbs.margin = -1;
        else {
            o->nbs.margin = atoi(v);
            if (o->nbs.margin < 0) { nbs_set_error("--nbs-margin must be >= 0 or auto"); return -1; }
        }
        return 1;
    }
    if (starts(arg, "--xi-floor=", &v)) {
        o->nbs.xi_floor = atof(v);
        return 1;
    }
    if (starts(arg, "--xi-norm=", &v)) {
        if (strcmp(v, "geo") == 0) o->nbs.norm_p = 0.0;
        else if (strcmp(v, "max") == 0) o->nbs.norm_p = 1e9;
        else o->nbs.norm_p = atof(v);
        return 1;
    }
    if (starts(arg, "--transform=", &v)) {
        if (strcmp(v, "spatial") == 0) o->nbs.transform = NBS_XFORM_SPATIAL;
        else if (strcmp(v, "dft") == 0) o->nbs.transform = NBS_XFORM_DFT;
        else { nbs_set_error("--transform must be spatial or dft"); return -1; }
        return 1;
    }
    if (starts(arg, "--band-decimate=", &v)) {
        if (strcmp(v, "auto") == 0 || strcmp(v, "on") == 0) o->nbs.decimate = 1;
        else if (strcmp(v, "off") == 0) o->nbs.decimate = 0;
        else { nbs_set_error("--band-decimate must be auto or off"); return -1; }
        return 1;
    }
    if (starts(arg, "--decim-eps=", &v)) {
        double e = atof(v);
        if (!(e > 0.0) || e >= 1.0) { nbs_set_error("--decim-eps must be in (0,1)"); return -1; }
        o->nbs.decim_eps = e;
        return 1;
    }
    if (starts(arg, "--weight-sym=", &v)) {
        if (strcmp(v, "on") == 0)  o->nbs.w_sym = 1;
        else if (strcmp(v, "off") == 0) o->nbs.w_sym = 0;
        else { nbs_set_error("--weight-sym must be on or off"); return -1; }
        return 1;
    }
    if (starts(arg, "--taper=", &v)) {
        if (strcmp(v, "off") == 0) o->nbs.bank.taper_mult = 0.0;
        else {
            double t = atof(v);
            if (!(t > 0.0)) { nbs_set_error("--taper must be positive or off"); return -1; }
            o->nbs.bank.taper_mult = t;
        }
        return 1;
    }
    return 0;
}

int nbs_opts_finalise(nbs_opts *o)
{
    int n, nw, i, j, k;

    /* The band-pass corners are calibrated against the discrete response at a
     * specific ppd, so the bank must always be rebuilt here -- after --ppd,
     * --rho and --bands have all been seen.  Rebuilding only when --peaks was
     * given would silently leave a stale calibration behind. */
    if (g_have_cutoffs) {
        if (o->nbs.bank.kind != NBS_BANK_CUMGAUSS) {
            nbs_set_error("--cutoffs only applies to the cumgauss bank");
            return 1;
        }
    } else {
        if (g_nfreq != o->nbs.bank.nbands && g_nfreq != o->nbs.bank.nbands - 1) {
            if (g_nfreq == 5 && o->nbs.bank.nbands != 6) {
                nbs_set_error("--bands requires --peaks (or --cutoffs for cumgauss)");
                return 1;
            }
            o->nbs.bank.nbands = g_nfreq + 1;
        }
        if (o->nbs.bank.kind == NBS_BANK_CUMGAUSS) {
            if (nbs_bank_set_from_peaks(&o->nbs.bank, g_freq) != 0) return 1;
        }
        if (nbs_bank_set(&o->nbs.bank, g_freq, g_nfreq, g_first_bp, g_last_bp) != 0)
            return 1;
    }

    n  = o->nbs.bank.nbands;

    /* Re-seed the exponents for the (possibly changed) band count. */
    for (i = 0; i < NBS_MAX_W; ++i) o->nbs.w[i] = 2.0 / (double)(n * n);
    nw = nbs_w_count_ex(n, o->nbs.w_sym);

    if (g_have_wfile) {
        if (nbs_load_weights_ex(g_wfile, o->nbs.w, n, o->nbs.w_sym) != 0) return 1;
    }
    if (g_nwtmp > 0) {
        if (g_nwtmp != nw) { nbs_set_error("--weights: wrong number of values"); return 1; }
        k = 0;
        for (i = 0; i < n; ++i)
            for (j = (o->nbs.w_sym ? i : 0); j < n; ++j)
                o->nbs.w[nbs_w_index_ex(i, j, n, o->nbs.w_sym)] = g_wtmp[k++];
    }
    return 0;
}

/* Read an exponent file.
 *
 * The file records the band count it was fitted at and whether the matrix
 * was constrained to be symmetric, and neither need match the bank in use.
 * A bank configured at a shorter viewing distance drops every channel above
 * the display Nyquist rate, so a vector fitted with six channels has to be
 * readable by a five-channel bank; the surviving channels are the same
 * channels at the same cycles per degree, so the mapping is by (i,j) and
 * not by position in the file.  A symmetric file can also fill an
 * asymmetric bank, by mirroring.  Files written before these fields existed
 * are assumed to match. */
int nbs_load_weights_ex(const char *path, double *w, int nbands, int sym)
{
    FILE *f = fopen(path, "r");
    int i, j, k = 0, file_nb = 0, file_sym = -1, file_nw;
    double tmp[NBS_MAX_W];
    char line[4096];

    if (f == NULL) { nbs_set_error("cannot open weight file"); return 1; }

    /* Pass 1: read the header. */
    while (fgets(line, (int)sizeof(line), f) != NULL) {
        char *p = line, *q;
        while (*p == ' ' || *p == '\t') p++;
        if (*p != '#') break;
        q = strstr(p, "bands");
        if (q != NULL && q > p) {
            char *r = q - 1;
            while (r > p && (*r == ' ' || *r == '\t')) r--;
            if (*r >= '0' && *r <= '9') {
                while (r > p && *(r - 1) >= '0' && *(r - 1) <= '9') r--;
                file_nb = atoi(r);
            }
        }
        if (strstr(p, "symmetric") != NULL) file_sym = 1;
        if (strstr(p, "asymmetric") != NULL) file_sym = 0;
    }
    if (file_nb <= 0 || file_nb > NBS_MAX_BANDS) file_nb = nbands;
    if (file_sym < 0) file_sym = sym;
    file_nw = nbs_w_count_ex(file_nb, file_sym);

    rewind(f);
    while (k < file_nw && fgets(line, (int)sizeof(line), f) != NULL) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0') continue;
        while (k < file_nw) {
            char *end;
            double v = strtod(p, &end);
            if (end == p) break;
            tmp[k++] = v;
            p = end;
        }
    }
    fclose(f);
    if (k != file_nw) { nbs_set_error("weight file: wrong number of values"); return 1; }
    if (file_nb < nbands) {
        nbs_set_error("weight file has fewer bands than the configured bank");
        return 1;
    }
    if (!file_sym && sym) {
        nbs_set_error("cannot read an asymmetric exponent file into a "
                      "symmetric bank -- w_ij and w_ji disagree");
        return 1;
    }
    for (i = 0; i < nbands; ++i)
        for (j = 0; j < nbands; ++j)
            w[nbs_w_index_ex(i, j, nbands, sym)] =
                tmp[nbs_w_index_ex(i, j, file_nb, file_sym)];
    return 0;
}

int nbs_load_weights(const char *path, double *w, int nbands)
{
    return nbs_load_weights_ex(path, w, nbands, 1);
}

int nbs_save_weights(const char *path, const nbs_nbssim_cfg *cfg)
{
    FILE *f = fopen(path, "w");
    int i, j, n = cfg->bank.nbands;
    if (f == NULL) { nbs_set_error("cannot write weight file"); return 1; }
    if (cfg->w_sym)
        fprintf(f, "# nbssim exponents, %d bands, symmetric, "
                   "upper triangle (w_ij == w_ji)\n", n);
    else
        fprintf(f, "# nbssim exponents, %d bands, asymmetric, full matrix "
                   "(row i = horizontal channel, column j = vertical)\n", n);
    for (i = 0; i < n; ++i) {
        for (j = (cfg->w_sym ? i : 0); j < n; ++j)
            fprintf(f, "%.10g%s", cfg->w[nbs_w_index_ex(i, j, n, cfg->w_sym)],
                    (j + 1 < n) ? " " : "");
        fprintf(f, "\n");
    }
    fclose(f);
    return 0;
}

void nbs_opts_print_help(void)
{
    printf("Common options:\n");
    printf("  --ppd=F              pixels per degree of visual angle (default 30)\n");
    printf("  --bank=KIND          dog (default) | cumgauss | loggabor\n");
    printf("  --peaks=a,b,c,...    band frequencies in cycles/degree.  The first\n");
    printf("                       is the low-pass corner, the middle ones are\n");
    printf("                       band-pass peaks, the last is the high-pass\n");
    printf("                       corner.  Give N-1 values and the high-pass\n");
    printf("                       corner is placed sqrt(2) above the last peak.\n");
    printf("                       Default: 0.8,1.7,2.8,4,8\n");
    printf("  --bands=N            number of 1-D bands (2..%d)\n", NBS_MAX_BANDS);
    printf("  --first=lowpass|bandpass    role of band 0 (default lowpass)\n");
    printf("  --last=highpass|bandpass    role of band N-1 (default highpass)\n");
    printf("  --rho=a,b,...        DoG sharpness per band, > 1 (default 3)\n");
    printf("  --bw-oct=a,b,...     half-height bandwidth per band, in octaves.\n");
    printf("                       For the DoG bank this overrides --rho\n");
    printf("  --band-norm=peak|none  scale each band to unit peak gain\n");
    printf("  --preset=wilson      six band-pass channels at 0.8/1.7/2.8/4/8/16\n");
    printf("                       cpd with 2.5..1.25 octave bandwidths, after\n");
    printf("                       Wilson, McFarlane & Phillips (1983)\n");
    printf("  --cutoffs=a,b,...    cumgauss nested cut-offs, cpd\n");
    printf("  --weights=...        N(N+1)/2 symmetric exponents, upper triangle\n");
    printf("  --weight-file=PATH   read the exponents from a text file\n");
    printf("  --kc=a,b,...         stabilising constants k, C=(k*L)^2, indexed by\n");
    printf("                       the coarser 1-D band (default 0.01, then 0.03)\n");
    printf("  --win=N              statistics window size, odd (default 11).\n");
    printf("                       1 means no window: xi is formed from the\n");
    printf("                       sub-band samples themselves, which is the\n");
    printf("                       natural choice once the filterbank has\n");
    printf("                       already band-limited them\n");
    printf("  --nbs-win=N          as --win, but for MB-SSIM only, leaving\n");
    printf("                       SSIM and MS-SSIM at their published 11x11\n");
    printf("  --win-sigma=F        window Gaussian sigma (default 1.5)\n");
    printf("  --radius-mult=F      analysis filter support in sigmas (default 4)\n");
    printf("  --luma=601|709|mean  luma extraction (default 601, MATLAB rgb2gray)\n");
    printf("  --ssim-downsample=auto|off|N   SSIM pre-downsampling (default auto)\n");
    printf("  --nbs-downsample=auto|off|N    NB-SSIM pre-downsampling (default off)\n");
    printf("  --levels=N           MS-SSIM scales (default 5)\n");
    printf("  --nbs-margin=N|auto  pixels NB-SSIM discards at each edge.\n");
    printf("                       auto (default) evaluates only where every\n");
    printf("                       filter fits; 0 evaluates everywhere and so\n");
    printf("                       uses the same support as SSIM\n");
    printf("  --xi-floor=F         clamp on xi before exponentiation (1e-6)\n");
    printf("  --xi-norm=P          power mean of order P normalises the ratio:\n"
           "                       1 arithmetic (SSIM, default), geo|0 geometric\n"
           "                       (a correlation coefficient), max the larger energy\n");
    printf("  --transform=spatial|dft   how the sub-bands are formed.\n");
    printf("                       spatial (default) convolves with truncated\n");
    printf("                       FIR kernels; dft samples the transfer\n");
    printf("                       functions exactly on the DFT grid and needs\n");
    printf("                       no margin\n");
    printf("  --band-decimate=auto|off  dft only: carry each sub-band at its\n");
    printf("                       own critical rate (default auto), so the\n");
    printf("                       statistics window spans a fixed number of\n");
    printf("                       cycles of each channel rather than a fixed\n");
    printf("                       number of pixels.  off keeps every band at\n");
    printf("                       full resolution\n");
    printf("  --decim-eps=F        relative amplitude treated as zero when\n");
    printf("                       choosing a band's critical rate (0.02)\n");
    printf("  --weight-sym=on|off  tie w_ij to w_ji (default on).  i indexes\n");
    printf("                       horizontal and j vertical frequency, so\n");
    printf("                       the tie asserts that the two cardinal\n");
    printf("                       axes count equally and makes the diagonal\n");
    printf("                       entries -- the oblique sub-bands -- the\n");
    printf("                       only ones without a partner.  off frees\n");
    printf("                       all N^2 exponents\n");
    printf("  --taper=F|off        give log-Gabor bands compact support by\n");
    printf("                       tapering to zero at F*bandwidth octaves\n");
    printf("                       from the centre (default off)\n");
}
