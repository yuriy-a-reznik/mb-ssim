/* image.c -- 8-bit image loading (PNM, BMP) and luma extraction.
 * PNG lives in png.c.  Strict C89.
 */

#include "nbs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* from png.c */
int nbs_png_load(FILE *f, nbs_image8 *out);

void nbs_image8_free(nbs_image8 *im)
{
    if (im == NULL) return;
    free(im->px);
    im->px = NULL;
    im->w = im->h = im->nch = 0;
}

/* ------------------------------------------------------------------ */
/* PNM                                                                 */
/* ------------------------------------------------------------------ */

static int pnm_getc_skip(FILE *f)
{
    int c;
    for (;;) {
        c = fgetc(f);
        if (c == '#') {
            while (c != EOF && c != '\n') c = fgetc(f);
            continue;
        }
        return c;
    }
}

static int pnm_read_int(FILE *f, int *out)
{
    int c, v = 0, got = 0;
    do {
        c = pnm_getc_skip(f);
        if (c == EOF) return 1;
    } while (c == ' ' || c == '\t' || c == '\n' || c == '\r');
    while (c >= '0' && c <= '9') {
        v = v * 10 + (c - '0');
        got = 1;
        c = pnm_getc_skip(f);
        if (c == EOF) break;
    }
    if (!got) return 1;
    *out = v;
    return 0;
}

static int pnm_load(FILE *f, int type, nbs_image8 *out)
{
    int w, h, maxv, nch;
    size_t n;

    if (pnm_read_int(f, &w) || pnm_read_int(f, &h) || pnm_read_int(f, &maxv)) {
        nbs_set_error("PNM: bad header");
        return 1;
    }
    if (w <= 0 || h <= 0 || maxv <= 0 || maxv > 255) {
        nbs_set_error("PNM: unsupported dimensions or maxval (only 8-bit supported)");
        return 1;
    }
    nch = (type == 6 || type == 3) ? 3 : 1;
    n = (size_t)w * (size_t)h * (size_t)nch;
    out->px = (unsigned char *)nbs_malloc(n);
    if (out->px == NULL) return 1;
    out->w = w; out->h = h; out->nch = nch;

    if (type == 5 || type == 6) {
        if (fread(out->px, 1, n, f) != n) {
            nbs_set_error("PNM: truncated raster");
            nbs_image8_free(out);
            return 1;
        }
    } else {
        size_t i;
        for (i = 0; i < n; ++i) {
            int v;
            if (pnm_read_int(f, &v)) {
                nbs_set_error("PNM: truncated ASCII raster");
                nbs_image8_free(out);
                return 1;
            }
            out->px[i] = (unsigned char)(v < 0 ? 0 : (v > 255 ? 255 : v));
        }
    }
    if (maxv != 255) {
        size_t i;
        for (i = 0; i < n; ++i)
            out->px[i] = (unsigned char)((out->px[i] * 255 + maxv / 2) / maxv);
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* BMP                                                                 */
/* ------------------------------------------------------------------ */

static unsigned long rd_u32(const unsigned char *p)
{
    return (unsigned long)p[0] | ((unsigned long)p[1] << 8) |
           ((unsigned long)p[2] << 16) | ((unsigned long)p[3] << 24);
}

static unsigned int rd_u16(const unsigned char *p)
{
    return (unsigned int)p[0] | ((unsigned int)p[1] << 8);
}

static long rd_s32(const unsigned char *p)
{
    unsigned long v = rd_u32(p);
    if (v & 0x80000000UL) return -(long)(0xFFFFFFFFUL - v + 1UL);
    return (long)v;
}

static int bmp_load(FILE *f, nbs_image8 *out)
{
    unsigned char hdr[54];
    unsigned char pal[1024];
    unsigned long off, dibsize, compression, ncolors;
    long bw, bh;
    unsigned int bpp;
    int w, h, flip, y, x, rowbytes;
    unsigned char *row;

    rewind(f);
    if (fread(hdr, 1, 54, f) != 54) { nbs_set_error("BMP: short header"); return 1; }
    if (hdr[0] != 'B' || hdr[1] != 'M') { nbs_set_error("BMP: bad magic"); return 1; }

    off        = rd_u32(hdr + 10);
    dibsize    = rd_u32(hdr + 14);
    if (dibsize < 40) { nbs_set_error("BMP: unsupported DIB header (< 40 bytes)"); return 1; }
    bw         = rd_s32(hdr + 18);
    bh         = rd_s32(hdr + 22);
    bpp        = rd_u16(hdr + 28);
    compression= rd_u32(hdr + 30);
    ncolors    = rd_u32(hdr + 46);

    if (compression != 0) { nbs_set_error("BMP: only uncompressed (BI_RGB) supported"); return 1; }
    if (bw <= 0) { nbs_set_error("BMP: bad width"); return 1; }
    flip = 1;
    if (bh < 0) { bh = -bh; flip = 0; }
    if (bh == 0) { nbs_set_error("BMP: bad height"); return 1; }
    w = (int)bw; h = (int)bh;

    if (bpp == 8) {
        unsigned long i, need;
        if (ncolors == 0) ncolors = 256;
        if (ncolors > 256) ncolors = 256;
        if (fseek(f, (long)(14 + dibsize), SEEK_SET) != 0) { nbs_set_error("BMP: seek failed"); return 1; }
        need = ncolors * 4;
        memset(pal, 0, sizeof(pal));
        if (fread(pal, 1, (size_t)need, f) != (size_t)need) { nbs_set_error("BMP: short palette"); return 1; }
        (void)i;
    } else if (bpp != 24 && bpp != 32) {
        nbs_set_error("BMP: only 8/24/32 bpp supported");
        return 1;
    }

    rowbytes = ((w * (int)bpp + 31) / 32) * 4;
    row = (unsigned char *)nbs_malloc((size_t)rowbytes);
    if (row == NULL) return 1;

    out->px = (unsigned char *)nbs_malloc((size_t)w * (size_t)h * 3u);
    if (out->px == NULL) { free(row); return 1; }
    out->w = w; out->h = h; out->nch = 3;

    if (fseek(f, (long)off, SEEK_SET) != 0) {
        nbs_set_error("BMP: seek to pixel data failed");
        free(row); nbs_image8_free(out); return 1;
    }

    for (y = 0; y < h; ++y) {
        int dy = flip ? (h - 1 - y) : y;
        unsigned char *dst = out->px + (size_t)dy * (size_t)w * 3u;
        if (fread(row, 1, (size_t)rowbytes, f) != (size_t)rowbytes) {
            nbs_set_error("BMP: truncated pixel data");
            free(row); nbs_image8_free(out); return 1;
        }
        if (bpp == 24) {
            for (x = 0; x < w; ++x) {
                dst[3 * x + 0] = row[3 * x + 2];
                dst[3 * x + 1] = row[3 * x + 1];
                dst[3 * x + 2] = row[3 * x + 0];
            }
        } else if (bpp == 32) {
            for (x = 0; x < w; ++x) {
                dst[3 * x + 0] = row[4 * x + 2];
                dst[3 * x + 1] = row[4 * x + 1];
                dst[3 * x + 2] = row[4 * x + 0];
            }
        } else {
            for (x = 0; x < w; ++x) {
                unsigned int idx = row[x];
                dst[3 * x + 0] = pal[4 * idx + 2];
                dst[3 * x + 1] = pal[4 * idx + 1];
                dst[3 * x + 2] = pal[4 * idx + 0];
            }
        }
    }
    free(row);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Dispatch                                                            */
/* ------------------------------------------------------------------ */

int nbs_image8_load(const char *path, nbs_image8 *out)
{
    FILE *f;
    unsigned char magic[8];
    int rc;

    if (out == NULL) return 1;
    out->px = NULL; out->w = out->h = out->nch = 0;

    f = fopen(path, "rb");
    if (f == NULL) {
        char buf[400];
        strcpy(buf, "cannot open file: ");
        strncat(buf, path, sizeof(buf) - strlen(buf) - 1);
        nbs_set_error(buf);
        return 1;
    }
    if (fread(magic, 1, 8, f) != 8) {
        nbs_set_error("file too short to identify");
        fclose(f);
        return 1;
    }

    if (magic[0] == 0x89 && magic[1] == 'P' && magic[2] == 'N' && magic[3] == 'G') {
        rewind(f);
        rc = nbs_png_load(f, out);
    } else if (magic[0] == 'B' && magic[1] == 'M') {
        rc = bmp_load(f, out);
    } else if (magic[0] == 'P' && magic[1] >= '2' && magic[1] <= '6') {
        int type = magic[1] - '0';
        if (type == 4) { nbs_set_error("PBM (P4) not supported"); fclose(f); return 1; }
        if (fseek(f, 2, SEEK_SET) != 0) { nbs_set_error("PNM: seek failed"); fclose(f); return 1; }
        rc = pnm_load(f, type, out);
    } else {
        nbs_set_error("unrecognised image format (expected PNG, BMP or PNM)");
        rc = 1;
    }
    fclose(f);
    return rc;
}

int nbs_image8_to_luma(const nbs_image8 *im, nbs_luma_mode mode, nbs_plane *out)
{
    int x, y;
    double cr, cg, cb;

    if (im == NULL || im->px == NULL) { nbs_set_error("luma: null image"); return 1; }
    if (nbs_plane_alloc(out, im->w, im->h) != 0) return 1;

    switch (mode) {
    case NBS_LUMA_BT709: cr = 0.2126; cg = 0.7152; cb = 0.0722; break;
    case NBS_LUMA_MEAN:  cr = cg = cb = 1.0 / 3.0;              break;
    case NBS_LUMA_BT601:
    default:             cr = 0.299;  cg = 0.587;  cb = 0.114;  break;
    }

    if (im->nch == 1) {
        for (y = 0; y < im->h; ++y)
            for (x = 0; x < im->w; ++x)
                NBS_AT(out, x, y) = (double)im->px[(size_t)y * (size_t)im->w + (size_t)x];
    } else {
        for (y = 0; y < im->h; ++y) {
            const unsigned char *s = im->px + (size_t)y * (size_t)im->w * 3u;
            for (x = 0; x < im->w; ++x)
                NBS_AT(out, x, y) = cr * (double)s[3 * x + 0]
                                  + cg * (double)s[3 * x + 1]
                                  + cb * (double)s[3 * x + 2];
        }
    }
    return 0;
}

int nbs_load_luma(const char *path, nbs_luma_mode mode, nbs_plane *out)
{
    nbs_image8 im;
    int rc;
    if (nbs_image8_load(path, &im) != 0) return 1;
    rc = nbs_image8_to_luma(&im, mode, out);
    nbs_image8_free(&im);
    return rc;
}
