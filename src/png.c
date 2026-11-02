/* png.c -- self-contained PNG reader with a built-in DEFLATE decoder.
 *
 * Supports: bit depths 1/2/4/8/16, colour types 0 (grey), 2 (truecolour),
 * 3 (palette), 4 (grey+alpha), 6 (truecolour+alpha), non-interlaced.
 * 16-bit samples are truncated to their high byte.  Alpha is discarded
 * (the datasets we target are opaque).  Adam7 interlacing is rejected.
 *
 * Strict C89.  No zlib dependency.
 */

#include "nbs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================== */
/* DEFLATE                                                             */
/* ================================================================== */

typedef struct {
    const unsigned char *src;
    size_t               len;
    size_t               pos;
    unsigned long        bitbuf;
    int                  bitcnt;
    int                  err;
} bitreader;

static void br_init(bitreader *b, const unsigned char *src, size_t len)
{
    b->src = src; b->len = len; b->pos = 0;
    b->bitbuf = 0; b->bitcnt = 0; b->err = 0;
}

static int br_bit(bitreader *b)
{
    int v;
    if (b->bitcnt == 0) {
        if (b->pos >= b->len) { b->err = 1; return 0; }
        b->bitbuf = b->src[b->pos++];
        b->bitcnt = 8;
    }
    v = (int)(b->bitbuf & 1UL);
    b->bitbuf >>= 1;
    b->bitcnt--;
    return v;
}

static int br_bits(bitreader *b, int n)
{
    int v = 0, i;
    for (i = 0; i < n; ++i) v |= br_bit(b) << i;
    return v;
}

/* Canonical Huffman table. */
#define HUF_MAXBITS 15
typedef struct {
    int count[HUF_MAXBITS + 1];
    int symbol[288 + 32];
    int nsym;
} huftab;

static int huf_build(huftab *h, const unsigned char *lengths, int n)
{
    int i, len;
    int offs[HUF_MAXBITS + 1];

    for (i = 0; i <= HUF_MAXBITS; ++i) h->count[i] = 0;
    for (i = 0; i < n; ++i) h->count[lengths[i]]++;
    h->count[0] = 0;

    offs[0] = 0; offs[1] = 0;
    for (len = 1; len < HUF_MAXBITS; ++len)
        offs[len + 1] = offs[len] + h->count[len];
    for (i = 0; i < n; ++i)
        if (lengths[i] != 0) h->symbol[offs[lengths[i]]++] = i;
    h->nsym = n;
    return 0;
}

static int huf_decode(bitreader *b, const huftab *h)
{
    int code = 0, first = 0, index = 0, len;
    for (len = 1; len <= HUF_MAXBITS; ++len) {
        code |= br_bit(b);
        if (b->err) return -1;
        if (code - first < h->count[len])
            return h->symbol[index + (code - first)];
        index += h->count[len];
        first += h->count[len];
        first <<= 1;
        code <<= 1;
    }
    return -1;
}

/* Output buffer that grows. */
typedef struct {
    unsigned char *d;
    size_t         n;
    size_t         cap;
    int            err;
} obuf;

static void ob_push(obuf *o, unsigned char c)
{
    if (o->err) return;
    if (o->n == o->cap) {
        size_t nc = o->cap ? o->cap * 2 : 65536;
        unsigned char *p = (unsigned char *)realloc(o->d, nc);
        if (p == NULL) { o->err = 1; return; }
        o->d = p; o->cap = nc;
    }
    o->d[o->n++] = c;
}

static const int len_base[29] = {
    3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258
};
static const int len_extra[29] = {
    0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0
};
static const int dist_base[30] = {
    1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,
    1025,1537,2049,3073,4097,6145,8193,12289,16385,24577
};
static const int dist_extra[30] = {
    0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13
};

static int inflate_block(bitreader *b, obuf *o, const huftab *lt, const huftab *dt)
{
    for (;;) {
        int sym = huf_decode(b, lt);
        if (sym < 0) return 1;
        if (sym < 256) {
            ob_push(o, (unsigned char)sym);
            if (o->err) return 1;
        } else if (sym == 256) {
            return 0;
        } else {
            int i, length, dsym, dist;
            sym -= 257;
            if (sym >= 29) return 1;
            length = len_base[sym] + br_bits(b, len_extra[sym]);
            dsym = huf_decode(b, dt);
            if (dsym < 0 || dsym >= 30) return 1;
            dist = dist_base[dsym] + br_bits(b, dist_extra[dsym]);
            if (b->err) return 1;
            if ((size_t)dist > o->n) return 1;
            for (i = 0; i < length; ++i) {
                ob_push(o, o->d[o->n - (size_t)dist]);
                if (o->err) return 1;
            }
        }
    }
}

static void build_fixed(huftab *lt, huftab *dt)
{
    unsigned char l[288], d[30];
    int i;
    for (i = 0;   i < 144; ++i) l[i] = 8;
    for (i = 144; i < 256; ++i) l[i] = 9;
    for (i = 256; i < 280; ++i) l[i] = 7;
    for (i = 280; i < 288; ++i) l[i] = 8;
    for (i = 0;   i < 30;  ++i) d[i] = 5;
    huf_build(lt, l, 288);
    huf_build(dt, d, 30);
}

static int build_dynamic(bitreader *b, huftab *lt, huftab *dt)
{
    static const int ord[19] = {16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15};
    unsigned char clen[19];
    unsigned char lengths[288 + 32];
    huftab ct;
    int hlit, hdist, hclen, i, n;

    hlit  = br_bits(b, 5) + 257;
    hdist = br_bits(b, 5) + 1;
    hclen = br_bits(b, 4) + 4;
    if (b->err || hlit > 288 || hdist > 32) return 1;

    for (i = 0; i < 19; ++i) clen[i] = 0;
    for (i = 0; i < hclen; ++i) clen[ord[i]] = (unsigned char)br_bits(b, 3);
    if (b->err) return 1;
    huf_build(&ct, clen, 19);

    n = 0;
    while (n < hlit + hdist) {
        int sym = huf_decode(b, &ct);
        if (sym < 0) return 1;
        if (sym < 16) {
            lengths[n++] = (unsigned char)sym;
        } else if (sym == 16) {
            int rep, prev;
            if (n == 0) return 1;
            prev = lengths[n - 1];
            rep = 3 + br_bits(b, 2);
            while (rep-- > 0 && n < hlit + hdist) lengths[n++] = (unsigned char)prev;
        } else if (sym == 17) {
            int rep = 3 + br_bits(b, 3);
            while (rep-- > 0 && n < hlit + hdist) lengths[n++] = 0;
        } else {
            int rep = 11 + br_bits(b, 7);
            while (rep-- > 0 && n < hlit + hdist) lengths[n++] = 0;
        }
        if (b->err) return 1;
    }
    huf_build(lt, lengths, hlit);
    huf_build(dt, lengths + hlit, hdist);
    return 0;
}

/* zlib stream -> raw bytes.  Caller frees *out. */
static int zlib_inflate(const unsigned char *src, size_t len,
                        unsigned char **out, size_t *outlen)
{
    bitreader b;
    obuf o;
    int final, type;

    if (len < 2) { nbs_set_error("PNG: zlib stream too short"); return 1; }
    if ((src[0] & 0x0f) != 8) { nbs_set_error("PNG: unsupported zlib method"); return 1; }
    if (src[1] & 0x20) { nbs_set_error("PNG: preset dictionary not supported"); return 1; }

    o.d = NULL; o.n = 0; o.cap = 0; o.err = 0;
    br_init(&b, src + 2, len - 2);

    do {
        final = br_bit(&b);
        type  = br_bits(&b, 2);
        if (b.err) { free(o.d); nbs_set_error("PNG: truncated deflate stream"); return 1; }
        if (type == 0) {
            int lenv, i;
            b.bitcnt = 0; b.bitbuf = 0;
            if (b.pos + 4 > b.len) { free(o.d); nbs_set_error("PNG: bad stored block"); return 1; }
            lenv = (int)b.src[b.pos] | ((int)b.src[b.pos + 1] << 8);
            b.pos += 4;
            if (b.pos + (size_t)lenv > b.len) { free(o.d); nbs_set_error("PNG: bad stored block length"); return 1; }
            for (i = 0; i < lenv; ++i) ob_push(&o, b.src[b.pos++]);
            if (o.err) { free(o.d); nbs_set_error("PNG: out of memory"); return 1; }
        } else if (type == 1 || type == 2) {
            huftab lt, dt;
            if (type == 1) build_fixed(&lt, &dt);
            else if (build_dynamic(&b, &lt, &dt) != 0) {
                free(o.d); nbs_set_error("PNG: bad dynamic huffman block"); return 1;
            }
            if (inflate_block(&b, &o, &lt, &dt) != 0) {
                free(o.d); nbs_set_error("PNG: deflate decode error"); return 1;
            }
        } else {
            free(o.d); nbs_set_error("PNG: reserved deflate block type"); return 1;
        }
    } while (!final);

    *out = o.d;
    *outlen = o.n;
    return 0;
}

/* ================================================================== */
/* PNG                                                                 */
/* ================================================================== */

static unsigned long be32(const unsigned char *p)
{
    return ((unsigned long)p[0] << 24) | ((unsigned long)p[1] << 16) |
           ((unsigned long)p[2] << 8)  |  (unsigned long)p[3];
}

static int paeth(int a, int b, int c)
{
    int p  = a + b - c;
    int pa = p > a ? p - a : a - p;
    int pb = p > b ? p - b : b - p;
    int pc = p > c ? p - c : c - p;
    if (pa <= pb && pa <= pc) return a;
    if (pb <= pc) return b;
    return c;
}

int nbs_png_load(FILE *f, nbs_image8 *out)
{
    unsigned char sig[8];
    unsigned char *idat = NULL, *raw = NULL;
    unsigned char pal[256 * 3];
    size_t idat_n = 0, idat_cap = 0, raw_n = 0;
    unsigned long w = 0, h = 0;
    int depth = 0, ctype = 0, interlace = 0;
    int nch_in, bpp_bytes, rowbytes;
    int have_ihdr = 0, npal = 0;
    unsigned long y;
    int x, rc = 1;

    rewind(f);
    if (fread(sig, 1, 8, f) != 8) { nbs_set_error("PNG: short file"); return 1; }

    for (;;) {
        unsigned char ch[8];
        unsigned long clen;
        if (fread(ch, 1, 8, f) != 8) { nbs_set_error("PNG: truncated chunk header"); goto done; }
        clen = be32(ch);
        if (clen > 0x7FFFFFFFUL) { nbs_set_error("PNG: absurd chunk length"); goto done; }

        if (memcmp(ch + 4, "IHDR", 4) == 0) {
            unsigned char b[13];
            if (clen != 13 || fread(b, 1, 13, f) != 13) { nbs_set_error("PNG: bad IHDR"); goto done; }
            w = be32(b); h = be32(b + 4);
            depth = b[8]; ctype = b[9]; interlace = b[12];
            have_ihdr = 1;
            if (w == 0 || h == 0) { nbs_set_error("PNG: zero dimensions"); goto done; }
            if (interlace != 0) { nbs_set_error("PNG: Adam7 interlacing not supported"); goto done; }
            if (depth != 1 && depth != 2 && depth != 4 && depth != 8 && depth != 16) {
                nbs_set_error("PNG: unsupported bit depth"); goto done;
            }
            if (ctype != 0 && ctype != 2 && ctype != 3 && ctype != 4 && ctype != 6) {
                nbs_set_error("PNG: unsupported colour type"); goto done;
            }
        } else if (memcmp(ch + 4, "PLTE", 4) == 0) {
            if (clen > 256 * 3 || (clen % 3) != 0) { nbs_set_error("PNG: bad PLTE"); goto done; }
            if (fread(pal, 1, (size_t)clen, f) != (size_t)clen) { nbs_set_error("PNG: short PLTE"); goto done; }
            npal = (int)(clen / 3);
        } else if (memcmp(ch + 4, "IDAT", 4) == 0) {
            if (idat_n + clen > idat_cap) {
                size_t nc = idat_cap ? idat_cap : 65536;
                unsigned char *p;
                while (nc < idat_n + clen) nc *= 2;
                p = (unsigned char *)realloc(idat, nc);
                if (p == NULL) { nbs_set_error("PNG: out of memory"); goto done; }
                idat = p; idat_cap = nc;
            }
            if (fread(idat + idat_n, 1, (size_t)clen, f) != (size_t)clen) {
                nbs_set_error("PNG: short IDAT"); goto done;
            }
            idat_n += clen;
        } else if (memcmp(ch + 4, "IEND", 4) == 0) {
            break;
        } else {
            if (fseek(f, (long)clen, SEEK_CUR) != 0) { nbs_set_error("PNG: seek failed"); goto done; }
        }
        if (fseek(f, 4, SEEK_CUR) != 0) { nbs_set_error("PNG: seek past CRC failed"); goto done; }
    }

    if (!have_ihdr || idat_n == 0) { nbs_set_error("PNG: missing IHDR or IDAT"); goto done; }

    switch (ctype) {
    case 0: nch_in = 1; break;
    case 2: nch_in = 3; break;
    case 3: nch_in = 1; break;
    case 4: nch_in = 2; break;
    default: nch_in = 4; break;
    }
    if (ctype == 3 && npal == 0) { nbs_set_error("PNG: palette image without PLTE"); goto done; }

    bpp_bytes = (nch_in * depth + 7) / 8;
    rowbytes  = (int)(((unsigned long)nch_in * (unsigned long)depth * w + 7UL) / 8UL);

    if (zlib_inflate(idat, idat_n, &raw, &raw_n) != 0) goto done;
    if (raw_n < (size_t)(rowbytes + 1) * (size_t)h) {
        nbs_set_error("PNG: decompressed data too short"); goto done;
    }

    /* Un-filter in place, producing a contiguous rowbytes*h buffer. */
    {
        unsigned char *cur, *prv = NULL;
        for (y = 0; y < h; ++y) {
            unsigned char ft = raw[(size_t)y * (size_t)(rowbytes + 1)];
            cur = raw + (size_t)y * (size_t)(rowbytes + 1) + 1;
            switch (ft) {
            case 0: break;
            case 1:
                for (x = bpp_bytes; x < rowbytes; ++x)
                    cur[x] = (unsigned char)(cur[x] + cur[x - bpp_bytes]);
                break;
            case 2:
                if (prv) for (x = 0; x < rowbytes; ++x)
                    cur[x] = (unsigned char)(cur[x] + prv[x]);
                break;
            case 3:
                for (x = 0; x < rowbytes; ++x) {
                    int a = (x >= bpp_bytes) ? cur[x - bpp_bytes] : 0;
                    int b = prv ? prv[x] : 0;
                    cur[x] = (unsigned char)(cur[x] + ((a + b) >> 1));
                }
                break;
            case 4:
                for (x = 0; x < rowbytes; ++x) {
                    int a = (x >= bpp_bytes) ? cur[x - bpp_bytes] : 0;
                    int b = prv ? prv[x] : 0;
                    int c = (prv && x >= bpp_bytes) ? prv[x - bpp_bytes] : 0;
                    cur[x] = (unsigned char)(cur[x] + paeth(a, b, c));
                }
                break;
            default:
                nbs_set_error("PNG: unknown row filter"); goto done;
            }
            prv = cur;
        }
    }

    /* Expand to 8-bit RGB (or 8-bit grey for ctype 0/4). */
    {
        int nch_out = (ctype == 0 || ctype == 4) ? 1 : 3;
        out->px = (unsigned char *)nbs_malloc((size_t)w * (size_t)h * (size_t)nch_out);
        if (out->px == NULL) goto done;
        out->w = (int)w; out->h = (int)h; out->nch = nch_out;

        for (y = 0; y < h; ++y) {
            const unsigned char *src = raw + (size_t)y * (size_t)(rowbytes + 1) + 1;
            unsigned char *dst = out->px + (size_t)y * (size_t)w * (size_t)nch_out;
            for (x = 0; x < (int)w; ++x) {
                int s[4];
                int c;
                if (depth == 8) {
                    for (c = 0; c < nch_in; ++c) s[c] = src[x * nch_in + c];
                } else if (depth == 16) {
                    for (c = 0; c < nch_in; ++c) s[c] = src[(x * nch_in + c) * 2];
                } else {
                    /* sub-byte, only meaningful for ctype 0 and 3 (nch_in == 1) */
                    int per = 8 / depth;
                    int bitpos = (x % per) * depth;
                    int v = src[x / per];
                    int mask = (1 << depth) - 1;
                    s[0] = (v >> (8 - depth - bitpos)) & mask;
                    if (ctype == 0) s[0] = s[0] * 255 / mask;
                }
                if (ctype == 3) {
                    int idx = s[0];
                    if (idx >= npal) idx = npal - 1;
                    dst[3 * x + 0] = pal[3 * idx + 0];
                    dst[3 * x + 1] = pal[3 * idx + 1];
                    dst[3 * x + 2] = pal[3 * idx + 2];
                } else if (nch_out == 1) {
                    dst[x] = (unsigned char)s[0];
                } else {
                    dst[3 * x + 0] = (unsigned char)s[0];
                    dst[3 * x + 1] = (unsigned char)s[1];
                    dst[3 * x + 2] = (unsigned char)s[2];
                }
            }
        }
    }
    rc = 0;

done:
    free(idat);
    free(raw);
    if (rc != 0) nbs_image8_free(out);
    return rc;
}
