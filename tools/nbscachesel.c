/* nbscachesel.c -- select entries of a log-xi cache by reference group.
 *
 *   nbscachesel IN.bin OUT.bin --keep=1,2     keep entries whose group is listed
 *   nbscachesel IN.bin OUT.bin --drop=3       drop entries whose group is listed
 *   nbscachesel IN.bin --list                 print the groups and entry counts
 *
 * Lets leave-one-content-out fits reuse one tensor per condition instead of
 * rebuilding a tensor per fold.  Strict C89.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nbs_cache.h"

static int in_list(const char *list, int g)
{
    const char *p = list;
    while (*p) {
        char *end;
        long v = strtol(p, &end, 10);
        if (end == p) break;
        if ((int)v == g) return 1;
        p = end;
        while (*p == ',' || *p == ' ') ++p;
    }
    return 0;
}

int main(int argc, char **argv)
{
    nbs_logxi_cache c, o;
    const char *keep = NULL, *drop = NULL;
    int list = 0, i, k;
    long total = 0;
    if (argc < 3) { fprintf(stderr, "usage: nbscachesel IN.bin OUT.bin --keep=g,g | --drop=g,g | IN.bin --list\n"); return 2; }
    for (i = 2; i < argc; ++i) {
        if (strncmp(argv[i], "--keep=", 7) == 0) keep = argv[i] + 7;
        else if (strncmp(argv[i], "--drop=", 7) == 0) drop = argv[i] + 7;
        else if (strcmp(argv[i], "--list") == 0) list = 1;
    }
    memset(&c, 0, sizeof c);
    if (nbs_cache_load(&c, argv[1]) != 0) { fprintf(stderr, "%s\n", nbs_last_error()); return 1; }
    if (list) {
        int counts[1024]; memset(counts, 0, sizeof counts);
        for (i = 0; i < c.n; ++i) if (c.group[i] >= 0 && c.group[i] < 1024) counts[c.group[i]]++;
        printf("%d entries, %d exponents, ppd %.4f, %d bands\n", c.n, c.nw, c.ppd, c.nbands);
        for (i = 0; i < 1024; ++i) if (counts[i]) printf("group %d: %d entries\n", i, counts[i]);
        nbs_cache_free(&c); return 0;
    }
    if (!keep && !drop) { fprintf(stderr, "need --keep or --drop\n"); return 2; }
    o = c;                                   /* header fields copied */
    o.n = 0;
    o.off  = (long *)malloc(sizeof(long) * (c.n + 1));
    o.np   = (int *)malloc(sizeof(int) * (c.n + 1));
    o.subj = (double *)malloc(sizeof(double) * (c.n + 1));
    o.group = (int *)malloc(sizeof(int) * (c.n + 1));
    for (i = 0; i < c.n; ++i) total += (long)c.np[i] * c.nw;
    o.L = (float *)malloc(sizeof(float) * (total > 0 ? total : 1));
    total = 0;
    for (i = 0; i < c.n; ++i) {
        int sel = keep ? in_list(keep, c.group[i]) : !in_list(drop, c.group[i]);
        if (!sel) continue;
        k = o.n++;
        o.np[k] = c.np[i]; o.subj[k] = c.subj[i]; o.group[k] = c.group[i];
        o.off[k] = total;
        memcpy(o.L + total, c.L + c.off[i], sizeof(float) * (size_t)c.np[i] * c.nw);
        total += (long)c.np[i] * c.nw;
    }
    if (nbs_cache_save(&o, argv[2]) != 0) { fprintf(stderr, "%s\n", nbs_last_error()); return 1; }
    printf("%s: kept %d of %d entries\n", argv[2], o.n, c.n);
    return 0;
}
