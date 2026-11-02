/* data.c -- dataset index file parsing.  Strict C89. */

#include "nbs.h"
#include "nbs_data.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *dupstr(const char *s)
{
    size_t n = strlen(s) + 1;
    char *p = (char *)nbs_malloc(n);
    if (p) memcpy(p, s, n);
    return p;
}

static void rstrip(char *s)
{
    size_t n = strlen(s);
    while (n > 0 && (s[n - 1] == '\n' || s[n - 1] == '\r' ||
                     s[n - 1] == ' '  || s[n - 1] == '\t')) {
        s[--n] = '\0';
    }
}

static char *lstrip(char *s)
{
    while (*s == ' ' || *s == '\t') s++;
    return s;
}

void nbs_path_join(char *buf, int bufn, const char *root, const char *rel)
{
    if (rel[0] == '/' || root == NULL || root[0] == '\0') {
        strncpy(buf, rel, (size_t)bufn - 1);
        buf[bufn - 1] = '\0';
        return;
    }
    strncpy(buf, root, (size_t)bufn - 1);
    buf[bufn - 1] = '\0';
    if (buf[0] != '\0' && buf[strlen(buf) - 1] != '/')
        strncat(buf, "/", (size_t)bufn - strlen(buf) - 1);
    strncat(buf, rel, (size_t)bufn - strlen(buf) - 1);
}

static int group_of(nbs_dataset *ds, const char *refid, int upto)
{
    int i;
    for (i = 0; i < upto; ++i)
        if (strcmp(ds->e[i].refid, refid) == 0) return ds->e[i].group;
    return -1;
}

void nbs_dataset_free(nbs_dataset *ds)
{
    int i;
    if (ds == NULL || ds->e == NULL) return;
    for (i = 0; i < ds->n; ++i) {
        free(ds->e[i].ref);
        free(ds->e[i].dist);
        free(ds->e[i].refid);
    }
    free(ds->e);
    ds->e = NULL;
    ds->n = 0;
}

int nbs_dataset_load(const char *path, nbs_dataset *ds)
{
    FILE *f;
    char line[4096];
    int cap = 256, i;
    char dirbuf[1024];
    const char *slash;

    memset(ds, 0, sizeof(*ds));
    ds->ppd = 30.0;
    ds->lower_better = 1;
    strcpy(ds->name, "dataset");
    strcpy(ds->score_name, "score");

    /* Default root = directory containing the index file. */
    slash = strrchr(path, '/');
    if (slash != NULL) {
        size_t n = (size_t)(slash - path);
        if (n >= sizeof(dirbuf)) n = sizeof(dirbuf) - 1;
        memcpy(dirbuf, path, n);
        dirbuf[n] = '\0';
    } else {
        strcpy(dirbuf, ".");
    }
    memcpy(ds->root, dirbuf, strlen(dirbuf) + 1);

    f = fopen(path, "r");
    if (f == NULL) { nbs_set_error("cannot open dataset index"); return 1; }

    ds->e = (nbs_entry *)nbs_malloc((size_t)cap * sizeof(nbs_entry));
    if (ds->e == NULL) { fclose(f); return 1; }

    while (fgets(line, (int)sizeof(line), f) != NULL) {
        char *p = lstrip(line);
        rstrip(p);
        if (p[0] == '\0') continue;
        if (p[0] == '#') {
            char *q = lstrip(p + 1);
            if (strncmp(q, "name:", 5) == 0) {
                strncpy(ds->name, lstrip(q + 5), sizeof(ds->name) - 1);
            } else if (strncmp(q, "root:", 5) == 0) {
                char *r = lstrip(q + 5);
                if (r[0] == '/') {
                    strncpy(ds->root, r, sizeof(ds->root) - 1);
                } else {
                    nbs_path_join(ds->root, (int)sizeof(ds->root), dirbuf, r);
                }
            } else if (strncmp(q, "score:", 6) == 0) {
                strncpy(ds->score_name, lstrip(q + 6), sizeof(ds->score_name) - 1);
            } else if (strncmp(q, "lower_better:", 13) == 0) {
                ds->lower_better = atoi(lstrip(q + 13));
            } else if (strncmp(q, "ppd:", 4) == 0) {
                ds->ppd = atof(lstrip(q + 4));
            }
            continue;
        }
        {
            char *f1, *f2, *f3, *f4, *c;
            f1 = p;
            c = strchr(f1, ','); if (c == NULL) continue; *c = '\0'; f2 = c + 1;
            c = strchr(f2, ','); if (c == NULL) continue; *c = '\0'; f3 = c + 1;
            c = strchr(f3, ',');
            if (c != NULL) { *c = '\0'; f4 = c + 1; } else { f4 = f1; }

            if (ds->n == cap) {
                nbs_entry *ne;
                cap *= 2;
                ne = (nbs_entry *)realloc(ds->e, (size_t)cap * sizeof(nbs_entry));
                if (ne == NULL) { nbs_set_error("out of memory"); fclose(f); return 1; }
                ds->e = ne;
            }
            ds->e[ds->n].ref   = dupstr(f1);
            ds->e[ds->n].dist  = dupstr(f2);
            ds->e[ds->n].score = atof(f3);
            ds->e[ds->n].refid = dupstr(f4);
            ds->e[ds->n].group = -1;
            if (ds->e[ds->n].ref == NULL || ds->e[ds->n].dist == NULL ||
                ds->e[ds->n].refid == NULL) { fclose(f); return 1; }
            ds->n++;
        }
    }
    fclose(f);

    for (i = 0; i < ds->n; ++i) {
        int g = group_of(ds, ds->e[i].refid, i);
        if (g < 0) g = ds->ngroups++;
        ds->e[i].group = g;
    }

    if (ds->n == 0) { nbs_set_error("dataset index contains no entries"); return 1; }
    return 0;
}
