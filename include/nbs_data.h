/* nbs_data.h -- dataset index files.
 *
 * An index file is plain text:
 *
 *   # nbsdataset 1
 *   # name: LIVE-R2
 *   # root: /abs/or/relative/path        (optional, default: index file dir)
 *   # score: dmos
 *   # lower_better: 1                    (1 if a larger score means worse)
 *   # ppd: 30.0                          (pixels per degree of visual angle)
 *   ref_path,dist_path,score,ref_id
 *   refimgs/bikes.bmp,jp2k/img1.bmp,63.9,bikes
 *
 * Fields are comma separated.  `ref_id` groups images by source content and
 * is used for content-disjoint cross-validation.
 *
 * Strict C89.
 */

#ifndef NBS_DATA_H
#define NBS_DATA_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char  *ref;
    char  *dist;
    double score;
    char  *refid;
    int    group;     /* dense index of refid */
} nbs_entry;

typedef struct {
    char       name[128];
    char       score_name[64];
    char       root[1024];
    int        lower_better;
    double     ppd;
    nbs_entry *e;
    int        n;
    int        ngroups;
} nbs_dataset;

int  nbs_dataset_load(const char *path, nbs_dataset *ds);
void nbs_dataset_free(nbs_dataset *ds);

/* Compose root + "/" + rel into buf (size bufn).  Absolute rel is copied. */
void nbs_path_join(char *buf, int bufn, const char *root, const char *rel);

#ifdef __cplusplus
}
#endif

#endif /* NBS_DATA_H */
