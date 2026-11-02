#!/usr/bin/env python3
"""Turn a subjective image quality database into an nbssim index file.

    python3 scripts/make_index.py live   /path/to/databaserelease2  -o data/live.csv
    python3 scripts/make_index.py tid2013 /path/to/tid2013          -o data/tid2013.csv
    python3 scripts/make_index.py csiq   /path/to/CSIQ              -o data/csiq.csv

The index format is documented in include/nbs_data.h.  Scores are written in
each database's own units and orientation; the `lower_better` header records
which way is which, and nbsbench fits a monotone logistic either way, so no
rescaling is applied here.

Pixels per degree
-----------------
None of these three databases publishes a ppd figure, and NB-SSIM's band
frequencies are specified in cycles per degree, so ppd has to be supplied.
The defaults below are derived, with the reasoning recorded in docs/DATASETS.md.
Override with --ppd and check the sensitivity before trusting any conclusion
that depends on it.
"""
import argparse
import os
import re
import sys

# CID:IQ ppd comes from its own RunMetrics.m (dpi 92, 50/100 cm).
CIDIQ_PPD = {50: 31.612, 100: 63.223}
DEFAULT_PPD = {'live': 30.0, 'live_r1': 30.0, 'tid2013': 32.0, 'csiq': 50.0,
               'cidiq': CIDIQ_PPD[50]}

# LIVE Release 2: the 982 score slots are one concatenated index in this order.
LIVE_BLOCKS = [('jp2k', 227), ('jpeg', 233), ('wn', 174),
               ('gblur', 174), ('fastfading', 174)]


def write_index(path, name, score, lower_better, ppd, root, rows):
    with open(path, 'w') as f:
        f.write('# nbsdataset 1\n')
        f.write('# name: %s\n' % name)
        f.write('# root: %s\n' % root)
        f.write('# score: %s\n' % score)
        f.write('# lower_better: %d\n' % (1 if lower_better else 0))
        f.write('# ppd: %g\n' % ppd)
        for ref, dis, s, g in rows:
            f.write('%s,%s,%.6f,%s\n' % (ref, dis, s, g))
    print('wrote %s: %d pairs, %d reference groups'
          % (path, len(rows), len(set(r[3] for r in rows))))


def build_live(root, args):
    from scipy.io import loadmat

    realigned = os.path.join(root, 'dmos_realigned.mat')
    plain = os.path.join(root, 'dmos.mat')
    if os.path.exists(realigned):
        m = loadmat(realigned)
        key = 'dmos_new' if 'dmos_new' in m else 'dmos'
        dmos = m[key].ravel()
        orgs = m['orgs'].ravel()
        which = 'dmos_realigned.mat (%s)' % key
    elif os.path.exists(plain):
        m = loadmat(plain)
        dmos = m['dmos'].ravel()
        orgs = m['orgs'].ravel()
        which = 'dmos.mat'
        print('note: dmos_realigned.mat not found, falling back to dmos.mat.\n'
              '      The realigned scores are the ones the literature uses; get\n'
              '      them from live.ece.utexas.edu/research/quality/release2/',
              file=sys.stderr)
    else:
        sys.exit('neither dmos_realigned.mat nor dmos.mat found in %s' % root)

    refm = loadmat(os.path.join(root, 'refnames_all.mat'))['refnames_all'].ravel()
    refnames = [str(x[0]) for x in refm]

    total = sum(n for _, n in LIVE_BLOCKS)
    if len(dmos) != total or len(refnames) != total:
        sys.exit('expected %d entries, got dmos=%d refnames=%d'
                 % (total, len(dmos), len(refnames)))

    rows = []
    i = 0
    for folder, count in LIVE_BLOCKS:
        for k in range(1, count + 1):
            if orgs[i] == 0:
                ref = 'refimgs/' + refnames[i]
                dis = '%s/img%d.bmp' % (folder, k)
                rows.append((ref, dis, float(dmos[i]),
                             os.path.splitext(refnames[i])[0]))
            i += 1
    print('LIVE: read scores from %s; %d distorted of %d slots'
          % (which, len(rows), total))
    return ('LIVE-R2', 'dmos', True, rows)


def build_tid2013(root, args):
    mos_path = None
    for cand in ('mos_with_names.txt', 'MOS_WITH_NAMES.TXT'):
        p = os.path.join(root, cand)
        if os.path.exists(p):
            mos_path = p
            break
    if mos_path is None:
        sys.exit('mos_with_names.txt not found in %s' % root)

    refdir = None
    for cand in ('reference_images', 'reference_images'.upper()):
        if os.path.isdir(os.path.join(root, cand)):
            refdir = cand
            break
    distdir = None
    for cand in ('distorted_images', 'distorted_images'.upper()):
        if os.path.isdir(os.path.join(root, cand)):
            distdir = cand
            break
    if refdir is None or distdir is None:
        sys.exit('reference_images/ or distorted_images/ not found in %s' % root)

    # TID2013 ships mixed-case names: mos_with_names.txt lists "i01_01_2.bmp"
    # while the file on disk may be "I01_01_2.BMP".  Harmless on Windows,
    # fatal on a case-sensitive filesystem, so resolve every name against a
    # case-folded index of what is actually there.
    refindex = {}
    for fn in os.listdir(os.path.join(root, refdir)):
        m = re.match(r'^[iI](\d\d)\.(bmp|BMP)$', fn)
        if m:
            refindex[m.group(1)] = fn

    distindex = {}
    for fn in os.listdir(os.path.join(root, distdir)):
        distindex[fn.lower()] = fn

    rows = []
    missing = 0
    with open(mos_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) != 2:
                continue
            score, fn = float(parts[0]), parts[1]
            m = re.match(r'^[iI](\d\d)_(\d\d)_(\d)\.(bmp|BMP)$', fn)
            if not m:
                continue
            rid = m.group(1)
            if rid not in refindex:
                continue
            if args.types and int(m.group(2)) not in args.types:
                continue
            real = distindex.get(fn.lower())
            if real is None:
                missing += 1
                continue
            rows.append(('%s/%s' % (refdir, refindex[rid]),
                         '%s/%s' % (distdir, real), score, 'i' + rid))
    if missing:
        print('warning: %d entries in %s had no matching image file'
              % (missing, os.path.basename(mos_path)), file=sys.stderr)
    if not rows:
        sys.exit('no usable rows parsed from %s' % mos_path)
    return ('TID2013', 'mos', False, rows)


def build_csiq(root, args):
    from openpyxl import load_workbook

    xls = None
    for cand in ('csiq.DMOS.xlsx', 'CSIQ.DMOS.xlsx', 'csiq.dmos.xlsx'):
        p = os.path.join(root, cand)
        if os.path.exists(p):
            xls = p
            break
    if xls is None:
        sys.exit('csiq.DMOS.xlsx not found in %s' % root)

    # The distribution unpacks in at least two layouts: src_imgs/ + dst_imgs/<TYPE>/,
    # or references at the top level beside bare <TYPE>/ folders.  Rather than
    # guess, index every PNG under root by its lower-cased basename and resolve
    # names from the spreadsheet against that.  This also absorbs the
    # inconsistent casing of the distortion token in the file names.
    byname = {}
    for dirpath, _dirs, files in os.walk(root):
        for fn in files:
            if fn.lower().endswith('.png'):
                rel = os.path.relpath(os.path.join(dirpath, fn), root)
                byname.setdefault(fn.lower(), rel.replace(os.sep, '/'))

    wb = load_workbook(xls, data_only=True)
    sheet = 'all_by_image' if 'all_by_image' in wb.sheetnames else wb.sheetnames[0]
    ws = wb[sheet]
    raw = list(ws.values)

    hdr_i = None
    for i, r in enumerate(raw[:12]):
        cells = [str(c).strip().lower() if c is not None else '' for c in r]
        if 'dmos' in cells and 'image' in cells:
            hdr_i = i
            break
    if hdr_i is None:
        sys.exit("could not locate the header row (need 'image' and 'dmos') "
                 "in sheet '%s'" % sheet)
    hdr = [str(c).strip().lower() if c is not None else '' for c in raw[hdr_i]]
    need = ('image', 'dst_type', 'dst_lev', 'dmos')
    col = dict((k, hdr.index(k)) for k in need if k in hdr)
    if set(need) - set(col):
        sys.exit('sheet %s is missing columns: %s'
                 % (sheet, sorted(set(need) - set(col))))

    # Spreadsheet distortion names vs the token used in file names.
    alias = {'noise': ['awgn', 'noise'],
             'jpeg 2000': ['jpeg2000', 'jp2k'],
             'jpeg2000': ['jpeg2000', 'jp2k'],
             'jpeg': ['jpeg'],
             'blur': ['blur'],
             'fnoise': ['fnoise'],
             'contrast': ['contrast']}

    rows, missing_img, missing_ref = [], 0, 0
    for r in raw[hdr_i + 1:]:
        if r is None or r[col['image']] is None or r[col['dmos']] is None:
            continue
        img = str(r[col['image']]).strip()
        typ = str(r[col['dst_type']]).strip().lower()
        try:
            lev = str(int(float(str(r[col['dst_lev']]).strip())))
            dm = float(r[col['dmos']])
        except (TypeError, ValueError):
            continue

        refrel = byname.get((img + '.png').lower())
        if refrel is None:
            missing_ref += 1
            continue
        distrel = None
        for tok in alias.get(typ, [typ.replace(' ', '')]):
            distrel = byname.get(('%s.%s.%s.png' % (img, tok, lev)).lower())
            if distrel is not None:
                break
        if distrel is None:
            missing_img += 1
            continue
        rows.append((refrel, distrel, dm, img))

    if missing_ref:
        print('warning: %d rows had no reference image' % missing_ref, file=sys.stderr)
    if missing_img:
        print('warning: %d rows had no distorted image' % missing_img, file=sys.stderr)
    if not rows:
        sys.exit('no usable rows built from %s' % xls)
    return ('CSIQ', 'dmos', True, rows)


def build_live_r1(root, args):
    """LIVE Release 1 (2003): the JPEG and JPEG2000 studies, reported as MOS.

    This is the dataset the 2004 SSIM paper was validated on: 175 JPEG and
    169 JPEG2000 distorted images over 29 references.

    Two structural cautions, not cosmetic ones:
      * Each study ran in two sessions with separate score files, and those
        sessions were never realigned onto a common scale -- that realignment
        is precisely what Release 2 adds.  Pooling sessions mixes two scales.
        --session=1|2 keeps them apart; --session=0 pools and warns.
      * Entries with bit rate 0 are lossless copies of the source, used as
        in-study anchors rather than distorted images, and are excluded.
    """
    from scipy.io import loadmat

    studies = []
    for sub, tag, info in (('jpeg2000 database release', 'jp2k', 'jp2info.txt'),
                           ('jpeg database release', 'jpeg', 'jpeginfo.txt')):
        d = os.path.join(root, sub, 'database release')
        if os.path.isdir(d):
            studies.append((d, tag, info))
    if not studies:
        sys.exit('no "<jpeg|jpeg2000> database release/database release" '
                 'directory found under %s' % root)

    want = getattr(args, 'session', 0) or 0
    rows = []
    for d, tag, info in studies:
        ipath = os.path.join(d, info)
        if not os.path.exists(ipath):
            sys.exit('missing %s' % ipath)
        srcof = {}
        for line in open(ipath):
            f = line.split()
            if len(f) >= 2:
                srcof[f[1].strip().lower()] = f[0].strip()

        base = 0
        for sess in (1, 2):
            mp = os.path.join(d, 'scores%d.mat' % sess)
            if not os.path.exists(mp):
                continue
            m = loadmat(mp)
            mmt = m['mmt'].ravel()
            br = m['br'].ravel()
            n = len(mmt)
            if want in (0, sess):
                rel = os.path.relpath(d, root).replace(os.sep, '/')
                for k in range(n):
                    if br[k] == 0:
                        continue
                    dist = 'img%d.bmp' % (base + k + 1)
                    src = srcof.get(dist.lower())
                    if src is None:
                        continue
                    rows.append(('%s/%s' % (rel, src),
                                 '%s/%s' % (rel, dist),
                                 float(mmt[k]),
                                 os.path.splitext(src)[0]))
            base += n

    if want == 0:
        print('warning: pooling sessions whose scores were never realigned '
              'onto a common scale; prefer --session=1 or 2, or Release 2',
              file=sys.stderr)
    if not rows:
        sys.exit('no usable rows built from %s' % root)
    return ('LIVE-R1', 'mos', False, rows)


def build_cidiq(root, args):
    """CID:IQ (Colourlab Image Database: Image Quality), 23 refs x 6 distortions
    x 5 levels = 690 images at 800x800, scored at TWO viewing distances.

    The 690 MOS entries are ordered exactly as the database's own RunMetrics.m
    walks the tree: the six distortion folders in alphabetical order, and within
    each, the reproduction files in alphabetical order.  The reference for a
    reproduction is the first seven characters of its name.

    Viewing geometry is not guessed here -- RunMetrics.m states dpi = 92 and
    viewing distances of 50 and 100 cm, giving 31.61 and 63.22 pixels per
    degree.  Those are the only measured ppd figures in this whole project.
    """
    from scipy.io import loadmat

    dist_cm = getattr(args, 'distance', 50) or 50
    if dist_cm not in (50, 100):
        sys.exit('--distance must be 50 or 100 (CID:IQ was run at both)')

    mp = os.path.join(root, 'MOS%d.mat' % dist_cm)
    if not os.path.exists(mp):
        sys.exit('missing %s' % mp)
    key = 'MOS%d' % dist_cm
    mos = loadmat(mp)[key].ravel()

    repro = os.path.join(root, 'Images', 'Reproduction')
    orig = os.path.join(root, 'Images', 'Original')
    if not os.path.isdir(repro) or not os.path.isdir(orig):
        sys.exit('expected Images/Original and Images/Reproduction under %s' % root)

    folders = sorted(d for d in os.listdir(repro)
                     if os.path.isdir(os.path.join(repro, d)))
    rows = []
    k = 0
    for fol in folders:
        for fn in sorted(f for f in os.listdir(os.path.join(repro, fol))
                         if f.lower().endswith('.bmp')):
            if k >= len(mos):
                break
            ref = fn[:7] + '.bmp'
            if not os.path.exists(os.path.join(orig, ref)):
                k += 1
                continue
            rows.append(('Images/Original/%s' % ref,
                         'Images/Reproduction/%s/%s' % (fol, fn),
                         float(mos[k]), os.path.splitext(ref)[0]))
            k += 1
    if k != len(mos):
        print('warning: walked %d files but MOS has %d entries -- ordering may '
              'not match the database convention' % (k, len(mos)), file=sys.stderr)
    if not rows:
        sys.exit('no usable rows built from %s' % root)
    # MOS here is a quality rating on a 1..9 category scale: higher is better.
    return ('CID:IQ-%dcm' % dist_cm, 'mos', False, rows)


BUILDERS = {'live': build_live, 'live_r1': build_live_r1,
            'tid2013': build_tid2013, 'csiq': build_csiq,
            'cidiq': build_cidiq}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('database', choices=sorted(BUILDERS))
    ap.add_argument('root', help='directory containing the unpacked database')
    ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--ppd', type=float, default=None)
    ap.add_argument('--distance', type=int, default=50,
                    help='CID:IQ only: viewing distance in cm, 50 or 100')
    ap.add_argument('--session', type=int, default=0,
                    help='LIVE R1 only: 1 or 2 to keep a single session, '
                         '0 (default) to pool them')
    ap.add_argument('--types', type=lambda s: [int(v) for v in s.split(',')],
                    default=None,
                    help='TID2013 only: keep just these distortion type numbers, '
                         'e.g. 8,10,11 for blur/JPEG/JPEG2000')
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    if not os.path.isdir(root):
        sys.exit('not a directory: %s' % root)

    name, score, lower_better, rows = BUILDERS[args.database](root, args)
    ppd = args.ppd
    if ppd is None:
        ppd = (CIDIQ_PPD[args.distance] if args.database == 'cidiq'
               else DEFAULT_PPD[args.database])
    outdir = os.path.dirname(os.path.abspath(args.out))
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    write_index(args.out, name, score, lower_better, ppd, root, rows)


if __name__ == '__main__':
    main()
