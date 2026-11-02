#!/usr/bin/env python3
"""Build a synthetic full-reference IQA set for exercising the pipeline.

This exists so the testbench, the tuner and the cross-validation machinery can
be verified end to end before real subjective data is available.  The "scores"
it produces come from a masked, CSF-weighted error model -- a crude stand-in
for a human observer.  Any correlation numbers obtained on this set say
something about the software, and nothing about the metric's agreement with
human judgement.  Do not report them as results.

Usage:
    python3 scripts/make_synth.py --out data/synth --refs 16 --levels 5
"""
import argparse
import os

import numpy as np
from scipy.fft import dctn, idctn
from scipy.ndimage import gaussian_filter, zoom


# ---------------------------------------------------------------- images

def fractal(rng, h, w, beta=1.6):
    """1/f^beta noise field, roughly natural-image-like in its spectrum."""
    fy = np.fft.fftfreq(h)[:, None]
    fx = np.fft.fftfreq(w)[None, :]
    r = np.sqrt(fy ** 2 + fx ** 2)
    r[0, 0] = 1.0
    amp = 1.0 / (r ** beta)
    amp[0, 0] = 0.0
    ph = rng.uniform(0, 2 * np.pi, (h, w))
    f = amp * np.exp(1j * ph)
    img = np.real(np.fft.ifft2(f))
    img -= img.min()
    if img.max() > 0:
        img /= img.max()
    return img


def make_reference(rng, h, w):
    img = 0.65 * fractal(rng, h, w, rng.uniform(1.3, 2.0))
    # A few hard edges and flat regions, so blocking and ringing have
    # something to attach to.
    for _ in range(rng.integers(3, 8)):
        y0, x0 = rng.integers(0, h - 8), rng.integers(0, w - 8)
        dy, dx = rng.integers(20, max(21, h // 3)), rng.integers(20, max(21, w // 3))
        img[y0:y0 + dy, x0:x0 + dx] += rng.uniform(-0.35, 0.35)
    # Oriented texture patch.
    yy, xx = np.mgrid[0:h, 0:w]
    th = rng.uniform(0, np.pi)
    fq = rng.uniform(0.15, 0.45)
    tex = 0.12 * np.sin(2 * np.pi * fq * (xx * np.cos(th) + yy * np.sin(th)))
    m = np.zeros((h, w))
    y0, x0 = rng.integers(0, h // 2), rng.integers(0, w // 2)
    m[y0:y0 + h // 3, x0:x0 + w // 3] = 1.0
    img += tex * gaussian_filter(m, 8)
    img = gaussian_filter(img, 0.6)
    img -= img.min()
    img /= max(img.max(), 1e-9)
    return np.round(20 + 215 * img)


# ------------------------------------------------------------ distortions

def blockdct(img, q):
    h, w = img.shape
    H, W = (h // 8) * 8, (w // 8) * 8
    out = img.copy()
    a = img[:H, :W].reshape(H // 8, 8, W // 8, 8).transpose(0, 2, 1, 3)
    c = dctn(a, axes=(2, 3), norm='ortho')
    qm = q * (1.0 + 0.9 * (np.add.outer(np.arange(8), np.arange(8))))
    c = np.round(c / qm) * qm
    b = idctn(c, axes=(2, 3), norm='ortho')
    out[:H, :W] = b.transpose(0, 2, 1, 3).reshape(H, W)
    return out


def distort(rng, ref, kind, level):
    t = (level + 1) / 5.0
    if kind == 'blur':
        return gaussian_filter(ref, 0.4 + 3.2 * t ** 1.5, mode='reflect')
    if kind == 'noise':
        return ref + (2.0 + 34.0 * t ** 1.6) * rng.standard_normal(ref.shape)
    if kind == 'block':
        return blockdct(ref, 4.0 + 70.0 * t ** 1.8)
    if kind == 'scale':
        f = 1.0 / (1.0 + 2.2 * t)
        small = zoom(ref, f, order=1)
        back = zoom(small, (ref.shape[0] / small.shape[0],
                            ref.shape[1] / small.shape[1]), order=1)
        return back[:ref.shape[0], :ref.shape[1]]
    if kind == 'contrast':
        m = ref.mean()
        return m + (ref - m) * (1.0 - 0.55 * t)
    raise ValueError(kind)


# --------------------------------------------------------- observer model

def proxy_dmos(ref, dis, ppd):
    """A masked multi-channel error model, deliberately unlike SSIM.

    Error is decomposed into three radial bands, each divided by a local
    masking energy taken from the reference, weighted by a Barten-like CSF at
    the band centre, and pooled with a Minkowski exponent.  Returns a DMOS-like
    number where larger means worse.
    """
    err = dis - ref
    total = 0.0
    for fc in (1.0, 4.0, 12.0):
        sigma = ppd / (2 * np.pi * fc)
        lo = gaussian_filter(err, sigma * 1.4, mode='reflect')
        hi = gaussian_filter(err, sigma / 1.4, mode='reflect')
        band = hi - lo
        rlo = gaussian_filter(ref, sigma * 1.4, mode='reflect')
        rhi = gaussian_filter(ref, sigma / 1.4, mode='reflect')
        rband = rhi - rlo
        mask = gaussian_filter(np.abs(rband), 6.0, mode='reflect')
        # Barten-ish contrast sensitivity, normalised to peak 1 near 4 cpd.
        csf = (fc ** 0.9) * np.exp(-0.25 * fc) / (4.0 ** 0.9 * np.exp(-1.0))
        norm = np.abs(band) / (1.0 + 0.55 * mask)
        total += float(np.mean((csf * norm) ** 3.0))
    d = total ** (1.0 / 3.0)
    return 100.0 * (1.0 - 1.0 / (1.0 + (d / 1.6) ** 1.25))


# ---------------------------------------------------------------- driver

def write_pgm(path, arr):
    a = np.clip(np.round(arr), 0, 255).astype(np.uint8)
    with open(path, 'wb') as f:
        f.write(b'P5\n%d %d\n255\n' % (a.shape[1], a.shape[0]))
        f.write(a.tobytes())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='data/synth')
    ap.add_argument('--refs', type=int, default=16)
    ap.add_argument('--levels', type=int, default=5)
    ap.add_argument('--width', type=int, default=512)
    ap.add_argument('--height', type=int, default=384)
    ap.add_argument('--ppd', type=float, default=30.0)
    ap.add_argument('--seed', type=int, default=11)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    os.makedirs(os.path.join(args.out, 'ref'), exist_ok=True)
    os.makedirs(os.path.join(args.out, 'dist'), exist_ok=True)

    kinds = ['blur', 'noise', 'block', 'scale', 'contrast']
    rows = []
    for r in range(args.refs):
        ref = make_reference(rng, args.height, args.width)
        rname = 'ref/r%02d.pgm' % r
        write_pgm(os.path.join(args.out, rname), ref)
        refq = np.clip(np.round(ref), 0, 255)
        for kind in kinds:
            for lv in range(args.levels):
                dis = np.clip(np.round(distort(rng, ref, kind, lv)), 0, 255)
                dname = 'dist/r%02d_%s_%d.pgm' % (r, kind, lv)
                write_pgm(os.path.join(args.out, dname), dis)
                s = proxy_dmos(refq, dis, args.ppd)
                s += rng.normal(0, 2.0)          # simulated observer noise
                rows.append((rname, dname, s, 'r%02d' % r))
        print('reference %d/%d' % (r + 1, args.refs), end='\r')
    print()

    idx = os.path.join(args.out, 'index.csv')
    with open(idx, 'w') as f:
        f.write('# nbsdataset 1\n')
        f.write('# name: SYNTHETIC (proxy scores -- not human data)\n')
        f.write('# score: proxy_dmos\n')
        f.write('# lower_better: 1\n')
        f.write('# ppd: %g\n' % args.ppd)
        for a, b, s, g in rows:
            f.write('%s,%s,%.4f,%s\n' % (a, b, s, g))
    print('wrote %s (%d pairs, %d references)' % (idx, len(rows), args.refs))


if __name__ == '__main__':
    main()
