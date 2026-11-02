#!/usr/bin/env python3
"""Independent numpy reference for SSIM and MS-SSIM, used to cross-check the
C implementation.  This is a direct transcription of Z. Wang's ssim_index.m
and ms_ssim.m, written without reference to the C code.

Usage:  python3 scripts/ref_check.py [--bin ./bin/nbsmetric]
"""
import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy.ndimage import correlate1d


def gauss_win(n=11, sigma=1.5):
    r = (n - 1) // 2
    x = np.arange(-r, r + 1, dtype=np.float64)
    w = np.exp(-(x ** 2) / (2.0 * sigma ** 2))
    return w / w.sum()


def filter_valid(img, w):
    """Separable 'valid' correlation, as filter2(win, img, 'valid')."""
    r = (len(w) - 1) // 2
    t = correlate1d(img, w, axis=1, mode='constant')
    t = correlate1d(t, w, axis=0, mode='constant')
    return t[r:img.shape[0] - r, r:img.shape[1] - r]


def box_downsample(img, f):
    """ones(f)/f^2 with 'symmetric' padding, then img[::f, ::f]."""
    lo = -((f - 1) // 2)
    hi = f // 2
    H, W = img.shape
    acc = np.zeros_like(img)
    for dy in range(lo, hi + 1):
        yy = np.clip(np.arange(H) + dy, None, None)
        yy = np.array([reflect(v, H) for v in np.arange(H) + dy])
        for dx in range(lo, hi + 1):
            xx = np.array([reflect(v, W) for v in np.arange(W) + dx])
            acc += img[np.ix_(yy, xx)]
    acc /= f * f
    return acc[::f, ::f]


def reflect(i, n):
    while i < 0 or i >= n:
        if i < 0:
            i = -i - 1
        if i >= n:
            i = 2 * n - i - 1
    return i


def ssim_index(x, y, K1=0.01, K2=0.03, L=255.0, win=11, sigma=1.5):
    w = gauss_win(win, sigma)
    C1, C2 = (K1 * L) ** 2, (K2 * L) ** 2
    mx = filter_valid(x, w)
    my = filter_valid(y, w)
    mxx, myy, mxy = mx * mx, my * my, mx * my
    sxx = filter_valid(x * x, w) - mxx
    syy = filter_valid(y * y, w) - myy
    sxy = filter_valid(x * y, w) - mxy
    lum = (2 * mxy + C1) / (mxx + myy + C1)
    cs = (2 * sxy + C2) / (sxx + syy + C2)
    return float((lum * cs).mean()), float(cs.mean())


def ssim_auto(x, y):
    f = max(1, int(round(min(x.shape) / 256.0)))
    if f > 1:
        x = box_downsample(x, f)
        y = box_downsample(y, f)
    return ssim_index(x, y)[0]


def ms_ssim(x, y, weights=(0.0448, 0.2856, 0.3001, 0.2363, 0.1333)):
    mssim, mcs = [], []
    a, b = x, y
    for i in range(len(weights)):
        s, c = ssim_index(a, b)
        mssim.append(s)
        mcs.append(c)
        if i + 1 < len(weights):
            a = box_downsample(a, 2)
            b = box_downsample(b, 2)
    out = 1.0
    for i in range(len(weights) - 1):
        out *= mcs[i] ** weights[i]
    out *= mssim[-1] ** weights[-1]
    return float(out)


def write_pgm(path, arr):
    a = np.clip(np.round(arr), 0, 255).astype(np.uint8)
    with open(path, 'wb') as f:
        f.write(b'P5\n%d %d\n255\n' % (a.shape[1], a.shape[0]))
        f.write(a.tobytes())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bin', default='./bin/nbsmetric')
    args = ap.parse_args()

    rng = np.random.default_rng(7)
    tmp = tempfile.mkdtemp()
    worst_ssim = worst_ms = 0.0
    print('%-28s %10s %10s %10s %10s' %
          ('case', 'C SSIM', 'py SSIM', 'C MS-SSIM', 'py MS-SSIM'))

    cases = []
    for (H, W) in [(256, 256), (384, 512), (200, 173)]:
        yy, xx = np.mgrid[0:H, 0:W]
        base = (128 + 60 * np.sin(xx * 0.11) * np.cos(yy * 0.07)
                + 25 * np.sin((xx + 2 * yy) * 0.31)
                + 30 * rng.standard_normal((H, W)))
        base = np.clip(base, 0, 255)
        base = np.round(base)
        for amp in (0.0, 4.0, 16.0, 48.0):
            d = np.clip(np.round(base + amp * rng.standard_normal((H, W))), 0, 255)
            cases.append(('%dx%d noise=%.0f' % (W, H, amp), base, d))
        from scipy.ndimage import gaussian_filter
        for s in (1.0, 3.0):
            d = np.round(np.clip(gaussian_filter(base, s, mode='reflect'), 0, 255))
            cases.append(('%dx%d blur=%.1f' % (W, H, s), base, d))

    for name, a, b in cases:
        pa = os.path.join(tmp, 'a.pgm')
        pb = os.path.join(tmp, 'b.pgm')
        write_pgm(pa, a)
        write_pgm(pb, b)
        out = subprocess.run([args.bin, '--quiet', pa, pb],
                             capture_output=True, text=True)
        if out.returncode != 0:
            print('FAILED to run: %s' % out.stderr.strip())
            return 1
        c_ssim, c_ms, _ = [float(v) for v in out.stdout.split()]
        p_ssim = ssim_auto(a, b)
        try:
            p_ms = ms_ssim(a, b)
        except Exception:
            p_ms = float('nan')
        worst_ssim = max(worst_ssim, abs(c_ssim - p_ssim))
        if p_ms == p_ms and c_ms >= 0:
            worst_ms = max(worst_ms, abs(c_ms - p_ms))
        print('%-28s %10.7f %10.7f %10.7f %10.7f' %
              (name, c_ssim, p_ssim, c_ms, p_ms))

    print('\nmax |C - numpy|  SSIM: %.3e   MS-SSIM: %.3e' % (worst_ssim, worst_ms))
    ok = worst_ssim < 1e-10 and worst_ms < 1e-10
    print('AGREEMENT: %s' % ('OK' if ok else 'MISMATCH'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
