# Subjective databases

The SSIM paper (Wang et al., IEEE TIP 2004, §IV-B) evaluated on 344 JPEG and
JPEG2000 images from what became **LIVE Release 2**; the MS-SSIM paper (Wang,
Simoncelli & Bovik, Asilomar 2003) used the same LIVE material plus its own
set. The paper's primary database is CID:IQ, the one that rates the same images at
two viewing distances; LIVE Release 2 is the second, as the database SSIM and
MS-SSIM were designed against. TID2013 and CSIQ are described for completeness.

None of the three can be downloaded from an automated environment: LIVE is
password-gated behind a request form, TID2013 ships as a 913 MB `.rar`, and
CSIQ's original host is dead. They have to be fetched by hand.

---

## CID:IQ (Colourlab Image Database: Image Quality)

**Get it**: from the Norwegian Colour and Visual Computing Laboratory
(NTNU), <https://www.ntnu.edu/colourlab/software>; cite Liu, Pedersen &
Hardeberg, *CID:IQ -- a new image quality database*, ICISP 2014, LNCS 8509.

**Contents**: 23 reference images (`Images/Original/final01.bmp` ...
`final23.bmp`, 800x800, 24-bit BMP) and 690 distorted versions under
`Images/Reproduction/<1_JPEG2000_Compression ... 6_DeltaE_Gamut_Mapping>/`,
six distortion types at five levels each, named `finalNN_dD_lL.bmp`.

**Scores**: `MOS50.mat` and `MOS100.mat` -- the same 690 images rated by
the same 17 observers at 50 cm and at 100 cm. **MOS, higher is better**,
1--5 scale. The 690 entries are ordered as the database's own `RunMetrics.m`
walks the tree (distortion folders alphabetically, files alphabetically
within each), which is what `scripts/make_index.py cidiq` reproduces.

**Viewing geometry -> ppd = 31.612 (50 cm) and 63.223 (100 cm).** These are
the only *measured* geometries in this project: `RunMetrics.m` fixes the
display at 92 dpi and computes `ppd = dpi x distance_in x tan(1 deg)`; eq. 7
of the paper gives the same values to three decimals. The image subtends
about 25 deg at 50 cm and 12.6 deg at 100 cm. Pass `--distance=50|100` to
`make_index.py`; the ppd is written into the index header.

**Why it matters here**: to our knowledge the only public full-reference
*image* database that rates identical stimuli at two viewing distances,
hence the direct test of a viewing-aware metric (Section 3 of the paper).

---

## BC-KU Multi-Screen Dataset (video, three devices)

**Get it**: <https://github.com/NabajeetBarman/Multiscreen-Dataset> (media on
the linked Google Drive folder). Cite Barman, Reznik & Martini, *A subjective
dataset for multi-screen video streaming applications*, QoMEX 2023.

**Contents**: three 1080p60 HEVC sources (`source_sequences/s0N_1080p.mp4`,
5 s), 27 test renditions (`test_sequences/s0N_{L,M,H}_{540p,720p,1080p}.mp4`)
and `Results.xlsx` with MOS per device (columns M, T, TV) from 26 observers.

**Viewing geometry**: 55" 3840x2160 TV at 3H (ppd 113.1); 12.4" 2800x1752
tablet at 18 in (ppd 83.7); 6.1" 2400x1080 phone at 12-14 in (ppd 97.9 at
13 in). Renditions are compared on the 1080p grid after bicubic upscaling,
so the *effective* ppd is divided by the video scale factor: TV 56.55,
tablet 57.4 (video 2800 px wide / 1920), phone 97.9 (1080p shown 1:1).

---

## LIVE Image Quality Assessment Database, Release 2

**Get it**

1. Request the archive password (the form is live; the password arrives by
   email): <https://docs.google.com/forms/d/e/1FAIpQLSfxXFQRRlQDSLNQbefFqTAN2FTwWTDvn00wgNALZ73Np8O8Uw/viewform>
   — linked from <https://live.ece.utexas.edu/research/quality/subjective.htm>
2. Download `databaserelease2.zip` (~700 MB) from
   <https://utexas.box.com/v/databaserelease2> and unlock it.
   The old direct URL under `live.ece.utexas.edu` now 404s; Box is the only
   working route.
3. **Separately** fetch the realigned scores, which are *not* in the zip:
   <https://live.ece.utexas.edu/research/quality/release2/dmos_realigned.mat>

**Contents**: 29 references (`refimgs/`), 779 distorted images across
`jp2k/ jpeg/ wn/ gblur/ fastfading/`, 24-bit BMP, mostly 768×512.

**Scores**: `dmos_realigned.mat` holds `dmos_new` (1×982) and `orgs` (1×982);
`refnames_all.mat` holds the reference filename for each slot. **DMOS —
higher is worse**, nominally 0–100. `orgs[i] == 1` marks a hidden-reference
slot with a meaningless score; filtering those leaves the 779 distorted
images. The 982 slots are one concatenated index in a fixed order:

| slots | folder | count | distorted after filtering |
|---|---|---|---|
| 1–227 | `jp2k` | 227 | 175 |
| 228–460 | `jpeg` | 233 | 169 |
| 461–634 | `wn` | 174 | 145 |
| 635–808 | `gblur` | 174 | 145 |
| 809–982 | `fastfading` | 174 | 145 |

`scripts/make_index.py live` implements exactly this, preferring
`dmos_realigned.mat` and warning if it falls back to `dmos.mat`.

**Viewing geometry → ppd = 30 (default).** The companion paper (Sheikh, Sabir
& Bovik, IEEE TIP 15(11):3440–3451, 2006, §II) states: *"The display monitors
were all 21-inch CRT monitors displaying at a resolution of 1024 × 768
pixels"* and *"Subjects viewed the monitors from an approximate viewing
distance of 2-2.5 screen heights."* No ppd is published. Deriving it from 768
vertical pixels over the full screen height:

- at 2.0 screen heights → 2·atan(0.25) = 28.1° vertical → **≈ 27 px/deg**
- at 2.5 screen heights → 2·atan(0.20) = 22.6° vertical → **≈ 34 px/deg**

so ≈ 27–34 ppd. The default is 30. A 768×512 image then subtends roughly
26° × 17°. **This is a derived figure with about ±15 % uncertainty** — worth a
sensitivity sweep (`--ppd=27`, `--ppd=34`) before trusting any conclusion that
turns on the band placement.

**Licence**: <https://live.ece.utexas.edu/research/quality/copyright2.htm>.
Publications must acknowledge LIVE and CPS at UT Austin and cite the database,
Sheikh/Sabir/Bovik 2006, and Wang et al. 2004. Nine images from Visual
Delights Inc. may not be used outside the database without that company's
permission.

---

## TID2013

**Get it**: <https://www.ponomarenko.info/tid2013.htm> →
`https://www.ponomarenko.info/tid2013/tid2013.rar` (~913 MB, no registration).

**Contents**: 25 references (`reference_images/I01.BMP` … `I25.BMP`), 3000
distorted (`distorted_images/iXX_YY_Z.bmp`, 25 refs × 24 distortion types × 5
levels), 24-bit BMP, 512×384.

**Scores**: `mos_with_names.txt`, one `"<MOS> <filename>"` per line. **MOS —
higher is better**, range 0–9. (`mos.txt` is the same values without names,
`mos_std.txt` the per-image standard deviations.)

**Viewing geometry → ppd = 32 (default, weakly supported).** The TID2013
papers state only that the Ukrainian sessions followed ITU-R BT.500-13 and
that images were displayed at 512×384 so three could be shown side by side.
**201 of 985 sessions were run over the Internet on uncontrolled displays**,
so no single geometry applies. Winkler's survey (IEEE JSTSP 6(6), 2012,
Table II) lists TID as 19″ LCD at *varying* distance. A 19″ 1280×1024 panel
(0.294 mm pitch) at 60 cm would give ≈ 36 ppd; the 32 default is a
conservative round number. **Treat ppd as a free parameter on this database
and report the sensitivity** rather than presenting a single number.

TID2013 spans 24 distortion types, many well outside what SSIM was designed
for (colour saturation changes, mean shift, non-eccentricity pattern noise).
`--types=8,10,11` restricts to blur / JPEG / JPEG2000 if you want a comparison
closer to the original SSIM paper's scope.

**Citation**: Ponomarenko et al., "Image database TID2013: Peculiarities,
results and perspectives", *Signal Processing: Image Communication* 30:57–77,
2015. No explicit licence text is published beyond the request to cite.

---

## CSIQ

**Get it** (easiest of the three — plain unauthenticated GETs, ~371 MB total).
Eric Larson's page at SMU: <https://s2.smu.edu/~eclarson/csiq.html>

- `https://s2.smu.edu/~eclarson/csiq/src_imgs.zip` (14 MB)
- `https://s2.smu.edu/~eclarson/csiq/dst_imgs.zip` (357 MB)
- `https://s2.smu.edu/~eclarson/csiq/csiq.DMOS.xlsx` (208 KB)

The Shizuoka University mirror (<http://vision.eng.shizuoka.ac.jp/mod/page/view.php?id=23>)
carries the same files. The original `vision.okstate.edu` host is gone.

**Contents**: 30 references (`src_imgs/*.png`), 866 distorted
(`dst_imgs/{AWGN,BLUR,contrast,fnoise,JPEG,jpeg2000}/<ref>.<TYPE>.<level>.png`),
8-bit RGB PNG, 512×512. Folder name casing is inconsistent in the
distribution; `make_index.py csiq` matches case-insensitively.

**Scores**: `csiq.DMOS.xlsx`, sheet `all_by_image`, header on row 4, columns
`image`, `dst_type`, `dst_lev`, `dmos`. **DMOS — higher is worse**, range
[0, 1].

**Viewing geometry → ppd = 50 (default, estimated).** The MAD paper (Larson &
Chandler, *J. Electronic Imaging* 19(1):011006, 2010) describes the setup only
qualitatively — four calibrated LCDs side by side at equal viewing distance.
The one numeric figure available is Winkler's survey: LCD at **80 cm**. Pixel
pitch is undocumented, so ppd is not rigorously derivable; at 80 cm a 24″
1920×1200 panel gives ≈ 52 ppd, a 19″ 1280×1024 ≈ 47 ppd. Hence ≈ 50, and a
512×512 image subtending ≈ 10°. **This is an estimate, not a sourced value.**

Note that CSIQ's ppd is roughly 1.6× LIVE's, which matters directly: the same
cpd channel set maps to different pixel scales on the two databases. That is
the point of specifying bands in cpd, but it does mean a single set of
exponents tuned on one database is not automatically right for the other —
worth testing explicitly.

**Citation**: E. C. Larson and D. M. Chandler, "Most apparent distortion:
full-reference image quality assessment and the role of strategy", *J.
Electronic Imaging* 19(1), 011006, 2010. Copyright (c) 2013 Eric Larson;
contact `eclarson@lyle.smu.edu` for usage permissions.

---

## A note on ppd

NB-SSIM's channels are specified in cycles per degree, so ppd is a genuine
model input, not a cosmetic setting. None of the three databases publishes it.
Everything above is either derived from a stated monitor and viewing distance
(LIVE), estimated from a viewing distance plus an assumed pixel pitch (CSIQ),
or essentially unavailable (TID2013).

The honest treatment is to sweep it. `nbsbench --ppd=...` overrides the value
in the index file, so a sweep is a shell loop; if the ranking of NB-SSIM
against SSIM and MS-SSIM is stable across the plausible range, the conclusion
survives the uncertainty. If it is not, that is itself the finding.
