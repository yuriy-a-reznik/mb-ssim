#!/bin/sh
# Band-count sweep on CID:IQ.  First and last channel peaks are Wilson's
# (0.8 and 16 c/deg); the intermediate peaks are spaced log-uniformly
# between them, and each bandwidth is Wilson's value interpolated at the new
# peak, scaled by 6/N so that the bank still covers the range.
cd "$(dirname "$0")/.."

# ---------------- N = 3 ----------------
B3="--bank=loggabor --bands=3 --peaks=0.8000,3.5777,16.0000 --first=bandpass --last=bandpass --bw=5.0000,3.5251,2.5000 --band-norm=peak"
C3="$B3 --transform=dft --taper=1.0 --weight-sym=off --stride=16 --progress --maxit=1 --restarts=1 --folds=0"
./bin/nbstune $C3 --ppd=31.612 --cache-save=runs/n3_50cm.bin  --out=/dev/null data/cidiq_50cm.csv  >runs/n3_c50.log  2>&1
echo n3_c50 >> runs/sweep.state
./bin/nbstune $C3 --ppd=63.223 --cache-save=runs/n3_100cm.bin --out=/dev/null data/cidiq_100cm.csv >runs/n3_c100.log 2>&1
echo n3_c100 >> runs/sweep.state
./bin/nbsjoint --label=50cm --cond=runs/n3_50cm.bin --label=100cm --cond=runs/n3_100cm.bin \
  --combine=mean --objective=plcc --folds=5 --restarts=2 --basin=6 --maxit=3000 \
  --eval-stride=8 --out=runs/w_n3.txt > runs/fit_n3.log 2>&1
echo n3_fit >> runs/sweep.state
rm -f runs/n3_50cm.bin runs/n3_100cm.bin

# ---------------- N = 4 ----------------
B4="--bank=loggabor --bands=4 --peaks=0.8000,2.1715,5.8945,16.0000 --first=bandpass --last=bandpass --bw=3.7500,3.0792,2.3822,1.8750 --band-norm=peak"
C4="$B4 --transform=dft --taper=1.0 --weight-sym=off --stride=16 --progress --maxit=1 --restarts=1 --folds=0"
./bin/nbstune $C4 --ppd=31.612 --cache-save=runs/n4_50cm.bin  --out=/dev/null data/cidiq_50cm.csv  >runs/n4_c50.log  2>&1
echo n4_c50 >> runs/sweep.state
./bin/nbstune $C4 --ppd=63.223 --cache-save=runs/n4_100cm.bin --out=/dev/null data/cidiq_100cm.csv >runs/n4_c100.log 2>&1
echo n4_c100 >> runs/sweep.state
./bin/nbsjoint --label=50cm --cond=runs/n4_50cm.bin --label=100cm --cond=runs/n4_100cm.bin \
  --combine=mean --objective=plcc --folds=5 --restarts=2 --basin=6 --maxit=3000 \
  --eval-stride=8 --out=runs/w_n4.txt > runs/fit_n4.log 2>&1
echo n4_fit >> runs/sweep.state
rm -f runs/n4_50cm.bin runs/n4_100cm.bin

# ---------------- N = 5 ----------------
B5="--bank=loggabor --bands=5 --peaks=0.8000,1.6918,3.5777,7.5659,16.0000 --first=bandpass --last=bandpass --bw=3.0000,2.6423,2.1151,1.8193,1.5000 --band-norm=peak"
C5="$B5 --transform=dft --taper=1.0 --weight-sym=off --stride=16 --progress --maxit=1 --restarts=1 --folds=0"
./bin/nbstune $C5 --ppd=31.612 --cache-save=runs/n5_50cm.bin  --out=/dev/null data/cidiq_50cm.csv  >runs/n5_c50.log  2>&1
echo n5_c50 >> runs/sweep.state
./bin/nbstune $C5 --ppd=63.223 --cache-save=runs/n5_100cm.bin --out=/dev/null data/cidiq_100cm.csv >runs/n5_c100.log 2>&1
echo n5_c100 >> runs/sweep.state
./bin/nbsjoint --label=50cm --cond=runs/n5_50cm.bin --label=100cm --cond=runs/n5_100cm.bin \
  --combine=mean --objective=plcc --folds=5 --restarts=2 --basin=6 --maxit=3000 \
  --eval-stride=8 --out=runs/w_n5.txt > runs/fit_n5.log 2>&1
echo n5_fit >> runs/sweep.state
rm -f runs/n5_50cm.bin runs/n5_100cm.bin

# ---------------- N = 6 ----------------
B6="--bank=loggabor --bands=6 --peaks=0.8000,1.4565,2.6516,4.8273,8.7885,16.0000 --first=bandpass --last=bandpass --bw=2.5000,2.2615,1.9327,1.6458,1.4661,1.2500 --band-norm=peak"
C6="$B6 --transform=dft --taper=1.0 --weight-sym=off --stride=16 --progress --maxit=1 --restarts=1 --folds=0"
./bin/nbstune $C6 --ppd=31.612 --cache-save=runs/n6_50cm.bin  --out=/dev/null data/cidiq_50cm.csv  >runs/n6_c50.log  2>&1
echo n6_c50 >> runs/sweep.state
./bin/nbstune $C6 --ppd=63.223 --cache-save=runs/n6_100cm.bin --out=/dev/null data/cidiq_100cm.csv >runs/n6_c100.log 2>&1
echo n6_c100 >> runs/sweep.state
./bin/nbsjoint --label=50cm --cond=runs/n6_50cm.bin --label=100cm --cond=runs/n6_100cm.bin \
  --combine=mean --objective=plcc --folds=5 --restarts=2 --basin=6 --maxit=3000 \
  --eval-stride=8 --out=runs/w_n6.txt > runs/fit_n6.log 2>&1
echo n6_fit >> runs/sweep.state
rm -f runs/n6_50cm.bin runs/n6_100cm.bin

