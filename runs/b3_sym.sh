#!/bin/sh
# Three-channel bank with the cardinal-isotropy constraint w_ij = w_ji,
# so the 3x3 sub-band matrix carries six free weights instead of nine.
cd "$(dirname "$0")/.."
B="--bank=loggabor --bands=3 --peaks=1.1662,3.3466,11.3137 --first=bandpass --last=bandpass --bw=4.7,3.6,2.75 --band-norm=peak"
C="$B --transform=dft --taper=1.0 --weight-sym=on --stride=16 --progress --maxit=1 --restarts=1 --folds=0"
./bin/nbstune $C --ppd=31.612 --cache-save=runs/b3s_50cm.bin  --out=/dev/null data/cidiq_50cm.csv  >runs/b3sc50.log  2>&1
echo c50 >> runs/b3s.state
./bin/nbstune $C --ppd=63.223 --cache-save=runs/b3s_100cm.bin --out=/dev/null data/cidiq_100cm.csv >runs/b3sc100.log 2>&1
echo c100 >> runs/b3s.state
./bin/nbsjoint --label=50cm --cond=runs/b3s_50cm.bin --label=100cm --cond=runs/b3s_100cm.bin \
  --combine=mean --objective=plcc --folds=5 --restarts=2 --basin=6 --maxit=3000 \
  --eval-stride=8 --out=runs/w_b3_sym.txt > runs/fit_b3_sym.log 2>&1
echo fit >> runs/b3s.state
