#!/bin/sh
cd "$(dirname "$0")/.."
B="--bank=loggabor --bands=3 --peaks=1.1662,3.3466,11.3137 --first=bandpass --last=bandpass --bw=4.7,3.6,2.75 --band-norm=peak"
C="$B --transform=dft --taper=1.0 --weight-sym=off --stride=16 --progress --maxit=1 --restarts=1 --folds=0"

./bin/nbstune $C --ppd=31.612 --cache-save=runs/b3_50cm.bin  --out=/dev/null data/cidiq_50cm.csv  >runs/b3c50.log  2>runs/b3c50.err
echo c50 >> runs/b3.state
./bin/nbstune $C --ppd=63.223 --cache-save=runs/b3_100cm.bin --out=/dev/null data/cidiq_100cm.csv >runs/b3c100.log 2>runs/b3c100.err
echo c100 >> runs/b3.state

# identical protocol to the 6-band headline fit (fit_win.sh)
./bin/nbsjoint --label=50cm --cond=runs/b3_50cm.bin --label=100cm --cond=runs/b3_100cm.bin \
  --combine=mean --objective=plcc --folds=5 --restarts=2 --basin=6 --maxit=3000 \
  --eval-stride=8 --out=runs/w_b3_mean.txt > runs/fit_b3.log 2>&1
echo fit >> runs/b3.state

# heavier search, to show the 3-band result is not search-limited
./bin/nbsjoint --label=50cm --cond=runs/b3_50cm.bin --label=100cm --cond=runs/b3_100cm.bin \
  --combine=mean --objective=plcc --folds=5 --restarts=5 --basin=14 --maxit=4000 \
  --eval-stride=8 --out=runs/w_b3_heavy.txt > runs/fit_b3_heavy.log 2>&1
echo fitheavy >> runs/b3.state
