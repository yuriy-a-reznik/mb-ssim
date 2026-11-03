#!/bin/sh
cd "$(dirname "$0")/.."
B="--bank=loggabor --bands=3 --peaks=1.1662,3.3466,11.3137 --first=bandpass --last=bandpass --bw=4.7,3.6,2.75 --band-norm=peak"
C="$B --transform=dft --taper=1.0 --weight-sym=off --stride=16 --progress --maxit=1 --restarts=1 --folds=0"
./bin/nbstune $C --ppd=31.612 --cache-save=runs/b3_50cm.bin  --out=/dev/null data/cidiq_50cm.csv  >runs/b3c50.log  2>runs/b3c50.err
echo c50 >> runs/b3.state
