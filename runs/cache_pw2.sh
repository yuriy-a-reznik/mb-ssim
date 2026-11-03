#!/bin/sh
cd "$(dirname "$0")/.."
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --nbs-win=1 --stride=16 --progress --maxit=1 --restarts=1 --folds=0"
./bin/nbstune $C --ppd=63.223 --cache-save=runs/pw_100cm.bin --out=/dev/null data/cidiq_100cm.csv >runs/pw100.log 2>runs/pw100.err
echo pw100 >> runs/pw.state
