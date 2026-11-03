#!/bin/sh
cd "$(dirname "$0")/.."
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --nbs-win=1 --stride=16 --progress --maxit=1 --restarts=1 --folds=0"
./bin/nbstune $C --ppd=31.612 --cache-save=runs/pw_50cm.bin  --out=/dev/null data/cidiq_50cm.csv  >runs/pw50.log  2>runs/pw50.err
echo pw50 >> runs/pw.state
