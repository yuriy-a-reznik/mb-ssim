#!/bin/sh
cd "$(dirname "$0")/.."
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --stride=16 --progress --maxit=1 --restarts=1 --folds=0"
./bin/nbstune $C --ppd=63.223 --cache-save=runs/cache_100cm.bin --out=/tmp/w100.txt data/cidiq_100cm.csv >runs/cache100.log 2>runs/cache100.err
echo done100 >> runs/cache.state
