#!/bin/sh
cd "$(dirname "$0")/.."
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --stride=16 --progress --maxit=1 --restarts=1 --folds=0"
./bin/nbstune $C --ppd=31.612 --cache-save=runs/cache_50cm.bin  --out=/tmp/w50.txt  data/cidiq_50cm.csv  >runs/cache50.log  2>runs/cache50.err
echo done50 >> runs/cache.state
