#!/bin/sh
cd "$(dirname "$0")/.."
C="--preset=wilson --transform=dft --taper=1.0 --band-decimate=off --stride=16 --progress --maxit=1 --restarts=1 --folds=0"
./bin/nbstune $C --ppd=63.223 --cache-save=runs/cache_100cm_nodec.bin --out=/tmp/wnd1.txt data/cidiq_100cm.csv >runs/cnd100.log 2>runs/cnd100.err
echo nodec100 >> runs/cache.state
./bin/nbstune $C --ppd=31.612 --cache-save=runs/cache_50cm_nodec.bin  --out=/tmp/wnd2.txt data/cidiq_50cm.csv  >runs/cnd50.log  2>runs/cnd50.err
echo nodec50 >> runs/cache.state
