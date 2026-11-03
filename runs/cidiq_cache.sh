#!/bin/sh
set -e
cd "$(dirname "$0")/.."
C="--preset=wilson --transform=dft --taper=1.0 --stride=16 --progress --maxit=1 --restarts=1 --folds=0"
[ -f runs/cache_50cm.bin ]  || ./bin/nbstune $C --ppd=31.612 --cache-save=runs/cache_50cm.bin  --out=/tmp/w50.txt  data/cidiq_50cm.csv  >runs/cache50.log  2>runs/cache50.err &
P1=$!
[ -f runs/cache_100cm.bin ] || ./bin/nbstune $C --ppd=63.223 --cache-save=runs/cache_100cm.bin --out=/tmp/w100.txt data/cidiq_100cm.csv >runs/cache100.log 2>runs/cache100.err &
P2=$!
wait $P1 $P2
