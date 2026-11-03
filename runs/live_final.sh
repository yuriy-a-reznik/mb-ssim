#!/bin/sh
cd "$(dirname "$0")/.."
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --stride=12 --progress --maxit=1 --restarts=1 --folds=0"
[ -f runs/live_full80.bin ] || ./bin/nbstune $C --ppd=80 \
    --cache-save=runs/live_full80.bin --out=/dev/null \
    data/live_r2_realigned.csv >runs/live_cfull.log 2>runs/live_cfull.err
echo cachefull >> runs/live3.state
./bin/nbsjoint --label=LIVE --cond=runs/live_full80.bin --combine=mean --objective=plcc \
    --folds=5 --restarts=2 --basin=6 --maxit=3000 --eval-stride=8 \
    --out=runs/w_live80.txt > runs/live_final.log 2>&1
echo fitfull >> runs/live3.state
