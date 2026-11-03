#!/bin/sh
cd "$(dirname "$0")/.."
W=runs/w_free_mean.txt
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --weight-file=$W --metrics=nbssim --progress --resume"
./bin/nbsbench $C --ppd=63.223 --out=runs/final_100cm.csv data/cidiq_100cm.csv >runs/bf100.log 2>runs/bf100.err
echo bf100 >> runs/bench.state
