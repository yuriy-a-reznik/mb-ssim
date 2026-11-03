#!/bin/sh
cd "$(dirname "$0")/.."
W=runs/w_free_mean.txt
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --weight-file=$W --metrics=nbssim --progress --resume"
./bin/nbsbench $C --ppd=31.612 --out=runs/final_50cm.csv  data/cidiq_50cm.csv  >runs/bf50.log  2>runs/bf50.err
echo bf50 >> runs/bench.state
