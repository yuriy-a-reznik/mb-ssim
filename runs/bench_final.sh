#!/bin/sh
cd "$(dirname "$0")/.."
W=runs/w_joint_pooled.txt
C="--preset=wilson --transform=dft --taper=1.0 --weight-file=$W --metrics=nbssim --progress --resume"
./bin/nbsbench $C --ppd=31.612 --out=runs/final_50cm.csv  data/cidiq_50cm.csv  >runs/bf50.log  2>runs/bf50.err
echo bf50 >> runs/bench.state
./bin/nbsbench $C --ppd=63.223 --out=runs/final_100cm.csv data/cidiq_100cm.csv >runs/bf100.log 2>runs/bf100.err
echo bf100 >> runs/bench.state
