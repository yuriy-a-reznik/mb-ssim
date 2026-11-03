#!/bin/sh
cd "$(dirname "$0")/.."
W=results/cidiq/w_cidiq_joint.txt
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --weight-file=$W --metrics=nbssim --progress --resume"
./bin/nbsbench $C --ppd=30 --out=runs/liveU_30.csv data/live_r2_realigned.csv >runs/liveU_30.log 2>runs/liveU_30.err
echo "u 30" >> runs/liveU.state
