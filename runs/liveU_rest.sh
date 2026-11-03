#!/bin/sh
cd "$(dirname "$0")/.."
W=results/cidiq/w_cidiq_joint.txt
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --weight-file=$W --metrics=nbssim --progress --resume"
for P in 25 40; do
  ./bin/nbsbench $C --ppd=$P --out=runs/liveU_$P.csv data/live_r2_realigned.csv >runs/liveU_$P.log 2>runs/liveU_$P.err
  echo "u $P" >> runs/liveU.state
done
