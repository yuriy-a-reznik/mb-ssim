#!/bin/sh
cd "$(dirname "$0")/.."
W=results/cidiq/w_cidiq_joint.txt
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --weight-file=$W --metrics=nbssim --progress --resume"
for P in 110 80 55 40 30 25 20 15; do
  ./bin/nbsbench $C --ppd=$P --out=runs/liveU_$P.csv data/live_r2_realigned.csv \
     >runs/liveU2_$P.log 2>runs/liveU2_$P.err
  echo "v $P" >> runs/liveU.state
done
