#!/bin/sh
# LIVE Release 2, MB-SSIM with the CID:IQ exponents applied verbatim.
# Nothing is fitted on LIVE except the single scalar ppd.
cd "$(dirname "$0")/.."
W=results/cidiq/w_cidiq_joint.txt
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --weight-file=$W --metrics=nbssim --progress --resume"
for P in 15 20 25 30 40 55 80 110; do
  ./bin/nbsbench $C --ppd=$P --out=runs/liveU_$P.csv data/live_r2_realigned.csv \
     >runs/liveU_$P.log 2>runs/liveU_$P.err
  echo "u $P" >> runs/liveU.state
done
echo done >> runs/liveU.state
