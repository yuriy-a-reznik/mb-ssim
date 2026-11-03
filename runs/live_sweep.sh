#!/bin/sh
# ppd sweep on LIVE Release 2, which publishes no viewing geometry.
cd "$(dirname "$0")/.."
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --nbs-win=1 --stride=12 --subset=2 --progress --maxit=1 --restarts=1 --folds=0"
for P in 15 20 25 30 40 55; do
  [ -f runs/live_p$P.bin ] && continue
  ./bin/nbstune $C --ppd=$P --cache-save=runs/live_p$P.bin --out=/dev/null \
      data/live_r2_realigned.csv >runs/live_c$P.log 2>runs/live_c$P.err
  echo "cache $P" >> runs/live.state
done
for P in 15 20 25 30 40 55; do
  ./bin/nbsjoint --label=ppd$P --cond=runs/live_p$P.bin --combine=mean --objective=plcc \
      --folds=0 --restarts=2 --basin=3 --maxit=2000 --eval-stride=4 \
      --out=runs/live_w$P.txt > runs/live_f$P.log 2>&1
  echo "fit $P" >> runs/live.state
done
echo sweepdone >> runs/live.state
