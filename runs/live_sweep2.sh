#!/bin/sh
cd "$(dirname "$0")/.."
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --stride=12 --subset=2 --progress --maxit=1 --restarts=1 --folds=0"
for P in 80 110; do
  [ -f runs/live_p$P.bin ] && continue
  ./bin/nbstune $C --ppd=$P --cache-save=runs/live_p$P.bin --out=/dev/null \
      data/live_r2_realigned.csv >runs/live_c$P.log 2>runs/live_c$P.err
  echo "cache $P" >> runs/live2.state
done
# held-out sweep: cheap settings, three folds split by reference image
for P in 15 20 25 30 40 55 80 110; do
  ./bin/nbsjoint --label=ppd$P --cond=runs/live_p$P.bin --combine=mean --objective=plcc \
      --folds=3 --restarts=2 --basin=2 --maxit=1500 --eval-stride=8 \
      > runs/live_cv$P.log 2>&1
  echo "cv $P" >> runs/live2.state
done
echo sweep2done >> runs/live2.state
