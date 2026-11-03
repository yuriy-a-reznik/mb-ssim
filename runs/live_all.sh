#!/bin/sh
# LIVE Release 2 ppd sweep.  The database publishes no viewing geometry, so
# ppd is the one fitted parameter of the setup.  Windowed xi, matching the
# main experiments.
cd "$(dirname "$0")/.."
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --stride=12 --subset=2 --progress --maxit=1 --restarts=1 --folds=0"
for P in 15 20 25 30 40 55 80 110; do
  if [ ! -f runs/live_p$P.bin ]; then
    ./bin/nbstune $C --ppd=$P --cache-save=runs/live_p$P.bin --out=/dev/null \
        data/live_r2_realigned.csv >runs/live_c$P.log 2>runs/live_c$P.err
  fi
  echo "cache $P" >> runs/live2.state
done
for P in 15 20 25 30 40 55 80 110; do
  ./bin/nbsjoint --label=ppd$P --cond=runs/live_p$P.bin --combine=mean --objective=plcc \
      --folds=3 --restarts=2 --basin=2 --maxit=1500 --eval-stride=8 \
      > runs/live_cv$P.log 2>&1
  echo "cv $P" >> runs/live2.state
done
echo sweepdone >> runs/live2.state
