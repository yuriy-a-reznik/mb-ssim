#!/bin/bash
# caches and fast-budget joint fits for the normaliser experiment
cd "$(dirname "$0")/../.."
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --stride=16 --maxit=1 --restarts=1 --folds=0 --out=/dev/null"
FIT="--combine=mean --objective=plcc --folds=5 --restarts=1 --basin=2 --maxit=3000 --eval-stride=8"
for P in geo max 0.5 2; do
  D=runs/norm/p_$P; mkdir -p $D
  bin/nbstune $C --xi-norm=$P --ppd=31.612 --cache-save=$D/c50.bin  data/cidiq_50cm.csv  > $D/c50.log 2>&1 &
  bin/nbstune $C --xi-norm=$P --ppd=63.223 --cache-save=$D/c100.bin data/cidiq_100cm.csv > $D/c100.log 2>&1 &
  wait
  bin/nbsjoint --label=50cm --cond=$D/c50.bin --label=100cm --cond=$D/c100.bin $FIT --out=$D/w.txt > $D/fit.log 2>&1
  echo "done $P" >> runs/norm/progress.log
done
