#!/bin/bash
cd "$(dirname "$0")/../.."
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --stride=16 --maxit=1 --restarts=1 --folds=0 --out=/dev/null"
FIT="--combine=mean --objective=plcc --folds=5 --restarts=1 --basin=2 --maxit=3000 --eval-stride=8"
# name:kc-list
for V in "third:0.00333,0.01" "triple:0.03,0.09" "tenth:0.001,0.003" "csf:0.01,0.031,0.030,0.033,0.055,0.147"; do
  N=${V%%:*}; K=${V#*:}; D=runs/constants/$N; mkdir -p $D
  bin/nbstune $C --kc=$K --ppd=31.612 --cache-save=$D/c50.bin  data/cidiq_50cm.csv  > $D/c50.log 2>&1 &
  bin/nbstune $C --kc=$K --ppd=63.223 --cache-save=$D/c100.bin data/cidiq_100cm.csv > $D/c100.log 2>&1 &
  wait
  python3 scripts/norm/normeval.py $D > $D/eval.txt 2>/dev/null
  echo "eval $N" >> runs/constants/progress.log
done
for N in third triple tenth csf; do
  D=runs/constants/$N
  bin/nbsjoint --label=50cm --cond=$D/c50.bin --label=100cm --cond=$D/c100.bin $FIT --out=$D/w.txt > $D/fit.log 2>&1
  echo "fit $N" >> runs/constants/progress.log
done
