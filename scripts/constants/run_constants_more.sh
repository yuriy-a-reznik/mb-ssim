#!/bin/bash
cd "$(dirname "$0")/../.."
C="--preset=wilson --transform=dft --taper=1.0 --weight-sym=off --stride=16 --maxit=1 --restarts=1 --folds=0 --out=/dev/null"
FIT="--combine=mean --objective=plcc --folds=5 --restarts=1 --basin=2 --maxit=3000 --eval-stride=8"
for V in "half:0.005,0.015" "fifth:0.002,0.006" "csf_all:0.046,0.031,0.030,0.033,0.055,0.147"; do
  N=${V%%:*}; K=${V#*:}; D=runs/constants/$N; mkdir -p $D
  if [ ! -s $D/eval.txt ]; then
    bin/nbstune $C --kc=$K --ppd=31.612 --cache-save=$D/c50.bin  data/cidiq_50cm.csv  > $D/c50.log 2>&1 &
    bin/nbstune $C --kc=$K --ppd=63.223 --cache-save=$D/c100.bin data/cidiq_100cm.csv > $D/c100.log 2>&1 &
    wait
    python3 scripts/norm/normeval.py $D > $D/eval.txt 2>/dev/null
    echo "eval $N" >> runs/constants/progress.log
  fi
done
for N in tenth csf csf3 half; do
  D=runs/constants/$N
  if ! grep -q "mean held-out" $D/fit.log 2>/dev/null; then
    bin/nbsjoint --label=50cm --cond=$D/c50.bin --label=100cm --cond=$D/c100.bin $FIT --out=$D/w.txt > $D/fit.log 2>&1
    echo "fit $N" >> runs/constants/progress.log
  fi
done
