#!/bin/sh
cd "$(dirname "$0")/.."
./bin/nbsjoint --label=50cm --cond=runs/pw_50cm.bin --label=100cm --cond=runs/pw_100cm.bin \
  --combine=mean --objective=plcc --folds=5 --restarts=2 --basin=6 --maxit=3000 \
  --eval-stride=8 --out=runs/w_pw_mean.txt > runs/fit_pw.log 2>&1
echo pw >> runs/fit.state
