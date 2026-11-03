#!/bin/sh
cd "$(dirname "$0")/.."
./bin/nbsjoint --label=50cm --cond=runs/cache_50cm.bin --label=100cm --cond=runs/cache_100cm.bin \
  --combine=pooled --objective=plcc --folds=5 --restarts=2 --basin=6 --maxit=3000 \
  --eval-stride=8 --out=runs/w_free_pooled.txt > runs/fit_pooled.log 2>&1
echo pooled >> runs/fit.state
