#!/bin/sh
cd "$(dirname "$0")/.."
./bin/nbsjoint --label=50cm --cond=runs/cache_50cm.bin --label=100cm --cond=runs/cache_100cm.bin \
  --combine=mean --objective=plcc --folds=0 --restarts=3 --basin=4 --maxit=1500 \
  --eval-stride=8 --per-cond --out=runs/w_joint_mean.txt > runs/joint_final.log 2>&1
echo finalfit >> runs/joint.state
