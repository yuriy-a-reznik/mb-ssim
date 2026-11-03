#!/bin/sh
# Window ablation: windowed xi against pointwise xi, identical search
# settings so the two rows are comparable.
cd "$(dirname "$0")/.."
A="--combine=mean --objective=plcc --folds=3 --restarts=2 --basin=2 --maxit=1500 --eval-stride=8"
./bin/nbsjoint --label=50cm --cond=runs/cache_50cm.bin --label=100cm --cond=runs/cache_100cm.bin \
   $A > runs/abl_win.log 2>&1
echo win >> runs/abl.state
./bin/nbsjoint --label=50cm --cond=runs/pw_50cm.bin --label=100cm --cond=runs/pw_100cm.bin \
   $A > runs/abl_pw.log 2>&1
echo pw >> runs/abl.state
