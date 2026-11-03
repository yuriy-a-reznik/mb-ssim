#!/bin/sh
cd "$(dirname "$0")/.."
./bin/nbsbench --metrics=ssim,msssim --progress --out=runs/live_base.csv \
    data/live_r2_realigned.csv > runs/live_base.log 2>runs/live_base.err
echo base >> runs/live2.state
