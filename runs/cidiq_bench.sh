#!/bin/sh
# CID:IQ benchmark, Wilson bank, DFT core, both viewing distances.
set -e
cd "$(dirname "$0")/.."
COMMON="--preset=wilson --transform=dft --taper=1.0 --progress --resume"
./bin/nbsbench $COMMON --ppd=31.612 --out=runs/cidiq_50cm.csv  data/cidiq_50cm.csv  > runs/cidiq_50cm.log  2>runs/cidiq_50cm.err &
P1=$!
./bin/nbsbench $COMMON --ppd=63.223 --out=runs/cidiq_100cm.csv data/cidiq_100cm.csv > runs/cidiq_100cm.log 2>runs/cidiq_100cm.err &
P2=$!
wait $P1 $P2
