#!/bin/sh
# Convenience driver: benchmark all three metrics on a dataset, sharded across
# CPUs, then tune the N-band exponents with content-disjoint cross-validation.
#
#   scripts/run_benchmark.sh data/live.csv results/live [extra nbsbench opts...]
#
# Everything after the second argument is passed to both nbsbench and nbstune,
# so e.g. --ppd=27 or --peaks=0.8,1.7,2.8,4,8 applies to both.

set -e

INDEX="$1"
OUTDIR="$2"
[ -n "$INDEX" ] && [ -n "$OUTDIR" ] || {
    echo "usage: $0 INDEX.csv OUTDIR [options...]" >&2
    exit 2
}
shift 2

JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)}"
mkdir -p "$OUTDIR"

echo "== benchmarking with $JOBS shards"
i=0
while [ "$i" -lt "$JOBS" ]; do
    ./bin/nbsbench --shard="$i/$JOBS" --out="$OUTDIR/scores.$i.csv" "$@" "$INDEX" \
        > "$OUTDIR/bench.$i.log" 2>&1 &
    i=$((i + 1))
done
wait

# Merge, keeping a single header.
{
    head -4 "$OUTDIR/scores.0.csv"
    for f in "$OUTDIR"/scores.*.csv; do
        tail -n +5 "$f"
    done
} > "$OUTDIR/scores.csv"

echo
echo "== combined result"
./bin/nbsbench --report="$OUTDIR/scores.csv" | tee "$OUTDIR/report.txt"

echo
echo "== tuning exponents (this is the slow part)"
./bin/nbstune --folds=5 --out="$OUTDIR/weights.txt" "$@" "$INDEX" \
    2>&1 | tee "$OUTDIR/tune.txt"

echo
echo "== re-benchmarking with the fitted exponents"
i=0
while [ "$i" -lt "$JOBS" ]; do
    ./bin/nbsbench --shard="$i/$JOBS" --weight-file="$OUTDIR/weights.txt" \
        --out="$OUTDIR/scores_tuned.$i.csv" "$@" "$INDEX" \
        > "$OUTDIR/bench_tuned.$i.log" 2>&1 &
    i=$((i + 1))
done
wait
{
    head -4 "$OUTDIR/scores_tuned.0.csv"
    for f in "$OUTDIR"/scores_tuned.*.csv; do
        tail -n +5 "$f"
    done
} > "$OUTDIR/scores_tuned.csv"
./bin/nbsbench --report="$OUTDIR/scores_tuned.csv" | tee "$OUTDIR/report_tuned.txt"

echo
echo "NOTE: the numbers above are an all-data fit.  The held-out figures in"
echo "      $OUTDIR/tune.txt are the ones to quote."
