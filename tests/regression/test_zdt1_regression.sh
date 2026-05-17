#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

echo "=== Compiling benchmark_fair.cpp ==="
g++-14 -std=c++23 -O2 -I include examples/benchmark_fair.cpp -o /tmp/benchmark_fair -lm

echo "=== Running benchmark (this may take a moment) ==="
OUTPUT=$(/tmp/benchmark_fair)

echo "$OUTPUT"

MEAN_IGD=$(echo "$OUTPUT" | grep "Mean IGD:" | awk '{print $3}')

if [ -z "$MEAN_IGD" ]; then
    echo "REGRESSION FAIL: Could not extract Mean IGD value"
    exit 1
fi

# Compare using bc
PASS=$(echo "$MEAN_IGD < 0.01" | bc -l)

if [ "$PASS" -eq 1 ]; then
    echo "REGRESSION PASS: IGD = $MEAN_IGD"
    exit 0
else
    echo "REGRESSION FAIL: IGD = $MEAN_IGD, expected < 0.01"
    exit 1
fi
