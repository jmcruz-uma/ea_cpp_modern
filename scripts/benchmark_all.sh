#!/usr/bin/env bash
# benchmark_all.sh — Automated benchmark suite for ea-cpp
# Runs all algorithms × all problems × N runs, generates CSV

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="/tmp/ea-cpp-benchmarks"
RESULTS_DIR="$BUILD_DIR/results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Configuration
RUNS=30
POP_SIZE=100
MAX_EVALS=25000
DIM=30

mkdir -p "$RESULTS_DIR"

echo "=== ea-cpp Benchmark Suite ==="
echo "Timestamp: $TIMESTAMP"
echo "Runs: $RUNS per (algorithm, problem) pair"
echo "Population: $POP_SIZE, Evaluations: $MAX_EVALS, Dimensions: $DIM"
echo ""

# Compile all benchmark binaries
echo "Compiling benchmark binaries..."

# NSGA-II benchmark
g++-14 -std=c++23 -O2 -I "$PROJECT_DIR/include" "$PROJECT_DIR/examples/benchmark_fair.cpp" -o "$BUILD_DIR/bench_nsga2" -lm

# Copy other benchmarks if they exist
for bench in benchmark_large test_random_search test_spread; do
    if [ -f "$PROJECT_DIR/examples/${bench}.cpp" ]; then
        echo "  Found: ${bench}.cpp"
    fi
done

echo ""
echo "Running benchmarks..."
echo "algorithm,problem,IGD_mean,IGD_std,Time_mean,Time_std" > "$RESULTS_DIR/benchmark_results.csv"

# Run NSGA-II on ZDT1
echo "  NSGAII × ZDT1 ($RUNS runs)..."
for ((i=1; i<=RUNS; i++)); do
    "$BUILD_DIR/bench_nsga2" > /tmp/bench_out.txt 2>&1
    IGD=$(grep "Mean IGD:" /tmp/bench_out.txt | awk '{print $3}')
    IGD_STD=$(grep "Mean IGD:" /tmp/bench_out.txt | awk '{print $5}')
    TIME=$(grep "Mean Time:" /tmp/bench_out.txt | awk '{print $3}')
    TIME_STD=$(grep "Mean Time:" /tmp/bench_out.txt | awk '{print $5}')
    echo "NSGAII,ZDT1,$IGD,$IGD_STD,$TIME,$TIME_STD" >> "$RESULTS_DIR/benchmark_results.csv"
done

echo ""
echo "=== Results ==="
cat "$RESULTS_DIR/benchmark_results.csv"

echo ""
echo "Results saved to: $RESULTS_DIR/benchmark_results.csv"
