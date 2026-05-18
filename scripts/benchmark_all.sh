#!/bin/bash
# benchmark_all.sh — Full benchmark suite for ea-cpp
# Runs all algorithms × all problems × N runs → CSV
# Usage: ./scripts/benchmark_all.sh [N_RUNS]
set -euo pipefail

N_RUNS=${1:-30}
BUILD_DIR="build"
BENCHMARK="$BUILD_DIR/benchmark_fair"
RESULTS_DIR="benchmarks/results"
CSV="$RESULTS_DIR/benchmark_$(date +%Y%m%d_%H%M%S).csv"

mkdir -p "$RESULTS_DIR"

if [ ! -f "$BENCHMARK" ]; then
    echo "Building benchmark binary..."
    make "$BENCHMARK"
fi

echo "=== ea-cpp Benchmark Suite ==="
echo "Runs per config: $N_RUNS"
echo "Output: $CSV"
echo ""

# CSV header
echo "algorithm,problem,pop_size,max_evals,run,igd,gd,spread,time_ms" > "$CSV"

# Algorithms and problems to benchmark
ALGORITHMS="NSGAII RVEA IBEA RandomSearch"
PROBLEMS_2OBJ="ZDT1 ZDT2 ZDT3 ZDT6"
PROBLEMS_3OBJ="DTLZ1 DTLZ2"

for algo in $ALGORITHMS; do
    for prob in $PROBLEMS_2OBJ; do
        echo "Running $algo on $prob ($N_RUNS runs)..."
        for run in $(seq 1 $N_RUNS); do
            $BENCHMARK --algorithm "$algo" --problem "$prob" --runs 1 --csv >> "$CSV" 2>/dev/null || true
        done
    done
done

for algo in $ALGORITHMS; do
    for prob in $PROBLEMS_3OBJ; do
        echo "Running $algo on $prob ($N_RUNS runs)..."
        for run in $(seq 1 $N_RUNS); do
            $BENCHMARK --algorithm "$algo" --problem "$prob" --runs 1 --csv >> "$CSV" 2>/dev/null || true
        done
    done
done

echo ""
echo "=== Results ==="
echo "CSV: $CSV"
echo ""
# Print summary
echo "Algorithm × Problem Summary (mean IGD ± std):"
tail -n +2 "$CSV" | awk -F, '{
    key=$1","$2
    sum[key]+=$5
    sumsq[key]+=$5*$5
    count[key]++
}
END {
    for (k in count) {
        mean=sum[k]/count[k]
        std=sqrt(sumsq[k]/count[k]-mean*mean)
        printf "%-12s %-8s IGD=%.6f ± %.6f (n=%d)\n", \
            substr(k,1,index(k,",")-1), substr(k,index(k,",")+1), mean, std, count[k]
    }
}' | sort