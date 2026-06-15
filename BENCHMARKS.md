# Benchmarks

Reproducible guide for all performance comparisons. Each section documents
prerequisites, exact commands, and how to interpret the results.

---

## Setup

### C++ compiler

The scripts default to `g++-14`. Override with the `CXX_BENCH` environment variable:
```bash
export CXX_BENCH=g++-14   # default
export CXX_BENCH=g++-16   # for C++26 migration experiments
```

### jMetal 7.4 (multi-objective comparison)

Clone and build once:
```bash
git clone https://github.com/jMetal/jMetal.git $HOME/jMetal
cd $HOME/jMetal && mvn clean install -DskipTests -q
```

The scripts default to `JMETAL_ROOT=$HOME/jMetal`. Override if you cloned elsewhere:
```bash
export JMETAL_ROOT=/path/to/your/jMetal
```

### ea-cpp-original (single-objective architectural comparison)

Clone and build once:
```bash
git clone https://github.com/Bio4Res/ea.git $HOME/ea_cpp_original
cd $HOME/ea_cpp_original
${CXX_BENCH:-g++-14} -O3 -DNDEBUG -I. ea_main.cpp -std=c++23 -march=native -o ea_main_gcc -lm
```

The scripts default to `ORIGINAL_ROOT=$HOME/ea_cpp_original`. Override if needed:
```bash
export ORIGINAL_ROOT=/path/to/your/ea_cpp_original
```

### Bare-metal Linux (recommended for publication results)

```bash
sudo cpupower frequency-set -g performance
echo 1 | sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo
```

Verify with `turbostat` or `watch sensors`. WSL2 results are indicative; bare-metal
results are the definitive numbers for the paper.

### Python (statistical analysis)

```bash
pip install scipy matplotlib
```

---

## Multi-objective comparison: ea-cpp-modern vs jMetal 7.4

JVM flags used: `-server -Xms512m -Xmx2g -XX:+UseG1GC -XX:+AlwaysPreTouch
-XX:+DisableExplicitGC -Djava.awt.headless=true`

### Summary (30 runs, WSL2, GCC 14)

| Algorithm | Median speedup | IGD parity |
|-----------|:--------------:|:----------:|
| NSGA-II   | 3.57×          | ✅ all ns   |
| SMPSO     | 4.28×          | ✅ all ns   |
| SPEA2     | 2.17×          | ✅ all ns   |
| IBEA      | 12.41×         | ✅ all ns   |
| MOEA/D-DE | 7.87×          | ✅ all ns   |
| NSGA-III  | 2.93×          | ✅ all ns   |
| SMS-EMOA  | 3.92×          | ✅ all ns   |
| GDE3      | 13.05×         | ✅ all ns   |
| AGE-MOEA  | 5.02×          | ✅ 3/4 ns   |
| MO-Cell   | 5.12×          | ✅ 3/4 ns   |
| PAES      | 4.09×          | ✅ all ns   |

### ZDT family — NSGA-II (pop=100, max\_evals=25000, d=30 except ZDT4 d=10)

| Problem | C++ median (ms) | Java median (ms) | Speedup | IGD p |
|---------|:---:|:---:|:---:|:---:|
| ZDT1 | 109 | 483 | 4.42× | 0.22 ns |
| ZDT2 | 110 | 630 | 5.71× | 0.82 ns |
| ZDT3 |  55 | 185 | 3.40× | 0.46 ns |
| ZDT4 |  79 | 283 | 3.60× | 0.27 ns |

### WFG family — NSGA-II (pop=100, max\_evals=25000, k=2 l=4 M=2 d=6)

Reference fronts: pre-computed CSVs from jMetal (~119–2600 points per problem).

| Problem | C++ median (ms) | Java median (ms) | Speedup | IGD p |
|---------|:---:|:---:|:---:|:---:|
| WFG1 | 125 | 417 | 3.33× | 0.009 **  |
| WFG2 | 123 | 392 | 3.19× | 0.70 ns   |
| WFG3 | 108 | 398 | 3.69× | 0.36 ns   |
| WFG4 | 109 | 413 | 3.80× | 0.61 ns   |
| WFG5 | 104 | 416 | 4.01× | 0.21 ns   |
| WFG6 | 113 | 412 | 3.65× | 0.82 ns   |
| WFG7 | 111 | 416 | 3.76× | 0.43 ns   |
| WFG8 | 134 | 443 | 3.30× | <0.0001 ***|
| WFG9 | 123 | 434 | 3.52× | 0.18 ns   |
| **median** | | | **3.65×** | 7/9 ns |

**WFG8 note** (p<0.0001, C++ IGD better — 0.060 vs 0.125): `bParam` over distance variables
uses `rSum` of previous variables. Java uses `float` internally; C++ uses `double`. The
precision difference changes the exponent (range 0.02–50), altering the fitness landscape.
C++ double gives NSGA-II a more differentiated surface → better convergence.

**WFG1 note** (p=0.009, Java IGD slightly better): WFG1's `bPoly(y, 0.02)` across all
variables creates a nearly flat landscape with high variance (σ=0.19 C++, σ=0.33 Java).
The difference is likely a seed-coverage effect; both systems are in the same IGD order
of magnitude.

### UF family — NSGA-II (pop=100, max\_evals=25000, d=30)

| Problem | C++ median (ms) | Java median (ms) | Speedup | IGD p |
|---------|:---:|:---:|:---:|:---:|
| UF1 | 161 | 510 | 3.16× | 0.57 ns |
| UF2 | 165 | 601 | 3.65× | 0.89 ns |
| UF3 | 185 | 579 | 3.14× | 0.33 ns |
| UF4 | 152 | 560 | 3.67× | 0.97 ns |
| UF5 | 181 | 616 | 3.41× | 0.58 ns |
| UF6 | 186 | 623 | 3.35× | 1.00 ns |
| UF7 | 164 | 590 | 3.61× | 0.58 ns |
| **mean** | | | **3.43×** | all ns |

**UF6 note**: C++ σ=98ms due to a single outlier (run 4, seed 45) in which NSGA-II found
both segments of the discontinuous Pareto front → more non-dominated solutions → more
expensive NDS pass. Without that run, σ=7ms. The outlier produces the best IGD — correct
algorithmic behavior.

### LZ09 family — NSGA-II / NSGA-III (pop=100/92, max\_evals=25000, d=30)

LZ09F1–F5, F7–F9: NSGA-II (pop=100).
LZ09F6: NSGA-III (divisions=12, pop=92) — 3-objective problem.

| Problem | C++ median (ms) | Java median (ms) | Speedup | IGD p |
|---------|:---:|:---:|:---:|:---:|
| LZ09F1 | 140 | 648 | 4.61× | 0.70 ns |
| LZ09F2 | 172 | 712 | 4.13× | 0.45 ns |
| LZ09F3 | 166 | 657 | 3.96× | 0.64 ns |
| LZ09F4 | 166 | 549 | 3.31× | 0.13 ns |
| LZ09F5 | 177 | 563 | 3.19× | 0.18 ns |
| LZ09F6 | 198 | 523 | 2.64× | 0.031 * |
| LZ09F7 | 150 | 492 | 3.29× | 0.84 ns |
| LZ09F8 | 161 | 491 | 3.04× | 0.98 ns |
| LZ09F9 | 191 | 564 | 2.96× | 0.92 ns |
| **mean** | | | **3.46×** | 8/9 ns |

**LZ09F6 note** (p=0.031, Java NSGA-III slightly better — 0.118 vs 0.140): high variance
in both (min ~0.10, max ~0.54); absolute difference is 0.022. Likely due to jMetal's more
mature niching/normalization for 3-objective problems.

### Reproducing the comparison

Build jMetal once, then run each script with `SKIP_JMETAL_BUILD=1`:
```bash
./scripts/run_nsga2_comparison.sh    results/nsga2    30

SKIP_JMETAL_BUILD=1 ./scripts/run_smpso_comparison.sh   results/smpso   30
SKIP_JMETAL_BUILD=1 ./scripts/run_spea2_comparison.sh   results/spea2   30
SKIP_JMETAL_BUILD=1 ./scripts/run_ibea_comparison.sh    results/ibea    30
SKIP_JMETAL_BUILD=1 ./scripts/run_moead_comparison.sh   results/moead   30
SKIP_JMETAL_BUILD=1 ./scripts/run_nsga3_comparison.sh   results/nsga3   30
SKIP_JMETAL_BUILD=1 ./scripts/run_smsemoa_comparison.sh results/smsemoa 30
SKIP_JMETAL_BUILD=1 ./scripts/run_gde3_comparison.sh    results/gde3    30
SKIP_JMETAL_BUILD=1 ./scripts/run_agemoea_comparison.sh results/agemoea 30
SKIP_JMETAL_BUILD=1 ./scripts/run_mocell_comparison.sh  results/mocell  30
SKIP_JMETAL_BUILD=1 ./scripts/run_paes_comparison.sh    results/paes    30
SKIP_JMETAL_BUILD=1 ./scripts/run_uf_comparison.sh      results/uf      30
SKIP_JMETAL_BUILD=1 ./scripts/run_lz09_comparison.sh    results/lz09    30
SKIP_JMETAL_BUILD=1 ./scripts/run_wfg_comparison.sh     results/wfg     30
```

Statistical analysis:
```bash
python3 scripts/analyze_mo_comparison.py results/nsga2
```

Each script generates `results/<name>/environment.md` with the exact hardware and
compiler configuration at time of execution, for reproducibility documentation.

---

## Single-objective comparison: ea-cpp-modern vs ea-cpp-original

This comparison isolates the **architectural** speedup (SoA + modern C++ vs AoS + virtual
dispatch) independent of the language speedup (C++ vs JVM).

Both systems implement the same (µ,λ)-ES on Rastrigin d=20 with identical parameters:

| Parameter | Value |
|---|---|
| µ (parents) | 100 |
| λ (offspring) | 99 |
| Selection | Tournament k=2 |
| Crossover | BLX-α, prob=0.9, α=0.1 |
| Mutation | Gaussian, rate=0.6321, step=1.0 |
| Replacement | Comma |
| Problem | Rastrigin d=20, [−5.12, 5.12] |
| RNG | mt19937\_64, seed=1..NUMRUNS |

`ea-cpp-original` retains the AoS layout and virtual dispatch of its Java source
([Bio4Res/ea](https://github.com/Bio4Res/ea)), compiled with the same C++23 flags.

### Run

```bash
# Quick test (~20 s)
./scripts/run_so_comparison.sh /tmp/test_so 1000000 5

# Paper run (~15 min)
./scripts/run_so_comparison.sh results/so_paper 10000000 30

python3 scripts/analyze_so_comparison.py results/so_paper --csv
```

Override `ORIGINAL_ROOT` if ea-cpp-original is not at `$HOME/ea_cpp_original`:
```bash
ORIGINAL_ROOT=/path/to/ea_cpp_original \
    ./scripts/run_so_comparison.sh results/so_paper 10000000 30
```

### Reference results (30 runs × 10M evals, WSL2)

| Metric | ea-cpp-original (AoS) | ea-cpp-modern (SoA) | Ratio |
|---|---|---|---|
| Median time (s) | 7.84 | 6.74 | **1.16×** |
| Median fitness | 0.0 | 0.0 | = |
| hits≈0 (of 30) | 30/30 | 30/30 | = |

> **WSL2 note**: the mean is unstable due to thermal throttling (one outlier at 14.2s
> was observed). Use **median** as the robust estimator. A 60 s cooldown between systems
> mitigates but does not eliminate the thermal effect.
>
> **Bare-metal results are pending**: the SoA advantage is expected to be more pronounced
> outside WSL2.
