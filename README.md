# ea-cpp-modern

High-performance multi-objective evolutionary algorithm framework in **C++23**, designed as a faithful re-implementation of [jMetal 7.4](https://github.com/jMetal/jMetal) with a zero-overhead C++ architecture.

[![C++23](https://img.shields.io/badge/C%2B%2B-23-blue.svg)](https://en.cppreference.com/w/cpp/23)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Speedup vs jMetal 7.4

30 independent runs, ZDT1–4 (where applicable), WSL2 Ubuntu 24.04, GCC 14, `-O3 -march=native`.
Wilcoxon signed-rank test on IGD: all results statistically non-significant (p > 0.05) unless noted.

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

See [BENCHMARKS.md](BENCHMARKS.md) for full results across ZDT, WFG, UF, LZ09, and DTLZ families.

## Requirements

| Dependency | Version | Purpose |
|---|---|---|
| GCC | 14+ (C++23) | Framework compilation |
| Java JDK | 17+ | jMetal comparison (optional) |
| Maven | 3.6+ | jMetal build (optional) |
| Python 3 | 3.8+ | Statistical analysis (optional) |
| scipy | any | Wilcoxon test (optional) |

Install on Ubuntu/Debian:
```bash
sudo apt-get install g++-14 default-jdk maven python3-scipy
```

## Quick Start

```bash
git clone https://github.com/OWNER/ea-cpp-modern.git
cd ea-cpp-modern
make
make test
```

Run NSGA-II on ZDT1–4 (5 runs, seed 42):
```bash
./build/nsga2_runner 5 42
```

Output is CSV: `problem,run,seed,time_ms,igd`.

## Reproducing the jMetal comparison

Clone and build jMetal once:
```bash
git clone https://github.com/jMetal/jMetal.git $HOME/jMetal
cd $HOME/jMetal && mvn clean install -DskipTests -q
```

Then run any comparison script from the repo root. The variable `JMETAL_ROOT` defaults to
`$HOME/jMetal`; override if you cloned jMetal elsewhere:

```bash
# NSGA-II vs jMetal, 30 runs
./scripts/run_nsga2_comparison.sh results/nsga2 30

# All algorithms — commands from BENCHMARKS.md
SKIP_JMETAL_BUILD=1 ./scripts/run_smpso_comparison.sh  results/smpso  30
SKIP_JMETAL_BUILD=1 ./scripts/run_ibea_comparison.sh   results/ibea   30
# ... etc.

# Statistical analysis
python3 scripts/analyze_mo_comparison.py results/nsga2
```

The compiler defaults to `g++-14`. Override with `CXX_BENCH=g++-16` for GCC 16.

## Architecture

The framework is built around three design decisions that jointly account for the speedup:

**1. Structure of Arrays (SoA) population**
All genes, objectives, and auxiliary data live in flat contiguous arrays.
No per-individual heap allocation; evaluation loops are SIMD-vectorizable.

**2. Deducing `this` (C++23) instead of virtual/CRTP**
Operators use `void apply(this OpType& self, Population& pop, ...)`.
Zero vtable overhead, full inlining, no inheritance hierarchy.

**3. C++20 Concepts for compile-time contracts**
`Crossover`, `Mutation`, `Selection`, `Replacement`, `EvalFunctor` — the compiler
rejects mismatched types with a clear diagnostic, not a runtime error.

```
include/ea/
├── core/          — Population (SoA), Encoding, Concepts, Config
├── algorithm/     — 11 MO algorithms
├── operator/
│   ├── crossover/ — 21 operators (SBX, BLX-αβ, DE, PMX, …)
│   ├── mutation/  — 17 operators (Polynomial, BitFlip, Swap, …)
│   ├── selection/ — 8 operators (BinaryTournament, Ranking+Crowding, …)
│   └── replacement/ — NSGA-II/III, Mu+Lambda, MOEA/D
├── problem/       — ZDT (6), DTLZ (7), WFG (9), UF (10), LZ09 (9)
├── indicator/     — IGD, GD, Hypervolume, Spread, Epsilon, R2, Hausdorff
└── util/          — ReferencePoint, Aggregation, Random (Xoshiro256**)
```

See [ARCHITECTURE.md](ARCHITECTURE.md) for patterns, porting rules, and the AoS→SoA rationale.

## Usage example

```cpp
#include <ea/ea.hpp>

int main() {
    ea::ZDT1 problem(30);
    ea::Population pop(100, 30, 2);
    pop.lower_bounds = std::vector<double>(30, 0.0);
    pop.upper_bounds = std::vector<double>(30, 1.0);

    for (int i = 0; i < 100; ++i)
        for (int j = 0; j < 30; ++j)
            pop.gene(i, j) = ea::Random::instance().uniform(0.0, 1.0);

    ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> nsga;
    nsga.pop_size  = 100;
    nsga.max_evals = 25000;

    nsga.run(pop, [&problem](ea::Population& p, int idx) {
        problem.evaluate(p, idx);
    });
}
```

## Bare-metal validation

For publication-quality results, run on bare-metal Linux with the CPU governor set to
`performance` and Turbo Boost disabled:

```bash
sudo cpupower frequency-set -g performance
echo 1 | sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo
```

Then re-run all comparison scripts as shown in [BENCHMARKS.md](BENCHMARKS.md).

## License

MIT — see [LICENSE](LICENSE).

## Contact

- Issues: GitHub Issues
- Email: jmcruz@uma.es
