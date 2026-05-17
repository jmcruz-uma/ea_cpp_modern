# ea-cpp

High-performance evolutionary algorithm framework in **C++23** — jMetal operators and algorithms, C++-native architecture.

[![C++23](https://img.shields.io/badge/C%2B%2B-23-blue.svg)](https://en.cppreference.com/w/cpp/23)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![IGD](https://img.shields.io/badge/IGD-0.005%20(ZDT1)-brightgreen.svg)](benchmarks)

## Features

- **Zero-overhead architecture**: SoA Population, deducing this, concepts — no virtual/CRTP
- **jMetal parity**: Same operators and algorithms, same IGD, **12x faster**
- **76 headers**: Core, 21 crossover, 17 mutation, 8 selection, 4 replacement, 13 algorithms, 22 problems, 7+ indicators
- **C++23 modern**: Template composition, compile-time dispatch, clean APIs

## Quick Start

```bash
# Requirements: g++-14+ with C++23 support
git clone https://github.com/elCanosail/ea-cpp.git
cd ea-cpp

# Run benchmark (NSGA-II on ZDT1)
make benchmark

# Run all tests
make test

# Full benchmark suite (30 runs per config)
make benchmark-full
```

## Architecture

```
ea-cpp/
├── include/ea/
│   ├── core/          # Population (SoA), Encoding, Concepts, Comparators, Config
│   ├── algorithm/     # 13 algorithms (NSGA-II, NSGA-III, MOEA/D, SPEA2, RVEA, IBEA, ...)
│   ├── operator/
│   │   ├── crossover/ # 21 operators (SBX, BLX-αβ, DE, PMX, ERX, ...)
│   │   ├── mutation/  # 17 operators (Polynomial, BitFlip, Swap, Inversion, ...)
│   │   ├── selection/ # 8 operators (BinaryTournament, RankingAndCrowding, ...)
│   │   └── replacement/ # 4 operators (NSGAII, NSGAIII, MuPlusLambda, ...)
│   ├── problem/       # 22 benchmark problems (ZDT, DTLZ, WFG families)
│   ├── indicator/     # 7 quality indicators (Hypervolume, IGD, GD, Spread, Epsilon, R2, Hausdorff)
│   └── util/          # ReferencePoint, Aggregation, Random
├── tests/             # Unit tests + regression tests + CI
├── examples/          # Benchmarks + JSON config demos
├── scripts/           # Automated benchmark suite + energy estimation
├── paper.md           # Paper (markdown)
├── paper.tex          # Paper (LaTeX)
├── Makefile           # Build automation
└── README.md          # This file
```

## Benchmarks

NSGA-II on ZDT1 (100 pop, 25000 evals, 10 runs):

| Metric | ea-cpp (C++23) | jMetal (Java 21) | Ratio |
|--------|----------------|-------------------|-------|
| Mean IGD | 0.005 ± 0.000 | 0.005 ± 0.000 | Parity |
| Mean Time | 153 ms | 1831 ms | **12x faster** |
| Time stddev | ±10 ms | ±156 ms | **15x more stable** |

Extended benchmarks (5 runs):

| Algorithm | Problem | IGD | Time (ms) |
|-----------|---------|-----|-----------|
| NSGAII | ZDT1 | 0.00019 | 162 |
| NSGAII | ZDT2 | 0.00020 | 163 |
| NSGAII | ZDT3 | 0.00948 | 161 |
| NSGAII | ZDT6 | 0.00892 | 198 |
| RVEA | ZDT1 | 0.01060 | 178 |
| IBEA | ZDT1 | 0.20652 | 233 |
| RandomSearch | ZDT1 | 0.05339 | 6265 |

## Usage

### Basic Example

```cpp
#include <ea/ea.hpp>

int main() {
    // Create problem
    ea::ZDT1 problem(30);

    // Create population
    ea::Population pop(100, 30, 2);
    pop.lower_bounds = std::vector<double>(30, 0.0);
    pop.upper_bounds = std::vector<double>(30, 1.0);

    // Initialize randomly
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 30; ++j) {
            pop.gene(i, j) = ea::Random::instance().uniform(0.0, 1.0);
        }
    }

    // Configure and run NSGA-II
    ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> nsga;
    nsga.pop_size = 100;
    nsga.max_evals = 25000;
    nsga.crossover.distribution_index = 20.0;
    nsga.mutation.distribution_index = 20.0;

    nsga.run(pop, [&problem](ea::Population& p, int idx) {
        problem.evaluate(p, idx);
    });

    return 0;
}
```

### JSON Configuration

```json
{
  "algorithm": "NSGAII",
  "crossover": { "type": "sbx", "probability": 0.9, "distribution_index": 20.0 },
  "mutation": { "type": "polynomial", "rate": -1.0, "distribution_index": 20.0 },
  "pop_size": 100,
  "max_evals": 25000
}
```

See `examples/config_example.cpp` for full usage.

## Tests

```bash
# Compile and run tests
make test

# Expected output:
# ========================================
# ea-cpp Operator Unit Tests
# ========================================
# TEST: SBXCrossover ... PASS
# TEST: PolynomialMutation ... PASS
# TEST: BinaryTournamentSelection ... PASS
# Results: 3 passed, 0 failed
```

## Version History

| Version | Date | Features |
|---------|------|----------|
| v0.1.0 | 2026-05-16 | Framework base, 76 headers, jMetal parity |
| v0.2.0 | 2026-05-17 | Robustez, tests, Spread indicator, CI GitHub Actions |
| v0.3.0 | 2026-05-17 | RVEA, IBEA, RandomSearch algorithms |
| v0.4.0 | 2026-05-17 | Epsilon, R2, Hausdorff indicators |
| v0.5.0 | 2026-05-17 | Makefile, README, benchmark script |
| v0.6.0 | 2026-05-17 | JSON configuration without recompiling |
| **v1.0.0** | **2026-05-17** | **Paper-ready release** |

## Paper

See `paper.md` (markdown) or `paper.tex` (LaTeX) for the full paper.

**Key results:**
- Algorithmic parity with jMetal (IGD 0.005)
- 12× faster execution
- 15× more temporally stable
- 13 algorithms, 22 problems, 7+ indicators

## License

MIT

## Contact

- Repository: https://github.com/elCanosail/ea-cpp
- Issues: GitHub Issues
- Paper review: See `REVIEW.md`
