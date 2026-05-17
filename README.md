# ea-cpp

High-performance evolutionary algorithm framework in **C++23** — jMetal operators and algorithms, C++-native architecture.

[![C++23](https://img.shields.io/badge/C%2B%2B-23-blue.svg)](https://en.cppreference.com/w/cpp/23)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![IGD](https://img.shields.io/badge/IGD-0.005%20(ZDT1)-brightgreen.svg)](benchmarks)

## Features

- **Zero-overhead architecture**: SoA Population, deducing this, concepts — no virtual/CRTP
- **jMetal parity**: Same operators and algorithms, same IGD, **12x faster**
- **76 headers**: Core, 21 crossover, 17 mutation, 8 selection, 4 replacement, 10 algorithms, 22 problems, 4+ indicators
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
│   ├── core/          # Population (SoA), Encoding, Concepts, Comparators
│   ├── algorithm/     # NSGA-II, NSGA-III, MOEA/D, SPEA2, RVEA, IBEA, ...
│   ├── operator/
│   │   ├── crossover/ # SBX, BLX-αβ, DE, PMX, ERX, ...
│   │   ├── mutation/  # Polynomial, BitFlip, Swap, Inversion, ...
│   │   ├── selection/ # BinaryTournament, RankingAndCrowding, ...
│   │   └── replacement/ # NSGAII, NSGAIII, MuPlusLambda, ...
│   ├── problem/       # ZDT1-6, DTLZ1-7, WFG1-9
│   ├── indicator/     # Hypervolume, IGD, GD, Spread, Epsilon, R2, Hausdorff
│   └── util/          # ReferencePoint, Aggregation, Random
├── tests/             # Unit tests + regression tests
├── examples/          # Benchmarks + demos
└── scripts/           # Automated benchmark suite
```

## Benchmarks

NSGA-II on ZDT1 (100 pop, 25000 evals, 10 runs):

| Metric | ea-cpp (C++23) | jMetal (Java 21) | Ratio |
|--------|----------------|-------------------|-------|
| Mean IGD | 0.005 ± 0.000 | 0.005 ± 0.000 | Parity |
| Mean Time | 153 ms | 1831 ms | **12x faster** |
| Time stddev | ±10 ms | ±156 ms | **15x more stable** |

## Version History

| Version | Date | Features |
|---------|------|----------|
| v0.1.0 | 2026-05-16 | Framework base, 76 headers, jMetal parity |
| v0.2.0 | 2026-05-17 | Robustez, tests, Spread indicator, CI |
| v0.3.0 | 2026-05-17 | RVEA, IBEA, RandomSearch algorithms |
| v0.4.0 | 2026-05-17 | Epsilon, R2, Hausdorff indicators |

## License

MIT
