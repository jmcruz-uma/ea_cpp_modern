# ea-cpp

High-performance evolutionary algorithm framework in C++23.

Port of jMetal operators and algorithms with a C++-native architecture optimized for
runtime performance and energy efficiency, using:

- **Structure of Arrays (SoA)** for cache-friendly population storage
- **Deducing this (C++23)** for zero-overhead explicit object parameters
- **Concepts (C++20)** for compile-time contract enforcement
- **`std::variant` + `std::visit`** for runtime-configurable dispatch without vtables
- **CRTP eliminated** — deducing this provides the same static dispatch without inheritance boilerplate

## Architecture

```
Population (SoA)
  ├── genes[]      contiguo en memoria
  ├── objectives[] contiguo en memoria
  ├── constraints[]contiguo en memoria
  └── flags[]      bitfield: evaluated, feasible, etc.

Operators (Concepts + Deducing this)
  ├── Crossover: SBX, BLX-α, DE, PMX, Cycle, ...
  ├── Mutation: Polynomial, BitFlip, Swap, Inversion, ...
  ├── Selection: BinaryTournament, NaryTournament, Ranking, ...
  └── Replacement: MuPlusLambda, NSGAII, NSGAIII, ...

Algorithms (Template compositions)
  ├── NSGA-II, NSGA-III, MOEA/D, SPEA2, SMPSO, ...
  └── Configured via compile-time (templates) or runtime (JSON → variant)

Problems
  ├── ZDT1-6, DTLZ1-7, WFG1-9
  └── Batch evaluation interface for SIMD
```

## Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

## Requirements

- C++23 compiler (GCC 14+, Clang 18+, MSVC 2022 19.34+)
- CMake 3.25+
- Google Test (fetched automatically for tests)
- nlohmann/json 3.11.3+ (fetched automatically if not found)

## License

GPL v3 — fork of Bio4Res/ea-cpp with jMetal algorithm/operator ports.
jMetal code used as specification reference (MIT license).
