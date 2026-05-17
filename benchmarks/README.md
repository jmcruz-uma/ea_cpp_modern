# ea-cpp — Benchmarks

Architecture comparison: ea-cpp v1.0 vs Bio4Res/ea-cpp (original).

## Hardware

- **CPU:** AMD/Intel x86_64 (server-grade recommended for mission-critical)
- **Compiler:** Clang 18, `-std=c++23 -O3 -march=native`
- **Runs:** 30 (unless noted)

## Methodology

We compare two C++ implementations of the same algorithmic lineage:

| | Original (Bio4Res) | Current (v1.0) |
|---|---|---|
| **Architecture** | AoS: `Individual` + `unique_ptr<Genotype>` | SoA: `Population` (contiguous memory) |
| **Dispatch** | Virtual (runtime polymorphism) | Templates + Concepts (zero-overhead) |
| **Initialization** | Dynamic (static init order risk) | `constinit` (zero UB, compile-time) |
| **Config** | JSON (runtime) | Compile-time composition + Factory |
| **Multi-objective** | ❌ None | ✅ 10 algorithms (NSGA-II, NSGA-III, MOEA/D, ...) |
| **LOC** | 4,942 (55 files) | 76+ headers (header-only) |

## Results

### 1. Evaluation Throughput (25K individuals)

| Problem | Original (full EA loop) | Current (pure eval) | Speedup |
|---------|------------------------|---------------------|---------|
| Sphere(30) | ~94 ms | 0.8 ms | ~117x |
| Rastrigin(30) | ~105 ms | 15.0 ms | ~7x |
| ZDT1(30, 2obj) | N/A | 0.4 ms | — |

> **Note:** Original's ~100ms includes Individual heap allocation, virtual dispatch, and JSON parsing overhead — not just function evaluation. Our pure evaluation is 7-117x faster.

### 2. Operator Throughput (SBX Crossover + Polynomial Mutation)

| Dimension | 25K ops |
|-----------|---------|
| dim=30 | 0.9 ms |
| dim=100 | 2.8 ms |
| dim=500 | 14.1 ms |

Linear scaling with dimension. Zero virtual dispatch overhead.

### 3. NSGA-II End-to-End (25K evaluations, 30 runs)

| Problem | Time | IGD |
|---------|------|-----|
| ZDT1(30) | 295 ms | 0.00492 |
| ZDT2(30) | 287 ms | 0.00478 |
| ZDT3(30) | 295 ms | 0.10521 |
| ZDT4(30) | 294 ms | 0.00455 |
| ZDT6(30) | 232 ms | 0.04494 |
| DTLZ1(3obj) | 236 ms | — |
| DTLZ2(3obj) | 239 ms | — |

> **⚠ Original has NO multi-objective capability.** This comparison cannot be made against Bio4Res/ea-cpp.

### 4. Scalability — Population Size (ZDT1, 250 generations)

| Population | Time | Evaluations | Scaling |
|-----------|------|-------------|---------|
| 50 | 132 ms | 12,500 | 1.0x |
| 100 | 312 ms | 25,000 | 2.4x |
| 200 | 760 ms | 50,000 | 5.8x |
| 500 | 3,130 ms | 125,000 | 23.7x |
| 1,000 | 10,882 ms | 250,000 | 82.5x |

Scaling follows O(MN²) from non-dominated sorting — inherent to NSGA-II, not our implementation. For large populations, NSGA-III or MOEA/D scale better (O(MN) or O(N log N)).

### 5. (µ,λ)-ES Comparison (Same Algorithm, Mono-objective)

| Problem | Original | Current | Ratio |
|---------|----------|---------|-------|
| Sphere(30) | 94 ms | 151 ms | 0.62x |
| Rastrigin(30) | 105 ms | 161 ms | 0.65x |
| Rosenbrock(30) | 93 ms | 146 ms | 0.64x |
| Ackley(30) | 102 ms | 157 ms | 0.65x |

> Current is slower in mono-objective because our (µ,λ)-ES is an ad-hoc benchmark, not an optimized algorithm. The original's `Individual` with heap allocation is surprisingly cache-friendly at pop=100. **This is not our value proposition.**

## Mission-Critical Relevance

| Property | Original | Current (v1.0) |
|----------|----------|-----------------|
| **Undefined behavior** | Possible (static init order) | Eliminated (`constinit`) |
| **Heap allocation in hot path** | Per-individual `unique_ptr<Genotype>` | Zero (SoA, contiguous) |
| **Dispatch overhead** | Virtual (vtable lookup) | Zero (templates, inlined) |
| **Reproducibility** | Depends on JSON config + init order | Deterministic (compile-time composition) |
| **Audit surface** | 4,942 LOC (OOP, virtual, JSON) | Header-only, Concepts-validated |
| **Multi-objective** | None | NSGA-II, NSGA-III, MOEA/D, SPEA2, IBEA, ... |
| **Extensibility** | Fork + modify source | `EA_REGISTER_CROSSOVER("name", Type)` |

## Running the Benchmarks

```bash
# Build
clang++-18 -std=c++23 -O3 -march=native -I include/ benchmarks/bench_v1_vs_current.cpp -o bench_current
clang++-18 -std=c++23 -O3 -march=native -I include/ benchmarks/bench_arch_comparison.cpp -o bench_arch

# Run (default: 30 runs)
./bench_current --runs 30 --verbose
./bench_arch --runs 30

# Full benchmark (includes scalability sweep)
./bench_arch --runs 30 --full

# CSV output
./bench_current --csv
```