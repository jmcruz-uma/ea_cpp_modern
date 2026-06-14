# ARCHITECTURE.md — ea-cpp Reference

This document defines the architecture that ALL implementations must follow.

## Core Design Principles

1. **Structure of Arrays (SoA)** — Population stores genes, objectives, constraints in contiguous arrays. No per-individual heap allocations. Access via `pop.gene(i, j)`, `pop.objective(i, o)`.

2. **Deducing this (C++23)** — Operators use `this auto& self` as first parameter. NO virtual dispatch, NO CRTP, NO inheritance hierarchies.

3. **Concepts (C++20)** — All operators satisfy `ea::Crossover`, `ea::Mutation`, `ea::Selection`, or `ea::Replacement` concepts. The compiler enforces contracts at compile time.

4. **std::variant for runtime dispatch** — When algorithm is configured from JSON, operators are held in `std::variant<SBX, BLXAlpha, ...>`. `std::visit` dispatches without vtable.

5. **Batch-friendly** — Problems implement `evaluate(Population&, int)` and optionally `evaluate_batch(Population&, int start, int count)`.

## File Structure

```
include/ea/
├── core/
│   ├── encoding.hpp       — Encoding enum (Real, Binary, Permutation, Integer, Composite)
│   ├── population.hpp     — SoA Population struct
│   ├── concepts.hpp       — C++20 Concepts for operators, problems, algorithms
│   └── comparator.hpp     — Dominance comparison, crowding distance, fast NDS
├── operator/
│   ├── crossover/         — SBX, BLXAlpha, DE, PMX, Cycle, Uniform, etc.
│   ├── mutation/          — Polynomial, BitFlip, Swap, Inversion, Uniform, NonUniform, etc.
│   ├── selection/         — BinaryTournament, NaryTournament, Random, etc.
│   └── replacement/       — MuPlusLambda, NSGA2Replacement, etc.
├── algorithm/             — NSGAII, NSGAIII, MOEAD, etc.
├── problem/               — ZDT, DTLZ, WFG, etc.
├── indicator/             — Hypervolume, IGD, GD, Epsilon, Spread
└── util/
    └── random.hpp         — Thread-safe RNG with batch generation
```

## Operator Pattern (Deducing this)

```cpp
// Every crossover follows this pattern:
struct SBXCrossover {
    double distribution_index = 20.0;
    double crossover_probability = 0.9;

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Real; }

    void apply(this SBXCrossover& self, Population& pop,
               int parent_a, int parent_b, int child_start) {
        // Implementation using self.distribution_index, etc.
    }
};
```

## Selection Pattern

```cpp
struct BinaryTournamentSelection {
    void select(this BinaryTournamentSelection& self, Population& pop,
                std::vector<int>& mating_pool,
                const std::vector<int>& ranks,
                const std::vector<double>& crowding_dist) {
        // Fills mating_pool with selected parent indices
    }
};
```

## Problem Pattern

```cpp
struct ZDT1 {
    static constexpr int num_objectives() { return 2; }

    void evaluate(Population& pop, int idx) const {
        double g = g_function(pop, idx);
        pop.objective(idx, 0) = pop.gene(idx, 0);
        pop.objective(idx, 1) = g * (1.0 - std::sqrt(pop.gene(idx, 0) / g));
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};
```

## Algorithm Pattern (Template Composition)

```cpp
template<typename CX, typename MT, typename SEL>
struct NSGAII {
    CX crossover;
    MT mutation;
    SEL selection;
    int pop_size = 100;
    int max_evals = 25000;

    static constexpr std::string_view name() { return "NSGA-II"; }

    template<typename Problem>
    void run(this auto& self, Population& pop, Problem&& problem) { ... }
};
```

## Evolutionary History: From ea_cpp_original to ea_cpp_modern

### The Three-Tier Lineage

```
Bio4Res/ea  (Java, single-objective)   +   jMetal 7.4  (Java, multi-objective)
      ↓ C++ port, OO style                       ↓ semantic reference
ea_cpp_original  (C++23, AoS, virtual dispatch, energy-reduction focus)
      ↓ architecture redesign
ea_cpp_modern    (C++23, SoA, Concepts, deducing-this, this paper)
```

`ea_cpp_original` is a direct C++ port of the Java framework `Bio4Res/ea`
(https://github.com/Bio4Res/ea), explicitly aimed at reducing energy consumption via C++.
It compiles with the same flags (`-O3 -march=native -std=c++23`) but retains the OO
structure of its Java source: inheritance hierarchies, virtual dispatch, and
Array-of-Structs population layout.

### What ea_cpp_original Got Right

- Same algorithmic operators as the Java reference (SBX, PM, tournament, etc.)
- Same RNG philosophy (Xoshiro, singleton `Random::instance()`)
- Same problem definitions (Sphere, Rastrigin, ZDT, DTLZ families)
- Correct `Island + Topology + Migration` model
- C++23 compilation required (uses `std::format`, `std::execution::par_unseq`)

### The Architectural Gap

| Aspect | ea_cpp_original | ea_cpp_modern |
|---|---|---|
| Population layout | `vector<Individual>` (AoS) — heap per individual | `Population<>` SoA — flat contiguous arrays |
| Polymorphism | Pure virtual + Factory registry | `deducing this` + Concepts |
| Problem interface | `ObjectiveFunction` inheritance (vtable) | `Problem<T>` Concept — no vtable |
| Operator interface | `VariationOperator::apply_(vector<Individual>)` returns new `Individual` | `apply(Population&, int, int, int)` — in-place, zero allocation |
| Encoding | `Union{int, double}` inside `Genotype` | `gene_t<Encoding>` — compile-time type |
| Algorithm config | JSON → `EAConfiguration` → Factory → `Island` setup | `NSGAII<SBX, PM>` — template composition |

### Why AoS → SoA Is the Core Contribution

In ea_cpp_original, evaluating the objectives of N individuals requires
touching N heap-allocated `Genotype` objects with pointer indirection:

```
for Individual ind : population:             // vtable + pointer chase
    double fit = objective.evaluate(ind);    // virtual call
    // ind.genome → heap → ind.genome->genes → heap → [gene0, gene1, ...]
```

In ea_cpp_modern, the same loop is:

```cpp
for (int i = 0; i < pop_size; ++i)          // sequential, no indirection
    problem.evaluate(pop, i);               // inline, no virtual call
    // pop.genes[i * dim + j]               // cache-line hot, SIMD-ready
```

The transition from AoS to SoA turns objective evaluation and NDS into
SIMD-vectorizable loops. Combined with eliminating vtable dispatch, this
accounts for the majority of the speedup observed vs jMetal.

### Implication for the Paper

The speedup ea_cpp_modern achieves over jMetal (3–13×, depending on algorithm)
has two separable components:

1. **Language advantage** (C++ over JVM): no GC pauses, no boxing overhead,
   deterministic inlining. This would be present even in ea_cpp_original.

2. **Architecture advantage** (SoA + Concepts over AoS + virtual):
   cache-friendly layout, SIMD vectorization, zero-allocation operators.
   This is the contribution *beyond* a naive C++ port.

A 3-way benchmark (jMetal → ea_cpp_original → ea_cpp_modern) would isolate
these two effects precisely. As a lower bound, the language advantage alone
is estimated at ~2–3× (typical C++/JVM ratio for compute-bound code);
the remaining factor comes from the architectural redesign.

---

## Key Differences from jMetal

| Aspect | jMetal | ea-cpp |
|---|---|---|
| Population | `List<Solution<Double>>` (AoS) | SoA with contiguous arrays |
| Solution | Object with `List<Double>` | View into Population arrays |
| Polymorphism | Interfaces + virtual | Deducing this + Concepts |
| Runtime config | Builder pattern | std::variant + std::visit |
| Evaluation | Per-solution | Batch-friendly |
| Threading | ForkJoinPool | std::execution::par_unseq |
| Cache locality | Poor (pointer chasing) | Excellent (SoA) |

## Porting Rules

When porting from jMetal:
1. Read the Java implementation as SPECIFICATION, not as code to translate
2. Implement using SoA access patterns (pop.gene(), pop.objective())
3. Use deducing this, not virtual/CRTP
4. Satisfy the appropriate Concept
5. Add `encoding()` static method
6. Ensure bounds checking with std::clamp
7. Document the original jMetal source class for reference