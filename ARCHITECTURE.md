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