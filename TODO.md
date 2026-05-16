# TODO.md — ea-cpp Implementation Roadmap

## ✅ Phase 0: Foundation
- [x] CMakeLists.txt with C++23, GTest, benchmarks
- [x] Population SoA (population.hpp)
- [x] Encoding enum (encoding.hpp)
- [x] Concepts (concepts.hpp)
- [x] Comparator (comparator.hpp) — dominance, crowding distance, fast NDS
- [x] Random (random.hpp) — thread-safe, batch generation

## ✅ Phase 1: Core Operators
### Crossover — 21 total
- [x] SBX, BLX-α, BLX-αβ, Arithmetic, WholeArithmetic
- [x] Laplace, DE (6 variants), ParentCentric/UNDX
- [x] PMX, Cycle, ERX, OXD, PositionBased, Uniform
- [x] FuzzyRecombination, IntegerSBX, HUX, Composite
- [x] SinglePoint, TwoPoint, NPoint

### Mutation — 17 total
- [x] Polynomial, BitFlip, Swap, Inversion, Insert, Scramble
- [x] Uniform, NonUniform, PowerLaw, LevyFlight, SimpleRandom
- [x] CDG, Displacement, SimpleInversion, GroupedPolynomial
- [x] LinkedPolynomial, IntegerPolynomial

### Selection — 8 total
- [x] BinaryTournament, NaryTournament, Random
- [x] RankingAndCrowding, BestSolution, NaryRandom
- [x] DifferentialEvolution, SpatialSpreadDeviation

### Replacement — 4 total
- [x] MuPlusLambda, MuCommaLambda, NSGAII, NSGAIII

## ✅ Phase 2: Algorithms — 9 total
- [x] NSGA-II — tested with ZDT1 (IGD=0.004)
- [x] NSGA-III — Das-Dennis reference points, niching
- [x] MOEA/D — Tchebycheff, WeightedSum, PBI
- [x] SPEA2 — strength fitness, archive truncation
- [x] SMS-EMOA — hypervolume selection
- [x] GDE3 — DE with NSGA-II replacement
- [x] SMPSO — PSO with velocity constriction
- [x] AGE-MOEA — adaptive geometry estimation
- [x] PAES — (1+1)-ES with adaptive grid
- [x] MOCell — cellular GA with archive feedback

## ✅ Phase 3: Problems — 22 benchmarks
- [x] ZDT1-6 (complete suite)
- [x] DTLZ1-7 (complete suite)
- [x] WFG1-9 (complete suite with all transformations)

## ✅ Phase 4: Quality Indicators
- [x] Hypervolume (WFG — 2D efficient, 3D+ recursive)
- [x] IGD (Inverted Generational Distance)
- [x] GD (Generational Distance)
- [x] Spread (Δ)

## ✅ Phase 5: Utilities
- [x] ReferencePoint (Das-Dennis generation, two-layer)
- [x] Aggregation (Tchebycheff, WeightedSum, PBI)
- [x] AdaptiveGrid (for PAES)

## 🔲 Phase 6: Remaining from jMetal
### Algorithms
- [ ] RandomSearch (trivial — just random generation + evaluation)

### Selection
- [ ] BoltzmannSelection, TruncationSelection, StochasticUniversalSampling
- [ ] RankingAndPreferenceSelection, RankingAndDirScoreSelection

### Replacement
- [ ] RVEA, AGEMOEA, SMSEMOA specific replacements
- [ ] PairwiseReplacement, RandomReplacement, TournamentReplacement

## 🔲 Phase 7: Infrastructure
- [ ] JSON configuration → std::variant construction
- [ ] Island model (from original ea-cpp)
- [ ] Benchmarks vs jMetal (speed + energy)
- [ ] More quality indicators (R2, Epsilon, AverageHausdorff)
- [ ] Documentation