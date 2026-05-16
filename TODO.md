# TODO.md — ea-cpp Implementation Roadmap

## Phase 0: Foundation ✅ (current)
- [x] CMakeLists.txt with C++23, GTest, benchmarks
- [x] Population SoA (population.hpp)
- [x] Encoding enum (encoding.hpp)
- [x] Concepts (concepts.hpp) — Crossover, Mutation, Selection, Replacement, Problem, Algorithm
- [x] Comparator (comparator.hpp) — dominance, crowding distance, fast NDS
- [x] Random (random.hpp) — thread-safe, batch generation
- [x] SBXCrossover (crossover/sbx.hpp)
- [x] PolynomialMutation (mutation/polynomial.hpp)
- [x] BinaryTournamentSelection (selection/binary_tournament.hpp)
- [x] NSGA-II initial skeleton (algorithm/nsga2.hpp)

## Phase 1: Core Operators (in progress)
### Crossover (continuous)
- [x] SBXCrossover
- [ ] BLXAlphaCrossover — jMetal: BLXAlphaCrossover.java
- [ ] BLXAlphaBetaCrossover — jMetal: BLXAlphaBetaCrossover.java
- [ ] ArithmeticCrossover — jMetal: ArithmeticCrossover.java
- [ ] WholeArithmeticCrossover — jMetal: WholeArithmeticCrossover.java
- [ ] LaplaceCrossover — jMetal: LaplaceCrossover.java
- [ ] DifferentialEvolutionCrossover — jMetal: DifferentialEvolutionCrossover.java
- [ ] ParentCentricCrossover (UNDX) — jMetal: UnimodalNormalDistributionCrossover.java

### Crossover (permutation)
- [ ] PMXCrossover — jMetal: PMXCrossover.java
- [ ] CycleCrossover — jMetal: CycleCrossover.java
- [ ] EdgeRecombinationCrossover — jMetal: EdgeRecombinationCrossover.java
- [ ] OXDCrossover — jMetal: OXDCrossover.java
- [ ] PositionBasedCrossover — jMetal: PositionBasedCrossover.java
- [ ] NPointCrossover — jMetal: NPointCrossover.java
- [ ] TwoPointCrossover — jMetal: TwoPointCrossover.java
- [ ] UniformCrossover (binary/real) — jMetal: UniformCrossover.java

### Mutation
- [x] PolynomialMutation
- [ ] BitFlipMutation — jMetal: BitFlipMutation.java
- [ ] SwapMutation — jMetal: PermutationSwapMutation.java
- [ ] InversionMutation — jMetal: InversionMutation.java
- [ ] InsertMutation — jMetal: InsertMutation.java
- [ ] ScrambleMutation — jMetal: ScrambleMutation.java
- [ ] UniformMutation — jMetal: UniformMutation.java
- [ ] NonUniformMutation — jMetal: NonUniformMutation.java
- [ ] PowerLawMutation — jMetal: PowerLawMutation.java
- [ ] LevyFlightMutation — jMetal: LevyFlightMutation.java

### Selection
- [x] BinaryTournamentSelection
- [ ] NaryTournamentSelection — jMetal: NaryTournamentSelection.java
- [ ] RandomSelection — jMetal: RandomSelection.java
- [ ] RankingAndCrowdingSelection — used in NSGA-II replacement

### Replacement
- [ ] MuPlusLambdaReplacement — jMetal: MuPlusLambdaReplacement.java
- [ ] NSGA2Replacement — NSGA-II environmental selection
- [ ] NSGAIIIReplacement — NSGA-III reference point based

## Phase 2: Algorithms
- [ ] NSGA-II (clean implementation — fix current skeleton)
- [ ] NSGA-III
- [ ] MOEA/D
- [ ] SPEA2
- [ ] SMPSO
- [ ] AGE-MOEA
- [ ] GDE3
- [ ] SMS-EMOA
- [ ] RVEA

## Phase 3: Problems
- [ ] ZDT1-6
- [ ] DTLZ1-7
- [ ] WFG1-9
- [ ] SCH, FON, KUR, POL (common benchmarks)

## Phase 4: Quality Indicators
- [ ] Hypervolume (WFG algorithm)
- [ ] IGD (Inverted Generational Distance)
- [ ] IGD+
- [ ] GD (Generational Distance)
- [ ] Epsilon
- [ ] Spread

## Phase 5: Integration
- [ ] JSON configuration → std::variant construction
- [ ] Island model (from original ea-cpp)
- [ ] Benchmarks vs jMetal (speed + energy)
- [ ] Documentation