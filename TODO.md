# TODO.md — ea-cpp Implementation Roadmap

## ✅ Phase 0: Foundation
- [x] CMakeLists.txt with C++23, GTest, benchmarks
- [x] Population SoA (population.hpp)
- [x] Encoding enum (encoding.hpp)
- [x] Concepts (concepts.hpp)
- [x] Comparator (comparator.hpp) — dominance, crowding distance, fast NDS
- [x] Random (random.hpp) — thread-safe, batch generation
- [x] ARCHITECTURE.md

## ✅ Phase 1: Core Operators
### Crossover (continuous) — 13 total
- [x] SBXCrossover
- [x] BLXAlphaCrossover
- [x] BLXAlphaBetaCrossover
- [x] ArithmeticCrossover
- [x] WholeArithmeticCrossover
- [x] LaplaceCrossover
- [x] DifferentialEvolutionCrossover (with DE variants: RAND/BEST/RAND-TO-BEST + BIN/EXP)
- [x] ParentCentricCrossover (UNDX)

### Crossover (permutation) — 6 total
- [x] PMXCrossover
- [x] CycleCrossover
- [x] EdgeRecombinationCrossover
- [x] OXDCrossover
- [x] PositionBasedCrossover
- [x] UniformCrossover

### Mutation — 11 total
- [x] PolynomialMutation
- [x] BitFlipMutation
- [x] SwapMutation
- [x] InversionMutation
- [x] InsertMutation
- [x] ScrambleMutation
- [x] UniformMutation
- [x] NonUniformMutation
- [x] PowerLawMutation
- [x] LevyFlightMutation
- [x] SimpleRandomMutation

### Selection — 6 total
- [x] BinaryTournamentSelection
- [x] NaryTournamentSelection
- [x] RandomSelection
- [x] RankingAndCrowdingSelection
- [x] BestSolutionSelection
- [x] NaryRandomSelection

### Replacement — 3 total
- [x] MuPlusLambdaReplacement
- [x] MuCommaLambdaReplacement
- [x] NSGAIIReplacement (environmental selection)

## ✅ Phase 2: Algorithms
- [x] NSGA-II — functional, tested with ZDT1 (IGD=0.034)

## ✅ Phase 3: Problems
- [x] ZDT1-6 (complete suite)
- [x] DTLZ1-7 (complete suite)

## ✅ Phase 4: Quality Indicators
- [x] Hypervolume (2D efficient, 3D+ naive)
- [x] IGD (Inverted Generational Distance)
- [x] GD (Generational Distance)
- [x] Spread (Δ)

## 🔲 Phase 5: More Algorithms
- [ ] NSGA-III
- [ ] MOEA/D
- [ ] SPEA2
- [ ] SMPSO
- [ ] AGE-MOEA
- [ ] GDE3
- [ ] SMS-EMOA
- [ ] RVEA
- [ ] PAES
- [ ] MOCell
- [ ] RandomSearch

## 🔲 Phase 6: Remaining Operators from jMetal
### Crossover
- [ ] FuzzyRecombinationCrossover
- [ ] IntegerSBXCrossover
- [ ] HUXCrossover (binary)
- [ ] CompositeCrossover
- [ ] NullCrossover

### Mutation
- [ ] CDGMutation
- [ ] GroupedPolynomialMutation
- [ ] GroupedAndLinkedPolynomialMutation
- [ ] LinkedPolynomialMutation
- [ ] IntegerPolynomialMutation
- [ ] IntegerSimpleRandomMutation
- [ ] CharSequenceRandomMutation
- [ ] CompositeMutation
- [ ] NullMutation

### Selection
- [ ] DifferentialEvolutionSelection
- [ ] SpatialSpreadDeviationSelection
- [ ] RankingAndPreferenceSelection
- [ ] RankingAndDirScoreSelection
- [ ] BoltzmannSelection
- [ ] TruncationSelection
- [ ] StochasticUniversalSampling
- [ ] NeighborhoodSelection
- [ ] PopulationAndNeighborhoodSelection

### Replacement
- [ ] NSGAIIIReplacement (reference point based)
- [ ] AGEMOEAReplacement
- [ ] SMSEMOAReplacement
- [ ] RVEAReplacement
- [ ] PairwiseReplacement
- [ ] RandomReplacement
- [ ] TournamentReplacement
- [ ] SingleSolutionReplacement

## 🔲 Phase 7: Infrastructure
- [ ] JSON configuration → std::variant construction
- [ ] Island model (from original ea-cpp)
- [ ] Benchmarks vs jMetal (speed + energy)
- [ ] WFG benchmark problems
- [ ] More quality indicators (R2, Epsilon, AverageHausdorff)
- [ ] Documentation