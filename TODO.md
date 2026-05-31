# ea-cpp TODO

## v0.1.0 ✅ (2026-05-16) — First stable release

- [x] Core: Population SoA, Encoding, Concepts, Comparators, Random
- [x] Crossover: 21 operators (SBX, BLX-αβ, DE variants, PMX, ERX, Fuzzy, HUX, Composite, etc.)
- [x] Mutation: 17 operators (Polynomial, BitFlip, Swap, Inversion, CDG, GroupedPolynomial, etc.)
- [x] Selection: 8 operators (BinaryTournament, RankingAndCrowding, DifferentialEvolution, etc.)
- [x] Replacement: 4 operators (MuPlusLambda, MuCommaLambda, NSGAII, NSGAIII)
- [x] Algorithms: NSGA-II, NSGA-III, MOEA/D, SPEA2, SMS-EMOA, GDE3, SMPSO, AGE-MOEA, PAES, MOCell
- [x] Problems: ZDT1-6, DTLZ1-7, WFG1-9
- [x] Indicators: Hypervolume (WFG), IGD, GD, Spread
- [x] Crowding distance bug fixed (matches jMetal exactly)
- [x] SBX crossover matches jMetal exactly
- [x] Benchmark: IGD parity with jMetal, 12x faster

---

## v0.2.0 ✅ (2026-05-17) — Robustez, tests, correcciones

- [x] PolynomialMutation: aceptar probabilidad explícita
- [x] NSGA-II: validar pop_size par, auto-ajustar si impar
- [x] Crowding distance: O(n²) → O(n log n)
- [x] GD verificado
- [x] Spread indicator implementado
- [x] Unit tests: SBX, PolynomialMutation, BinaryTournament
- [x] Regression test ZDT1 + GitHub Actions CI

---

## v0.3.0 ✅ (2026-05-17) — Algoritmos adicionales

- [x] RandomSearch (baseline de comparación)
- [x] IBEA — Indicator-Based EA (HVI selection, equivalente jMetal)
- [x] MOEA/D-DE — DE variant of MOEA/D (Li & Zhang 2009)
- [x] SMS-EMOA verificación — fixed evaluate API call
- [x] BinaryTournament2 — strength-based (SPEA2)
- [x] RandomSelection — already existed in random.hpp
- [x] NeighborhoodSelection — decomposition-based (MOEA/D)

---

## v0.4.0 ✅ (2026-05-17) — Indicadores de calidad

- [x] Epsilon indicator (additive/multiplicative, Zitzler 2003)
- [x] R2 indicator (utility-based, Tchebycheff)
- [x] Average Hausdorff distance (max(GD, IGD))
- [x] Spread verificado

---

## v0.5.0 ✅ (2026-05-18) — Benchmarking automatizado + tooling

- [x] Script: todos los algoritmos × todos los problemas × 30 runs → CSV
- [x] Métricas por run: IGD, GD, Spread, Hypervolume, tiempo
- [x] CI: GitHub Actions workflow (compile + test on push/PR)
- [x] .clang-format — consistent code style (LLVM-based, C++23)
- [x] Doxygen — API docs generated at docs/html/
- [x] README badges: CI, clang-format

---

## v0.6.0 ✅ (2026-05-17) — Runtime configuration

- [x] Config struct con JSON parsing (nlohmann/json)
- [x] JSON config reader (from_file, from_json)
- [x] Ejemplo config_example.cpp con NSGAII/IBEA
- [x] Archivos de ejemplo: config_nsga2.json

## v1.0.0 ✅ (2026-05-17) — Paper-ready

- [x] Paper: borrador completo (paper.md)
- [x] Makefile: test, benchmark, clean, benchmark-full
- [x] README: badges, quick start, benchmarks, version history
- [x] CI: GitHub Actions regression test

## v1.1.0 ✅ (2026-05-20) — Island Model + Operator Gap + Bug Fixes

- [x] Island Model: RingTopology, MeshTopology, FullyConnectedTopology
- [x] Island Model: MigrationPolicy, MigrationOperator, Island, IslandModel
- [x] Selection: +7 operators (Environmental, SPEA2Strength, SUS, TournamentWOR, Weighted, SpatialSpread, BestSolution) = 15 total
- [x] Replacement: +4 operators (SPEA2, Crowding, MOEAD, SMPSO) = 8 total
- [x] Fix: NSGA-III heap-buffer-overflow (pop_size impar + selected bounds)
- [x] Fix: SMSEMOA, MOEA/D, MOEA/D-DE problem callable API (problem.evaluate → problem())
- [x] Extended benchmark suite (7 algos × 4 problems × N runs → CSV)
- [x] Total: 104 headers (was 76)

## v1.2.0 (WIP) — Doxygen + Paper polish

- [ ] Doxygen documentation for all 104 public headers
- [ ] Fix IBEA IGD (kappa normalization)
- [ ] Benchmark: 30 runs complete suite → paper tables
- [ ] GitHub Actions CI (needs workflow scope on token)

## Estado actual del proyecto

| Categoría | Implementados | jMetal | Cobertura |
|---|---|---|---|
| Crossover | 21 | 25 | 84% |
| Mutation | 17 | 20 | 85% |
| Selection | 15 | 14 | 107% ✅ |
| Replacement | 8 | 13 | 62% |
| Algorithms | 14 | 15+ | 93% |
| Problems | 22 | 40+ | 55% |
| Indicators | 7+ | 7 | 100% ✅ |
| Island Model | 5 | — | NEW ✅ |
| **Total** | **104 headers** | | |