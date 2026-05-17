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
- [x] RVEA — Reference Vector-guided EA (many-objective)
- [x] IBEA — Indicator-Based EA (Hypervolume/Epsilon selection)
- [ ] MOEA/D-DE — pendiente
- [ ] SMS-EMOA verificación — pendiente
- [ ] BinaryTournament2, RandomSelection, NeighborhoodSelection, NaryTournament

---

## v0.4.0 ✅ (2026-05-17) — Indicadores de calidad

- [x] Epsilon indicator (additive/multiplicative, Zitzler 2003)
- [x] R2 indicator (utility-based, Tchebycheff)
- [x] Average Hausdorff distance (max(GD, IGD))
- [x] Spread verificado

---

## v0.5.0 — Benchmarking automatizado

- [ ] Script: todos los algoritmos × todos los problemas × 30 runs → CSV
- [ ] Métricas por run: IGD, GD, Spread, Hypervolume, tiempo
- [ ] Comparación automática vs jMetal
- [ ] Tabla LaTeX/markdown generada
- [ ] CI: ejecutar benchmarks en push

## v0.6.0 — Runtime configuration

- [ ] Config struct con std::variant
- [ ] JSON config reader
- [ ] Factory: JSON → algorithm
- [ ] CLI para sweeps de parámetros

## v1.0.0 — Paper-ready

- [ ] Island model
- [ ] Consumo eléctrico (RAPL/perf)
- [ ] Escalabilidad: 1K, 5K, 10K individuos
- [ ] Problemas con restricciones
- [ ] Multi-objective >2 objetivos
- [ ] Documentación API (Doxygen)
- [ ] README completo
