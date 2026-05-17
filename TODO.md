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

## v0.2.0 — Robustez y corrección

- [x] PolynomialMutation: aceptar probabilidad explícita (ya era configurable, documentado)
- [x] NSGA-II: validar que pop_size sea par, auto-ajustar a par+1 si impar
- [x] Crowding distance: O(n²) → O(n log n) con unordered_map
- [x] Verificar GD — funciona correctamente (0 para match perfecto)
- [ ] Spread indicator — no existe, crearlo
- [ ] Tests unitarios para cada operador (crossover, mutation, selection)
- [ ] Test de regresión: IGD < 0.01 en ZDT1 con config estándar

## v0.3.0 — Algoritmos que faltan

- [ ] RandomSearch (baseline de comparación)
- [ ] RVEA (usa ReferencePoints existentes)
- [ ] IBEA (usa Hypervolume existente)
- [ ] MOEA/D-DE (variante con DifferentialEvolution)
- [ ] SMS-EMOA: verificar contra jMetal
- [ ] BinaryTournament2 (modificado, permite repetidos)
- [ ] RandomSelection (para RandomSearch)
- [ ] NeighborhoodSelection (para MOEA/D)
- [ ] NaryTournament (generalización de binary)

## v0.4.0 — Indicadores de calidad

- [ ] Epsilon indicator (additive, Zitzler 2003)
- [ ] R2 indicator (IGD-based, usa Tchebycheff existente)
- [ ] Average Hausdorff distance (max(GD, IGD))
- [ ] GD y Spread: verificación contra jMetal

## v0.5.0 — Benchmarking automatizado

- [ ] Script: todos los algoritmos × todos los problemas × 30 runs → CSV
- [ ] Métricas por run: IGD, GD, Spread, Hypervolume, tiempo
- [ ] Comparación automática vs jMetal (misma config, mismas evaluaciones)
- [ ] Tabla LaTeX/markdown generada lista para paper
- [ ] CI: ejecutar benchmarks en push a main (GitHub Actions)

## v0.6.0 — Runtime configuration

- [ ] Config struct con std::variant para crossover/mutation/selection
- [ ] JSON config reader (nlohmann/json o similar)
- [ ] Factory: JSON → algorithm configurado sin recompilar
- [ ] CLI para sweeps de parámetros
- [ ] Permite benchmarking masivo sin recompilar

## v1.0.0 — Paper-ready

- [ ] Island model (migración asíncrona, topología ring/torus)
- [ ] Consumo eléctrico: medir joules con RAPL/perf
- [ ] Escalabilidad: benchmarks 1K, 5K, 10K individuos
- [ ] Problemas con restricciones: DTLZ constrained, SRN, TNK
- [ ] Multi-objective >2 objetivos: DTLZ con 3 y 5 objetivos
- [ ] Documentación API completa (Doxygen)
- [ ] README con badges, instalación, ejemplos