# Guía de Revisión para jmcruz@uma.es

## Estado del proyecto

**ea-cpp v1.0.0 — Paper-ready release**
- Repo: https://github.com/elCanosail/ea-cpp
- Paper: `paper.md` (markdown legible) y `paper.tex` (LaTeX)

## Estructura del repo

```
ea-cpp/
├── include/ea/          # 76 headers, código fuente
│   ├── core/            # Population, Config
│   ├── algorithm/       # 13 algoritmos (NSGAII, IBEA, MOEAD, etc.)
│   ├── operator/        # Crossover, Mutation, Selection, Replacement
│   ├── problem/         # ZDT, DTLZ, WFG families
│   ├── indicator/       # IGD, GD, Spread, Hypervolume, Epsilon, R2, Hausdorff
│   └── util/            # ReferencePoint, Aggregation, Random
├── tests/               # Tests unitarios + regresión
├── examples/            # Benchmarks + config demos
├── scripts/             # Benchmark suite + energy estimation
├── paper.md             # Paper en markdown
├── paper.tex            # Paper en LaTeX
├── Makefile             # Build automation
└── README.md            # Documentación
```

## Cómo empezar

```bash
# Clonar
git clone https://github.com/elCanosail/ea-cpp.git
cd ea-cpp

# Requisitos: g++-14+ con C++23
make test          # Ejecutar tests
make benchmark     # Benchmark rápido (NSGAII ZDT1)
make benchmark-full # Suite completa (30 runs)
```

## Qué revisar

1. **Paper**: `paper.md` — estructura, resultados, referencias
2. **Código**: `include/ea/algorithm/` — implementaciones de algoritmos
3. **Benchmarks**: Ejecutar `scripts/benchmark_suite.cpp` en tu hardware
4. **Paridad jMetal**: Verificar IGD con misma configuración

## Notas técnicas

- **C++23 features**: concepts, deducing-this, SoA Population
- **Compilación**: `g++-14 -std=c++23 -O2`
- **Sin dependencias**: Solo STL + nlohmann/json (para config)

## Resultados de benchmark (WSL2, 30 runs, pop=100, max_evals=25000)

| Algoritmo  | ZDT1 speedup | ZDT2 speedup | ZDT3 speedup | ZDT4 speedup | Mediana |
|------------|-------------|-------------|-------------|-------------|---------|
| NSGA-II    | 3.57×       | —           | —           | —           | 3.57×   |
| SMPSO      | 3.87×       | —           | —           | —           | 3.87×   |
| SPEA2      | 2.17×       | —           | —           | —           | 2.17×   |
| IBEA       | 12.41×      | —           | —           | —           | 12.41×  |
| MOEA/D-DE  | 8.37×       | 8.01×       | 7.43×       | 7.65×       | 7.87×   |

Todos los IGD son estadísticamente equivalentes a jMetal (Mann-Whitney p ≥ 0.05).

## Contacto

- Issues: GitHub Issues
- Email: jmcruz@uma.es (revisor)
- Autor: Chema (elCano Research)
