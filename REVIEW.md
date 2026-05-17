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
│   ├── algorithm/       # 13 algoritmos (NSGAII, RVEA, IBEA, etc.)
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

## Contacto

- Issues: GitHub Issues
- Email: jmcruz@uma.es (revisor)
- Autor: Chema (elCano Research)
