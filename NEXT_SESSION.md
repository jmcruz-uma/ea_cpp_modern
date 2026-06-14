# Prompt de continuación — siguiente sesión ea_cpp_modern

## Contexto del proyecto

Paper de rendimiento: **ea_cpp_modern** (C++23, SoA, zero-overhead) vs **jMetal 7.4** (Java).
Repos:
- `/home/alumno/ea_cpp_modern` — nuestro framework C++23, branch main
- `/home/alumno/jMetal` — referencia canónica Java 7.4
- `/home/alumno/ea_cpp_jmetal` — fork v1.3 de un compañero; útil como punto de partida pero NO de confianza ciega (tiene bugs; siempre contrastar con jMetal)

**Compilador obligatorio**: `g++-14 -std=c++23 -O3 -march=native -DNDEBUG`
**JVM benchmark**: `-server -Xms512m -Xmx2g -XX:+UseG1GC -XX:+AlwaysPreTouch -XX:+DisableExplicitGC -Djava.awt.headless=true`

---

## Estado completo del proyecto

### Algoritmos MO — TODOS completados ✅

| Algoritmo | Speedup mediano | Wilcoxon IGD |
|-----------|:---------------:|:------------:|
| NSGA-II | 3.57× | todos ns |
| SMPSO | 3.87× | todos ns |
| SPEA2 | 2.17× | todos ns |
| IBEA | 12.41× | todos ns |
| MOEA/D-DE | 7.87× | todos ns |
| NSGA-III | 2.93× | todos ns |
| SMS-EMOA | 3.92× | todos ns |
| GDE3 | 13.05× | todos ns |
| AGE-MOEA | 5.02× | 3/4 ns |
| MO-Cell | 5.12× | 3/4 ns |
| PAES | 4.09× | C++ mejor ZDT1-2 p<0.05 |

### Problemas benchmark — estado

| Familia | Portado | Auditado | Benchmarked |
|---------|:-------:|:--------:|:-----------:|
| ZDT (6) | ✅ | parcial | ✅ ZDT1-4 (todos los algoritmos) |
| DTLZ (7) | ✅ | ❌ | ✅ DTLZ1-4 (NSGA-III) |
| WFG (9) | ✅ | ✅ | ✅ WFG1-9 (NSGA-II, m=2 k=2 l=4) |
| UF (10) | ✅ | ✅ | ✅ UF1-7 / ❌ UF8-10 |
| LZ09 (9) | ✅ | ✅ | ✅ LZ09F1-F9 (NSGA-II+NSGA-III) |

### UF1-7 — COMPLETO ✅ (sesión 2026-06-14)

NSGA-II, pop=100, max_evals=25000, 30 runs, 30 vars cada uno:

| Problema | C++ med (ms) | Java med (ms) | Speedup | IGD p-val |
|----------|:---:|:---:|:---:|:---:|
| UF1 | 161 | 510 | 3.16× | 0.57 ns |
| UF2 | 165 | 601 | 3.65× | 0.89 ns |
| UF3 | 185 | 579 | 3.14× | 0.33 ns |
| UF4 | 152 | 560 | 3.67× | 0.97 ns |
| UF5 | 181 | 616 | 3.41× | 0.58 ns |
| UF6 | 186 | 623 | 3.35× | 1.00 ns |
| UF7 | 164 | 590 | 3.61× | 0.58 ns |
| **media** | | | **3.43×** | todos ns |

Infraestructura: `examples/uf_runner.cpp`, `benchmarks/jmetal/JMetalUFBenchmark.java`,
`scripts/run_uf_comparison.sh`, `results/uf_comparison/`.

Nota UF6: σ C++=98ms por un único outlier (run 4 = seed 45) en el que NSGA-II encontró
ambos segmentos del frente discontinuo → más no dominadas → NDS más caro. Sin ese run,
σ=7ms (normal). El outlier tiene el MEJOR IGD — comportamiento algorítmico correcto.

### LZ09F1-F9 — COMPLETO ✅ (sesión 2026-06-14)

LZ09F1-F5,F7-F9: NSGA-II, pop=100, max_evals=25000, 30 runs.
LZ09F6: NSGA-III, divisions=12, pop=92, max_evals=25000, 30 runs.

| Problema | C++ med (ms) | Java med (ms) | Speedup | IGD p-val |
|----------|:---:|:---:|:---:|:---:|
| LZ09F1 | 140 | 648 | 4.61× | 0.70 ns |
| LZ09F2 | 172 | 712 | 4.13× | 0.45 ns |
| LZ09F3 | 166 | 657 | 3.96× | 0.64 ns |
| LZ09F4 | 166 | 549 | 3.31× | 0.13 ns |
| LZ09F5 | 177 | 563 | 3.19× | 0.18 ns |
| LZ09F6 | 198 | 523 | 2.64× | 0.031 * |
| LZ09F7 | 150 | 492 | 3.29× | 0.84 ns |
| LZ09F8 | 161 | 491 | 3.04× | 0.98 ns |
| LZ09F9 | 191 | 564 | 2.96× | 0.92 ns |
| **media** | | | **3.46×** | 8/9 ns |

Nota LZ09F6 (p=0.031): Java NSGA-III obtiene IGD mediana ligeramente mejor (0.118 vs 0.140).
Ambas distribuciones tienen alta varianza (min ~0.10, max ~0.54) — la diferencia absoluta es
pequeña (0.022). El speedup 2.64× es el más bajo de la familia LZ09. Causa probable:
implementación NSGA-III de jMetal más madura en niching/normalización para este problema 3-obj.
Para el paper: speedup real, calidad comparable (mismo orden de magnitud).

Infraestructura: `examples/lz09_runner.cpp`, `benchmarks/jmetal/JMetalLZ09Benchmark.java`,
`scripts/run_lz09_comparison.sh`, `results/lz09_comparison/`.

### WFG1-9 — COMPLETO ✅ (sesión 2026-06-14)

NSGA-II, pop=100, max_evals=25000, 30 runs, k=2 l=4 M=2 dim=6 (jMetal DefaultWFGSettings).
Frentes de referencia: CSVs precomputados de jMetal (~119-2600 pts por problema).

**Auditoría: sin bugs.** wfg.hpp es fiel a jMetal canónico en todos los problemas.
Únicas diferencias: (1) epsilon normalise 1e-10 C++ vs 1e-7 Java — irrelevante. (2) double C++ vs float Java — efecto en WFG8 (ver notas).

| Problema | C++ med (ms) | Java med (ms) | Speedup | IGD p-val |
|----------|:---:|:---:|:---:|:---:|
| WFG1 | 125 | 417 | 3.33× | 0.009 ** |
| WFG2 | 123 | 392 | 3.19× | 0.70 ns |
| WFG3 | 108 | 398 | 3.69× | 0.36 ns |
| WFG4 | 109 | 413 | 3.80× | 0.61 ns |
| WFG5 | 104 | 416 | 4.01× | 0.21 ns |
| WFG6 | 113 | 412 | 3.65× | 0.82 ns |
| WFG7 | 111 | 416 | 3.76× | 0.43 ns |
| WFG8 | 134 | 443 | 3.30× | <0.0001 *** |
| WFG9 | 123 | 434 | 3.52× | 0.18 ns |
| **mediana** | | | **3.65×** | 7/9 ns |

**Nota WFG1** (p=0.009, Java IGD mejor — med 0.565 vs C++ 0.845):
WFG1 tiene el bPoly(y, 0.02) sobre todas las variables → landscape muy plana y alta varianza de IGD
(σ=0.19 C++, σ=0.33 Java). La diferencia podría ser seed effect (C++ seeds 42..71, Java 45..74),
o efecto de float vs double en bPoly. Para el paper: reportar honestamente; speedup real 3.33×.

**Nota WFG8** (p<0.0001, C++ IGD mejor — med 0.060 vs Java 0.125; W=0, 30/30 wins):
bParam sobre variables de distancia usa rSum de variables previas. Con float en Java, la pequeña
imprecisión en el parámetro u del bParam cambia el exponente (rango 0.02–50), alterando el
landscape. C++ double produce función más diferenciada → mejor convergencia de NSGA-II.
Cat-2 finding: precisión afecta la calidad de solución, no solo la exactitud numérica.

Infraestructura: `examples/wfg_runner.cpp`, `benchmarks/jmetal/JMetalWFGBenchmark.java`,
`scripts/run_wfg_comparison.sh`, `results/wfg_comparison/`.

---

## PRÓXIMA TAREA: bare-metal validation + escritura del paper

**WFG1-9 completado en sesión 2026-06-14. Ver sección de resultados arriba.**

Todos los benchmarks de problemas están completos. La siguiente tarea es bare-metal validation.

---

## Pendiente (menor prioridad)

### UF8-10 (opcional)

- 3 objetivos, 30 vars; NSGA-III (pop=92, divisions=12).
- LZ09F6 ya cubre el caso 3-obj → baja prioridad.

### Modelo distribuido island

`taps_transport.hpp` de ea_cpp_jmetal es incorrecto. Ver memoria `feedback_distributed_transport.md`.
Pendiente de decisión de middleware. Baja prioridad para el paper.

### Bare-metal validation

Repetir todos los benchmarks fuera de WSL2. Comandos:
```bash
sudo cpupower frequency-set -g performance
echo 1 | sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo
SKIP_JMETAL_BUILD=1 ./scripts/run_uf_comparison.sh results/uf_baremetal 30
SKIP_JMETAL_BUILD=1 ./scripts/run_lz09_comparison.sh results/lz09_baremetal 30
# ... (igual para todos los algoritmos)
```
Bloquea la escritura del paper — hacer antes de redactar.

### Escritura del paper

Ver `project_paper_framing.md` en memoria: Cat-1/2/3, tabla speedups,
bugs reportables (GDE3 doble NDS, PAES density). Pendiente hasta bare-metal.

---

## Encuadre académico

**Enfoque A** (principal): speedup vs jMetal as-shipped.
**Enfoque B** (lower bound): speedup vs jMetal con bugs corregidos.
- GDE3: sin doble NDS → Java ~553ms → speedup ~10× en vez de 13×
- PAES: con O(1) incremental → Java ~45ms → speedup ~2× en vez de 4×

**Tres categorías de ventaja**:
- Cat-1: boxing, GC, JIT → aplica a todos los algoritmos
- Cat-2: C++ más fiel al paper original (PAES density, GDE3 NDS, MO-Cell coin_flip)
- Cat-3: optimización natural en C++ pero posible en Java (PAES O(1) grid)
