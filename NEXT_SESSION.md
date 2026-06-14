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
| WFG (9) | ✅ | ❌ | ❌ pendiente |
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

---

## PRÓXIMA TAREA: WFG1-9 (auditoría + benchmark)

### Contexto general WFG

El WFG (Walking Fish Group) es escalable en M objetivos y n_vars = k + l variables:
- k: parámetros de posición (relacionados con el PS)
- l: parámetros de distancia
- Bounds: variable i ∈ [0, 2*(i+1)]

Archivos clave:
- `include/ea/problem/wfg.hpp` (995 líneas) — portado, **SIN AUDITAR**
- `/home/alumno/jMetal/jmetal-problem/src/main/java/org/uma/jmetal/problem/multiobjective/wfg/WFG{1-9}.java` — referencia canónica
- `/home/alumno/ea_cpp_jmetal/include/ea/problem/wfg.hpp` — fork secundario (NO de confianza ciega)

### Configuración de benchmark (jMetal default)

jMetal usa `new WFG1()` → k=2, l=4, M=2, n_vars=6.
DefaultWFGSettings: `numberOfPositionParameters=2, numberOfDistanceParameters=4, numberOfObjectives=2`.

**ATENCIÓN**: ea_cpp_modern WFG usa defaults distintos: `WFG1(int m=3, int k=-1, int l=20)`.
Para el benchmark, instanciar con los parámetros jMetal: `WFG1(/*m=*/2, /*k=*/2, /*l=*/4)`.

Algoritmo benchmark: **NSGA-II** (pop=100, max_evals=25000, SBX+PM), igual que WFGStudy.java.

### Frentes de referencia

jMetal tiene CSVs precomputados para todos WFG1-9, 2D y 3D:
```
/home/alumno/jMetal/resources/referenceFrontsCSV/WFG{1-9}.2D.csv
```
Formato: `f1,f2` (comma-separated), ~1113 puntos por archivo.
**Usar estos CSVs como frente de referencia** — no hay frentes analíticos simples para WFG.

### Protocolo de auditoría

Orden obligatorio: jMetal canónico → ea_cpp_jmetal (solo como pista) → ea_cpp_modern.

Archivos a auditar en jMetal (un archivo por problema):
```
WFG.java       — pipeline evaluate(): normalize → t1..tN → calculate_x → shape
WFG1.java      — t1=bFlat+bPoly, t2=rSum, t3=rSum; shape: convex+mixed
WFG2.java      — t3=rNonsep; shape: convex+disconnected
WFG3.java      — t3=rNonsep; shape: linear (degenerate)
WFG4.java      — t1=sMulti, t2=rSum, t3=rSum; shape: spherical
WFG5.java      — t1=sDecept, ...; shape: spherical
WFG6.java      — t3=rNonsep; shape: spherical
WFG7.java      — t1=bParam, t2=rSum; shape: spherical
WFG8.java      — t1=bParam (backward); shape: spherical
WFG9.java      — t1=sDecept+bParam, t3=rNonsep; shape: spherical
```

Puntos especialmente críticos (fuente de bugs en otros portados):
1. `correct_to_01`: tolerancia epsilon (Java usa float, C++ usa double — revisar precisión)
2. `r_nonsep`: reduce con floor y mod (fácil de meter off-by-one)
3. `b_param`: usa aleatoriedad interna en WFG7/WFG8 (¿cómo la maneja ea_cpp_modern?)
4. `calculate_x`: fórmula x[i] = 2*(i+1)*y[i] (bounds-dependent)
5. `s` array: `s[i] = 2*(i+1)` — escalado de objetivos
6. Orden de transformaciones por problema: verificar t1→t2→t3 uno a uno

### Archivos a crear (tras auditoría)

1. `examples/wfg_runner.cpp` — NSGA-II sobre WFG1-9 (m=2, k=2, l=4)
2. `benchmarks/jmetal/JMetalWFGBenchmark.java` — contraparte Java
3. `scripts/run_wfg_comparison.sh` — patrón: run_uf_comparison.sh

### Java imports para WFG

```java
import org.uma.jmetal.problem.multiobjective.wfg.WFG1; // ... hasta WFG9
// Constructor: new WFG1() — usa DefaultWFGSettings: k=2, l=4, M=2
```

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
