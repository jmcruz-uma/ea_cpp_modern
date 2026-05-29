# Benchmarks

Guía reproducible de las comparativas del paper. Cada sección documenta
prerequisitos, comandos exactos e interpretación de resultados.

---

## Comparativa monobjetivo: ea_cpp_original vs ea_cpp_modern

Compara dos implementaciones del mismo (µ,λ)-ES sobre Rastrigin d=20:
- **ea_cpp_original** — AoS clásico, sin abstracciones modernas
- **ea_cpp_modern** — SoA + C++23, `Population<E>` templado

### Prerequisitos

**ea_cpp_original** — compilar una vez:
```bash
cd /home/alumno/ea_cpp_original
g++-14 -O3 -DNDEBUG -I. ea_main.cpp -std=c++23 -march=native -o ea_main_gcc -lm
```

**ea_cpp_modern** — compilar una vez:
```bash
cd /home/alumno/ea_cpp_modern
make
```

### Paso 1 — Ejecutar ambos sistemas

```bash
./scripts/run_so_comparison.sh OUTDIR MAXEVALS NUMRUNS
```

| Argumento  | Descripción                          | Valor paper | Valor rápido |
|------------|--------------------------------------|-------------|--------------|
| `OUTDIR`   | Directorio de salida                 | `results/so_paper` | `/tmp/test_so` |
| `MAXEVALS` | Evaluaciones por run                 | `10000000`  | `1000000`    |
| `NUMRUNS`  | Número de runs independientes        | `30`        | `5`          |

Ejemplo rápido (≈20 s):
```bash
./scripts/run_so_comparison.sh /tmp/test_so 1000000 5
```

Ejemplo paper (≈15 min):
```bash
./scripts/run_so_comparison.sh results/so_paper 10000000 30
```

El script genera en `OUTDIR`:
- `numeric.json`, `experiment.json` — configuración usada
- `rastrigin-20-1-stats.json` — salida de ea_cpp_original (formato original)
- `modern-rastrigin-20-1-stats.json` — salida de ea_cpp_modern (mismo formato)
- `original_stdout.txt`, `modern_stdout.txt` — salida de consola

### Paso 2 — Analizar resultados

```bash
python3 scripts/analyze_so_comparison.py OUTDIR [--csv] [--plot]
```

| Flag     | Efecto                                                        |
|----------|---------------------------------------------------------------|
| `--csv`  | Exporta `OUTDIR/comparison_results.csv` (una fila por run)   |
| `--plot` | Genera `OUTDIR/comparison_plot.png` (requiere matplotlib)    |

Ejemplo:
```bash
python3 scripts/analyze_so_comparison.py /tmp/test_so --csv --plot
```

El script imprime:
- Tiempo por run: media, mediana, Q1/Q3, speedup
- Fitness final: mediana, media, hits≈0 (convergencia al óptimo global)
- Fitness medio en hitos de evaluaciones (convergencia temporal)
- Test de Wilcoxon sobre fitness final (requiere scipy)

### Configuración del algoritmo

Ambos sistemas usan exactamente la misma configuración (generada por el script):

| Parámetro         | Valor  | Fuente                    |
|-------------------|--------|---------------------------|
| µ (padres)        | 100    | numeric.json → popsize    |
| λ (offspring)     | 99     | numeric.json → offspring  |
| Selección         | Torneo k=2 | numeric.json → selection |
| Cruce BLX-α       | prob=0.9, α=0.1 | numeric.json → variation[0] |
| Mutación Gaussian | rate=0.6321, step=1.0 | numeric.json → variation[1] |
| Reemplazo         | Comma (elitista si λ<µ) | numeric.json → replacement |
| Problema          | Rastrigin d=20, [-5.12, 5.12] | experiment.json |
| RNG               | mt19937_64, seed=1..NUMRUNS | ambos repos |

### Resultados de referencia (30 runs × 1M evals, WSL2)

| Métrica               | ea_cpp_original | ea_cpp_modern | Ratio  |
|-----------------------|-----------------|---------------|--------|
| Tiempo medio (s)      | 0.902           | 0.534         | **1.69×** |
| Tiempo mediana (s)    | 0.920           | 0.542         | **1.70×** |
| Fitness mediano       | 0.0             | 0.0           | =      |
| hits≈0 (de 30)        | 30/30           | 29/30         | ≈      |

> Nota: los tiempos en bare-metal Linux con RAPL diferirán de WSL2.
> Los resultados definitivos del paper deben medirse en bare-metal.

### Dependencias Python opcionales

```bash
pip install scipy      # test de Wilcoxon
pip install matplotlib # gráficas
```

---

## Próximas comparativas (pendientes)

- **Multiobjetivo SO vs MO**: impacto en rendimiento al añadir objetivos (NSGA-II sobre ZDT1–ZDT4)
- **ea_cpp_modern vs jMetal**: comparativa Java vs C++23 en NSGA-II, MOEA/D, SPEA2
- **Medición energética**: `perf stat -e power/energy-pkg/` en bare-metal Linux con RAPL
