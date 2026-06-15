# results/

This directory holds benchmark output produced by the scripts in `../scripts/`.

CSV files and logs are **not tracked in git** (see `.gitignore`). Each subdirectory
is created by its corresponding script.

## Directory layout

```
results/
├── nsga2/            ← run_nsga2_comparison.sh
├── smpso/            ← run_smpso_comparison.sh
├── spea2/            ← run_spea2_comparison.sh
├── ibea/             ← run_ibea_comparison.sh
├── moead/            ← run_moead_comparison.sh
├── nsga3/            ← run_nsga3_comparison.sh
├── smsemoa/          ← run_smsemoa_comparison.sh
├── gde3/             ← run_gde3_comparison.sh
├── agemoea/          ← run_agemoea_comparison.sh
├── mocell/           ← run_mocell_comparison.sh
├── paes/             ← run_paes_comparison.sh
├── uf/               ← run_uf_comparison.sh
├── lz09/             ← run_lz09_comparison.sh
├── wfg/              ← run_wfg_comparison.sh
└── so_paper/         ← run_so_comparison.sh
```

## Files generated per run

Each script writes into its output directory:

| File | Contents |
|---|---|
| `results_eacpp.csv` | C++ timings and IGD per run |
| `results_jmetal.csv` | Java timings and IGD per run |
| `environment.md` | Hardware, OS, compiler, JVM at time of execution |

## Reproducing all results

See [BENCHMARKS.md](../BENCHMARKS.md) for the full command sequence. Quick start:

```bash
# Build jMetal once
git clone https://github.com/jMetal/jMetal.git $HOME/jMetal
cd $HOME/jMetal && mvn clean install -DskipTests -q && cd -

# Run all MO comparisons (30 runs each)
./scripts/run_nsga2_comparison.sh results/nsga2 30
SKIP_JMETAL_BUILD=1 ./scripts/run_smpso_comparison.sh  results/smpso  30
# ... (see BENCHMARKS.md for all algorithms)
```
