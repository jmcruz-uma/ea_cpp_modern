#!/usr/bin/env bash
# =============================================================================
# run_so_comparison.sh
# Comparativa monobjetivo: ea_cpp_original vs ea_cpp_modern (MuLambdaES)
#
# Uso:
#   ./scripts/run_so_comparison.sh [OUTDIR] [MAXEVALS] [NUMRUNS]
#
# Argumentos opcionales (con valores por defecto):
#   OUTDIR    directorio de salida    (default: results/so_comparison)
#   MAXEVALS  evaluaciones por run   (default: 1000000)
#   NUMRUNS   número de runs         (default: 30)
#
# Ejemplo rápido (1 run, 100K evals):
#   ./scripts/run_so_comparison.sh /tmp/test 100000 1
#
# Ejemplo paper (30 runs, 10M evals):
#   ./scripts/run_so_comparison.sh results/so_paper 10000000 30
#
# Prerequisitos:
#   - ea_cpp_original compilado:  cd /home/alumno/ea_cpp_original && \
#       g++-14 -O3 -DNDEBUG -I. ea_main.cpp -std=c++23 -march=native -o ea_main_gcc -lm
#   - ea_cpp_modern compilado:    cd /home/alumno/ea_cpp_modern && make
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODERN_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ORIGINAL_ROOT="/home/alumno/ea_cpp_original"

OUTDIR="${1:-${MODERN_ROOT}/results/so_comparison}"
MAXEVALS="${2:-1000000}"
NUMRUNS="${3:-30}"

ORIGINAL_BIN="${ORIGINAL_ROOT}/ea_main_gcc"
MODERN_BIN="${MODERN_ROOT}/build/mu_lambda_runner"

# ── Validaciones ──────────────────────────────────────────────────────────────
if [[ ! -x "$ORIGINAL_BIN" ]]; then
    echo "ERROR: no encontrado $ORIGINAL_BIN"
    echo "Compila con: cd $ORIGINAL_ROOT && g++-14 -O3 -DNDEBUG -I. ea_main.cpp -std=c++23 -march=native -o ea_main_gcc -lm"
    exit 1
fi
if [[ ! -x "$MODERN_BIN" ]]; then
    echo "ERROR: no encontrado $MODERN_BIN"
    echo "Compila con: cd $MODERN_ROOT && make"
    exit 1
fi

mkdir -p "$OUTDIR"

# ── Generar JSONs de configuración ────────────────────────────────────────────
cat > "${OUTDIR}/numeric.json" << EOF
{
    "numruns" : ${NUMRUNS},
    "seed" : 1,
    "islands" : [{
        "numislands" : 1,
        "popsize" : 100,
        "offspring" : 99,
        "maxevals" : ${MAXEVALS},
        "initialization" : {"name": "randomvector"},
        "selection"      : {"name": "tournament", "parameters": ["2"]},
        "variation"      : [
            {"name": "blx",      "parameters": ["0.9", "0.1"]},
            {"name": "gaussian", "parameters": ["0.6321", "1"]}
        ],
        "replacement" : {"name": "comma"}
    }]
}
EOF

cat > "${OUTDIR}/experiment.json" << EOF
{
    "numruns" : ${NUMRUNS},
    "problem" : "rastrigin",
    "numvars" : 20,
    "range"   : 5.12,
    "numbits" : 0
}
EOF

echo "============================================================"
echo " Comparativa monobjetivo: ea_cpp_original vs ea_cpp_modern"
echo "============================================================"
echo " Problema  : Rastrigin d=20, rango=[-5.12, 5.12]"
echo " Algoritmo : (µ=100, λ=99)-ES, BLX-α (prob=0.9, α=0.1)"
echo "             Gaussian (rate=0.6321, step=1.0), torneo k=2"
echo " Budget    : ${MAXEVALS} evaluaciones × ${NUMRUNS} runs"
echo " Salida    : ${OUTDIR}"
echo "============================================================"
echo ""

# ── Ejecutar ea_cpp_original ─────────────────────────────────────────────────
echo "[1/2] Ejecutando ea_cpp_original..."
rm -f "${OUTDIR}/lastseed.txt"
pushd "${OUTDIR}" > /dev/null
START=$(date +%s%N)
"${ORIGINAL_BIN}" experiment.json > original_stdout.txt 2>&1
END=$(date +%s%N)
popd > /dev/null
ORIG_WALL=$(echo "scale=3; ($END - $START) / 1000000000" | bc)
# El original escribe rastrigin-20-1-stats.json en OUTDIR
echo "   Completado en ${ORIG_WALL}s → ${OUTDIR}/rastrigin-20-1-stats.json"
echo ""

# ── Ejecutar ea_cpp_modern ───────────────────────────────────────────────────
echo "[2/2] Ejecutando ea_cpp_modern (MuLambdaES)..."
rm -f "${OUTDIR}/lastseed.txt"
START=$(date +%s%N)
"${MODERN_BIN}" \
    "${OUTDIR}/numeric.json" \
    "${OUTDIR}/experiment.json" \
    "${OUTDIR}/modern-rastrigin-20-1-stats.json" \
    > "${OUTDIR}/modern_stdout.txt" 2>&1
END=$(date +%s%N)
MODERN_WALL=$(echo "scale=3; ($END - $START) / 1000000000" | bc)
echo "   Completado en ${MODERN_WALL}s → ${OUTDIR}/modern-rastrigin-20-1-stats.json"
echo ""

# ── Resumen rápido ────────────────────────────────────────────────────────────
echo "============================================================"
echo " Tiempos de pared (${NUMRUNS} runs × ${MAXEVALS} evals)"
echo "   ea_cpp_original : ${ORIG_WALL}s"
echo "   ea_cpp_modern   : ${MODERN_WALL}s"
SPEEDUP=$(echo "scale=3; ${ORIG_WALL} / ${MODERN_WALL}" | bc)
echo "   Speedup         : ${SPEEDUP}×"
echo "============================================================"
echo ""
echo "Análisis detallado:"
echo "   python3 ${SCRIPT_DIR}/analyze_so_comparison.py ${OUTDIR}"
