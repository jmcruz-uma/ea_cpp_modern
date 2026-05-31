#!/usr/bin/env bash
# =============================================================================
# run_nsga3_comparison.sh
# Comparativa many-objective: ea_cpp_modern (NSGA-III, C++) vs jMetal 7.4 (Java)
# Problemas: DTLZ1, DTLZ2, DTLZ3, DTLZ4 (3 objetivos)
#
# Uso:
#   ./scripts/run_nsga3_comparison.sh [OUTDIR] [NUMRUNS]
#
# Argumentos opcionales:
#   OUTDIR    directorio de salida    (default: results/nsga3_comparison)
#   NUMRUNS   número de runs medidos  (default: 30)
#
# Variables de entorno:
#   SKIP_JAVA=1              omite la parte jMetal (solo corre C++)
#   SKIP_JMETAL_BUILD=1      asume que jMetal ya está compilado
#
# Ejemplo rápido (3 runs, solo C++):
#   SKIP_JAVA=1 ./scripts/run_nsga3_comparison.sh /tmp/nsga3_test 3
#
# Ejemplo paper (30 runs):
#   ./scripts/run_nsga3_comparison.sh results/nsga3_paper 30
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODERN_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
JMETAL_ROOT="/home/alumno/jMetal"
BENCH_DIR="${MODERN_ROOT}/benchmarks/jmetal"

OUTDIR="${1:-${MODERN_ROOT}/results/nsga3_comparison}"
NUMRUNS="${2:-30}"

CPP_BIN="${MODERN_ROOT}/build/nsga3_runner"
JAVA_CSV="${OUTDIR}/results_jmetal.csv"
CPP_CSV="${OUTDIR}/results_eacpp.csv"

SKIP_JAVA="${SKIP_JAVA:-0}"

# ── Funciones auxiliares ──────────────────────────────────────────────────────

check_prereqs() {
    local missing=0
    if ! command -v g++-14 &>/dev/null; then
        echo "ERROR: g++-14 no encontrado. Instalar con: sudo apt-get install g++-14"
        missing=1
    fi
    if [[ "$SKIP_JAVA" == "0" ]]; then
        if ! command -v java &>/dev/null; then
            echo "ERROR: java no encontrado. Instalar con: sudo apt-get install default-jdk"
            missing=1
        fi
        if ! command -v mvn &>/dev/null; then
            echo "ERROR: mvn no encontrado. Instalar con: sudo apt-get install maven"
            missing=1
        fi
    fi
    if [[ "$missing" == "1" ]]; then exit 1; fi
}

build_cpp() {
    echo "[1/4] Compilando ea_cpp_modern (NSGA-III runner)..."
    mkdir -p "${MODERN_ROOT}/build"
    g++-14 -std=c++23 -O3 -march=native -DNDEBUG \
        -I "${MODERN_ROOT}/include" \
        "${MODERN_ROOT}/examples/nsga3_runner.cpp" \
        -o "${CPP_BIN}" -lm
    echo "      → ${CPP_BIN}"
}

build_jmetal() {
    if [[ "${SKIP_JMETAL_BUILD:-0}" == "1" ]]; then
        echo "[2/4] Omitiendo build de jMetal (SKIP_JMETAL_BUILD=1)."
        return
    fi
    echo "[2/4] Construyendo jMetal 7.4 con Maven (puede tardar varios minutos)..."
    cd "${JMETAL_ROOT}"
    mvn clean install -DskipTests -Dgpg.skip=true -q
    echo "      jMetal construido."
    cd "${MODERN_ROOT}"
}

build_java_benchmark() {
    echo "[3/4] Compilando JMetalNSGA3Benchmark.java..."

    local CORE_JAR ALGO_JAR PROB_JAR
    CORE_JAR=$(ls "${JMETAL_ROOT}/jmetal-core/target/jmetal-core-"*.jar 2>/dev/null | grep -v 'sources\|javadoc\|dependencies' | head -1)
    ALGO_JAR=$(ls "${JMETAL_ROOT}/jmetal-algorithm/target/jmetal-algorithm-"*.jar 2>/dev/null | grep -v 'sources\|javadoc\|dependencies' | head -1)
    PROB_JAR=$(ls "${JMETAL_ROOT}/jmetal-problem/target/jmetal-problem-"*.jar 2>/dev/null | grep -v 'sources\|javadoc\|dependencies' | head -1)
    local COMMONS_JAR
    COMMONS_JAR=$(find "${HOME}/.m2/repository/org/apache/commons/commons-lang3" -name "commons-lang3-*.jar" 2>/dev/null | sort -V | tail -1)

    if [[ -z "$CORE_JAR" || -z "$ALGO_JAR" || -z "$PROB_JAR" ]]; then
        echo "ERROR: no se encontraron los jars de jMetal."
        echo "Ejecuta primero: cd ${JMETAL_ROOT} && mvn clean install -DskipTests"
        exit 1
    fi
    if [[ -z "$COMMONS_JAR" ]]; then
        echo "ERROR: commons-lang3 no encontrado en ~/.m2"
        exit 1
    fi

    JMETAL_CP="${CORE_JAR}:${ALGO_JAR}:${PROB_JAR}:${COMMONS_JAR}"
    export JMETAL_CP

    javac -cp "${JMETAL_CP}" "${BENCH_DIR}/JMetalNSGA3Benchmark.java" -d "${OUTDIR}"
    echo "      → ${OUTDIR}/JMetalNSGA3Benchmark.class"
}

measure_idle_temp() {
    THERMAL_FILE="/sys/class/thermal/thermal_zone0/temp"
    if [[ -r "${THERMAL_FILE}" ]]; then
        IDLE_TEMP_mC=$(cat "${THERMAL_FILE}")
        echo "Temperatura idle: $((IDLE_TEMP_mC/1000))°C"
        export IDLE_TEMP_mC THERMAL_FILE
    fi
}

thermal_cooldown() {
    echo ""
    echo "   Esperando enfriamiento de CPU..."
    if [[ -v IDLE_TEMP_mC && -r "${THERMAL_FILE}" ]]; then
        local TARGET_mC=$(( IDLE_TEMP_mC + 3000 ))
        local WAITED=0
        while true; do
            local NOW=$(cat "${THERMAL_FILE}")
            local TEMP_C=$(( NOW / 1000 ))
            local TEMP_D=$(( (NOW % 1000) / 100 ))
            echo -ne "   ${WAITED}s — ${TEMP_C}.${TEMP_D}°C  (objetivo ≤$(( TARGET_mC/1000 ))°C)\r"
            if (( NOW <= TARGET_mC )); then
                echo ""
                echo "   CPU a ${TEMP_C}.${TEMP_D}°C — dentro del margen."
                break
            fi
            sleep 5
            WAITED=$(( WAITED + 5 ))
            if (( WAITED >= 300 )); then
                echo ""
                echo "   AVISO: 300s sin alcanzar objetivo. Continuando (${TEMP_C}°C)."
                break
            fi
        done
    else
        echo "   Sensor térmico no disponible — sleep 60s."
        sleep 60
    fi
    echo ""
}

# ── Main ──────────────────────────────────────────────────────────────────────

echo "============================================================"
echo " Comparativa many-objective: ea_cpp_modern vs jMetal 7.4"
echo "============================================================"
echo " Algoritmo : NSGA-III (Das-Dennis, 3 objetivos)"
echo " Problemas : DTLZ1 (d=7), DTLZ2 (d=12), DTLZ3 (d=12), DTLZ4 (d=12)"
echo " Config    : divisions=12 → pop=92, max_evals≈25000"
echo "             SBX(p=0.9, η=20), PolMut(p=1/dim, η=20)"
echo " Runs      : ${NUMRUNS}"
echo " Salida    : ${OUTDIR}"
echo "============================================================"
echo ""

check_prereqs
mkdir -p "${OUTDIR}"
measure_idle_temp

# ── Capturar entorno de ejecución ─────────────────────────────────────────────
{
    echo "# Entorno de benchmark NSGA-III"
    echo ""
    echo "## Fecha"
    date -u +%Y-%m-%dT%H:%M:%SZ
    echo ""
    echo "## Hardware"
    grep 'model name' /proc/cpuinfo | head -1 || echo "CPU: desconocida"
    grep MemTotal /proc/meminfo      || echo "RAM: desconocida"
    echo ""
    echo "## Sistema operativo"
    uname -a
    echo ""
    echo "## Compilador C++"
    g++-14 --version | head -1
    echo "Flags: -std=c++23 -O3 -march=native -DNDEBUG"
    echo ""
    echo "## Java"
    java -version 2>&1 | head -1 || echo "Java: no instalado"
    echo ""
    echo "## Parámetros del benchmark"
    echo "Algoritmo    : NSGA-III"
    echo "Objetivos    : 3"
    echo "Divisions    : 12 (Das-Dennis) → 91 ref points → pop_size=92"
    echo "max_evals C++: 25000"
    echo "maxIter Java : 272 (≈ 25024 evals totales)"
    echo "Crossover    : SBX p=0.9 η_c=20"
    echo "Mutation     : PolynomialMutation p=1/dim η_m=20"
    echo "Problemas    : DTLZ1(d=7), DTLZ2(d=12), DTLZ3(d=12), DTLZ4(d=12)"
    echo "Runs         : ${NUMRUNS}"
    echo "Seed base    : 42"
} > "${OUTDIR}/benchmark_env.md"
echo "      → ${OUTDIR}/benchmark_env.md"
echo ""

build_cpp
echo ""

if [[ "$SKIP_JAVA" == "0" ]]; then
    build_jmetal
    echo ""
    build_java_benchmark
    echo ""
fi

# ── C++ benchmark ─────────────────────────────────────────────────────────────
echo "[4/4] Corriendo benchmark C++..."
echo "      → ${CPP_CSV}"
"${CPP_BIN}" "${NUMRUNS}" 42 > "${CPP_CSV}"
echo "      C++ completado."

# ── Java benchmark ────────────────────────────────────────────────────────────
if [[ "$SKIP_JAVA" == "0" ]]; then
    thermal_cooldown

    echo "[5/4] Corriendo benchmark jMetal..."
    echo "      → ${JAVA_CSV}"
    java -server -Xms512m -Xmx2g \
         -XX:+UseG1GC -XX:+AlwaysPreTouch -XX:+DisableExplicitGC \
         -Djava.awt.headless=true \
         -cp "${JMETAL_CP}:${OUTDIR}" \
         JMetalNSGA3Benchmark "${JAVA_CSV}" "${NUMRUNS}" 42
    echo "      jMetal completado."
fi

# ── Análisis estadístico ──────────────────────────────────────────────────────
if [[ "$SKIP_JAVA" == "0" && -f "${CPP_CSV}" && -f "${JAVA_CSV}" ]]; then
    echo ""
    echo "[6/4] Análisis estadístico (Mann-Whitney U + speedup)..."
    if command -v python3 &>/dev/null && python3 -c "import scipy" &>/dev/null 2>&1; then
        python3 "${SCRIPT_DIR}/analyze_mo_comparison.py" \
            "${CPP_CSV}" "${JAVA_CSV}" 2>&1 | tee "${OUTDIR}/analysis.txt"
    else
        echo "      AVISO: scipy no disponible — análisis omitido."
        echo "      Instalar: pip3 install scipy pandas"
    fi
fi

# ── Resumen rápido C++ ────────────────────────────────────────────────────────
echo ""
echo "── Resumen C++ (${CPP_CSV}) ──────────────────────────────────────────────"
if command -v python3 &>/dev/null; then
    python3 -c "
import csv, statistics
from collections import defaultdict
data = defaultdict(lambda: {'times':[], 'igds':[]})
try:
    with open('${CPP_CSV}') as f:
        for row in csv.DictReader(f):
            p = row['problem']
            data[p]['times'].append(float(row['time_ms']))
            data[p]['igds'].append(float(row['igd']))
    for p, d in sorted(data.items()):
        tm = statistics.median(d['times'])
        im = statistics.median(d['igds'])
        print(f'  {p:<8}: median time={tm:8.1f}ms  median IGD={im:.6f}')
except Exception as e:
    print(f'  (error: {e})')
"
fi

echo ""
echo "Benchmark NSGA-III completado. Resultados en: ${OUTDIR}"
