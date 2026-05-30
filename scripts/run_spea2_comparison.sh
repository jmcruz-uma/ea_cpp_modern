#!/usr/bin/env bash
# =============================================================================
# run_spea2_comparison.sh
# Comparativa multiobjetivo: ea_cpp_modern (SPEA2, C++) vs jMetal 7.4 (Java)
# Problemas: ZDT1, ZDT2, ZDT3, ZDT4
#
# Uso:
#   ./scripts/run_spea2_comparison.sh [OUTDIR] [NUMRUNS]
#
# Argumentos opcionales:
#   OUTDIR    directorio de salida    (default: results/spea2_comparison)
#   NUMRUNS   número de runs medidos  (default: 30)
#
# Variables de entorno:
#   SKIP_JAVA=1              omite la parte jMetal (solo corre C++)
#   SKIP_JMETAL_BUILD=1      asume que jMetal ya está compilado
#
# Ejemplo rápido (1 run, solo C++):
#   SKIP_JAVA=1 ./scripts/run_spea2_comparison.sh /tmp/spea2_test 1
#
# Ejemplo paper (30 runs):
#   ./scripts/run_spea2_comparison.sh results/spea2_paper 30
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODERN_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
JMETAL_ROOT="/home/alumno/jMetal"
BENCH_DIR="${MODERN_ROOT}/benchmarks/jmetal"

OUTDIR="${1:-${MODERN_ROOT}/results/spea2_comparison}"
NUMRUNS="${2:-30}"

CPP_BIN="${MODERN_ROOT}/build/spea2_runner"
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
    echo "[1/4] Compilando ea_cpp_modern (SPEA2 runner)..."
    mkdir -p "${MODERN_ROOT}/build"
    g++-14 -std=c++23 -O3 -march=native -DNDEBUG \
        -I "${MODERN_ROOT}/include" \
        "${MODERN_ROOT}/examples/spea2_runner.cpp" \
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
    echo "[3/4] Compilando JMetalSPEA2Benchmark.java..."

    local cp=""
    for jar in "${JMETAL_ROOT}"/*/target/jmetal-*.jar; do
        [[ "$jar" == *"-sources.jar"              ]] && continue
        [[ "$jar" == *"-javadoc.jar"              ]] && continue
        [[ "$jar" == *"-jar-with-dependencies.jar" ]] && continue
        [[ -f "$jar" ]] && cp="${cp}:${jar}"
    done
    cp="${cp:1}"

    if [[ -z "$cp" ]]; then
        echo "ERROR: no se encontraron los jars de jMetal."
        echo "Ejecuta primero: cd ${JMETAL_ROOT} && mvn clean install -DskipTests"
        exit 1
    fi

    JMETAL_CP="$cp"
    export JMETAL_CP

    javac -cp "${JMETAL_CP}" "${BENCH_DIR}/JMetalSPEA2Benchmark.java" -d "${OUTDIR}"
    echo "      → ${OUTDIR}/JMetalSPEA2Benchmark.class"
}

# ── Temperatura idle de referencia ────────────────────────────────────────────
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
echo " Comparativa multiobjetivo: ea_cpp_modern vs jMetal 7.4"
echo "============================================================"
echo " Algoritmo : SPEA2"
echo " Problemas : ZDT1 (d=30), ZDT2 (d=30), ZDT3 (d=30), ZDT4 (d=10)"
echo " Config    : pop=100, max_evals=25000, archive=100"
echo "             SBX pc=0.9 eta_c=20 | PM pm=1/D eta_m=20"
echo "             Torneo: BinaryTournament (Dominance)"
echo " Runs      : ${NUMRUNS}"
echo " Salida    : ${OUTDIR}"
echo "============================================================"
echo ""

check_prereqs
mkdir -p "${OUTDIR}"
measure_idle_temp

# ── Capturar entorno de ejecución ─────────────────────────────────────────────
{
    echo "# Entorno de benchmark SPEA2"
    echo ""
    echo "## Fecha"
    date -u +%Y-%m-%dT%H:%M:%SZ
    echo ""
    echo "## Hardware"
    grep 'model name' /proc/cpuinfo | head -1 || echo "CPU: desconocida"
    grep MemTotal /proc/meminfo || echo "RAM: desconocida"
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
    echo "Algoritmo    : SPEA2"
    echo "Problemas    : ZDT1 (d=30), ZDT2 (d=30), ZDT3 (d=30), ZDT4 (d=10)"
    echo "pop_size     : 100"
    echo "archive_size : 100"
    echo "max_evals    : 25000"
    echo "SBX          : pc=0.9, eta_c=20"
    echo "PM           : pm=1/D, eta_m=20"
    echo "Runs         : ${NUMRUNS}"
    echo "Base seed    : 42"
} > "${OUTDIR}/environment.md"
echo "Entorno documentado → ${OUTDIR}/environment.md"

build_cpp

if [[ "$SKIP_JAVA" == "0" ]]; then
    build_jmetal
    build_java_benchmark
fi

# ── Ejecutar ea_cpp_modern ────────────────────────────────────────────────────
echo ""
echo "── Ejecutando ea_cpp_modern (C++) ──────────────────────────────────"
START=$(date +%s%N)
"${CPP_BIN}" "${NUMRUNS}" 42 > "${CPP_CSV}"
END=$(date +%s%N)
CPP_WALL=$(echo "scale=2; ($END - $START) / 1000000000" | bc)
echo "   Completado en ${CPP_WALL}s → ${CPP_CSV}"

if [[ "$SKIP_JAVA" == "0" ]]; then
    thermal_cooldown

    # ── Ejecutar jMetal ───────────────────────────────────────────────────────
    echo "── Ejecutando jMetal 7.4 (Java) ────────────────────────────────────"
    START=$(date +%s%N)
    java -server -Xms512m -Xmx2g \
        -cp "${JMETAL_CP}:${OUTDIR}" \
        JMetalSPEA2Benchmark "${JAVA_CSV}" "${NUMRUNS}" 42
    END=$(date +%s%N)
    JAVA_WALL=$(echo "scale=2; ($END - $START) / 1000000000" | bc)
    echo "   Completado en ${JAVA_WALL}s → ${JAVA_CSV}"
fi

# ── Resumen ───────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " Tiempos totales (${NUMRUNS} runs × 4 problemas)"
echo "   ea_cpp_modern : ${CPP_WALL}s"
if [[ "$SKIP_JAVA" == "0" ]]; then
    echo "   jMetal 7.4    : ${JAVA_WALL}s  (incluye warmup JVM)"
fi
echo "============================================================"
echo ""
echo "Análisis detallado:"
echo "   python3 ${SCRIPT_DIR}/analyze_mo_comparison.py ${OUTDIR}"
