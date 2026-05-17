#!/bin/bash
# Energy estimation for ea-cpp benchmarks
# Since RAPL is not available on cloud VPS, we estimate based on CPU load

# AMD EPYC Rome TDP (estimated for cloud instance)
TDP_WATTS=120
UTILIZATION=0.30  # Estimated from CPU load during benchmark
RUN_TIME_MS=$1

if [ -z "$RUN_TIME_MS" ]; then
    echo "Usage: $0 <benchmark_time_ms>"
    echo "Example: $0 153"
    exit 1
fi

# Convert ms to seconds
RUN_TIME_S=$(echo "$RUN_TIME_MS / 1000" | bc -l)

# Energy = TDP * utilization * time
ENERGY_J=$(echo "$TDP_WATTS * $UTILIZATION * $RUN_TIME_S" | bc -l)
ENERGY_KJ=$(echo "$ENERGY_J / 1000" | bc -l)

echo "Energy Estimation"
echo "================="
echo "TDP: ${TDP_WATTS}W"
echo "Estimated utilization: ${UTILIZATION}"
echo "Run time: ${RUN_TIME_MS}ms = ${RUN_TIME_S}s"
echo "Energy: ${ENERGY_J} Joules = ${ENERGY_KJ} kJ"
echo ""
echo "For comparison (jMetal, 1831ms):"
JMETAL_ENERGY=$(echo "$TDP_WATTS * $UTILIZATION * 1.831" | bc -l)
echo "Estimated: ${JMETAL_ENERGY} Joules"
SAVINGS=$(echo "$JMETAL_ENERGY / $ENERGY_J" | bc -l)
echo "ea-cpp savings: ${SAVINGS}x"
