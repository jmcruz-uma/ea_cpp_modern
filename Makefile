.PHONY: all test clean benchmark mu-lambda-runner nsga2-runner

CXX ?= g++-14
CXXFLAGS = -std=c++23 -O3 -march=native -DNDEBUG -I include
LDFLAGS = -lm

# Directories
BUILD_DIR = build
EXAMPLES_DIR = examples
TEST_DIR = tests/unit

# Targets
BENCHMARK_FAIR    = $(BUILD_DIR)/benchmark_fair
BENCHMARK_LARGE   = $(BUILD_DIR)/benchmark_large
MU_LAMBDA_RUNNER  = $(BUILD_DIR)/mu_lambda_runner
NSGA2_RUNNER      = $(BUILD_DIR)/nsga2_runner
TEST_COMPILE      = $(BUILD_DIR)/test_compile
TEST_OPERATORS    = $(BUILD_DIR)/test_operators
TEST_ISLAND       = $(BUILD_DIR)/test_island_model
TEST_NEW_OPERATORS    = $(BUILD_DIR)/test_new_operators
TEST_SEED_MANAGER     = $(BUILD_DIR)/test_seed_manager
TEST_STATS_COLLECTOR  = $(BUILD_DIR)/test_stats_collector
TEST_CONCEPTS         = $(BUILD_DIR)/test_concepts

all: $(BUILD_DIR) $(BENCHMARK_FAIR) $(MU_LAMBDA_RUNNER) $(NSGA2_RUNNER) $(TEST_COMPILE)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BENCHMARK_FAIR): $(EXAMPLES_DIR)/benchmark_fair.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(BENCHMARK_LARGE): $(EXAMPLES_DIR)/benchmark_large.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(MU_LAMBDA_RUNNER): $(EXAMPLES_DIR)/mu_lambda_runner.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(NSGA2_RUNNER): $(EXAMPLES_DIR)/nsga2_runner.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(TEST_COMPILE): $(TEST_DIR)/test_compile.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(TEST_OPERATORS): $(TEST_DIR)/test_operators.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(TEST_ISLAND): tests/integration/test_island_model.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(TEST_NEW_OPERATORS): $(TEST_DIR)/test_new_operators.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(TEST_SEED_MANAGER): $(TEST_DIR)/test_seed_manager.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(TEST_STATS_COLLECTOR): $(TEST_DIR)/test_stats_collector.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(TEST_CONCEPTS): $(TEST_DIR)/test_concepts.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# Quick test
test: $(TEST_COMPILE) $(TEST_OPERATORS) $(TEST_ISLAND) $(TEST_NEW_OPERATORS) $(TEST_SEED_MANAGER) $(TEST_STATS_COLLECTOR) $(TEST_CONCEPTS)
	@echo "=== Running tests ==="
	$(TEST_COMPILE)
	$(TEST_OPERATORS)
	$(TEST_ISLAND)
	$(TEST_NEW_OPERATORS)
	$(TEST_SEED_MANAGER)
	$(TEST_STATS_COLLECTOR)
	$(TEST_CONCEPTS)

# (µ,λ)-ES runner — reads numeric.json + experiment.json, writes stats JSON
mu-lambda-runner: $(MU_LAMBDA_RUNNER)
	$(MU_LAMBDA_RUNNER)

# NSGA-II MO runner — outputs CSV for MO comparison benchmark
nsga2-runner: $(NSGA2_RUNNER)

# Benchmark
benchmark: $(BENCHMARK_FAIR)
	@echo "=== Running benchmark ==="
	$(BENCHMARK_FAIR)

# Full benchmark suite (30 runs)
benchmark-full: $(BENCHMARK_FAIR)
	@echo "=== Running full benchmark suite ==="
	@bash scripts/benchmark_all.sh

clean:
	rm -rf $(BUILD_DIR)
