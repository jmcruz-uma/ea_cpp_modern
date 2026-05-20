.PHONY: all test clean benchmark

CXX = g++-14
CXXFLAGS = -std=c++23 -O2 -I include
LDFLAGS = -lm

# Directories
BUILD_DIR = build
EXAMPLES_DIR = examples
TEST_DIR = tests/unit

# Targets
BENCHMARK_FAIR = $(BUILD_DIR)/benchmark_fair
BENCHMARK_LARGE = $(BUILD_DIR)/benchmark_large
TEST_COMPILE = $(BUILD_DIR)/test_compile
TEST_OPERATORS = $(BUILD_DIR)/test_operators

TEST_ISLAND = $(BUILD_DIR)/test_island_model
TEST_NEW_OPERATORS = $(BUILD_DIR)/test_new_operators

all: $(BUILD_DIR) $(BENCHMARK_FAIR) $(TEST_COMPILE)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BENCHMARK_FAIR): $(EXAMPLES_DIR)/benchmark_fair.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(BENCHMARK_LARGE): $(EXAMPLES_DIR)/benchmark_large.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(TEST_COMPILE): $(TEST_DIR)/test_compile.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(TEST_OPERATORS): $(TEST_DIR)/test_operators.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(TEST_ISLAND): tests/integration/test_island_model.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(TEST_NEW_OPERATORS): $(TEST_DIR)/test_new_operators.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# Quick test
test: $(TEST_COMPILE) $(TEST_OPERATORS) $(TEST_ISLAND) $(TEST_NEW_OPERATORS)
	@echo "=== Running tests ==="
	$(TEST_COMPILE)
	$(TEST_OPERATORS)
	$(TEST_ISLAND)
	$(TEST_NEW_OPERATORS)

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
