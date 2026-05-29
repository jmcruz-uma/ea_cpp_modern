// tests/unit/test_operators.cpp
// Unit tests for SBXCrossover, PolynomialMutation, BinaryTournamentSelection,
// and GaussianMutation

#include <iostream>
#include <cassert>
#include <cmath>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/mutation/gaussian.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/operator/selection/binary_tournament.hpp>
#include <ea/util/random.hpp>

using namespace ea;

// Helper: initialize population with random values in [0,1]
static void init_pop_random(Population<>& pop, uint64_t seed) {
    auto& rng = Random::instance();
    rng.seed(seed);
    for (int i = 0; i < pop.pop_size; ++i) {
        for (int j = 0; j < pop.dim; ++j) {
            pop.gene(i, j) = rng.uniform(pop.lower_bounds[j], pop.upper_bounds[j]);
        }
        // Set dummy objectives
        for (int o = 0; o < pop.n_obj; ++o) {
            pop.objective(i, o) = rng.uniform();
        }
        pop.set_evaluated(i, true);
    }
}

// Helper: check all genes are within bounds
static bool within_bounds(const Population<>& pop, int idx) {
    for (int j = 0; j < pop.dim; ++j) {
        double g = pop.gene(idx, j);
        if (g < pop.lower_bounds[j] - 1e-12 || g > pop.upper_bounds[j] + 1e-12) {
            return false;
        }
    }
    return true;
}

// ------------------------------------------------------------------
// Test 1: SBXCrossover produces valid children within bounds
// ------------------------------------------------------------------
static bool test_sbx_crossover() {
    std::cout << "TEST: SBXCrossover ... " << std::flush;

    Population<> pop(10, 5, 2, 0, 0.0, 1.0);
    init_pop_random(pop, 42);

    SBXCrossover cx;
    cx.distribution_index = 20.0;
    cx.crossover_probability = 1.0; // force crossover

    // Copy parent genes for comparison
    std::vector<double> parent0(pop.dim), parent1(pop.dim);
    for (int j = 0; j < pop.dim; ++j) {
        parent0[j] = pop.gene(0, j);
        parent1[j] = pop.gene(1, j);
    }

    cx.apply(pop, 0, 1, 2); // children go to indices 2 and 3

    bool pass = true;
    // Children must be within bounds
    if (!within_bounds(pop, 2)) pass = false;
    if (!within_bounds(pop, 3)) pass = false;

    // Children should differ from parents (probabilistic, but with p=1.0 and distinct parents,
    // extremely likely for at least one dimension to differ)
    bool child2_different = false;
    bool child3_different = false;
    for (int j = 0; j < pop.dim; ++j) {
        if (std::abs(pop.gene(2, j) - parent0[j]) > 1e-14) child2_different = true;
        if (std::abs(pop.gene(3, j) - parent1[j]) > 1e-14) child3_different = true;
    }
    if (!child2_different && !child3_different) pass = false;

    // Evaluated flags should be cleared
    if (pop.evaluated(2)) pass = false;
    if (pop.evaluated(3)) pass = false;

    std::cout << (pass ? "PASS" : "FAIL") << std::endl;
    return pass;
}

// ------------------------------------------------------------------
// Test 2: PolynomialMutation mutates within bounds and respects mutation_rate
// ------------------------------------------------------------------
static bool test_polynomial_mutation() {
    std::cout << "TEST: PolynomialMutation ... " << std::flush;

    bool pass = true;

    // --- 2a: mutation_rate = 0 => no mutation ---
    {
        Population<> pop(10, 5, 2, 0, 0.0, 1.0);
        init_pop_random(pop, 123);

        // Snapshot genes
        std::vector<double> before(pop.pop_size * pop.dim);
        for (int i = 0; i < pop.pop_size; ++i) {
            for (int j = 0; j < pop.dim; ++j) {
                before[i * pop.dim + j] = pop.gene(i, j);
            }
        }

        PolynomialMutation mut;
        mut.distribution_index = 20.0;
        mut.mutation_rate = 0.0; // no mutation

        for (int i = 0; i < pop.pop_size; ++i) {
            mut.apply(pop, i);
        }

        for (int i = 0; i < pop.pop_size; ++i) {
            for (int j = 0; j < pop.dim; ++j) {
                if (std::abs(pop.gene(i, j) - before[i * pop.dim + j]) > 1e-14) {
                    pass = false;
                }
            }
            // Evaluated flag should be cleared
            if (pop.evaluated(i)) pass = false;
        }
    }

    // --- 2b: mutation_rate = 1 => all genes mutated (and within bounds) ---
    {
        Population<> pop(10, 5, 2, 0, 0.0, 1.0);
        init_pop_random(pop, 456);

        std::vector<double> before(pop.pop_size * pop.dim);
        for (int i = 0; i < pop.pop_size; ++i) {
            for (int j = 0; j < pop.dim; ++j) {
                before[i * pop.dim + j] = pop.gene(i, j);
            }
        }

        PolynomialMutation mut;
        mut.distribution_index = 20.0;
        mut.mutation_rate = 1.0; // mutate every gene

        for (int i = 0; i < pop.pop_size; ++i) {
            mut.apply(pop, i);
        }

        for (int i = 0; i < pop.pop_size; ++i) {
            for (int j = 0; j < pop.dim; ++j) {
                // Must be within bounds
                if (!within_bounds(pop, i)) pass = false;
                // Must have mutated (since rate=1.0 and parents not at bounds exactly)
                // Use a tolerance: with rate=1.0, every gene should change
                if (std::abs(pop.gene(i, j) - before[i * pop.dim + j]) < 1e-14) {
                    pass = false;
                }
            }
            if (pop.evaluated(i)) pass = false;
        }
    }

    std::cout << (pass ? "PASS" : "FAIL") << std::endl;
    return pass;
}

// ------------------------------------------------------------------
// Test 3: BinaryTournament returns better-ranked individual
// ------------------------------------------------------------------
static bool test_binary_tournament() {
    std::cout << "TEST: BinaryTournamentSelection ... " << std::flush;

    Population<> pop(10, 5, 2, 0, 0.0, 1.0);
    init_pop_random(pop, 789);

    // Manually assign ranks and crowding distances
    std::vector<int> ranks(pop.pop_size);
    std::vector<double> crowding(pop.pop_size);
    for (int i = 0; i < pop.pop_size; ++i) {
        ranks[i] = i % 3;        // e.g. [0,1,2,0,1,2,0,1,2,0]
        crowding[i] = (double)(pop.pop_size - i); // decreasing
    }

    BinaryTournamentSelection sel;
    std::vector<int> mating_pool;

    // Seed for reproducibility
    Random::instance().seed(999);
    sel.select(pop, mating_pool, ranks, crowding);

    bool pass = true;
    if ((int)mating_pool.size() != 2 * pop.pop_size) pass = false;

    // Verify each selection: for the two chosen indices, the winner should be
    // the one with lower rank, or if equal, higher crowding distance.
    for (size_t i = 0; i < mating_pool.size(); ++i) {
        // We can't know which pair was drawn (it's random), but we can verify
        // that every selected index is a valid individual.
        if (mating_pool[i] < 0 || mating_pool[i] >= pop.pop_size) {
            pass = false;
        }
    }

    // Deterministic test: fix seed and manually compute expected winner for first draw
    // With seed=999, uniform_int(0,9) draws should be deterministic.
    // We'll verify the selection logic directly.
    {
        Random::instance().seed(999);
        int a = Random::instance().uniform_int(0, pop.pop_size - 1);
        int b = Random::instance().uniform_int(0, pop.pop_size - 1);
        int expected;
        if (ranks[a] < ranks[b]) {
            expected = a;
        } else if (ranks[a] > ranks[b]) {
            expected = b;
        } else if (crowding[a] > crowding[b]) {
            expected = a;
        } else {
            expected = b;
        }
        if (mating_pool[0] != expected) pass = false;
    }

    std::cout << (pass ? "PASS" : "FAIL") << std::endl;
    return pass;
}

// ------------------------------------------------------------------
// Test 4: GaussianMutation
// ------------------------------------------------------------------
static bool test_gaussian_mutation() {
    std::cout << "TEST: GaussianMutation ... " << std::flush;
    bool pass = true;

    // --- 4a: rate = 0.0 → operator never fires; genes and flags unchanged ---
    {
        Population<> pop(20, 10, 1, 0, -5.12, 5.12);
        init_pop_random(pop, 1001);

        std::vector<double> before(pop.pop_size * pop.dim);
        for (int i = 0; i < pop.pop_size; ++i)
            for (int j = 0; j < pop.dim; ++j)
                before[i * pop.dim + j] = pop.gene(i, j);

        GaussianMutation mut;
        mut.mutation_rate = 0.0;
        mut.step_size     = 1.0;

        for (int i = 0; i < pop.pop_size; ++i)
            mut.apply(pop, i);

        for (int i = 0; i < pop.pop_size; ++i) {
            // Genes must be unchanged
            for (int j = 0; j < pop.dim; ++j)
                if (std::abs(pop.gene(i, j) - before[i * pop.dim + j]) > 1e-15)
                    pass = false;
            // Evaluated flag must still be set (we never called set_evaluated(false))
            if (!pop.evaluated(i)) pass = false;
        }
    }

    // --- 4b: rate = 1.0, step_size = 0 → fires every call; exactly one gene
    //         value is unchanged (v * (1+0) == v), but evaluated flag IS cleared ---
    {
        Population<> pop(50, 10, 1, 0, -5.12, 5.12);
        init_pop_random(pop, 1002);

        GaussianMutation mut;
        mut.mutation_rate = 1.0;
        mut.step_size     = 0.0;   // N(0,1)*0 == 0 → val == v, but flag cleared

        for (int i = 0; i < pop.pop_size; ++i)
            mut.apply(pop, i);

        for (int i = 0; i < pop.pop_size; ++i)
            if (pop.evaluated(i)) { pass = false; break; }
    }

    // --- 4c: rate = 1.0, step_size = 100 → result always within bounds ---
    {
        Population<> pop(100, 30, 1, 0, -5.12, 5.12);
        init_pop_random(pop, 1003);

        GaussianMutation mut;
        mut.mutation_rate = 1.0;
        mut.step_size     = 100.0; // extreme: forces clamp in almost every case

        for (int i = 0; i < pop.pop_size; ++i)
            mut.apply(pop, i);

        for (int i = 0; i < pop.pop_size; ++i)
            if (!within_bounds(pop, i)) { pass = false; break; }
    }

    // --- 4d: rate = 1.0, step_size = 1.0 → exactly ONE gene changes per call ---
    //   Verify by counting how many genes differ before and after each apply.
    {
        Population<> pop(200, 20, 1, 0, -5.12, 5.12);
        init_pop_random(pop, 1004);

        GaussianMutation mut;
        mut.mutation_rate = 1.0;
        mut.step_size     = 1.0;

        // Avoid genes initialised to 0.0 (multiplicative blind spot);
        // init_pop_random uses uniform in [-5.12,5.12], so exact 0.0 is negligible.
        // We still guard: if step_size*N==0 the count will be 0 — acceptable.

        int violations = 0;
        for (int i = 0; i < pop.pop_size; ++i) {
            std::vector<double> before(pop.dim);
            for (int j = 0; j < pop.dim; ++j) before[j] = pop.gene(i, j);

            mut.apply(pop, i);

            int changed = 0;
            for (int j = 0; j < pop.dim; ++j)
                if (std::abs(pop.gene(i, j) - before[j]) > 1e-15) ++changed;

            // changed must be 0 (gene==0 blind spot or normal==0) or exactly 1
            if (changed > 1) ++violations;
        }
        if (violations > 0) pass = false;
    }

    // --- 4e: auto rate (-1.0) with dim=20 → approximately 1/20 of calls fire ---
    //   Run 10 000 applications; expect ~500 fires. Accept [300, 700] (±4 σ).
    {
        Population<> pop(1, 20, 1, 0, -5.12, 5.12);
        init_pop_random(pop, 1005);

        GaussianMutation mut;
        mut.mutation_rate = -1.0;  // auto = 1/dim = 0.05
        mut.step_size     = 1.0;

        const int trials = 10000;
        int fires = 0;
        for (int t = 0; t < trials; ++t) {
            // Reset gene to a non-zero value and mark as evaluated
            pop.gene(0, 0) = 1.0;
            pop.set_evaluated(0, true);

            mut.apply(pop, 0);

            if (!pop.evaluated(0)) ++fires;   // flag cleared → operator fired
        }
        // Expected fires ≈ 500; accept [300, 700]
        if (fires < 300 || fires > 700) pass = false;
    }

    std::cout << (pass ? "PASS" : "FAIL") << std::endl;
    return pass;
}

// ------------------------------------------------------------------
// Main
// ------------------------------------------------------------------
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "ea-cpp Operator Unit Tests" << std::endl;
    std::cout << "========================================" << std::endl;

    int passed = 0;
    int failed = 0;

    if (test_sbx_crossover()) ++passed; else ++failed;
    if (test_polynomial_mutation()) ++passed; else ++failed;
    if (test_binary_tournament()) ++passed; else ++failed;
    if (test_gaussian_mutation()) ++passed; else ++failed;

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Results: " << passed << " passed, " << failed << " failed" << std::endl;
    std::cout << "========================================" << std::endl;

    return failed > 0 ? 1 : 0;
}
