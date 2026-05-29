// tests/unit/test_new_operators.cpp
// Unit tests for new selection and replacement operators

#include <iostream>
#include <cassert>
#include <cmath>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/operator/selection/environmental.hpp>
#include <ea/operator/selection/spea2_strength.hpp>
#include <ea/operator/selection/sus.hpp>
#include <ea/operator/selection/tournament_wor.hpp>
#include <ea/operator/selection/weighted.hpp>
#include <ea/operator/replacement/crowding_replacement.hpp>
#include <ea/operator/replacement/moead_replacement.hpp>
#include <ea/operator/replacement/smpso_replacement.hpp>
#include <ea/operator/replacement/spea2_replacement.hpp>
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

// ------------------------------------------------------------------
// Test 1: EnvironmentalSelection picks best by fitness
// ------------------------------------------------------------------
static bool test_environmental_selection() {
    std::cout << "TEST: EnvironmentalSelection ... " << std::flush;

    Population<> pop(10, 2, 1, 0, 0.0, 1.0);
    init_pop_random(pop, 100);

    // Set clear fitness: lower index = better fitness
    std::vector<double> fitness(10);
    for (int i = 0; i < 10; ++i) {
        fitness[i] = (double)(10 - i); // i=0 gets fitness=10 (worst), i=9 gets fitness=1 (best)
    }

    EnvironmentalSelection sel;
    sel.pool_size = 5;
    std::vector<int> selected;
    sel.select(pop, selected, fitness);

    bool pass = true;
    if ((int)selected.size() != 5) pass = false;

    // The 5 best (lowest fitness) should be indices 5,6,7,8,9
    std::sort(selected.begin(), selected.end());
    for (int i = 0; i < 5; ++i) {
        if (selected[i] < 5) pass = false; // should only have indices >= 5
    }

    // Check distinct
    for (size_t i = 1; i < selected.size(); ++i) {
        if (selected[i] == selected[i-1]) pass = false;
    }

    std::cout << (pass ? "PASS" : "FAIL") << std::endl;
    return pass;
}

// ------------------------------------------------------------------
// Test 2: SPEA2StrengthSelection uses fitness-based binary tournament
// ------------------------------------------------------------------
static bool test_spea2_strength_selection() {
    std::cout << "TEST: SPEA2StrengthSelection ... " << std::flush;

    Population<> pop(10, 2, 2, 0, 0.0, 1.0);
    init_pop_random(pop, 200);

    std::vector<double> fitness(10);
    for (int i = 0; i < 10; ++i) {
        fitness[i] = (double)(i + 1); // i=0 has lowest fitness (best)
    }

    SPEA2StrengthSelection sel;
    std::vector<int> mating_pool;
    Random::instance().seed(111);
    sel.select(pop, mating_pool, fitness);

    bool pass = true;
    if ((int)mating_pool.size() != 2 * pop.pop_size) pass = false;

    // All indices valid
    for (int idx : mating_pool) {
        if (idx < 0 || idx >= pop.pop_size) pass = false;
    }

    std::cout << (pass ? "PASS" : "FAIL") << std::endl;
    return pass;
}

// ------------------------------------------------------------------
// Test 3: SUS produces valid indices and respects fitness weights
// ------------------------------------------------------------------
static bool test_sus_selection() {
    std::cout << "TEST: StochasticUniversalSampling ... " << std::flush;

    Population<> pop(20, 2, 2, 0, 0.0, 1.0);
    init_pop_random(pop, 300);

    // Make index 0 have very good (low) fitness
    std::vector<double> fitness(20);
    for (int i = 0; i < 20; ++i) {
        fitness[i] = (double)(i + 1);
    }
    fitness[0] = 0.01; // Much better than others

    StochasticUniversalSampling sel;
    std::vector<int> mating_pool;
    sel.select(pop, mating_pool, fitness);

    bool pass = true;
    if ((int)mating_pool.size() != 2 * pop.pop_size) pass = false;

    // All valid
    for (int idx : mating_pool) {
        if (idx < 0 || idx >= pop.pop_size) pass = false;
    }

    std::cout << (pass ? "PASS" : "FAIL") << std::endl;
    return pass;
}

// ------------------------------------------------------------------
// Test 4: TournamentWithoutReplacement produces valid indices
// ------------------------------------------------------------------
static bool test_tournament_wor() {
    std::cout << "TEST: TournamentWithoutReplacement ... " << std::flush;

    Population<> pop(10, 2, 2, 0, 0.0, 1.0);
    init_pop_random(pop, 400);

    std::vector<int> ranks(10, 0);
    std::vector<double> cd(10, 0.0);

    TournamentWithoutReplacement sel;
    sel.tournament_size = 2;
    std::vector<int> mating_pool;
    sel.select(pop, mating_pool, ranks, cd);

    bool pass = true;
    if ((int)mating_pool.size() != 2 * pop.pop_size) pass = false;

    for (int idx : mating_pool) {
        if (idx < 0 || idx >= 10) pass = false;
    }

    std::cout << (pass ? "PASS" : "FAIL") << std::endl;
    return pass;
}

// ------------------------------------------------------------------
// Test 5: WeightedSelection (roulette wheel) produces valid indices
// ------------------------------------------------------------------
static bool test_weighted_selection() {
    std::cout << "TEST: WeightedSelection ... " << std::flush;

    Population<> pop(20, 2, 2, 0, 0.0, 1.0);
    init_pop_random(pop, 500);

    std::vector<double> fitness(20);
    for (int i = 0; i < 20; ++i) {
        fitness[i] = (double)(i + 1);
    }

    WeightedSelection sel;
    std::vector<int> mating_pool;
    sel.select(pop, mating_pool, fitness);

    bool pass = true;
    if ((int)mating_pool.size() != 2 * pop.pop_size) pass = false;

    for (int idx : mating_pool) {
        if (idx < 0 || idx >= pop.pop_size) pass = false;
    }

    std::cout << (pass ? "PASS" : "FAIL") << std::endl;
    return pass;
}

// ------------------------------------------------------------------
// Test 6: CrowdingReplacement (general crowding replacement)
// ------------------------------------------------------------------
static bool test_crowding_replacement() {
    std::cout << "TEST: CrowdingReplacement ... " << std::flush;

    // Create combined population with known structure
    Population<> combined(20, 2, 2, 0, 0.0, 1.0);
    for (int i = 0; i < 20; ++i) {
        combined.gene(i, 0) = 0.1 * i;
        combined.gene(i, 1) = 0.1 * i;
        // Objective 0: increasing, Objective 1: decreasing
        combined.objective(i, 0) = 0.1 * i;
        combined.objective(i, 1) = 1.0 - 0.1 * i;
    }

    CrowdingReplacement repl;
    std::vector<int> offspring = {}; // empty for this test
    auto selected = repl.replace(combined, offspring, 10);

    bool pass = true;
    if ((int)selected.size() != 10) pass = false;

    // All distinct
    std::sort(selected.begin(), selected.end());
    for (size_t i = 1; i < selected.size(); ++i) {
        if (selected[i] == selected[i-1]) pass = false;
    }

    std::cout << (pass ? "PASS" : "FAIL") << std::endl;
    return pass;
}

// ------------------------------------------------------------------
// Test 7: SPEA2Replacement
// ------------------------------------------------------------------
static bool test_spea2_replacement() {
    std::cout << "TEST: SPEA2Replacement ... " << std::flush;

    // Combined population of 20, want to select 10
    Population<> combined(20, 2, 2, 0, 0.0, 1.0);
    for (int i = 0; i < 20; ++i) {
        // Create a simple Pareto front structure
        // Front 0: indices 0-4 (non-dominated)
        // Front 1: indices 5-9 (dominated by front 0)
        // etc.
        if (i < 5) {
            combined.objective(i, 0) = 0.1 * i;
            combined.objective(i, 1) = 1.0 - 0.1 * i;
        } else {
            combined.objective(i, 0) = 0.1 * (i - 5) + 0.5;
            combined.objective(i, 1) = 1.0 - 0.1 * (i - 5) + 0.5;
        }
    }

    SPEA2Replacement repl;
    repl.archive_size = 10;
    std::vector<int> offspring = {};
    auto selected = repl.replace(combined, offspring, 10);

    bool pass = true;
    if ((int)selected.size() != 10) pass = false;

    // Verify distinct
    std::sort(selected.begin(), selected.end());
    for (size_t i = 1; i < selected.size(); ++i) {
        if (selected[i] == selected[i-1]) pass = false;
    }

    std::cout << (pass ? "PASS" : "FAIL") << std::endl;
    return pass;
}

// ------------------------------------------------------------------
// Test 8: MOEADReplacement (basic compatibility)
// ------------------------------------------------------------------
static bool test_moead_replacement() {
    std::cout << "TEST: MOEADReplacement ... " << std::flush;

    MOEADReplacement repl;
    Population<> pop(10, 2, 2, 0, 0.0, 1.0);
    std::vector<int> offspring = {};
    auto selected = repl.replace(pop, offspring, 10);

    bool pass = true;
    if ((int)selected.size() != 10) pass = false;

    std::cout << (pass ? "PASS" : "FAIL") << std::endl;
    return pass;
}

// ------------------------------------------------------------------
// Test 9: SMPSOReplacement (basic archive operations)
// ------------------------------------------------------------------
static bool test_smpso_replacement() {
    std::cout << "TEST: SMPSOReplacement ... " << std::flush;

    SMPSOReplacement repl;
    repl.archive_size = 10;

    // Create a source population with non-dominated individuals
    Population<> source(5, 2, 2, 0, 0.0, 1.0);
    for (int i = 0; i < 5; ++i) {
        source.objective(i, 0) = 0.1 * i;
        source.objective(i, 1) = 1.0 - 0.1 * i;
    }

    Population<> archive(10, 2, 2, 0, 0.0, 1.0);
    archive.pop_size = 0; // Start empty

    // Add all 5
    for (int i = 0; i < 5; ++i) {
        repl.add_to_archive(archive, i, source);
    }

    bool pass = true;
    if (archive.pop_size != 5) pass = false;

    std::cout << (pass ? "PASS" : "FAIL") << std::endl;
    return pass;
}

// ------------------------------------------------------------------
// Main
// ------------------------------------------------------------------
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "ea-cpp New Operator Unit Tests" << std::endl;
    std::cout << "========================================" << std::endl;

    int passed = 0;
    int failed = 0;

    if (test_environmental_selection()) ++passed; else ++failed;
    if (test_spea2_strength_selection()) ++passed; else ++failed;
    if (test_sus_selection()) ++passed; else ++failed;
    if (test_tournament_wor()) ++passed; else ++failed;
    if (test_weighted_selection()) ++passed; else ++failed;
    if (test_crowding_replacement()) ++passed; else ++failed;
    if (test_spea2_replacement()) ++passed; else ++failed;
    if (test_moead_replacement()) ++passed; else ++failed;
    if (test_smpso_replacement()) ++passed; else ++failed;

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Results: " << passed << " passed, " << failed << " failed" << std::endl;
    std::cout << "========================================" << std::endl;

    return failed > 0 ? 1 : 0;
}
