// tests/unit/test_stats_collector.cpp
// Unit tests for ea::StatsCollector
//
// Tests cover:
//   A. JSON structure matches ea_cpp_original format
//   B. VarianceDiversity formula (exact replica of Bio4Res/ea VarianceDiversity.h)
//   C. Improvement tracking fires only when best strictly improves
//   D. Multiple runs accumulate correctly (seed, time, snapshots independent)
//   E. close_run() guards against double-close

#include <ea/util/stats_collector.hpp>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>

#include <cmath>
#include <iostream>
#include <string>

using ea::Population;
using ea::StatsCollector;
using ea::Encoding;

// ------------------------------------------------------------------
// Helpers
// ------------------------------------------------------------------

static Population<> make_pop(int mu, int dim, double obj_val,
                            double gene_val = 0.5) {
    Population<> pop(mu, dim, 1, 0, 0.0, 1.0);
    for (int i = 0; i < mu; ++i) {
        for (int j = 0; j < dim; ++j) pop.gene(i, j) = gene_val;
        pop.objective(i, 0) = obj_val;
        pop.set_evaluated(i, true);
    }
    return pop;
}

static bool approx(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) <= tol;
}

// ------------------------------------------------------------------
// Test A: JSON structure
// ------------------------------------------------------------------
static bool test_json_structure() {
    std::cout << "TEST: JSON structure ... " << std::flush;

    StatsCollector sc;
    sc.new_run(42);
    auto pop = make_pop(4, 2, 10.0);
    sc.take(pop, 100);
    sc.close_run();

    auto j = sc.to_json();

    bool pass = true;
    // top-level: array of 1 run
    pass &= j.is_array() && j.size() == 1;

    auto& run = j[0];
    pass &= run["run"]  == 0;
    pass &= run["seed"] == 42;
    pass &= run.contains("time");
    pass &= run["rundata"].is_array() && run["rundata"].size() == 1;

    auto& island = run["rundata"][0];
    pass &= island.contains("idata");
    pass &= island.contains("isols");

    auto& idata = island["idata"];
    pass &= idata["evals"].size() == 1;
    pass &= idata["best"].size()  == 1;
    pass &= idata["mean"].size()  == 1;
    pass &= idata["diversity"].size() == 1;
    pass &= static_cast<uint64_t>(idata["evals"][0]) == 100;
    pass &= approx(static_cast<double>(idata["best"][0]), 10.0);

    auto& isols = island["isols"];
    pass &= isols["evals"].size()   == 1;  // one improvement (first snapshot)
    pass &= isols["fitness"].size() == 1;
    pass &= isols["genome"].size()  == 1;
    pass &= isols["genome"][0].size() == 2;  // dim=2

    std::cout << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

// ------------------------------------------------------------------
// Test B: VarianceDiversity formula
//
// Manually verify for a population where we know the answer:
//   mu=4, dim=1, genes = {0, 1, 2, 3}
//   K = 0
//   sum  = 0+1+2+3 = 6
//   sum2 = 0+1+4+9 = 14
//   var  = (14 - 36/4) / 4 = (14 - 9) / 4 = 5/4 = 1.25
//   sigma_0 = sqrt(1.25) ≈ 1.11803...
//   diversity = sigma_0 / 1 = 1.11803...
// ------------------------------------------------------------------
static bool test_diversity_formula() {
    std::cout << "TEST: VarianceDiversity formula ... " << std::flush;

    Population<> pop(4, 1, 1, 0, 0.0, 10.0);
    double genes[] = {0.0, 1.0, 2.0, 3.0};
    for (int i = 0; i < 4; ++i) {
        pop.gene(i, 0) = genes[i];
        pop.objective(i, 0) = genes[i];
        pop.set_evaluated(i, true);
    }

    StatsCollector sc;
    sc.new_run(1);
    sc.take(pop, 4);
    sc.close_run();

    const double expected_div = std::sqrt(1.25); // 1.11803...
    const double got = sc.runs[0].snapshots[0].diversity;
    bool pass = approx(got, expected_div, 1e-9);

    std::cout << (pass ? "PASS" : "FAIL")
              << "  (got=" << got << " expected=" << expected_div << ")\n";
    return pass;
}

// ------------------------------------------------------------------
// Test C: Improvement tracking
//
// Sequence: obj=10 → 10 → 8 → 8 → 5
// Improvements should fire at evals: 100, 300, 500 (indices 0, 2, 4)
// Not at 200 and 400 (same or higher fitness)
// ------------------------------------------------------------------
static bool test_improvement_tracking() {
    std::cout << "TEST: Improvement tracking ... " << std::flush;

    StatsCollector sc;
    sc.new_run(1);

    double obj_seq[] = {10.0, 10.0, 8.0, 8.0, 5.0};
    for (int i = 0; i < 5; ++i) {
        auto pop = make_pop(2, 3, obj_seq[i], 0.1 * (i + 1));
        sc.take(pop, static_cast<uint64_t>((i + 1) * 100));
    }
    sc.close_run();

    const auto& imps = sc.runs[0].improvements;
    bool pass = true;
    pass &= imps.size() == 3;
    if (imps.size() == 3) {
        pass &= imps[0].evals == 100 && approx(imps[0].fitness, 10.0);
        pass &= imps[1].evals == 300 && approx(imps[1].fitness, 8.0);
        pass &= imps[2].evals == 500 && approx(imps[2].fitness, 5.0);
        // genome of improvement 2 (evals=300, gene_val=0.3) has 3 dims all ≈ 0.3
        pass &= approx(imps[1].genome[0], 0.3, 1e-9);
    }

    std::cout << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

// ------------------------------------------------------------------
// Test D: Multiple runs accumulate independently
// ------------------------------------------------------------------
static bool test_multiple_runs() {
    std::cout << "TEST: Multiple runs ... " << std::flush;

    StatsCollector sc;

    for (int r = 0; r < 3; ++r) {
        sc.new_run(static_cast<uint64_t>(r + 1));
        for (int g = 0; g < 4; ++g) {
            auto pop = make_pop(5, 2, 100.0 - g * 10.0);
            sc.take(pop, static_cast<uint64_t>((g + 1) * 100));
        }
        sc.close_run();
    }

    bool pass = true;
    pass &= sc.num_runs() == 3;

    auto j = sc.to_json();
    pass &= j.size() == 3;

    for (int r = 0; r < 3; ++r) {
        pass &= j[r]["run"]  == r;
        pass &= j[r]["seed"] == static_cast<uint64_t>(r + 1);
        pass &= j[r]["rundata"][0]["idata"]["evals"].size() == 4;
    }

    // Seeds should be independent
    pass &= static_cast<uint64_t>(j[0]["seed"]) == 1;
    pass &= static_cast<uint64_t>(j[1]["seed"]) == 2;
    pass &= static_cast<uint64_t>(j[2]["seed"]) == 3;

    std::cout << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

// ------------------------------------------------------------------
// Test E: close_run guard — calling close_run() twice is a no-op
// ------------------------------------------------------------------
static bool test_double_close_guard() {
    std::cout << "TEST: double close_run guard ... " << std::flush;

    StatsCollector sc;
    sc.new_run(7);
    auto pop = make_pop(3, 2, 5.0);
    sc.take(pop, 100);
    sc.close_run();
    sc.close_run();  // must not add a duplicate run

    bool pass = sc.num_runs() == 1;
    std::cout << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

// ------------------------------------------------------------------
// Test F: best and mean values
// ------------------------------------------------------------------
static bool test_best_and_mean() {
    std::cout << "TEST: best and mean values ... " << std::flush;

    // Population<>: objectives = {1, 2, 3, 4}  → best=1, mean=2.5
    Population<> pop(4, 1, 1, 0, 0.0, 10.0);
    double objs[] = {1.0, 2.0, 3.0, 4.0};
    for (int i = 0; i < 4; ++i) {
        pop.objective(i, 0) = objs[i];
        pop.gene(i, 0) = objs[i];
        pop.set_evaluated(i, true);
    }

    StatsCollector sc;
    sc.new_run(1);
    sc.take(pop, 4);
    sc.close_run();

    const auto& s = sc.runs[0].snapshots[0];
    bool pass = approx(s.best, 1.0) && approx(s.mean, 2.5);

    std::cout << (pass ? "PASS" : "FAIL")
              << "  (best=" << s.best << " mean=" << s.mean << ")\n";
    return pass;
}

// ------------------------------------------------------------------
// Main
// ------------------------------------------------------------------
int main() {
    std::cout << "========================================\n";
    std::cout << "ea-cpp StatsCollector Unit Tests\n";
    std::cout << "========================================\n";

    int passed = 0, failed = 0;
    auto run = [&](bool r) { if (r) ++passed; else ++failed; };

    run(test_json_structure());
    run(test_diversity_formula());
    run(test_improvement_tracking());
    run(test_multiple_runs());
    run(test_double_close_guard());
    run(test_best_and_mean());

    std::cout << "----------------------------------------\n";
    std::cout << "Results: " << passed << " passed, " << failed << " failed\n";
    std::cout << "========================================\n";
    return failed > 0 ? 1 : 0;
}
