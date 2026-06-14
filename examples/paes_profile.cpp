/// @file paes_profile.cpp
/// @brief PAES component-level profiler — mide tiempo en cada parte del algoritmo.
///
/// Compila:
///   g++-14 -std=c++23 -O3 -march=native -DNDEBUG -I include examples/paes_profile.cpp -o build/paes_profile -lm
///
/// Uso:
///   ./build/paes_profile        # ZDT1, dim=30, 10 runs

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/problem/zdt.hpp>
#include <ea/util/adaptive_grid.hpp>
#include <ea/util/random.hpp>

using Clock = std::chrono::high_resolution_clock;
using Ns    = std::chrono::duration<double, std::nano>;

struct ProfCounters {
    double t_mutation_ns  = 0;
    double t_eval_ns      = 0;
    double t_dominance_ns = 0;
    double t_archive_ns   = 0;  // archive.add() — now includes internal density when needed
    int    n_iters        = 0;
    int    n_nondom       = 0;
    int    n_added        = 0;
    int    n_incremental  = 0;  // times O(1) incremental path was taken
    int    n_rebuild      = 0;  // times O(N) full rebuild was triggered inside add()
};

template <typename Problem>
ProfCounters run_profiled(Problem& prob, int dim, int max_evals, int archive_size,
                          int bisections, uint64_t seed) {
    ea::Random::instance().seed(seed);
    const int n_obj = Problem::num_objectives();

    ea::Population<> pop(1, dim, n_obj);
    pop.lower_bounds.assign(prob.lower_bounds().begin(), prob.lower_bounds().end());
    pop.upper_bounds.assign(prob.upper_bounds().begin(), prob.upper_bounds().end());

    ea::PolynomialMutation pm;
    pm.distribution_index = 20.0;
    pm.mutation_rate = 1.0 / dim;

    ProfCounters cnt;

    auto& rng = ea::Random::instance();
    for (int j = 0; j < dim; ++j)
        pop.gene(0, j) = rng.uniform(pop.lower_bounds[j], pop.upper_bounds[j]);
    prob.evaluate(pop, 0);
    pop.set_evaluated(0, true);
    int evals = 1;

    ea::AdaptiveGridArchive archive(archive_size, bisections, n_obj);
    {
        std::vector<double> ig(dim), io(n_obj);
        for (int j = 0; j < dim; ++j) ig[j] = pop.gene(0, j);
        for (int o = 0; o < n_obj; ++o) io[o] = pop.objective(0, o);
        archive.add(ig, io);
    }

    ea::Population<> current(1, dim, n_obj, 0);
    current.lower_bounds = pop.lower_bounds;
    current.upper_bounds = pop.upper_bounds;
    for (int j = 0; j < dim; ++j) current.gene(0, j) = pop.gene(0, j);
    for (int o = 0; o < n_obj; ++o) current.objective(0, o) = pop.objective(0, o);
    current.flags[0] = pop.flags[0];

    ea::Population<> mutant(1, dim, n_obj, 0);
    mutant.lower_bounds = pop.lower_bounds;
    mutant.upper_bounds = pop.upper_bounds;

    std::vector<double> mut_genes(dim), mut_obj(n_obj), curr_obj(n_obj);

    while (evals < max_evals) {
        // ── 1. Mutación ──
        auto t0 = Clock::now();
        for (int j = 0; j < dim; ++j) mutant.gene(0, j) = current.gene(0, j);
        mutant.set_evaluated(0, false);
        pm.apply(mutant, 0);
        auto t1 = Clock::now();
        cnt.t_mutation_ns += Ns(t1 - t0).count();

        // ── 2. Evaluación ──
        auto t2 = Clock::now();
        prob.evaluate(mutant, 0);
        mutant.set_evaluated(0, true);
        evals++;
        auto t3 = Clock::now();
        cnt.t_eval_ns += Ns(t3 - t2).count();

        for (int j = 0; j < dim; ++j) mut_genes[j] = mutant.gene(0, j);
        for (int o = 0; o < n_obj; ++o) mut_obj[o]  = mutant.objective(0, o);
        for (int o = 0; o < n_obj; ++o) curr_obj[o] = current.objective(0, o);

        // ── 3. Dominancia current vs mutant ──
        auto t4 = Clock::now();
        bool a_dom = false, b_dom = false;
        for (int o = 0; o < n_obj; ++o) {
            if (mut_obj[o] < curr_obj[o]) a_dom = true;
            else if (curr_obj[o] < mut_obj[o]) b_dom = true;
        }
        ea::Dominance dom;
        if      (a_dom && !b_dom) dom = ea::Dominance::Dominates;
        else if (b_dom && !a_dom) dom = ea::Dominance::Dominated;
        else                      dom = ea::Dominance::Equal;
        auto t5 = Clock::now();
        cnt.t_dominance_ns += Ns(t5 - t4).count();

        if (dom == ea::Dominance::Equal) cnt.n_nondom++;

        // ── 4. archive.add() — density now maintained INSIDE add() ──
        //    O(1) in steady-state (bounds stable); O(N) only when bounds expand.
        auto t6 = Clock::now();
        bool added = false;
        if (dom == ea::Dominance::Dominates || dom == ea::Dominance::Equal) {
            added = archive.add(mut_genes, mut_obj);
        }
        auto t7 = Clock::now();
        cnt.t_archive_ns += Ns(t7 - t6).count();
        if (added) cnt.n_added++;

        // ── 5. Replacement decision using fresh grid_occupancy_ ──
        if (dom == ea::Dominance::Equal && added) {
            // density_is_fresh() guaranteed true; no explicit density call needed
            int curr_loc = archive.find_grid_location_for(curr_obj);
            int mut_loc  = archive.find_grid_location_for(mut_obj);
            if (archive.grid_occupancy_[curr_loc] > archive.grid_occupancy_[mut_loc]) {
                for (int j = 0; j < dim; ++j) current.gene(0, j) = mutant.gene(0, j);
                for (int o = 0; o < n_obj; ++o) current.objective(0, o) = mutant.objective(0, o);
                current.flags[0] = mutant.flags[0];
            }
        } else if (dom == ea::Dominance::Dominates && added) {
            for (int j = 0; j < dim; ++j) current.gene(0, j) = mutant.gene(0, j);
            for (int o = 0; o < n_obj; ++o) current.objective(0, o) = mutant.objective(0, o);
            current.flags[0] = mutant.flags[0];
        }

        cnt.n_iters++;
    }

    // Capture archive path counters
    cnt.n_incremental = archive.n_incremental_;
    cnt.n_rebuild     = archive.n_rebuild_;

    return cnt;
}

int main() {
    const int MAX_EVALS    = 25000;
    const int ARCHIVE_SIZE = 100;
    const int BISECTIONS   = 5;
    const int RUNS         = 10;
    const int DIM          = 30;

    std::cout << "=== PAES Component Profiler (ZDT1, dim=" << DIM
              << ", max_evals=" << MAX_EVALS << ", " << RUNS << " runs) ===\n\n";

    ea::ZDT1 prob(DIM);

    ProfCounters total{};
    double total_wall_ms = 0.0;

    for (int r = 0; r < RUNS; ++r) {
        auto tw0 = Clock::now();
        auto cnt = run_profiled(prob, DIM, MAX_EVALS, ARCHIVE_SIZE, BISECTIONS, 42 + r);
        auto tw1 = Clock::now();

        double wall_ms = Ns(tw1 - tw0).count() / 1e6;
        total_wall_ms += wall_ms;

        total.t_mutation_ns  += cnt.t_mutation_ns;
        total.t_eval_ns      += cnt.t_eval_ns;
        total.t_dominance_ns += cnt.t_dominance_ns;
        total.t_archive_ns   += cnt.t_archive_ns;
        total.n_iters        += cnt.n_iters;
        total.n_nondom       += cnt.n_nondom;
        total.n_added        += cnt.n_added;
        total.n_incremental  += cnt.n_incremental;
        total.n_rebuild      += cnt.n_rebuild;
    }

    double wall_ms = total_wall_ms / RUNS;
    double mut_ms  = total.t_mutation_ns  / RUNS / 1e6;
    double eval_ms = total.t_eval_ns      / RUNS / 1e6;
    double dom_ms  = total.t_dominance_ns / RUNS / 1e6;
    double arc_ms  = total.t_archive_ns   / RUNS / 1e6;
    double other_ms = wall_ms - mut_ms - eval_ms - dom_ms - arc_ms;

    int n_iters       = total.n_iters       / RUNS;
    int n_nondom      = total.n_nondom      / RUNS;
    int n_added       = total.n_added       / RUNS;
    int n_incremental = total.n_incremental / RUNS;
    int n_rebuild     = total.n_rebuild     / RUNS;

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Tiempo total (wall):               " << wall_ms  << " ms  (100.0%)\n";
    std::cout << "  Mutación (copy+PM):              " << mut_ms   << " ms  ("
              << std::setprecision(1) << 100*mut_ms/wall_ms  << "%)\n";
    std::cout << "  Evaluación ZDT1:                 " << eval_ms  << " ms  ("
              << 100*eval_ms/wall_ms  << "%)\n";
    std::cout << "  Dominancia current vs mutant:    " << dom_ms   << " ms  ("
              << 100*dom_ms/wall_ms   << "%)\n";
    std::cout << "  archive.add() [incl. density]:   " << arc_ms   << " ms  ("
              << 100*arc_ms/wall_ms   << "%)\n";
    std::cout << "    — O(1) incr path: " << n_incremental << " calls\n";
    std::cout << "    — O(N) rebuild:   " << n_rebuild     << " calls\n";
    std::cout << "  Overhead timing/otro:            " << other_ms << " ms  ("
              << 100*other_ms/wall_ms << "%)\n";

    std::cout << "\nEstadísticas de iteración (media por run):\n";
    std::cout << "  Iteraciones totales:  " << n_iters  << "\n";
    std::cout << "  Casos no dominados:   " << n_nondom << " ("
              << std::setprecision(1) << 100.0*n_nondom/n_iters << "% de iters)\n";
    std::cout << "  archive.add() true:   " << n_added  << "\n";

    double arc_per_call_us = (total.t_archive_ns / total.n_iters) / 1000.0;
    double arc_per_add_us  = (total.t_archive_ns / (total.n_added > 0 ? total.n_added : 1)) / 1000.0;
    std::cout << "\nCostes por llamada:\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  archive.add() por iteración:     " << arc_per_call_us << " µs\n";
    std::cout << "  archive.add() por add() efectivo:" << arc_per_add_us  << " µs\n";

    return 0;
}
