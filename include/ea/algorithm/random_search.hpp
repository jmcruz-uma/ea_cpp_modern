#pragma once
/// @file random_search.hpp
/// @brief Random Search — the simplest baseline multi-objective algorithm.
/// Generates random solutions within bounds, evaluates them, and keeps the
/// non-dominated front as the final result.
///
/// This is useful as a baseline for comparison: any "real" EA should
/// outperform RandomSearch on standard benchmarks.

#include <algorithm>
#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <iostream>
#include <string_view>
#include <vector>

namespace ea {

/// Random Search algorithm — baseline that generates random solutions and
/// keeps only the non-dominated ones.
///
/// Usage:
///   RandomSearch rs;
///   rs.pop_size = 100;
///   rs.max_evals = 25000;
///   rs.run(pop, problem);
struct RandomSearch {
    int pop_size = 100;    ///< Number of random solutions to generate per iteration
    int max_evals = 25000; ///< Maximum number of fitness evaluations

    static constexpr std::string_view name() { return "Random Search"; }

    /// Run Random Search on the given population.
    /// @param pop Population with bounds set. Will be resized to pop_size.
    /// @param problem Callable: void(Population&, int) — evaluates individual's objectives
    template <typename Problem> void run(this auto& self, Population& pop, Problem&& problem) {
        if (self.max_evals < self.pop_size) {
            std::cerr << "[ea::RandomSearch] Warning: max_evals (" << self.max_evals
                      << ") must be >= pop_size (" << self.pop_size << "). Adjusting max_evals to "
                      << self.pop_size << ".\n";
            self.max_evals = self.pop_size;
        }

        const int dim = pop.dim;
        const int n_obj = pop.n_obj;
        const int n = self.pop_size;

        // Ensure population is sized correctly
        if (pop.pop_size != n)
            pop.resize(n);

        auto& rng = Random::instance();
        int evals = 0;

        // Allocate storage for all evaluated solutions (up to max_evals)
        Population archive(self.max_evals, dim, n_obj, pop.encoding, pop.n_const);
        archive.lower_bounds = pop.lower_bounds;
        archive.upper_bounds = pop.upper_bounds;

        // === Main loop: generate random solutions ===
        while (evals < self.max_evals) {
            // Generate a batch of random individuals
            int batch_size = std::min(n, self.max_evals - evals);

            for (int i = 0; i < batch_size; ++i) {
                for (int j = 0; j < dim; ++j) {
                    pop.gene(i, j) = rng.uniform(pop.lower_bounds[j], pop.upper_bounds[j]);
                }
                pop.set_evaluated(i, false);
            }

            // Evaluate the batch
            for (int i = 0; i < batch_size; ++i) {
                if (!pop.evaluated(i)) {
                    problem(pop, i);
                    pop.set_evaluated(i, true);
                }
            }

            // Copy evaluated solutions to archive
            for (int i = 0; i < batch_size; ++i) {
                int arch_idx = evals + i;
                for (int j = 0; j < dim; ++j) {
                    archive.gene(arch_idx, j) = pop.gene(i, j);
                }
                for (int o = 0; o < n_obj; ++o) {
                    archive.objective(arch_idx, o) = pop.objective(i, o);
                }
                archive.set_evaluated(arch_idx, true);
            }

            evals += batch_size;
        }

        archive.pop_size = evals;

        // === Extract non-dominated front ===
        auto fronts = fast_non_dominated_sort(archive);

        // Copy non-dominated front back to pop
        const auto& nd_front = fronts[0];
        int nd_size = static_cast<int>(nd_front.size());

        if (nd_size > n) {
            // Truncate to pop_size if non-dominated front is too large
            nd_size = n;
        }

        pop.resize(nd_size);
        for (int i = 0; i < nd_size; ++i) {
            int src = nd_front[i];
            for (int j = 0; j < dim; ++j) {
                pop.gene(i, j) = archive.gene(src, j);
            }
            for (int o = 0; o < n_obj; ++o) {
                pop.objective(i, o) = archive.objective(src, o);
            }
            pop.set_evaluated(i, true);
        }
    }
};

} // namespace ea