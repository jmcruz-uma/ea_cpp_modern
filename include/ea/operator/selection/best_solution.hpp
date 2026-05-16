#pragma once
/// @file best_solution.hpp
/// @brief BestSolutionSelection — returns the index of the best individual.
///
/// Reference: jMetal
///   jmetal-core/src/main/java/org/uma/jmetal/operator/selection/impl/BestSolutionSelection.java

#include <ea/core/population.hpp>
#include <vector>
#include <limits>

namespace ea {

/// BestSolutionSelection — returns the index of the best individual.
///
/// For single-objective: picks the individual with lowest objective[0].
/// For multi-objective: picks the individual with lowest rank.
/// The comparator logic is embedded; extend with a custom comparator if needed.
struct BestSolutionSelection {
    /// Find the index of the best individual.
    /// @param pop   Population
    /// @param ranks Rank of each individual (lower is better, used for multi-objective)
    /// @return Index of the best individual
    int select(this BestSolutionSelection& self, Population& pop,
               const std::vector<int>& ranks) {
        (void)self;
        if (pop.pop_size == 0) return -1;

        int best_idx = 0;

        if (pop.n_obj == 1) {
            // Single-objective: minimize objective[0]
            double best_val = pop.objective(0, 0);
            for (int i = 1; i < pop.pop_size; ++i) {
                double val = pop.objective(i, 0);
                if (val < best_val) {
                    best_val = val;
                    best_idx = i;
                }
            }
        } else {
            // Multi-objective: minimize rank
            int best_rank = ranks[0];
            for (int i = 1; i < pop.pop_size; ++i) {
                if (ranks[i] < best_rank) {
                    best_rank = ranks[i];
                    best_idx = i;
                }
            }
        }
        return best_idx;
    }
};

} // namespace ea
