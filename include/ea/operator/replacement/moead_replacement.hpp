#pragma once
/// @file moead_replacement.hpp
/// @brief MOEA/D Replacement — decomposition-based neighborhood update.
///
/// Reference: Zhang & Li, "MOEA/D: A Multiobjective Evolutionary Algorithm Based
/// on Decomposition", IEEE Trans. Evol. Comput., 2007.
///
/// MOEA/D replacement updates neighbor subproblems if the offspring improves
/// their scalarized fitness (using Tchebycheff or weighted sum).
///
/// This is a lightweight replacement operator for use with MOEA/D algorithms.
/// It updates the population in-place rather than returning indices.

#include <algorithm>
#include <cmath>
#include <ea/core/population.hpp>
#include <ea/util/aggregation.hpp>
#include <limits>
#include <vector>

namespace ea {

/// MOEA/D decomposition-based replacement operator.
/// Updates neighbors when offspring improves their scalar fitness.
struct MOEADReplacement {
    int max_replaced_solutions = 2;       ///< Max replacements per offspring
    AggregationType aggregation_type = AggregationType::Tchebycheff;
    bool normalize = false;               ///< Whether to normalize objectives

    /// Update neighborhood with offspring.
    /// @param pop               Population<> to update
    /// @param offspring         Offspring population (single individual at index 0)
    /// @param subproblem_id     Current subproblem index
    /// @param neighborhood      Neighborhood indices for subproblem
    /// @param lambda            Weight vectors for all subproblems
    /// @param ideal_point       Ideal point [n_obj]
    /// @param nadir_point       Nadir point [n_obj] (can be empty if not normalizing)
    /// @return Number of replacements made
    int update(this MOEADReplacement& self, Population<>& pop, const Population<>& offspring,
               int subproblem_id, const std::vector<int>& neighborhood,
               const std::vector<std::vector<double>>& lambda,
               const std::vector<double>& ideal_point,
               const std::vector<double>& nadir_point) {
        (void)subproblem_id;
        const int n_obj = pop.n_obj;
        int replaced = 0;

        AggregationFunction aggr;
        aggr.type = self.aggregation_type;
        aggr.normalize = self.normalize;

        for (int neighbor_idx : neighborhood) {
            // Current fitness of neighbor
            double f_old = aggr.compute(pop.objectives_ptr(neighbor_idx),
                                        lambda[neighbor_idx].data(),
                                        ideal_point.data(),
                                        nadir_point.empty() ? nullptr : nadir_point.data(),
                                        n_obj);

            // Offspring fitness with neighbor's weight
            double f_new = aggr.compute(offspring.objectives_ptr(0),
                                        lambda[neighbor_idx].data(),
                                        ideal_point.data(),
                                        nadir_point.empty() ? nullptr : nadir_point.data(),
                                        n_obj);

            if (f_new < f_old) {
                // Replace neighbor with offspring
                for (int j = 0; j < pop.dim; ++j) {
                    pop.gene(neighbor_idx, j) = offspring.gene(0, j);
                }
                for (int o = 0; o < n_obj; ++o) {
                    pop.objective(neighbor_idx, o) = offspring.objective(0, o);
                }
                pop.set_evaluated(neighbor_idx, true);
                replaced++;

                if (replaced >= self.max_replaced_solutions)
                    break;
            }
        }

        return replaced;
    }

    /// Compatibility: replace() signature for algorithm templates.
    /// Returns indices of all individuals (MOEA/D updates in-place).
    std::vector<int> replace(this MOEADReplacement& self, Population<>& pop,
                             const std::vector<int>& offspring_indices, int target_size) {
        (void)self;
        (void)offspring_indices;
        std::vector<int> selected(target_size);
        std::iota(selected.begin(), selected.end(), 0);
        return selected;
    }
};

} // namespace ea
