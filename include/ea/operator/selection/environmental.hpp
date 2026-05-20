#pragma once
/// @file environmental.hpp
/// @brief Environmental Selection (Truncation Selection) — keeps best N individuals.
///
/// Reference: jMetal
///   jmetal-core/src/main/java/org/uma/jmetal/operator/selection/impl/EnvironmentalSelection.java
///
/// General-purpose truncation selection: sorts by fitness and keeps the top N.
/// Used as a standalone selection operator or within replacement strategies.

#include <algorithm>
#include <ea/core/population.hpp>
#include <limits>
#include <vector>

namespace ea {

/// Environmental Selection — truncation selection keeping best N individuals.
///
/// Given a population and per-individual fitness values, selects `pool_size`
/// individuals with the lowest fitness (for minimization).
/// For single-objective: directly uses fitness values.
/// For multi-objective: uses a combined scalar fitness.
struct EnvironmentalSelection {
    int pool_size = 0; ///< Number of individuals to select

    /// Select `pool_size` individuals with lowest fitness.
    /// @param pop          Population
    /// @param mating_pool  Output: selected indices (resized)
    /// @param fitness      Fitness of each individual (lower = better)
    void select(this EnvironmentalSelection& self, Population& pop,
                std::vector<int>& mating_pool, const std::vector<double>& fitness) {
        const int n = pop.pop_size;
        mating_pool.clear();

        if (self.pool_size <= 0)
            return;
        if (self.pool_size > n)
            self.pool_size = n;

        // Create index-fitness pairs and sort by fitness ascending
        std::vector<std::pair<double, int>> fit_idx(n);
        for (int i = 0; i < n; ++i) {
            fit_idx[i] = {fitness[i], i};
        }
        std::sort(fit_idx.begin(), fit_idx.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        // Take top pool_size
        mating_pool.reserve(self.pool_size);
        for (int i = 0; i < self.pool_size; ++i) {
            mating_pool.push_back(fit_idx[i].second);
        }
    }

    /// Convenience: select best N by a single objective value.
    /// @param pop          Population
    /// @param mating_pool  Output: selected indices
    /// @param obj_idx      Objective index to use for ranking (default: 0)
    void select_by_objective(this EnvironmentalSelection& self, Population& pop,
                             std::vector<int>& mating_pool, int obj_idx = 0) {
        const int n = pop.pop_size;
        mating_pool.clear();

        if (self.pool_size <= 0)
            return;
        if (self.pool_size > n)
            self.pool_size = n;
        if (obj_idx < 0 || obj_idx >= pop.n_obj)
            obj_idx = 0;

        std::vector<std::pair<double, int>> fit_idx(n);
        for (int i = 0; i < n; ++i) {
            fit_idx[i] = {pop.objective(i, obj_idx), i};
        }
        std::sort(fit_idx.begin(), fit_idx.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        mating_pool.reserve(self.pool_size);
        for (int i = 0; i < self.pool_size; ++i) {
            mating_pool.push_back(fit_idx[i].second);
        }
    }
};

} // namespace ea
