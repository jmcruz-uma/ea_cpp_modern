#pragma once
/// @file differential_evolution.hpp
/// @brief Differential Evolution Selection — selects distinct random individuals for DE.
/// Reference: jMetal DifferentialEvolutionSelection.java
///
/// Selects `number_of_solutions` distinct individuals from the population,
/// optionally including the current solution as the last selected index.
/// All selected individuals are distinct from each other and from any excluded index.

#include <vector>
#include <unordered_set>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Differential Evolution Selection.
/// Selects N distinct random individuals from the population.
/// Designed to provide parent indices for DE crossover operators.
struct DifferentialEvolutionSelection {
    int number_of_solutions = 3;     ///< Number of distinct solutions to select
    bool select_current = false;      ///< If true, current index is appended as last

    /// Select distinct individuals and store their indices in `mating_pool`.
    /// @param pop          Population
    /// @param mating_pool  Output: selected indices (resized to number_of_solutions)
    /// @param current_idx  Index of current solution (excluded from random selection)
    void select(this DifferentialEvolutionSelection& self, Population& pop,
                std::vector<int>& mating_pool,
                int current_idx = -1) {
        auto& rng = Random::instance();
        int n = pop.pop_size;
        if (n <= 0) return;

        int needed = self.number_of_solutions;
        if (self.select_current) {
            needed = self.number_of_solutions - 1;
        }

        mating_pool.clear();
        mating_pool.reserve(self.number_of_solutions);

        std::unordered_set<int> selected;

        while (static_cast<int>(selected.size()) < needed) {
            int idx = rng.uniform_int(0, n - 1);
            if (idx == current_idx) continue;
            if (selected.insert(idx).second) {
                mating_pool.push_back(idx);
            }
        }

        if (self.select_current && current_idx >= 0 && current_idx < n) {
            mating_pool.push_back(current_idx);
        }
    }
};

} // namespace ea
