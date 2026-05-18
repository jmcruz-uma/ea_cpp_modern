#pragma once
/// @file neighborhood.hpp
/// @brief Neighborhood Selection — selects individuals from the neighborhood
///        of a given subproblem. Used in MOEA/D and decomposition-based algorithms.
/// @reference Zhang & Li, "MOEA/D: A Multiobjective Evolutionary Algorithm Based
///            on Decomposition", IEEE Trans. Evol. Comput., 2007.

#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Neighborhood Selection operator.
/// Selects parents from the T nearest neighbors of a subproblem,
/// or from the entire population with probability (1 - delta).
/// Core selection mechanism in MOEA/D.
struct NeighborhoodSelection {
    double neighborhood_prob = 0.9; ///< Probability of selecting from neighborhood (delta)

    /// Select a single parent from the neighborhood of a subproblem.
    /// @param neighborhood Neighborhood indices for the current subproblem [T]
    /// @param pop_size Total population size
    /// @return Index of selected parent
    int select_one(this NeighborhoodSelection& self,
                   const std::vector<int>& neighborhood, int pop_size) {
        auto& rng = Random::instance();
        if (rng.coin_flip(self.neighborhood_prob)) {
            int idx = rng.uniform_int(0, static_cast<int>(neighborhood.size()) - 1);
            return neighborhood[idx];
        } else {
            return rng.uniform_int(0, pop_size - 1);
        }
    }

    /// Select multiple distinct parents from the neighborhood.
    /// @param neighborhood Neighborhood indices for the current subproblem [T]
    /// @param pop_size Total population size
    /// @param n_select Number of distinct parents to select
    /// @return Vector of selected parent indices (all distinct)
    std::vector<int> select(this NeighborhoodSelection& self,
                            const std::vector<int>& neighborhood, int pop_size,
                            int n_select) {
        auto& rng = Random::instance();
        std::vector<int> selected;
        selected.reserve(n_select);

        while (static_cast<int>(selected.size()) < n_select) {
            int candidate = self.select_one(neighborhood, pop_size);

            // Ensure distinct parents
            bool duplicate = false;
            for (int s : selected) {
                if (s == candidate) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                selected.push_back(candidate);
            }
        }
        return selected;
    }
};

} // namespace ea