#pragma once
/// @file spatial_spread.hpp
/// @brief Spatial Spread Deviation (SSD) Selection — tournament selection using
///        ranking and spatial spread deviation comparator.
/// Reference: jMetal SpatialSpreadDeviationSelection.java
///        Alejandro Santiago <aurelio.santiago@upalt.edu.mx>
///
/// Performs N-ary tournament selection where individuals are compared by:
/// 1. Lower Pareto rank (better)
/// 2. Higher spatial spread deviation (better diversity within rank)

#include <algorithm>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <vector>

namespace ea {

/// Spatial Spread Deviation Selection.
/// N-ary tournament selection preferring lower rank, then higher SSD.
struct SpatialSpreadDeviationSelection {
    int number_of_tournaments = 2; ///< Tournament size (default = binary tournament)

    /// Select mating pool using SSD-based tournament selection.
    /// @param pop              Population<>
    /// @param mating_pool      Output: selected indices (resized to 2 * pop_size)
    /// @param ranks            Pareto rank of each individual (lower = better)
    /// @param spatial_spread   Spatial spread deviation of each individual (higher = better)
    void select(this SpatialSpreadDeviationSelection& self, Population<>& pop,
                std::vector<int>& mating_pool, const std::vector<int>& ranks,
                const std::vector<double>& spatial_spread) {
        auto& rng = Random::instance();
        int n = pop.pop_size;
        if (n <= 0)
            return;

        int pool_size = 2 * n;
        mating_pool.resize(pool_size);

        int tournament_size = std::max(2, self.number_of_tournaments);

        for (int i = 0; i < pool_size; ++i) {
            // Pick first candidate
            int winner = rng.uniform_int(0, n - 1);

            // Tournament: compare with other random candidates
            for (int t = 1; t < tournament_size; ++t) {
                int candidate = rng.uniform_int(0, n - 1);

                // Prefer lower rank, then higher spatial spread
                if (ranks[candidate] < ranks[winner]) {
                    winner = candidate;
                } else if (ranks[candidate] == ranks[winner]) {
                    if (spatial_spread[candidate] > spatial_spread[winner]) {
                        winner = candidate;
                    }
                }
            }

            mating_pool[i] = winner;
        }
    }
};

} // namespace ea
