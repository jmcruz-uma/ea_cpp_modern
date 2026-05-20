#pragma once
/// @file tournament_wor.hpp
/// @brief Tournament Selection Without Replacement — samples individuals without replacement.
///
/// Reference: jMetal
///   jmetal-core/src/main/java/org/uma/jmetal/operator/selection/impl/TournamentSelection.java
///
/// Unlike standard tournament selection (with replacement), each individual can be
/// selected at most once per generation. This increases diversity in the mating pool
/// and reduces duplicate selections.

#include <algorithm>
#include <cassert>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <numeric>
#include <vector>

namespace ea {

/// Tournament Selection Without Replacement.
/// Selects distinct individuals via tournament; no individual appears twice.
struct TournamentWithoutReplacement {
    int tournament_size = 2; ///< Number of individuals in each tournament

    /// Select mating pool without replacement.
    /// @param pop            Population
    /// @param mating_pool    Output: selected indices (resized to 2 * pop_size)
    /// @param ranks          Rank of each individual (lower = better)
    /// @param crowding_dist  Crowding distance of each individual (higher = better)
    void select(this TournamentWithoutReplacement& self, Population& pop,
                std::vector<int>& mating_pool, const std::vector<int>& ranks,
                const std::vector<double>& crowding_dist) {
        assert(self.tournament_size >= 1);
        assert(self.tournament_size <= pop.pop_size);

        auto& rng = Random::instance();
        const int n = pop.pop_size;
        const int pool_size = 2 * n;
        mating_pool.clear();
        mating_pool.reserve(pool_size);

        // All individual indices
        std::vector<int> all_indices(n);
        std::iota(all_indices.begin(), all_indices.end(), 0);

        int selected_count = 0;

        // Keep shuffling and running tournaments until we have enough
        while (selected_count < pool_size) {
            // Shuffle all indices
            std::shuffle(all_indices.begin(), all_indices.end(), rng.engine());

            // Run tournaments on consecutive groups
            for (int start = 0; start + self.tournament_size <= n && selected_count < pool_size;
                 start += self.tournament_size) {
                int best = all_indices[start];
                for (int t = 1; t < self.tournament_size; ++t) {
                    int candidate = all_indices[start + t];
                    // Prefer: lower rank, then higher crowding distance
                    if (ranks[candidate] < ranks[best]) {
                        best = candidate;
                    } else if (ranks[candidate] == ranks[best]) {
                        if (crowding_dist[candidate] > crowding_dist[best]) {
                            best = candidate;
                        }
                    }
                }

                mating_pool.push_back(best);
                ++selected_count;
            }
        }
    }
};

} // namespace ea
