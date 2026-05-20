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

        // Track which individuals have been selected
        std::vector<bool> used(n, false);
        int selected_count = 0;

        // Full passes through available individuals
        while (selected_count < pool_size) {
            // Gather available indices
            std::vector<int> available;
            available.reserve(n);
            for (int i = 0; i < n; ++i) {
                if (!used[i]) {
                    available.push_back(i);
                }
            }

            if (available.empty()) {
                // All used — reset and continue
                std::fill(used.begin(), used.end(), false);
                for (int i = 0; i < n; ++i) {
                    available.push_back(i);
                }
            }

            // Shuffle available indices
            std::shuffle(available.begin(), available.end(), rng.engine());

            // Process available indices in groups of tournament_size
            // For WOR, we may have leftover individuals if available.size() % tournament_size != 0
            // We handle complete groups and then reshuffle leftovers
            size_t max_start = available.size() / self.tournament_size * self.tournament_size;
            
            for (size_t start = 0; start < max_start && selected_count < pool_size;
                 start += self.tournament_size) {
                // Tournament among available[start .. start+tournament_size-1]
                int best = available[start];
                for (int t = 1; t < self.tournament_size; ++t) {
                    int candidate = available[start + t];
                    // Prefer: lower rank, then higher crowding distance
                    if (ranks[candidate] < ranks[best]) {
                        best = candidate;
                    } else if (ranks[candidate] > ranks[best]) {
                        // best stays
                    } else if (crowding_dist[candidate] > crowding_dist[best]) {
                        best = candidate;
                    }
                }

                mating_pool.push_back(best);
                used[best] = true;
                ++selected_count;
            }
        }
    }
};

} // namespace ea
