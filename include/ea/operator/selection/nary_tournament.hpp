#pragma once
/// @file nary_tournament.hpp
/// @brief N-ary Tournament Selection — selects best among N randomly chosen individuals.
///
/// Reference: jMetal
/// jmetal-core/src/main/java/org/uma/jmetal/operator/selection/impl/NaryTournamentSelection.java
///
/// The N-ary tournament works by picking `tournament_size` distinct individuals uniformly
/// at random and keeping the best according to a comparator. When `tournament_size == 2`,
/// this degenerates to binary tournament selection. Larger sizes increase selection pressure.

#include <algorithm>
#include <cassert>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <vector>

namespace ea {

/// N-ary Tournament Selection for dominance-based MOEAs.
/// Selects the best among N randomly chosen individuals based on a comparator.
/// For NSGA-II: prefer lower rank, then higher crowding distance.
struct NaryTournamentSelection {
    int tournament_size = 2; ///< Number of individuals competing in each tournament

    /// Select mating pool from population. Fills mating_pool with selected indices.
    /// @param pop Population
    /// @param mating_pool Output indices (will be resized to 2 * pop_size)
    /// @param ranks Rank of each individual (from NDS)
    /// @param crowding_dist Crowding distance of each individual
    void select(Population& pop, std::vector<int>& mating_pool, const std::vector<int>& ranks,
                const std::vector<double>& crowding_dist) {
        assert(tournament_size >= 1);
        assert(tournament_size <= pop.pop_size);

        auto& rng = Random::instance();
        mating_pool.resize(2 * pop.pop_size);

        for (int i = 0; i < 2 * pop.pop_size; ++i) {
            // Select tournament_size distinct individuals without replacement
            std::vector<int> indices(pop.pop_size);
            std::iota(indices.begin(), indices.end(), 0);
            std::shuffle(indices.begin(), indices.end(), rng.engine());

            int best = indices[0];
            for (int j = 1; j < tournament_size; ++j) {
                int candidate = indices[j];
                // Prefer: lower rank, then higher crowding distance
                if (ranks[candidate] < ranks[best]) {
                    best = candidate;
                } else if (ranks[candidate] > ranks[best]) {
                    // best stays
                } else if (crowding_dist[candidate] > crowding_dist[best]) {
                    best = candidate;
                }
            }
            mating_pool[i] = best;
        }
    }
};

} // namespace ea
