#pragma once
/// @file binary_tournament.hpp
/// @brief Binary Tournament Selection — standard selection for NSGA-II and many MOEAs.

#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Binary Tournament Selection for dominance-based MOEAs.
/// Selects the better of two randomly chosen individuals based on a comparator.
/// For NSGA-II: prefer lower rank, then higher crowding distance.
struct BinaryTournamentSelection {
    /// Select mating pool from population. Fills mating_pool with selected indices.
    /// @param pop Population<>
    /// @param mating_pool Output indices (will be resized to 2 * pop_size)
    /// @param ranks Rank of each individual (from NDS)
    /// @param crowding_dist Crowding distance of each individual
    void select(this BinaryTournamentSelection& self, Population<>& pop,
                std::vector<int>& mating_pool, const std::vector<int>& ranks,
                const std::vector<double>& crowding_dist) {
        auto& rng = Random::instance();
        mating_pool.resize(2 * pop.pop_size);

        for (int i = 0; i < 2 * pop.pop_size; ++i) {
            int a = rng.uniform_int(0, pop.pop_size - 1);
            int b = rng.uniform_int(0, pop.pop_size - 1);

            // Prefer: lower rank, then higher crowding distance
            if (ranks[a] < ranks[b]) {
                mating_pool[i] = a;
            } else if (ranks[a] > ranks[b]) {
                mating_pool[i] = b;
            } else if (crowding_dist[a] > crowding_dist[b]) {
                mating_pool[i] = a;
            } else {
                mating_pool[i] = b;
            }
        }
    }
};

} // namespace ea