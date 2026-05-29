#pragma once
/// @file spea2_strength.hpp
/// @brief SPEA2 Strength Selection — selects based on SPEA2 fitness (strength + density).
///
/// Reference: Zitzler, E., Laumanns, M., & Thiele, L. "SPEA2: Improving the
/// Strength Pareto Evolutionary Algorithm", TIK Report 103, 2001.
///
/// Uses strength (raw fitness) and density (k-th nearest neighbor) for selection.
/// Lower fitness is better. Binary tournament variant with fitness-based comparison.

#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <vector>

namespace ea {

/// SPEA2 Strength-based Selection.
/// Performs binary tournament selection using SPEA2 fitness values.
/// Lower fitness (strength + density) is better.
struct SPEA2StrengthSelection {
    /// Select mating pool from population using SPEA2 fitness values.
    /// @param pop          Population<>
    /// @param mating_pool  Output indices (resized to 2 * pop_size)
    /// @param fitness      SPEA2 fitness of each individual (lower = better)
    void select(this SPEA2StrengthSelection& self, Population<>& pop,
                std::vector<int>& mating_pool, const std::vector<double>& fitness) {
        (void)self;
        auto& rng = Random::instance();
        mating_pool.resize(2 * pop.pop_size);

        for (int i = 0; i < 2 * pop.pop_size; ++i) {
            int a = rng.uniform_int(0, pop.pop_size - 1);
            int b = rng.uniform_int(0, pop.pop_size - 1);

            // Lower fitness is better in SPEA2
            if (fitness[a] < fitness[b]) {
                mating_pool[i] = a;
            } else if (fitness[b] < fitness[a]) {
                mating_pool[i] = b;
            } else {
                // Tied — pick randomly
                mating_pool[i] = rng.coin_flip(0.5) ? a : b;
            }
        }
    }
};

} // namespace ea
