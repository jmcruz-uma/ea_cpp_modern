#pragma once
/// @file binary_tournament2.hpp
/// @brief Binary Tournament Selection 2 — picks the first individual if comparator says so,
///        otherwise random. Used in SPEA2 and other algorithms.
/// @par Reference Zitzler, Thiele, Laumanns "SPEA2: Improving the Strength Pareto Evolutionary
///            Algorithm", TIK-Report 103, 2001.

#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Binary Tournament Selection variant 2.
/// Unlike BinaryTournamentSelection (which uses rank + crowding),
/// this variant takes a fitness comparator and picks the better of two
/// random individuals. If tied, picks randomly.
/// Typically used with strength-based fitness (SPEA2).
struct BinaryTournament2Selection {
    /// Select mating pool from population using strength/fitness values.
    /// @param pop Population
    /// @param mating_pool Output indices (will be resized to 2 * pop_size)
    /// @param fitness Fitness value of each individual (lower = better)
    void select(this BinaryTournament2Selection& self, Population& pop,
                std::vector<int>& mating_pool, const std::vector<double>& fitness) {
        auto& rng = Random::instance();
        mating_pool.resize(2 * pop.pop_size);

        for (int i = 0; i < 2 * pop.pop_size; ++i) {
            int a = rng.uniform_int(0, pop.pop_size - 1);
            int b = rng.uniform_int(0, pop.pop_size - 1);

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