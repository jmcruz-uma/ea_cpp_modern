#pragma once
/// @file random.hpp
/// @brief Random Selection — selects individuals uniformly at random from the population.
///
/// Reference: jMetal
/// jmetal-core/src/main/java/org/uma/jmetal/operator/selection/impl/RandomSelection.java
///
/// This operator performs no quality-based filtering; every individual has equal probability
/// of being selected. Useful for baseline comparisons or when combined with other operators
/// that provide selective pressure elsewhere.

#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <vector>

namespace ea {

/// Random Selection — uniform random selection of mating pool indices.
struct RandomSelection {
    /// Select mating pool from population uniformly at random.
    /// @param pop Population
    /// @param mating_pool Output indices (will be resized to 2 * pop_size)
    void select(Population& pop, std::vector<int>& mating_pool) {
        auto& rng = Random::instance();
        mating_pool.resize(2 * pop.pop_size);

        for (int i = 0; i < 2 * pop.pop_size; ++i) {
            mating_pool[i] = rng.uniform_int(0, pop.pop_size - 1);
        }
    }
};

} // namespace ea
