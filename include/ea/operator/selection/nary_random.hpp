#pragma once
/// @file nary_random.hpp
/// @brief NaryRandomSelection — selects N individuals uniformly at random (with replacement).
///
/// Reference: jMetal
///   jmetal-core/src/main/java/org/uma/jmetal/operator/selection/impl/NaryRandomSelection.java

#include <algorithm>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <vector>

namespace ea {

/// NaryRandomSelection — selects N individuals uniformly at random.
///
/// Selects `n_solutions` distinct individuals from the population without replacement.
/// If `n_solutions` >= pop_size, returns a shuffled copy of all indices.
struct NaryRandomSelection {
    int n_solutions = 1; ///< Number of solutions to select

    /// Select N random individuals into `selected`.
    /// @param pop      Population
    /// @param selected Output: selected individual indices (resized to n_solutions)
    void select(this NaryRandomSelection& self, Population& pop, std::vector<int>& selected) {
        auto& rng = Random::instance();
        selected.resize(self.n_solutions);

        if (self.n_solutions >= pop.pop_size) {
            // Return all indices shuffled
            selected.resize(pop.pop_size);
            std::iota(selected.begin(), selected.end(), 0);
            std::shuffle(selected.begin(), selected.end(), rng.engine());
            return;
        }

        // Fisher-Yates shuffle of [0, pop_size) then take first n
        std::vector<int> pool(pop.pop_size);
        std::iota(pool.begin(), pool.end(), 0);
        std::shuffle(pool.begin(), pool.end(), rng.engine());
        std::copy(pool.begin(), pool.begin() + self.n_solutions, selected.begin());
    }
};

} // namespace ea
