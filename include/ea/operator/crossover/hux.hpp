#pragma once
/// @file hux.hpp
/// @brief Half Uniform Crossover (HUX) — for binary encodings.
/// Reference: jMetal HUXCrossover.java
///
/// HUX swaps exactly half of the differing bits between two parents.
/// For each bit position where parents differ, swap with probability 0.5.
/// This ensures on average 50% of differing bits are exchanged.

#include <algorithm>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <vector>

namespace ea {

/// Half Uniform Crossover (HUX) for binary encodings.
/// Swaps bits where parents differ with probability 0.5.
struct HUXCrossover {
    double crossover_probability = 0.9; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Binary; }

    void apply(this HUXCrossover& self, Population<>& pop, int parent_a, int parent_b,
               int child_start) {
        auto& rng = Random::instance();

        if (!rng.coin_flip(self.crossover_probability)) {
            // No crossover — copy parents directly
            for (int j = 0; j < pop.dim; ++j) {
                pop.gene(child_start, j) = pop.gene(parent_a, j);
                pop.gene(child_start + 1, j) = pop.gene(parent_b, j);
            }
            pop.set_evaluated(child_start, false);
            pop.set_evaluated(child_start + 1, false);
            return;
        }

        // Collect positions where parents differ
        std::vector<int> differing_positions;
        differing_positions.reserve(pop.dim);

        for (int j = 0; j < pop.dim; ++j) {
            double v1 = pop.gene(parent_a, j);
            double v2 = pop.gene(parent_b, j);
            // For binary: 0 vs 1
            bool b1 = v1 >= 0.5;
            bool b2 = v2 >= 0.5;
            if (b1 != b2) {
                differing_positions.push_back(j);
            }
        }

        // Shuffle differing positions
        std::shuffle(differing_positions.begin(), differing_positions.end(), rng.engine());

        // Swap exactly half of differing bits
        int num_to_swap = static_cast<int>(differing_positions.size()) / 2;
        std::vector<bool> swap_mask(pop.dim, false);
        for (int i = 0; i < num_to_swap; ++i) {
            swap_mask[differing_positions[i]] = true;
        }

        for (int j = 0; j < pop.dim; ++j) {
            double p1 = pop.gene(parent_a, j);
            double p2 = pop.gene(parent_b, j);

            if (swap_mask[j]) {
                pop.gene(child_start, j) = p2;
                pop.gene(child_start + 1, j) = p1;
            } else {
                pop.gene(child_start, j) = p1;
                pop.gene(child_start + 1, j) = p2;
            }
        }

        pop.set_evaluated(child_start, false);
        pop.set_evaluated(child_start + 1, false);
    }
};

} // namespace ea
