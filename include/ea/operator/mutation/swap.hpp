#pragma once
/// @file swap.hpp
/// @brief Swap Mutation — exchanges two positions in a permutation.
/// Reference: org.uma.jmetal.operator.mutation.impl.PermutationSwapMutation

#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Swap Mutation for permutation encodings.
/// Swaps the values at two randomly chosen distinct positions.
struct SwapMutation {
    double mutation_rate = -1.0; ///< Per-individual probability (-1 = always if valid)

    static constexpr Encoding encoding() { return Encoding::Permutation; }

    void apply(this SwapMutation& self, Population& pop, int idx) {
        if (pop.dim <= 1)
            return;

        auto& rng = Random::instance();
        if (self.mutation_rate > 0 && !rng.coin_flip(self.mutation_rate))
            return;

        int pos1 = rng.uniform_int(0, pop.dim - 1);
        int pos2 = rng.uniform_int(0, pop.dim - 1);

        while (pos1 == pos2) {
            pos2 = rng.uniform_int(0, pop.dim - 1);
        }

        std::swap(pop.gene(idx, pos1), pop.gene(idx, pos2));
        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
