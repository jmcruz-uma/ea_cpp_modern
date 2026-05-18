#pragma once
/// @file simple_inversion.hpp
/// @brief Simple Inversion Mutation — reverses the entire permutation with a probability.
/// Reference: jMetal SimpleInversionMutation.java
///
/// Unlike standard InversionMutation which reverses a subsequence, this operator
/// reverses the entire permutation when mutation is applied.

#include <algorithm>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Simple Inversion Mutation for permutation encodings.
/// Reverses the entire permutation (not just a subsequence).
struct SimpleInversionMutation {
    double mutation_rate = -1.0; ///< Per-individual probability (-1 = always if valid)

    static constexpr Encoding encoding() { return Encoding::Permutation; }

    void apply(this SimpleInversionMutation& self, Population& pop, int idx) {
        if (pop.dim <= 1)
            return;

        auto& rng = Random::instance();
        if (self.mutation_rate > 0 && !rng.coin_flip(self.mutation_rate))
            return;

        // Reverse entire permutation in-place
        int left = 0;
        int right = pop.dim - 1;
        while (left < right) {
            std::swap(pop.gene(idx, left), pop.gene(idx, right));
            ++left;
            --right;
        }

        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
