#pragma once
/// @file inversion.hpp
/// @brief Inversion Mutation — reverses a subsequence of a permutation.
/// Reference: org.uma.jmetal.operator.mutation.impl.InversionMutation

#include <algorithm>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Inversion Mutation for permutation encodings.
/// Selects two positions, reverses the subsequence between them (inclusive).
struct InversionMutation {
    double mutation_rate = -1.0; ///< Per-individual probability (-1 = always if valid)

    static constexpr Encoding encoding() { return Encoding::Permutation; }

    void apply(this InversionMutation& self, Population<>& pop, int idx) {
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

        if (pos1 > pos2) {
            std::swap(pos1, pos2);
        }

        // Reverse subsequence [pos1, pos2] in-place
        int left = pos1;
        int right = pos2;
        while (left < right) {
            std::swap(pop.gene(idx, left), pop.gene(idx, right));
            ++left;
            --right;
        }

        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
