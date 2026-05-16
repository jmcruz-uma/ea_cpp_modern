#pragma once
/// @file scramble.hpp
/// @brief Scramble Mutation — shuffles a subsequence of a permutation.
/// Reference: org.uma.jmetal.operator.mutation.impl.ScrambleMutation

#include <algorithm>
#include <vector>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Scramble Mutation for permutation encodings.
/// Selects two positions, extracts the subsequence, shuffles it, and puts it back.
struct ScrambleMutation {
    double mutation_rate = -1.0;    ///< Per-individual probability (-1 = always if valid)

    static constexpr Encoding encoding() { return Encoding::Permutation; }

    void apply(this ScrambleMutation& self, Population& pop, int idx) {
        if (pop.dim <= 1) return;

        auto& rng = Random::instance();
        if (self.mutation_rate > 0 && !rng.coin_flip(self.mutation_rate)) return;

        int pos1 = rng.uniform_int(0, pop.dim - 1);
        int pos2 = rng.uniform_int(0, pop.dim - 1);

        while (pos1 == pos2) {
            pos2 = rng.uniform_int(0, pop.dim - 1);
        }

        if (pos1 > pos2) {
            std::swap(pos1, pos2);
        }

        // Extract subsequence into a temporary vector and shuffle
        int length = pos2 - pos1 + 1;
        std::vector<double> sub(length);
        for (int i = 0; i < length; ++i) {
            sub[i] = pop.gene(idx, pos1 + i);
        }

        // Fisher-Yates shuffle using our RNG
        for (int i = length - 1; i > 0; --i) {
            int j = rng.uniform_int(0, i);
            std::swap(sub[i], sub[j]);
        }

        // Write back
        for (int i = 0; i < length; ++i) {
            pop.gene(idx, pos1 + i) = sub[i];
        }

        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
