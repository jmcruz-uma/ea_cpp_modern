#pragma once
/// @file displacement.hpp
/// @brief Displacement Mutation — removes a subsequence and inserts it elsewhere (permutation).
/// Reference: jMetal DisplacementMutation.java
///
/// Selects two positions, extracts the subsequence between them, removes it,
/// and reinserts it at a new random position.

#include <vector>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Displacement Mutation for permutation encodings.
struct DisplacementMutation {
    double mutation_rate = -1.0;    ///< Per-individual probability (-1 = always if valid)

    static constexpr Encoding encoding() { return Encoding::Permutation; }

    void apply(this DisplacementMutation& self, Population& pop, int idx) {
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

        // Extract subsequence [pos1, pos2]
        int sublen = pos2 - pos1 + 1;
        std::vector<double> sub(sublen);
        for (int i = 0; i < sublen; ++i) {
            sub[i] = pop.gene(idx, pos1 + i);
        }

        // Remove subsequence: shift remaining elements left
        int remaining = pop.dim - sublen;
        std::vector<double> temp(remaining);
        int t = 0;
        for (int i = 0; i < pop.dim; ++i) {
            if (i < pos1 || i > pos2) {
                temp[t++] = pop.gene(idx, i);
            }
        }

        // Choose new insertion position in the remaining sequence
        int new_pos = rng.uniform_int(0, remaining);

        // Rebuild: temp[0..new_pos-1] + sub + temp[new_pos..]
        int k = 0;
        for (int i = 0; i < new_pos; ++i) {
            pop.gene(idx, k++) = temp[i];
        }
        for (int i = 0; i < sublen; ++i) {
            pop.gene(idx, k++) = sub[i];
        }
        for (int i = new_pos; i < remaining; ++i) {
            pop.gene(idx, k++) = temp[i];
        }

        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
