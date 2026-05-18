#pragma once
/// @file insert.hpp
/// @brief Insert Mutation — removes an element and reinserts it elsewhere.
/// Reference: org.uma.jmetal.operator.mutation.impl.InsertMutation

#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <vector>

namespace ea {

/// Insert Mutation for permutation encodings.
/// Removes the element at pos1 and inserts it at pos2, shifting elements in between.
struct InsertMutation {
    double mutation_rate = -1.0; ///< Per-individual probability (-1 = always if valid)

    static constexpr Encoding encoding() { return Encoding::Permutation; }

    void apply(this InsertMutation& self, Population& pop, int idx) {
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

        double value = pop.gene(idx, pos1);

        if (pos1 < pos2) {
            // Shift elements left to fill gap at pos1
            for (int i = pos1; i < pos2; ++i) {
                pop.gene(idx, i) = pop.gene(idx, i + 1);
            }
            pop.gene(idx, pos2 - 1) = value;
        } else {
            // pos1 > pos2: shift elements right to fill gap at pos1
            for (int i = pos1; i > pos2; --i) {
                pop.gene(idx, i) = pop.gene(idx, i - 1);
            }
            pop.gene(idx, pos2) = value;
        }

        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
