#pragma once
/// @file bit_flip.hpp
/// @brief Bit-Flip Mutation — flips individual bits with given probability.
/// Reference: org.uma.jmetal.operator.mutation.impl.BitFlipMutation

#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Bit-Flip Mutation for binary encodings.
/// Each bit in each gene is flipped independently with probability mutation_rate.
/// Genes are stored as doubles but interpreted as bit vectors.
struct BitFlipMutation {
    double mutation_rate = -1.0; ///< Per-bit flip probability

    static constexpr Encoding encoding() { return Encoding::Binary; }

    void apply(this BitFlipMutation& self, Population<>& pop, int idx) {
        auto& rng = Random::instance();
        // For binary: flip each bit in each gene with probability
        // In SoA, each "gene" is a double that encodes multiple bits
        // We treat each gene as a binary value and flip it
        double rate = self.mutation_rate < 0 ? 1.0 / pop.dim : self.mutation_rate;

        for (int j = 0; j < pop.dim; ++j) {
            if (rng.coin_flip(rate)) {
                // Flip the bit: 0 -> 1, 1 -> 0
                double val = pop.gene(idx, j);
                pop.gene(idx, j) = (val < 0.5) ? 1.0 : 0.0;
            }
        }
        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
