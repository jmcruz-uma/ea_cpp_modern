#pragma once
/// @file uniform.hpp
/// @brief Uniform Crossover for real-valued and binary encodings.
/// Reference: jMetal UniformCrossover — https://github.com/jMetal/jMetal
/// C++23 with deducing this, SoA Population access.

#include <algorithm>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <vector>

namespace ea {

/// Uniform Crossover for real-valued and binary encodings.
/// For each gene, swaps values between parents with probability 0.5.
/// Produces two children from two parents.
struct UniformCrossover {
    double crossover_probability = 0.9; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Real; }
    // Note: Also works for Binary encoding at runtime

    /// Apply Uniform crossover. Produces 2 children starting at child_start.
    void apply(this UniformCrossover& self, Population& pop, int parent_a, int parent_b,
               int child_start) {
        auto& rng = Random::instance();
        int dim = pop.dim;

        for (int j = 0; j < dim; ++j) {
            double p1 = pop.gene(parent_a, j);
            double p2 = pop.gene(parent_b, j);

            if (!rng.coin_flip(self.crossover_probability)) {
                // No crossover — copy parents directly
                pop.gene(child_start, j) = p1;
                pop.gene(child_start + 1, j) = p2;
            } else {
                // With probability 0.5, swap the genes
                if (rng.coin_flip(0.5)) {
                    pop.gene(child_start, j) = p2;
                    pop.gene(child_start + 1, j) = p1;
                } else {
                    pop.gene(child_start, j) = p1;
                    pop.gene(child_start + 1, j) = p2;
                }
            }
        }

        pop.set_evaluated(child_start, false);
        pop.set_evaluated(child_start + 1, false);
    }
};

} // namespace ea
