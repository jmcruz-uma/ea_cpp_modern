#pragma once
/// @file cycle.hpp
/// @brief Cycle Crossover (CX) for permutation encodings.
/// Reference: jMetal CycleCrossover — https://github.com/jMetal/jMetal
/// C++23 with deducing this, SoA Population access.

#include <algorithm>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <vector>

namespace ea {

/// Cycle Crossover (CX) for permutation-based encodings.
/// Identifies cycles between two parent permutations and exchanges them
/// to produce offspring, preserving positional information.
/// Genes are integers stored as doubles in the Population array.
struct CycleCrossover {
    double crossover_probability = 1.0; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Permutation; }

    /// Apply Cycle crossover. Produces 2 children starting at child_start.
    void apply(this CycleCrossover& self, Population& pop, int parent_a, int parent_b,
               int child_start) {
        auto& rng = Random::instance();
        int dim = pop.dim;

        // Copy parents to children
        for (int j = 0; j < dim; ++j) {
            pop.gene(child_start, j) = pop.gene(parent_a, j);
            pop.gene(child_start + 1, j) = pop.gene(parent_b, j);
        }

        if (dim < 2 || !rng.coin_flip(self.crossover_probability)) {
            pop.set_evaluated(child_start, false);
            pop.set_evaluated(child_start + 1, false);
            return;
        }

        // Build reverse lookup: for parent_a, map value -> position
        std::vector<int> pos_in_a(dim);
        for (int j = 0; j < dim; ++j) {
            int val = static_cast<int>(pop.gene(parent_a, j));
            pos_in_a[val] = j;
        }

        // Track which positions belong to the cycle
        std::vector<bool> in_cycle(dim, false);

        // Find the cycle starting at a random position
        int start = rng.uniform_int(0, dim - 1);
        int idx = start;
        do {
            in_cycle[idx] = true;
            int val_b = static_cast<int>(pop.gene(parent_b, idx));
            idx = pos_in_a[val_b];
        } while (idx != start);

        // Fill offspring:
        // Child 0: cycle positions from parent_a, others from parent_b
        // Child 1: cycle positions from parent_b, others from parent_a
        for (int i = 0; i < dim; ++i) {
            if (in_cycle[i]) {
                // Child 0 keeps parent_a's value at cycle positions
                pop.gene(child_start, i) = pop.gene(parent_a, i);
                // Child 1 gets parent_b's value at cycle positions
                pop.gene(child_start + 1, i) = pop.gene(parent_b, i);
            } else {
                // Child 0 gets parent_b's value at non-cycle positions
                pop.gene(child_start, i) = pop.gene(parent_b, i);
                // Child 1 keeps parent_a's value at non-cycle positions
                pop.gene(child_start + 1, i) = pop.gene(parent_a, i);
            }
        }

        pop.set_evaluated(child_start, false);
        pop.set_evaluated(child_start + 1, false);
    }
};

} // namespace ea
