#pragma once
/// @file oxd.hpp
/// @brief Order Crossover with Duplicate Elimination (OXD) for permutation encodings.
/// Reference: jMetal OXDCrossover — https://github.com/jMetal/jMetal
/// C++23 with deducing this, SoA Population access.

#include <vector>
#include <algorithm>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Order Crossover with Duplicate Elimination (OXD) for permutation-based encodings.
/// Selects a subsequence from one parent and fills remaining positions from the other,
/// maintaining order and avoiding duplicates to produce valid permutations.
/// Genes are integers stored as doubles in the Population array.
struct OXDCrossover {
    double crossover_probability = 1.0; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Permutation; }

    /// Apply OXD crossover. Produces 2 children starting at child_start.
    void apply(this OXDCrossover& self, Population& pop,
               int parent_a, int parent_b, int child_start) {
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

        // Select two cutting points
        int cp1 = rng.uniform_int(0, dim - 2);
        int cp2 = rng.uniform_int(cp1 + 1, dim - 1);

        // Swap segments between parents for initial offspring
        for (int i = cp1; i <= cp2; ++i) {
            pop.gene(child_start, i) = pop.gene(parent_b, i);
            pop.gene(child_start + 1, i) = pop.gene(parent_a, i);
        }

        // Repair offspring
        repair_child(pop, parent_a, parent_b, child_start, cp1, cp2);
        repair_child(pop, parent_b, parent_a, child_start + 1, cp1, cp2);

        pop.set_evaluated(child_start, false);
        pop.set_evaluated(child_start + 1, false);
    }

private:
    /// Repair a child to ensure it contains a valid permutation.
    void repair_child(Population& pop, int parent, int other_parent,
                      int child_idx, int cp1, int cp2) {
        int dim = pop.dim;

        // Track which values are already in the segment
        std::vector<bool> in_segment(dim, false);
        for (int i = cp1; i <= cp2; ++i) {
            int val = static_cast<int>(pop.gene(child_idx, i));
            in_segment[val] = true;
        }

        // Fill remaining positions from other_parent, skipping duplicates
        int current_idx = (cp2 + 1) % dim;
        int other_idx = (cp2 + 1) % dim;

        while (current_idx != cp1) {
            int gene = static_cast<int>(pop.gene(other_parent, other_idx));
            if (!in_segment[gene]) {
                pop.gene(child_idx, current_idx) = static_cast<double>(gene);
                current_idx = (current_idx + 1) % dim;
            }
            other_idx = (other_idx + 1) % dim;
        }
    }
};

} // namespace ea
