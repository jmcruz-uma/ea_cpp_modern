#pragma once
/// @file single_point.hpp
/// @brief Single Point Crossover — for binary and general encodings.
/// Reference: jMetal SinglePointCrossover.java
///
/// Selects a single crossover point. Genes before the point come from one parent,
/// genes after from the other parent.

#include <vector>
#include <algorithm>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Single Point Crossover for binary and general encodings.
/// A random point is chosen; values before the point come from parent A,
/// values after from parent B (and vice versa for the second child).
struct SinglePointCrossover {
    double crossover_probability = 0.9; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Binary; }
    // Also works for Real, Integer at runtime

    void apply(this SinglePointCrossover& self, Population& pop,
               int parent_a, int parent_b, int child_start) {
        auto& rng = Random::instance();

        if (!rng.coin_flip(self.crossover_probability)) {
            for (int j = 0; j < pop.dim; ++j) {
                pop.gene(child_start, j) = pop.gene(parent_a, j);
                pop.gene(child_start + 1, j) = pop.gene(parent_b, j);
            }
            pop.set_evaluated(child_start, false);
            pop.set_evaluated(child_start + 1, false);
            return;
        }

        // Single crossover point in [0, dim-1]
        // If point == 0, child 1 gets all from parent_b, child 2 gets all from parent_a
        // If point == dim-1, standard single-point behavior
        int point = rng.uniform_int(0, pop.dim - 1);

        for (int j = 0; j <= point; ++j) {
            pop.gene(child_start, j) = pop.gene(parent_a, j);
            pop.gene(child_start + 1, j) = pop.gene(parent_b, j);
        }
        for (int j = point + 1; j < pop.dim; ++j) {
            pop.gene(child_start, j) = pop.gene(parent_b, j);
            pop.gene(child_start + 1, j) = pop.gene(parent_a, j);
        }

        pop.set_evaluated(child_start, false);
        pop.set_evaluated(child_start + 1, false);
    }
};

} // namespace ea
