#pragma once
/// @file two_point.hpp
/// @brief Two Point Crossover — selects two points and swaps the segment between them.
/// Reference: jMetal TwoPointCrossover.java (delegates to NPointCrossover with N=2)
///
/// Selects two distinct crossover points. The segment between them is swapped
/// between parents to create offspring.

#include <algorithm>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Two Point Crossover for binary and general encodings.
/// Swaps the segment between two randomly chosen points.
struct TwoPointCrossover {
    double crossover_probability = 0.9; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Binary; }
    // Also works for Real, Integer at runtime

    void apply(this TwoPointCrossover& self, Population<>& pop, int parent_a, int parent_b,
               int child_start) {
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

        int point1 = rng.uniform_int(0, pop.dim - 1);
        int point2 = rng.uniform_int(0, pop.dim - 1);

        while (point1 == point2) {
            point2 = rng.uniform_int(0, pop.dim - 1);
        }

        if (point1 > point2) {
            std::swap(point1, point2);
        }

        for (int j = 0; j < point1; ++j) {
            pop.gene(child_start, j) = pop.gene(parent_a, j);
            pop.gene(child_start + 1, j) = pop.gene(parent_b, j);
        }
        for (int j = point1; j <= point2; ++j) {
            pop.gene(child_start, j) = pop.gene(parent_b, j);
            pop.gene(child_start + 1, j) = pop.gene(parent_a, j);
        }
        for (int j = point2 + 1; j < pop.dim; ++j) {
            pop.gene(child_start, j) = pop.gene(parent_a, j);
            pop.gene(child_start + 1, j) = pop.gene(parent_b, j);
        }

        pop.set_evaluated(child_start, false);
        pop.set_evaluated(child_start + 1, false);
    }
};

} // namespace ea
