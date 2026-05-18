#pragma once
/// @file arithmetic.hpp
/// @brief Arithmetic Crossover for real-valued encodings.
/// Reference: jMetal ArithmeticCrossover.java
///
/// For each variable, a random weight alpha is generated and the offspring
/// are computed as:
///   child1 = alpha * p1 + (1 - alpha) * p2
///   child2 = (1 - alpha) * p1 + alpha * p2
///
/// Reference: Michalewicz, Z. (1996). Genetic Algorithms + Data Structures
/// = Evolution Programs. Springer-Verlag, Berlin.

#include <algorithm>
#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Arithmetic Crossover for real-valued encodings.
/// Produces two children as a weighted arithmetic mean of two parents.
struct ArithmeticCrossover {
    double crossover_probability = 0.9; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Real; }

    /// Apply arithmetic crossover. Produces 2 children starting at child_start.
    void apply(this auto& self, Population& pop, int parent_a, int parent_b, int child_start) {
        auto& rng = Random::instance();

        for (int j = 0; j < pop.dim; ++j) {
            double p1 = pop.gene(parent_a, j);
            double p2 = pop.gene(parent_b, j);

            // Copy parents as default
            pop.gene(child_start, j) = p1;
            pop.gene(child_start + 1, j) = p2;

            if (!rng.coin_flip(self.crossover_probability)) {
                continue;
            }

            double alpha = rng.uniform();

            double c1 = alpha * p1 + (1.0 - alpha) * p2;
            double c2 = (1.0 - alpha) * p1 + alpha * p2;

            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];

            // Clamp to bounds
            c1 = std::clamp(c1, lb, ub);
            c2 = std::clamp(c2, lb, ub);

            pop.gene(child_start, j) = c1;
            pop.gene(child_start + 1, j) = c2;
        }

        pop.set_evaluated(child_start, false);
        pop.set_evaluated(child_start + 1, false);
    }
};

} // namespace ea
