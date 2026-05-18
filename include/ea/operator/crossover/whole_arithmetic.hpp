#pragma once
/// @file whole_arithmetic.hpp
/// @brief Whole Arithmetic Crossover for real-valued encodings.
/// Reference: jMetal WholeArithmeticCrossover.java
///
/// A single random weight alpha is used for ALL variables (not per-variable).
/// Offspring are computed as:
///   child1 = alpha * p1 + (1 - alpha) * p2
///   child2 = alpha * p2 + (1 - alpha) * p1
///
/// Reference: Michalewicz, Z. (1996). Genetic Algorithms + Data Structures
/// = Evolution Programs. Springer-Verlag, Berlin.

#include <algorithm>
#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Whole Arithmetic Crossover for real-valued encodings.
/// Uses a single alpha for all variables, producing two children.
struct WholeArithmeticCrossover {
    double crossover_probability = 0.9; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Real; }

    /// Apply whole arithmetic crossover. Produces 2 children starting at child_start.
    /// Uses a single alpha for all decision variables.
    void apply(this auto& self, Population& pop, int parent_a, int parent_b, int child_start) {
        auto& rng = Random::instance();

        // Copy parents as default
        for (int j = 0; j < pop.dim; ++j) {
            pop.gene(child_start, j) = pop.gene(parent_a, j);
            pop.gene(child_start + 1, j) = pop.gene(parent_b, j);
        }

        if (!rng.coin_flip(self.crossover_probability)) {
            pop.set_evaluated(child_start, false);
            pop.set_evaluated(child_start + 1, false);
            return;
        }

        double alpha = rng.uniform();

        for (int j = 0; j < pop.dim; ++j) {
            double p1 = pop.gene(parent_a, j);
            double p2 = pop.gene(parent_b, j);

            double c1 = alpha * p1 + (1.0 - alpha) * p2;
            double c2 = alpha * p2 + (1.0 - alpha) * p1;

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
