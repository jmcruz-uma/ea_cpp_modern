#pragma once
/// @file laplace.hpp
/// @brief Laplace Crossover (LX) for real-valued encodings.
/// Reference: jMetal LaplaceCrossover.java
///
/// Uses the Laplace distribution to generate offspring. For each variable,
/// a beta value is generated from the Laplace distribution and combined
/// with parent values to produce offspring.
///
/// Reference: Deep, K., & Thakur, M. (2007). A new crossover operator for
/// real coded genetic algorithms. Applied Mathematics and Computation,
/// 188(1), 895-911.

#include <algorithm>
#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Laplace Crossover for real-valued encodings.
/// Uses Laplace distribution with configurable scale parameter.
struct LaplaceCrossover {
    double scale = 0.5;                        ///< Scale parameter (b > 0)
    double crossover_probability = 0.9;        ///< Probability of applying crossover
    static constexpr double EPSILON = 1.0e-14; ///< Minimum difference to apply crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Real; }

    /// Apply Laplace crossover. Produces 2 children starting at child_start.
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

            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];

            if (std::abs(p1 - p2) <= self.EPSILON) {
                continue;
            }

            double u = rng.uniform();
            double r = rng.uniform();

            // Generate beta from Laplace distribution
            double beta;
            if (r <= 0.5) {
                beta = p1 - self.scale * std::log(1.0 - u);
            } else {
                beta = p1 + self.scale * std::log(u);
            }

            // Calculate offspring values
            double c1 = beta;
            double c2 = p1 + p2 - beta;

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
