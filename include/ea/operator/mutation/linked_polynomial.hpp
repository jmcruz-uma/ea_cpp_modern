#pragma once
/// @file linked_polynomial.hpp
/// @brief Linked Polynomial Mutation — polynomial mutation with a shared random value
///        across all mutated variables.
/// Reference: jMetal LinkedPolynomialMutation.java
///        https://doi.org/10.1109/SSCI.2016.7850214
///
/// The key difference from standard polynomial mutation is that a single random value `rnd`
/// is drawn once and used for all variables that are mutated (the "linked" aspect).
/// This creates correlated perturbations across variables.

#include <algorithm>
#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Linked Polynomial Mutation for real-valued encodings.
/// Uses a single shared random value for all mutated variables, creating correlated perturbations.
struct LinkedPolynomialMutation {
    double distribution_index = 20.0; ///< Distribution index (η)
    double mutation_rate = -1.0;      ///< Per-gene probability (-1 = 1/dim auto)

    static constexpr Encoding encoding() { return Encoding::Real; }

    void apply(this LinkedPolynomialMutation& self, Population<>& pop, int idx) {
        auto& rng = Random::instance();
        double rate = self.mutation_rate < 0 ? 1.0 / pop.dim : self.mutation_rate;

        // Single shared random value for all mutations (the "linked" aspect)
        double shared_rnd = rng.uniform();
        double mut_pow = 1.0 / (self.distribution_index + 1.0);

        for (int j = 0; j < pop.dim; ++j) {
            if (!rng.coin_flip(rate))
                continue;

            double y = pop.gene(idx, j);
            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];

            if (lb >= ub)
                continue;

            double delta1 = (y - lb) / (ub - lb);
            double delta2 = (ub - y) / (ub - lb);
            double deltaq;

            if (shared_rnd <= 0.5) {
                double xy = 1.0 - delta1;
                double val = 2.0 * shared_rnd +
                             (1.0 - 2.0 * shared_rnd) * std::pow(xy, self.distribution_index + 1.0);
                deltaq = std::pow(val, mut_pow) - 1.0;
            } else {
                double xy = 1.0 - delta2;
                double val = 2.0 * (1.0 - shared_rnd) +
                             2.0 * (shared_rnd - 0.5) * std::pow(xy, self.distribution_index + 1.0);
                deltaq = 1.0 - std::pow(val, mut_pow);
            }

            y += deltaq * (ub - lb);
            pop.gene(idx, j) = std::clamp(y, lb, ub);
        }

        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
