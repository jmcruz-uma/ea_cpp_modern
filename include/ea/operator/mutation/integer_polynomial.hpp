#pragma once
/// @file integer_polynomial.hpp
/// @brief Integer Polynomial Mutation — polynomial mutation for integer encodings.
/// Reference: jMetal IntegerPolynomialMutation.java
///
/// Applies polynomial mutation and then rounds the result to the nearest integer.
/// Bounds are treated as integer values and the result is clamped and rounded.

#include <cmath>
#include <algorithm>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Integer Polynomial Mutation for integer encodings.
/// Mutates using polynomial distribution, then rounds to nearest integer.
struct IntegerPolynomialMutation {
    double distribution_index = 20.0; ///< Distribution index (η)
    double mutation_rate = -1.0;       ///< Per-gene probability (-1 = 1/dim auto)

    static constexpr Encoding encoding() { return Encoding::Integer; }

    void apply(this IntegerPolynomialMutation& self, Population& pop, int idx) {
        auto& rng = Random::instance();
        double rate = self.mutation_rate < 0 ? 1.0 / pop.dim : self.mutation_rate;
        double mut_pow = 1.0 / (self.distribution_index + 1.0);

        for (int j = 0; j < pop.dim; ++j) {
            if (!rng.coin_flip(rate)) continue;

            // Work with integer values
            double yd = pop.gene(idx, j);
            double yl = pop.lower_bounds[j];
            double yu = pop.upper_bounds[j];

            if (yl >= yu) continue;

            double delta1 = (yd - yl) / (yu - yl);
            double delta2 = (yu - yd) / (yu - yl);
            double rnd = rng.uniform();
            double deltaq;

            if (rnd <= 0.5) {
                double xy = 1.0 - delta1;
                double val = 2.0 * rnd + (1.0 - 2.0 * rnd) * std::pow(xy, self.distribution_index + 1.0);
                deltaq = std::pow(val, mut_pow) - 1.0;
            } else {
                double xy = 1.0 - delta2;
                double val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * std::pow(xy, self.distribution_index + 1.0);
                deltaq = 1.0 - std::pow(val, mut_pow);
            }

            yd += deltaq * (yu - yl);
            yd = std::clamp(yd, yl, yu);

            // Round to nearest integer
            int mutated = static_cast<int>(std::round(yd));
            // Clamp again after rounding to be safe
            int lb_int = static_cast<int>(std::ceil(yl));
            int ub_int = static_cast<int>(std::floor(yu));
            mutated = std::clamp(mutated, lb_int, ub_int);

            pop.gene(idx, j) = static_cast<double>(mutated);
        }

        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
