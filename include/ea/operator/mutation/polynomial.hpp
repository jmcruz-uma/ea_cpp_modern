#pragma once
/// @file polynomial.hpp
/// @brief Polynomial Mutation — standard companion to SBX crossover.
/// Reference: Deb, K. & Goyal, M. "A combined genetic adaptive search (GeneAS)
/// for engineering design"

#include <cmath>
#include <algorithm>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Polynomial Mutation for real-valued encodings.
/// Perturbs each gene with probability mutation_rate using a polynomial distribution.
struct PolynomialMutation {
    double distribution_index = 20.0; ///< Distribution index (η): higher = smaller perturbation
    double mutation_rate = -1.0;       ///< Per-gene probability (-1 = 1/dim auto)

    static constexpr Encoding encoding() { return Encoding::Real; }

    void apply(this PolynomialMutation& self, Population& pop, int idx) {
        auto& rng = Random::instance();
        double rate = self.mutation_rate < 0 ? 1.0 / pop.dim : self.mutation_rate;

        for (int j = 0; j < pop.dim; ++j) {
            if (!rng.coin_flip(rate)) continue;

            double y = pop.gene(idx, j);
            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];
            double delta = ub - lb;

            if (delta < 1e-14) continue;

            double u = rng.uniform();
            double deltaq;

            if (u < 0.5) {
                deltaq = std::pow(2.0 * u, 1.0 / (self.distribution_index + 1.0)) - 1.0;
            } else {
                deltaq = 1.0 - std::pow(2.0 * (1.0 - u), 1.0 / (self.distribution_index + 1.0));
            }

            y += deltaq * delta;
            pop.gene(idx, j) = std::clamp(y, lb, ub);
        }
        pop.set_evaluated(idx, false);
    }
};

} // namespace ea