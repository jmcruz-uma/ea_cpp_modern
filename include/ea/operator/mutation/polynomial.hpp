#pragma once
/// @file polynomial.hpp
/// @brief Polynomial Mutation — faithful to jMetal 7.4 PolynomialMutation.java.
///
/// Uses the bound-aware formula from the original NSGA-II C code (KANGAL):
///   delta1 = (y - lb) / (ub - lb)
///   delta2 = (ub - y) / (ub - lb)
///   if rnd ≤ 0.5: xy=1-delta1, val=2*rnd+(1-2*rnd)*xy^(η+1), Δq=val^(1/(η+1))-1
///   else:          xy=1-delta2, val=2*(1-rnd)+2*(rnd-0.5)*xy^(η+1), Δq=1-val^(1/(η+1))
///   y' = y + Δq*(ub-lb), clamped to [lb, ub]

#include <algorithm>
#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Polynomial Mutation for real-valued encodings.
/// Perturbs each gene independently with probability mutation_rate.
struct PolynomialMutation {
    double distribution_index = 20.0; ///< η_m: higher → smaller perturbation
    double mutation_rate      = -1.0; ///< Per-gene probability; -1 → auto: 1/dim

    static constexpr Encoding encoding() { return Encoding::Real; }

    void apply(this PolynomialMutation& self, Population<>& pop, int idx) {
        auto& rng  = Random::instance();
        double rate = self.mutation_rate < 0.0 ? 1.0 / pop.dim : self.mutation_rate;

        for (int j = 0; j < pop.dim; ++j) {
            if (!rng.coin_flip(rate))
                continue;

            double y  = pop.gene(idx, j);
            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];

            if (lb == ub) {
                pop.gene(idx, j) = lb;
                continue;
            }

            double delta1  = (y - lb) / (ub - lb);
            double delta2  = (ub - y) / (ub - lb);
            double rnd     = rng.uniform();
            double mut_pow = 1.0 / (self.distribution_index + 1.0);
            double deltaq;

            if (rnd <= 0.5) {
                double xy  = 1.0 - delta1;
                double val = 2.0 * rnd + (1.0 - 2.0 * rnd) * std::pow(xy, self.distribution_index + 1.0);
                deltaq = std::pow(val, mut_pow) - 1.0;
            } else {
                double xy  = 1.0 - delta2;
                double val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * std::pow(xy, self.distribution_index + 1.0);
                deltaq = 1.0 - std::pow(val, mut_pow);
            }

            y += deltaq * (ub - lb);
            pop.gene(idx, j) = std::clamp(y, lb, ub);
        }
        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
