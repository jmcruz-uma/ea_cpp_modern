#pragma once
/// @file power_law.hpp
/// @brief Power-Law Mutation — heavy-tailed perturbations for exploration/exploitation.
/// Reference: org.uma.jmetal.operator.mutation.impl.PowerLawMutation

#include <algorithm>
#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Power-Law Mutation for real-valued encodings.
/// Produces mostly small perturbations with occasional large jumps.
/// delta controls distribution shape: <1 = uniform-ish, 1 = balanced, >1 = heavy-tailed.
struct PowerLawMutation {
    double delta = 1.0;          ///< Power-law exponent (> 0)
    double mutation_rate = -1.0; ///< Per-gene probability (-1 = 1/dim auto)

    static constexpr Encoding encoding() { return Encoding::Real; }

    void apply(this PowerLawMutation& self, Population& pop, int idx) {
        auto& rng = Random::instance();
        double rate = self.mutation_rate < 0 ? 1.0 / pop.dim : self.mutation_rate;

        for (int j = 0; j < pop.dim; ++j) {
            if (!rng.coin_flip(rate))
                continue;

            double y = pop.gene(idx, j);
            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];
            double range = ub - lb;

            if (range < 1e-14)
                continue;

            double rnd = rng.uniform();
            // Clamp to avoid division by zero / extreme values
            rnd = std::clamp(rnd, 1e-10, 1.0 - 1e-10);

            double temp_delta = std::pow(rnd, -self.delta);
            double deltaq = 0.5 * (rnd - 0.5) * (1.0 - temp_delta);

            y += deltaq * range;
            pop.gene(idx, j) = std::clamp(y, lb, ub);
        }
        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
