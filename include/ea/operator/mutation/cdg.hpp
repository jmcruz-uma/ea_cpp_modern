#pragma once
/// @file cdg.hpp
/// @brief Constrained Differential Growth (CDG) Mutation — real-valued mutation using
///        a differential growth perturbation scheme.
/// Reference: jMetal CDGMutation.java
///        https://github.com/jMetal/jMetal
///
/// The mutation adds a perturbation deltaq = 0.5 * (rnd - 0.5) * (1 - rnd^(-delta))
/// where delta controls the magnitude of the perturbation.

#include <cmath>
#include <algorithm>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Constrained Differential Growth Mutation for real-valued encodings.
/// Uses a power-law perturbation that is larger near the boundaries.
struct CDGMutation {
    double mutation_rate = -1.0;   ///< Per-gene probability (-1 = 1/dim auto)
    double delta = 0.5;             ///< Perturbation scaling factor

    static constexpr Encoding encoding() { return Encoding::Real; }

    void apply(this CDGMutation& self, Population& pop, int idx) {
        auto& rng = Random::instance();
        double rate = self.mutation_rate < 0 ? 1.0 / pop.dim : self.mutation_rate;

        for (int j = 0; j < pop.dim; ++j) {
            if (!rng.coin_flip(rate)) continue;

            double y = pop.gene(idx, j);
            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];
            double range = ub - lb;

            if (range < 1e-14) continue;

            double rnd = rng.uniform();
            double temp_delta = std::pow(rnd, -self.delta);
            double deltaq = 0.5 * (rnd - 0.5) * (1.0 - temp_delta);

            y += deltaq * range;
            pop.gene(idx, j) = std::clamp(y, lb, ub);
        }
        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
