#pragma once
/// @file uniform.hpp
/// @brief Uniform Mutation — perturbs each gene with uniform distribution.
/// Reference: org.uma.jmetal.operator.mutation.impl.UniformMutation

#include <algorithm>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Uniform Mutation for real-valued encodings.
/// Replaces selected genes with a value uniformly perturbed around the current value.
struct UniformMutation {
    double perturbation = 0.5;      ///< Half-width of uniform perturbation (delta)
    double mutation_rate = -1.0;    ///< Per-gene probability (-1 = 1/dim auto)

    static constexpr Encoding encoding() { return Encoding::Real; }

    void apply(this UniformMutation& self, Population& pop, int idx) {
        auto& rng = Random::instance();
        double rate = self.mutation_rate < 0 ? 1.0 / pop.dim : self.mutation_rate;

        for (int j = 0; j < pop.dim; ++j) {
            if (!rng.coin_flip(rate)) continue;

            double y = pop.gene(idx, j);
            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];
            double delta = ub - lb;

            if (delta < 1e-14) continue;

            double rand = rng.uniform();
            double tmp = (rand - 0.5) * self.perturbation * delta;
            y += tmp;

            pop.gene(idx, j) = std::clamp(y, lb, ub);
        }
        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
