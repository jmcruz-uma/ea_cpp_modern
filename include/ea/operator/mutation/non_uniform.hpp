#pragma once
/// @file non_uniform.hpp
/// @brief Non-Uniform Mutation — perturbation magnitude decreases over iterations.
/// Reference: org.uma.jmetal.operator.mutation.impl.NonUniformMutation

#include <cmath>
#include <algorithm>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Non-Uniform Mutation for real-valued encodings.
/// Perturbation shrinks as iterations progress, favoring fine local search later.
struct NonUniformMutation {
    double perturbation = 0.1;      ///< Mutation perturbation parameter (b)
    double mutation_rate = -1.0;    ///< Per-gene probability (-1 = 1/dim auto)
    int max_iterations = 100;       ///< Maximum expected iterations
    int current_iteration = 0;      ///< Current iteration (set by algorithm)

    static constexpr Encoding encoding() { return Encoding::Real; }

    void apply(this NonUniformMutation& self, Population& pop, int idx) {
        auto& rng = Random::instance();
        double rate = self.mutation_rate < 0 ? 1.0 / pop.dim : self.mutation_rate;

        for (int j = 0; j < pop.dim; ++j) {
            if (!rng.coin_flip(rate)) continue;

            double y = pop.gene(idx, j);
            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];

            double rand = rng.uniform();
            double tmp;

            if (rand <= 0.5) {
                tmp = self.delta(ub - y);
            } else {
                tmp = self.delta(lb - y);
            }
            y += tmp;

            pop.gene(idx, j) = std::clamp(y, lb, ub);
        }
        pop.set_evaluated(idx, false);
    }

private:
    /// Calculates the delta value: y * (1 - rand^((1 - t/T)^b))
    double delta(this NonUniformMutation& self, double y) {
        auto& rng = Random::instance();
        double rand = rng.uniform();
        double ratio = static_cast<double>(self.current_iteration) / static_cast<double>(self.max_iterations);
        double exponent = std::pow(1.0 - ratio, self.perturbation);
        return y * (1.0 - std::pow(rand, exponent));
    }
};

} // namespace ea
