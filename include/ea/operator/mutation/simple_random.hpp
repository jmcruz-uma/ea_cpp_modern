#pragma once
/// @file simple_random.hpp
/// @brief Simple Random Mutation — replaces gene with uniform random value in bounds.
/// Reference: org.uma.jmetal.operator.mutation.impl.SimpleRandomMutation

#include <algorithm>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Simple Random Mutation for real-valued encodings.
/// Replaces selected genes with a uniformly random value within [lower, upper].
struct SimpleRandomMutation {
    double mutation_rate = -1.0; ///< Per-gene probability (-1 = 1/dim auto)

    static constexpr Encoding encoding() { return Encoding::Real; }

    void apply(this SimpleRandomMutation& self, Population& pop, int idx) {
        auto& rng = Random::instance();
        double rate = self.mutation_rate < 0 ? 1.0 / pop.dim : self.mutation_rate;

        for (int j = 0; j < pop.dim; ++j) {
            if (!rng.coin_flip(rate))
                continue;

            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];

            pop.gene(idx, j) = rng.uniform(lb, ub);
        }
        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
