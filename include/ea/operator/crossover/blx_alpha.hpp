#pragma once
/// @file blx_alpha.hpp
/// @brief BLX-Alpha Crossover for real-valued encodings.
/// Reference: jMetal BLXAlphaCrossover.java
///
/// Creates offspring from an interval that extends beyond the parent values
/// by a factor alpha. For each variable, offspring are drawn uniformly from
/// [min - alpha*d, max + alpha*d] where d = max - min and min/max are the
/// parent values.
///
/// Reference: Eshelman, L. J., & Schaffer, J. D. (1993). Real-coded genetic
/// algorithms and interval-schemata. Foundations of Genetic Algorithms 2.

#include <algorithm>
#include <cmath>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// BLX-Alpha Crossover for real-valued encodings.
/// Produces two children from two parents.
struct BLXAlphaCrossover {
    double alpha = 0.5;               ///< Exploration extension factor
    double crossover_probability = 0.9; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Real; }

    /// Apply BLX-alpha crossover. Produces 2 children starting at child_start.
    void apply(this auto& self, Population& pop,
               int parent_a, int parent_b, int child_start) {
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

            double min_val = std::min(p1, p2);
            double max_val = std::max(p1, p2);
            double range = max_val - min_val;

            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];

            double min_range = min_val - range * self.alpha;
            double max_range = max_val + range * self.alpha;

            double c1 = min_range + rng.uniform() * (max_range - min_range);
            double c2 = min_range + rng.uniform() * (max_range - min_range);

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
