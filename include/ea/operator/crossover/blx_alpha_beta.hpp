#pragma once
/// @file blx_alpha_beta.hpp
/// @brief BLX-Alpha-Beta Crossover for real-valued encodings.
/// Reference: jMetal BLXAlphaBetaCrossover.java
///
/// Extension of BLX-alpha where the exploration below the smaller parent
/// is controlled by alpha and above the larger parent by beta.
/// Offspring are drawn uniformly from [min - alpha*d, max + beta*d].
///
/// Reference: Eshelman, L. J., & Schaffer, J. D. (1993). Real-coded genetic
/// algorithms and interval-schemata. Foundations of Genetic Algorithms 2.

#include <algorithm>
#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// BLX-Alpha-Beta Crossover for real-valued encodings.
/// Produces two children from two parents with asymmetric exploration.
struct BLXAlphaBetaCrossover {
    double alpha = 0.5;                 ///< Exploration factor below smaller parent
    double beta = 0.5;                  ///< Exploration factor above larger parent
    double crossover_probability = 0.9; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Real; }

    /// Apply BLX-alpha-beta crossover. Produces 2 children starting at child_start.
    void apply(this auto& self, Population<>& pop, int parent_a, int parent_b, int child_start) {
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

            // Ensure p1 <= p2 for consistent range calculation
            double min_val = std::min(p1, p2);
            double max_val = std::max(p1, p2);
            double d = max_val - min_val;

            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];

            double c_min = min_val - self.alpha * d;
            double c_max = max_val + self.beta * d;

            double c1 = c_min + rng.uniform() * (c_max - c_min);
            double c2 = c_min + rng.uniform() * (c_max - c_min);

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
