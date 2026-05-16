#pragma once
/// @file sbx.hpp
/// @brief Simulated Binary Crossover (SBX) — the most widely used crossover in MOEAs.
/// Reference: Deb, K., et al. "A fast and elitist multiobjective genetic algorithm: NSGA-II"
/// C++23 with deducing this.

#include <cmath>
#include <algorithm>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Simulated Binary Crossover (SBX) for real-valued encodings.
/// Produces two children from two parents, mimicking the probability distribution
/// of single-point crossover in binary strings.
struct SBXCrossover {
    double distribution_index = 20.0;  ///< Distribution index (η): higher = closer to parents
    double crossover_probability = 0.9; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Real; }

    /// Apply SBX crossover. Produces 2 children starting at child_start.
    void apply(this SBXCrossover& self, Population& pop,
               int parent_a, int parent_b, int child_start) {
        auto& rng = Random::instance();

        for (int j = 0; j < pop.dim; ++j) {
            double p1 = pop.gene(parent_a, j);
            double p2 = pop.gene(parent_b, j);

            if (std::abs(p1 - p2) < 1e-14 || !rng.coin_flip(self.crossover_probability)) {
                // No crossover — copy parents to children
                pop.gene(child_start, j) = p1;
                pop.gene(child_start + 1, j) = p2;
                continue;
            }

            double y1 = std::min(p1, p2);
            double y2 = std::max(p1, p2);
            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];

            double u = rng.uniform();
            double beta;

            if (u <= 0.5) {
                beta = std::pow(2.0 * u, 1.0 / (self.distribution_index + 1.0));
            } else {
                beta = std::pow(1.0 / (2.0 * (1.0 - u)), 1.0 / (self.distribution_index + 1.0));
            }

            double c1 = 0.5 * ((1.0 + beta) * y1 + (1.0 - beta) * y2);
            double c2 = 0.5 * ((1.0 - beta) * y1 + (1.0 + beta) * y2);

            // Bound handling
            c1 = std::clamp(c1, lb, ub);
            c2 = std::clamp(c2, lb, ub);

            // Swap with probability 0.5 to maintain symmetry
            if (rng.coin_flip()) std::swap(c1, c2);

            pop.gene(child_start, j) = c1;
            pop.gene(child_start + 1, j) = c2;
        }

        pop.set_evaluated(child_start, false);
        pop.set_evaluated(child_start + 1, false);
    }
};

} // namespace ea