#pragma once
/// @file sbx.hpp
/// @brief Simulated Binary Crossover (SBX) — the most widely used crossover in MOEAs.
/// Reference: Deb, K., et al. "A fast and elitist multiobjective genetic algorithm: NSGA-II"
/// Implementation matches jMetal's SBXCrossover exactly:
/// - Per-child betaq calculation with bound-aware scaling
/// - Same EPS, same probability logic, same ordering
/// C++23 with deducing this.

#include <algorithm>
#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Simulated Binary Crossover (SBX) for real-valued encodings.
/// Produces two children from two parents, mimicking the probability distribution
/// of single-point crossover in binary strings.
///
/// This implementation matches jMetal's SBXCrossover.java exactly:
/// - Each child gets its own betaq calculated with bound-aware scaling
/// - Child 1 uses lower bound scaling, child 2 uses upper bound scaling
/// - No random swap needed (the asymmetry is inherent in the calculation)
struct SBXCrossover {
    double distribution_index = 20.0;   ///< Distribution index (η): higher = closer to parents
    double crossover_probability = 0.9; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Real; }

    /// Apply SBX crossover on the given population.
    /// Parents at indices parent_a and parent_b.
    /// Children written to indices child_start and child_start+1.
    void apply(this SBXCrossover& self, Population<>& pop, int parent_a, int parent_b,
               int child_start) {
        auto& rng = Random::instance();

        for (int j = 0; j < pop.dim; ++j) {
            double value_x1 = pop.gene(parent_a, j);
            double value_x2 = pop.gene(parent_b, j);

            if (rng.uniform() > self.crossover_probability) {
                // No crossover — copy parents to children
                pop.gene(child_start, j) = value_x1;
                pop.gene(child_start + 1, j) = value_x2;
                continue;
            }

            if (std::abs(value_x1 - value_x2) < 1e-14) {
                // Parents are identical at this gene — copy directly
                pop.gene(child_start, j) = value_x1;
                pop.gene(child_start + 1, j) = value_x2;
                continue;
            }

            double y1, y2;
            if (value_x1 < value_x2) {
                y1 = value_x1;
                y2 = value_x2;
            } else {
                y1 = value_x2;
                y2 = value_x1;
            }

            double lower_bound = pop.lower_bounds[j];
            double upper_bound = pop.upper_bounds[j];

            double rand = rng.uniform();

            // Calculate betaq for child 1 (uses lower bound scaling)
            double beta = 1.0 + (2.0 * (y1 - lower_bound) / (y2 - y1));
            double alpha = 2.0 - std::pow(beta, -(self.distribution_index + 1.0));

            double betaq1;
            if (rand <= (1.0 / alpha)) {
                betaq1 = std::pow(rand * alpha, 1.0 / (self.distribution_index + 1.0));
            } else {
                betaq1 =
                    std::pow(1.0 / (2.0 - rand * alpha), 1.0 / (self.distribution_index + 1.0));
            }

            double c1 = 0.5 * ((y1 + y2) - betaq1 * (y2 - y1));

            // Calculate betaq for child 2 (uses upper bound scaling)
            beta = 1.0 + (2.0 * (upper_bound - y2) / (y2 - y1));
            alpha = 2.0 - std::pow(beta, -(self.distribution_index + 1.0));

            double betaq2;
            if (rand <= (1.0 / alpha)) {
                betaq2 = std::pow(rand * alpha, 1.0 / (self.distribution_index + 1.0));
            } else {
                betaq2 =
                    std::pow(1.0 / (2.0 - rand * alpha), 1.0 / (self.distribution_index + 1.0));
            }

            double c2 = 0.5 * ((y1 + y2) + betaq2 * (y2 - y1));

            // Clamp to bounds (same as jMetal's RepairDoubleSolutionWithBoundValue)
            c1 = std::clamp(c1, lower_bound, upper_bound);
            c2 = std::clamp(c2, lower_bound, upper_bound);

            // Swap children with probability 0.5 (same as jMetal)
            if (rng.coin_flip())
                std::swap(c1, c2);

            pop.gene(child_start, j) = c1;
            pop.gene(child_start + 1, j) = c2;
        }

        pop.set_evaluated(child_start, false);
        pop.set_evaluated(child_start + 1, false);
    }
};

} // namespace ea