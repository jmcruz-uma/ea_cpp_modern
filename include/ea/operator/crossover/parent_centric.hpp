#pragma once
/// @file parent_centric.hpp
/// @brief Parent-Centric Crossover (UNDX — Unimodal Normal Distribution Crossover)
/// for real-valued encodings.
/// Reference: jMetal UnimodalNormalDistributionCrossover.java
///
/// A multi-parent crossover that generates offspring based on the normal
/// distribution defined by three parent solutions. Particularly effective
/// for continuous optimization as it preserves population statistics.
///
/// Reference: Ono, I., & Kobayashi, S. (1999). A real-coded genetic algorithm
/// for function optimization using unimodal normal distribution crossover.
/// Proceedings of the 7th International Conference on Genetic Algorithms.

#include <algorithm>
#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Parent-Centric Crossover (UNDX) for real-valued encodings.
/// Uses 3 parents to define a distribution from which 2 offspring are sampled.
struct ParentCentricCrossover {
    double zeta = 0.5;                            ///< Spread along the line connecting parents
    double eta = 0.35;                            ///< Spread in orthogonal direction
    double crossover_probability = 0.9;           ///< Probability of applying crossover
    static constexpr double MIN_DISTANCE = 1e-10; ///< Minimum distance to avoid division by zero

    static constexpr int arity() { return 3; }
    static constexpr Encoding encoding() { return Encoding::Real; }

    /// Apply UNDX crossover. Produces 2 children starting at child_start.
    /// parent_a and parent_b define the primary axis; parent_c provides
    /// the orthogonal direction.
    void apply(this auto& self, Population<>& pop, int parent_a, int parent_b, int parent_c,
               int child_start) {
        auto& rng = Random::instance();

        // Copy parents as default
        for (int j = 0; j < pop.dim; ++j) {
            pop.gene(child_start, j) = pop.gene(parent_a, j);
            pop.gene(child_start + 1, j) = pop.gene(parent_b, j);
        }

        if (!rng.coin_flip(self.crossover_probability)) {
            pop.set_evaluated(child_start, false);
            pop.set_evaluated(child_start + 1, false);
            return;
        }

        // Calculate center of mass between parent_a and parent_b
        std::vector<double> center(pop.dim);
        for (int j = 0; j < pop.dim; ++j) {
            center[j] = (pop.gene(parent_a, j) + pop.gene(parent_b, j)) * 0.5;
        }

        // Calculate difference vector and distance between parent_a and parent_b
        std::vector<double> diff(pop.dim);
        double distance_sq = 0.0;
        for (int j = 0; j < pop.dim; ++j) {
            diff[j] = pop.gene(parent_b, j) - pop.gene(parent_a, j);
            distance_sq += diff[j] * diff[j];
        }
        double distance = std::sqrt(distance_sq);

        // If parents are too close, return copies
        if (distance < self.MIN_DISTANCE) {
            pop.set_evaluated(child_start, false);
            pop.set_evaluated(child_start + 1, false);
            return;
        }

        for (int j = 0; j < pop.dim; ++j) {
            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];

            // Generate values along the line connecting the parents
            double alpha = rng.normal(0.0, self.zeta * distance);

            // Generate values in the orthogonal direction
            double beta = 0.0;
            for (int k = 0; k < 2; ++k) {
                beta += rng.normal(0.0, self.eta * distance);
            }

            // Calculate the orthogonal component from parent_c
            double orthogonal = (pop.gene(parent_c, j) - center[j]) / distance;

            // Create the new values
            double c1 = center[j] + alpha * diff[j] / distance + beta * orthogonal;
            double c2 = center[j] - alpha * diff[j] / distance - beta * orthogonal;

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
