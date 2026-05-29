#pragma once
/// @file sus.hpp
/// @brief Stochastic Universal Sampling (SUS) — deterministic version of roulette wheel.
///
/// Reference: jMetal
///   jmetal-core/src/main/java/org/uma/jmetal/operator/selection/impl/StochasticUniversalSampling.java
///
/// SUS uses a single random starting point and selects individuals at regular
/// intervals based on cumulative fitness. Reduces variance compared to standard
/// roulette wheel sampling, ensuring better representation of high-fitness individuals.
///
/// For minimization problems, fitness is inverted: selection probability ∝ 1/fitness.

#include <algorithm>
#include <cmath>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <numeric>
#include <vector>

namespace ea {

/// Stochastic Universal Sampling Selection.
/// Deterministic alternative to roulette wheel that maintains diversity.
struct StochasticUniversalSampling {
    /// Select mating pool using SUS.
    /// @param pop          Population<>
    /// @param mating_pool  Output indices (resized to 2 * pop_size)
    /// @param fitness      Fitness of each individual (lower = better for minimization)
    void select(this StochasticUniversalSampling& self, Population<>& pop,
                std::vector<int>& mating_pool, const std::vector<double>& fitness) {
        (void)self;
        auto& rng = Random::instance();
        const int n = pop.pop_size;
        mating_pool.resize(2 * n);

        // Convert minimization fitness to maximization selection probability
        // Invert and shift to positive: w_i = max_fitness - fitness_i + epsilon
        double max_fit = *std::max_element(fitness.begin(), fitness.end());
        double min_fit = *std::min_element(fitness.begin(), fitness.end());
        double shift = (max_fit == min_fit) ? 1.0 : 0.0;
        double range = max_fit - min_fit + shift;

        std::vector<double> weights(n);
        double total_weight = 0.0;
        for (int i = 0; i < n; ++i) {
            // Higher weight for lower fitness (minimization)
            weights[i] = (max_fit - fitness[i]) + 1e-12;
            total_weight += weights[i];
        }

        // Build cumulative distribution
        std::vector<double> cumulative(n);
        double acc = 0.0;
        for (int i = 0; i < n; ++i) {
            acc += weights[i] / total_weight;
            cumulative[i] = acc;
        }
        cumulative[n - 1] = 1.0; // Ensure exact 1.0 at end

        // SUS: single random start, evenly spaced pointers
        double step = 1.0 / (2 * n);
        double start = rng.uniform(0.0, step);

        for (int i = 0; i < 2 * n; ++i) {
            double pointer = start + i * step;
            if (pointer > 1.0)
                pointer -= 1.0;

            // Find first cumulative >= pointer
            auto it = std::lower_bound(cumulative.begin(), cumulative.end(), pointer);
            int selected = static_cast<int>(it - cumulative.begin());
            if (selected >= n)
                selected = n - 1;

            mating_pool[i] = selected;
        }
    }
};

} // namespace ea
