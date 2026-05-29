#pragma once
/// @file weighted.hpp
/// @brief Weighted Selection (Roulette Wheel) — fitness-proportionate selection.
///
/// Reference: jMetal
///   jmetal-core/src/main/java/org/uma/jmetal/operator/selection/impl/RouletteWheelSelection.java
///
/// Standard roulette-wheel selection where probability of selection is proportional
/// to fitness. For minimization, fitness is inverted so lower fitness → higher weight.
///
/// Note: For MOEA, typically used with scalarized fitness (aggregation function output).

#include <algorithm>
#include <cmath>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <numeric>
#include <vector>

namespace ea {

/// Weighted / Roulette-Wheel Selection.
/// Probability proportional to inverted fitness (for minimization).
struct WeightedSelection {
    /// Select mating pool using fitness-proportionate selection.
    /// @param pop          Population<>
    /// @param mating_pool  Output indices (resized to 2 * pop_size)
    /// @param fitness      Fitness of each individual (lower = better)
    void select(this WeightedSelection& self, Population<>& pop,
                std::vector<int>& mating_pool, const std::vector<double>& fitness) {
        (void)self;
        auto& rng = Random::instance();
        const int n = pop.pop_size;
        mating_pool.resize(2 * n);

        // Invert fitness for minimization → maximization weights
        double max_fit = *std::max_element(fitness.begin(), fitness.end());
        double min_fit = *std::min_element(fitness.begin(), fitness.end());
        double shift = (max_fit == min_fit) ? 1.0 : 0.0;

        std::vector<double> weights(n);
        double total_weight = 0.0;
        for (int i = 0; i < n; ++i) {
            weights[i] = (max_fit - fitness[i]) + shift + 1e-12;
            total_weight += weights[i];
        }

        // Build cumulative distribution
        std::vector<double> cumulative(n);
        double acc = 0.0;
        for (int i = 0; i < n; ++i) {
            acc += weights[i] / total_weight;
            cumulative[i] = acc;
        }
        cumulative[n - 1] = 1.0;

        // Spin the wheel 2N times
        for (int i = 0; i < 2 * n; ++i) {
            double spin = rng.uniform(0.0, 1.0);
            auto it = std::lower_bound(cumulative.begin(), cumulative.end(), spin);
            int selected = static_cast<int>(it - cumulative.begin());
            if (selected >= n)
                selected = n - 1;
            mating_pool[i] = selected;
        }
    }
};

} // namespace ea
