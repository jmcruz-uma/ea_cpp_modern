#pragma once
/// @file fuzzy_recombination.hpp
/// @brief Fuzzy Recombination (FR) Crossover — real-valued crossover using fuzzy connectives.
/// Reference: Herrera, F., Lozano, M., & Verdegay, J. L. (1998).
///        "Tackling real-coded genetic algorithms: Operators and tools for behavioural analysis"
///        Artificial Intelligence Review, 12(4), 265-319.
///
/// Uses a fuzzy membership function to blend parent values, with alpha controlling
/// the spread of the membership function.

#include <algorithm>
#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Fuzzy Recombination Crossover for real-valued encodings.
/// Produces two children from two parents using fuzzy connectives.
struct FuzzyRecombinationCrossover {
    double crossover_probability = 0.9; ///< Probability of applying crossover
    double alpha = 1.0;                 ///< Spread of fuzzy membership function

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Real; }

    void apply(this FuzzyRecombinationCrossover& self, Population& pop, int parent_a, int parent_b,
               int child_start) {
        auto& rng = Random::instance();
        constexpr double EPSILON = 1e-10;

        for (int j = 0; j < pop.dim; ++j) {
            double p1 = pop.gene(parent_a, j);
            double p2 = pop.gene(parent_b, j);

            if (std::abs(p1 - p2) < EPSILON) {
                pop.gene(child_start, j) = p1;
                pop.gene(child_start + 1, j) = p2;
                continue;
            }

            if (!rng.coin_flip(self.crossover_probability)) {
                pop.gene(child_start, j) = p1;
                pop.gene(child_start + 1, j) = p2;
                continue;
            }

            double min_val = std::min(p1, p2);
            double max_val = std::max(p1, p2);
            double range = max_val - min_val;
            double spread = self.alpha * range;

            double beta = rng.uniform();
            double c1 = min_val + beta * spread;
            double c2 = max_val - (1.0 - beta) * spread;

            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];

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
