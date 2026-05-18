#pragma once
/// @file integer_sbx.hpp
/// @brief Integer Simulated Binary Crossover (SBX) — SBX for integer encodings.
/// Reference: jMetal IntegerSBXCrossover.java
///
/// Applies standard SBX crossover and rounds the result to the nearest integer.
/// The algorithm is identical to real-valued SBX but with integer bounds and rounding.

#include <algorithm>
#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Integer SBX Crossover for integer encodings.
/// Produces two children from two parents using SBX, then rounds to integers.
struct IntegerSBXCrossover {
    double distribution_index = 20.0;   ///< Distribution index (η)
    double crossover_probability = 0.9; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Integer; }

    void apply(this IntegerSBXCrossover& self, Population& pop, int parent_a, int parent_b,
               int child_start) {
        auto& rng = Random::instance();
        constexpr double EPS = 1.0e-14;

        for (int j = 0; j < pop.dim; ++j) {
            double p1 = pop.gene(parent_a, j);
            double p2 = pop.gene(parent_b, j);

            if (std::abs(p1 - p2) < EPS) {
                pop.gene(child_start, j) = p1;
                pop.gene(child_start + 1, j) = p2;
                continue;
            }

            if (!rng.coin_flip(self.crossover_probability)) {
                pop.gene(child_start, j) = p1;
                pop.gene(child_start + 1, j) = p2;
                continue;
            }

            double y1 = std::min(p1, p2);
            double y2 = std::max(p1, p2);
            double yl = pop.lower_bounds[j];
            double yu = pop.upper_bounds[j];
            double rand_val = rng.uniform();

            // Compute beta for child 1 (lower side)
            double beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1));
            double alpha_val = 2.0 - std::pow(beta, -(self.distribution_index + 1.0));
            double betaq;

            if (rand_val <= (1.0 / alpha_val)) {
                betaq = std::pow(rand_val * alpha_val, 1.0 / (self.distribution_index + 1.0));
            } else {
                betaq = std::pow(1.0 / (2.0 - rand_val * alpha_val),
                                 1.0 / (self.distribution_index + 1.0));
            }

            double c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));

            // Compute beta for child 2 (upper side)
            beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
            alpha_val = 2.0 - std::pow(beta, -(self.distribution_index + 1.0));

            if (rand_val <= (1.0 / alpha_val)) {
                betaq = std::pow(rand_val * alpha_val, 1.0 / (self.distribution_index + 1.0));
            } else {
                betaq = std::pow(1.0 / (2.0 - rand_val * alpha_val),
                                 1.0 / (self.distribution_index + 1.0));
            }

            double c2 = 0.5 * (y1 + y2 + betaq * (y2 - y1));

            // Clamp
            c1 = std::clamp(c1, yl, yu);
            c2 = std::clamp(c2, yl, yu);

            // Round to integers
            int ic1 = static_cast<int>(std::round(c1));
            int ic2 = static_cast<int>(std::round(c2));

            int lb_int = static_cast<int>(std::ceil(yl));
            int ub_int = static_cast<int>(std::floor(yu));
            ic1 = std::clamp(ic1, lb_int, ub_int);
            ic2 = std::clamp(ic2, lb_int, ub_int);

            // Swap with probability 0.5
            if (rng.coin_flip(0.5)) {
                std::swap(ic1, ic2);
            }

            pop.gene(child_start, j) = static_cast<double>(ic1);
            pop.gene(child_start + 1, j) = static_cast<double>(ic2);
        }

        pop.set_evaluated(child_start, false);
        pop.set_evaluated(child_start + 1, false);
    }
};

} // namespace ea
