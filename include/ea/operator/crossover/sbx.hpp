#pragma once
/// @file sbx.hpp
/// @brief Simulated Binary Crossover (SBX) — faithful to jMetal 7.4 SBXCrossover.java.
///
/// Protocol (matching jMetal exactly):
///   1. One uniform draw per pair: if > crossover_probability, return (offspring keep parent copies).
///   2. Per gene: one uniform draw against 0.5 threshold.
///      - If ≤ 0.5 AND |x1-x2| > EPS: compute c1/c2 via betaq, then 50% final swap.
///      - If ≤ 0.5 AND |x1-x2| ≤ EPS: no change (offspring already hold parent copies).
///      - If > 0.5: swap parent values between child_start and child_start+1.

#include <algorithm>
#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Simulated Binary Crossover (SBX) for real-valued encodings.
/// Caller must copy parent genes into child_start / child_start+1 before calling apply()
/// so that the "no crossover" and "identical genes" cases are handled by fallback.
struct SBXCrossover {
    double distribution_index    = 20.0; ///< η_c: higher → children closer to parents
    double crossover_probability = 0.9;  ///< Probability of applying crossover to the pair

    static constexpr int      arity()    { return 2; }
    static constexpr Encoding encoding() { return Encoding::Real; }

    void apply(this SBXCrossover& self, Population<>& pop, int parent_a, int parent_b,
               int child_start) {
        auto& rng = Random::instance();

        // One draw per pair — if fails, offspring already contain parent copies
        if (rng.uniform() > self.crossover_probability)
            return;

        for (int j = 0; j < pop.dim; ++j) {
            double value_x1 = pop.gene(parent_a, j);
            double value_x2 = pop.gene(parent_b, j);

            if (rng.coin_flip()) {
                // 50% path: apply SBX for this gene if parents differ
                if (std::abs(value_x1 - value_x2) > 1e-14) {
                    double y1, y2;
                    if (value_x1 < value_x2) { y1 = value_x1; y2 = value_x2; }
                    else                      { y1 = value_x2; y2 = value_x1; }

                    double lb   = pop.lower_bounds[j];
                    double ub   = pop.upper_bounds[j];
                    double rand = rng.uniform();

                    // betaq for child 1 (lower-bound scaling)
                    double beta  = 1.0 + (2.0 * (y1 - lb) / (y2 - y1));
                    double alpha = 2.0 - std::pow(beta, -(self.distribution_index + 1.0));
                    double betaq;
                    if (rand <= 1.0 / alpha)
                        betaq = std::pow(rand * alpha, 1.0 / (self.distribution_index + 1.0));
                    else
                        betaq = std::pow(1.0 / (2.0 - rand * alpha),
                                         1.0 / (self.distribution_index + 1.0));
                    double c1 = 0.5 * (y1 + y2 - betaq * (y2 - y1));

                    // betaq for child 2 (upper-bound scaling)
                    beta  = 1.0 + (2.0 * (ub - y2) / (y2 - y1));
                    alpha = 2.0 - std::pow(beta, -(self.distribution_index + 1.0));
                    if (rand <= 1.0 / alpha)
                        betaq = std::pow(rand * alpha, 1.0 / (self.distribution_index + 1.0));
                    else
                        betaq = std::pow(1.0 / (2.0 - rand * alpha),
                                         1.0 / (self.distribution_index + 1.0));
                    double c2 = 0.5 * (y1 + y2 + betaq * (y2 - y1));

                    c1 = std::clamp(c1, lb, ub);
                    c2 = std::clamp(c2, lb, ub);

                    if (rng.coin_flip()) std::swap(c1, c2);
                    pop.gene(child_start,     j) = c1;
                    pop.gene(child_start + 1, j) = c2;
                }
                // else: genes identical — offspring keep parent copies, nothing to write
            } else {
                // 50% path: swap parent values between the two offspring
                pop.gene(child_start,     j) = value_x2;
                pop.gene(child_start + 1, j) = value_x1;
            }
        }

        pop.set_evaluated(child_start,     false);
        pop.set_evaluated(child_start + 1, false);
    }
};

} // namespace ea
