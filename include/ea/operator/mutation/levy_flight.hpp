#pragma once
/// @file levy_flight.hpp
/// @brief Lévy Flight Mutation — heavy-tailed steps for global exploration.
/// Reference: org.uma.jmetal.operator.mutation.impl.LevyFlightMutation

#include <cmath>
#include <algorithm>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Lévy Flight Mutation for real-valued encodings.
/// Uses the Mantegna algorithm for heavy-tailed step generation.
/// Mostly small steps with occasional large jumps — good for escaping local optima.
struct LevyFlightMutation {
    double beta = 1.5;             ///< Lévy index (1 < beta <= 2)
    double step_size = 0.01;       ///< Step scaling factor
    double mutation_rate = -1.0;   ///< Per-gene probability (-1 = 1/dim auto)

    static constexpr Encoding encoding() { return Encoding::Real; }

    void apply(this LevyFlightMutation& self, Population& pop, int idx) {
        auto& rng = Random::instance();
        double rate = self.mutation_rate < 0 ? 1.0 / pop.dim : self.mutation_rate;

        for (int j = 0; j < pop.dim; ++j) {
            if (!rng.coin_flip(rate)) continue;

            double y = pop.gene(idx, j);
            double lb = pop.lower_bounds[j];
            double ub = pop.upper_bounds[j];
            double range = ub - lb;

            if (range < 1e-14) continue;

            double levy_step = self.generate_levy_step(rng);
            double perturbation = levy_step * self.step_size * range;

            y += perturbation;
            pop.gene(idx, j) = std::clamp(y, lb, ub);
        }
        pop.set_evaluated(idx, false);
    }

private:
    /// Generate a Lévy-distributed step using the Mantegna algorithm.
    double generate_levy_step(this LevyFlightMutation& self, Random& rng) {
        // Calculate sigma_u for the Mantegna algorithm
        double beta = self.beta;
        double gamma_1pb = std::tgamma(1.0 + beta);
        double sin_term = std::sin(M_PI * beta / 2.0);
        double gamma_half = std::tgamma((1.0 + beta) / 2.0);
        double denom = gamma_half * beta * std::pow(2.0, (beta - 1.0) / 2.0);
        double sigma_u = std::pow(gamma_1pb * sin_term / denom, 1.0 / beta);

        // u ~ Normal(0, sigma_u^2)
        double u = rng.normal(0.0, sigma_u);
        // v ~ Normal(0, 1)
        double v = rng.normal(0.0, 1.0);

        // Lévy step = u / |v|^(1/beta)
        return u / std::pow(std::abs(v), 1.0 / beta);
    }
};

} // namespace ea
