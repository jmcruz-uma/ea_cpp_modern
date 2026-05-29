#pragma once
/// @file gaussian.hpp
/// @brief Gaussian Mutation for real-valued encodings.
///
/// Replicates the operator from Bio4Res/ea (ea_cpp_original):
///   GaussianMutation.h — RegisteredExtension<VariationOperator, GaussianMutation, "GAUSSIAN">
///
/// Semantics (faithful to the original):
///   - With probability mutation_rate, the operator fires for the individual.
///   - When it fires, ONE gene is chosen uniformly at random and perturbed
///     with *multiplicative* Gaussian noise:
///       gene[j] = gene[j] * (1 + step_size * N(0, 1))
///   - Result is clamped to [lower_bound[j], upper_bound[j]].
///
/// Known limitation (inherited from the original):
///   If gene[j] == 0.0 the multiplicative perturbation produces 0.0 again.
///   For problems whose optimal lies at zero (e.g. Sphere), convergence near
///   the optimum requires that other genes are non-zero or that a additive
///   fallback is used. We replicate the original behaviour to allow fair
///   comparison; do not use this operator as a general-purpose mutation.

#include <algorithm>
#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Gaussian Mutation for real-valued encodings.
///
/// Matches the Bio4Res/ea original operator exactly:
///   - Fires with per-individual probability @p mutation_rate.
///   - Perturbs a single randomly-chosen gene multiplicatively.
///
/// Typical configuration matching ea_cpp_original's numeric.json:
/// @code
///   GaussianMutation mut;
///   mut.mutation_rate = 0.6321;   // pars[0] in original
///   mut.step_size     = 1.0;      // pars[1] in original
/// @endcode
struct GaussianMutation {
    /// Probability that the operator fires for a given individual.
    /// -1.0 → automatic: 1/dim (one expected application per individual).
    double mutation_rate = -1.0;

    /// σ scaling factor for the Gaussian perturbation.
    /// gene' = gene × (1 + step_size × N(0,1)).
    double step_size = 1.0;

    static constexpr Encoding encoding() { return Encoding::Real; }

    void apply(this GaussianMutation& self, Population<>& pop, int idx) {
        auto& rng = Random::instance();

        const double rate = self.mutation_rate < 0.0
                                ? 1.0 / static_cast<double>(pop.dim)
                                : self.mutation_rate;

        if (!rng.coin_flip(rate))
            return;

        // Choose one gene uniformly at random — matches original apply_()
        const int j  = rng.uniform_int(0, pop.dim - 1);
        const double v  = pop.gene(idx, j);
        const double lb = pop.lower_bounds[j];
        const double ub = pop.upper_bounds[j];

        // Multiplicative perturbation — matches original:
        //   val = v * (1.0 + stepsize * nrandom())
        const double val = v * (1.0 + self.step_size * rng.normal(0.0, 1.0));
        pop.gene(idx, j) = std::clamp(val, lb, ub);
        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
