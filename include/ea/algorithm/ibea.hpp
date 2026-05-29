#pragma once
/// @file ibea.hpp
/// @brief IBEA (Indicator-Based Evolutionary Algorithm)
///
/// Reference: E. Zitzler and S. Kunzli, "Indicator-Based Selection in Multiobjective Search"
/// PPSN 2004.
///
/// IBEA uses a binary quality indicator (typically Hypervolume or Epsilon) for selection
/// instead of Pareto dominance. This allows it to handle many objectives better than
/// dominance-based algorithms.

#include <algorithm>
#include <cmath>
#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <ea/indicator/hypervolume.hpp>
#include <ea/util/random.hpp>
#include <limits>
#include <vector>

namespace ea {

/// IBEA — Indicator-Based Evolutionary Algorithm.
///
/// Usage:
///   IBEA<ea::SBXCrossover, ea::PolynomialMutation> ibea;
///   ibea.pop_size = 100;
///   ibea.max_evals = 25000;
///   ibea.run(pop, problem);
template <typename CX, typename MT> struct IBEA {
    CX crossover;
    MT mutation;

    int pop_size = 100;
    int max_evals = 25000;
    double kappa = 0.05; ///< Scaling factor for indicator values

    static constexpr std::string_view name() { return "IBEA"; }

    void run(this auto& self, Population<>& pop, auto&& problem) {
        const int n = self.pop_size;
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;

        if (pop.pop_size != n)
            pop.resize(n);

        // === Evaluate initial population ===
        for (int i = 0; i < n; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }
        int evals = n;

        auto& rng = Random::instance();

        while (evals < self.max_evals) {
            // === 1. Generate offspring ===
            Population<> offspring(n, dim, n_obj);
            offspring.lower_bounds = pop.lower_bounds;
            offspring.upper_bounds = pop.upper_bounds;

            for (int i = 0; i < n; i += 2) {
                int p1 = rng.uniform_int(0, n - 1);
                int p2 = rng.uniform_int(0, n - 1);

                for (int j = 0; j < dim; ++j) {
                    offspring.gene(i, j) = pop.gene(p1, j);
                    if (i + 1 < n)
                        offspring.gene(i + 1, j) = pop.gene(p2, j);
                }

                self.crossover.apply(offspring, i, (i + 1 < n) ? i + 1 : i, i);
                self.mutation.apply(offspring, i);
                if (i + 1 < n)
                    self.mutation.apply(offspring, i + 1);
            }

            // Evaluate offspring
            for (int i = 0; i < n; ++i) {
                if (!offspring.evaluated(i)) {
                    problem(offspring, i);
                    offspring.set_evaluated(i, true);
                    evals++;
                    if (evals >= self.max_evals)
                        break;
                }
            }

            // === 2. Combine parent + offspring ===
            Population<> combined(2 * n, dim, n_obj);
            combined.lower_bounds = pop.lower_bounds;
            combined.upper_bounds = pop.upper_bounds;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < dim; ++j) {
                    combined.gene(i, j) = pop.gene(i, j);
                    combined.gene(n + i, j) = offspring.gene(i, j);
                }
                for (int o = 0; o < n_obj; ++o) {
                    combined.objective(i, o) = pop.objective(i, o);
                    combined.objective(n + i, o) = offspring.objective(i, o);
                }
            }
            combined.pop_size = 2 * n;

            // === 3. Compute fitness based on epsilon indicator ===
            // Binary epsilon indicator: I_epsilon+(a, b) = min {epsilon | a dominates b + epsilon}
            // We use additive epsilon
            std::vector<double> fitness(2 * n, 0.0);

            for (int i = 0; i < 2 * n; ++i) {
                double sum_indicator = 0.0;
                for (int j = 0; j < 2 * n; ++j) {
                    if (i == j)
                        continue;

                    // Compute epsilon indicator value I(i, j)
                    double epsilon = std::numeric_limits<double>::lowest();
                    for (int o = 0; o < n_obj; ++o) {
                        double val = combined.objective(i, o) - combined.objective(j, o);
                        if (val > epsilon)
                            epsilon = val;
                    }

                    // Indicator value: -exp(-I(i,j) / kappa)
                    // This gives higher values when i is better than j
                    double indicator_value = -std::exp(-epsilon / self.kappa);
                    sum_indicator += indicator_value;
                }
                fitness[i] = sum_indicator;
            }

            // === 4. Environmental selection: remove worst until n remain ===
            std::vector<bool> removed(2 * n, false);
            int remaining = 2 * n;

            while (remaining > n) {
                // Find individual with lowest fitness (worst)
                double worst_fitness = std::numeric_limits<double>::infinity();
                int worst_idx = -1;

                for (int i = 0; i < 2 * n; ++i) {
                    if (!removed[i] && fitness[i] < worst_fitness) {
                        worst_fitness = fitness[i];
                        worst_idx = i;
                    }
                }

                if (worst_idx < 0)
                    break;

                // Remove worst individual
                removed[worst_idx] = true;
                remaining--;

                // Update fitness of remaining individuals
                for (int i = 0; i < 2 * n; ++i) {
                    if (removed[i] || i == worst_idx)
                        continue;

                    // Recompute indicator contribution from i to worst
                    double epsilon = std::numeric_limits<double>::lowest();
                    for (int o = 0; o < n_obj; ++o) {
                        double val = combined.objective(i, o) - combined.objective(worst_idx, o);
                        if (val > epsilon)
                            epsilon = val;
                    }
                    double indicator_value = -std::exp(-epsilon / self.kappa);
                    fitness[i] -= indicator_value;
                }
            }

            // === 5. Copy selected back to population ===
            int idx = 0;
            for (int i = 0; i < 2 * n; ++i) {
                if (!removed[i]) {
                    for (int j = 0; j < dim; ++j) {
                        pop.gene(idx, j) = combined.gene(i, j);
                    }
                    for (int o = 0; o < n_obj; ++o) {
                        pop.objective(idx, o) = combined.objective(i, o);
                    }
                    pop.set_evaluated(idx, true);
                    idx++;
                }
            }
        }
    }
};

} // namespace ea
