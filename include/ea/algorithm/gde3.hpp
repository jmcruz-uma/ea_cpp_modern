#pragma once
/// @file gde3.hpp
/// @brief GDE3 (Generalized Differential Evolution 3)
/// Reference: Kukkonen, S. & Lampinen, J. "GDE3: A third step of differential
/// evolution", EMO 2005.
///
/// Differential evolution variant for multi-objective optimization.
/// Uses DE/rand/1/bin crossover and non-dominated sorting + crowding for selection.

#include <algorithm>
#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <ea/operator/replacement/nsga2_replacement.hpp>
#include <ea/util/random.hpp>
#include <string_view>
#include <vector>

namespace ea {

/// GDE3 — Generalized Differential Evolution 3.
/// Uses DE crossover operator (in ea::DECrossover) with NSGA-II style environmental selection.
template <typename CX, typename MT> struct GDE3 {
    CX crossover; // Should be DECrossover or similar
    MT mutation;

    int pop_size = 100;
    int max_evals = 25000;

    static constexpr std::string_view name() { return "GDE3"; }

    /// Run GDE3.
    template <typename Problem> void run(this auto& self, Population& pop, Problem&& problem) {
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;
        const int N = self.pop_size;

        if (pop.pop_size != N)
            pop.resize(N);

        // Evaluate initial population
        for (int i = 0; i < N; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }
        int evals = N;

        Population offspring(N, dim, n_obj, pop.encoding, pop.n_const);
        offspring.lower_bounds = pop.lower_bounds;
        offspring.upper_bounds = pop.upper_bounds;

        Population combined(2 * N, dim, n_obj, pop.encoding, pop.n_const);
        combined.lower_bounds = pop.lower_bounds;
        combined.upper_bounds = pop.upper_bounds;

        NSGAIIReplacement replacer;
        auto& rng = Random::instance();

        while (evals < self.max_evals) {
            // Create trial vectors (offspring) using DE crossover + mutation
            for (int i = 0; i < N; ++i) {
                // Copy current individual as base
                for (int j = 0; j < dim; ++j) {
                    offspring.gene(i, j) = pop.gene(i, j);
                }
            }

            // Apply DE crossover: crossover works on pop, creates children at indices
            // In GDE3, each individual generates one trial vector
            for (int i = 0; i < N; i += 2) {
                if (i + 1 < N) {
                    self.crossover.apply(pop, i, i + 1, i);
                }
            }

            // Apply mutation
            for (int i = 0; i < N; ++i) {
                self.mutation.apply(offspring, i);
            }

            // Evaluate offspring
            for (int i = 0; i < N; ++i) {
                if (!offspring.evaluated(i)) {
                    problem(offspring, i);
                    offspring.set_evaluated(i, true);
                    evals++;
                    if (evals >= self.max_evals)
                        break;
                }
            }

            // Combine parent + offspring
            for (int i = 0; i < N; ++i) {
                std::copy_n(pop.genes_ptr(i), dim, combined.genes_ptr(i));
                std::copy_n(pop.objectives_ptr(i), n_obj, combined.objectives_ptr(i));
                combined.flags[i] = pop.flags[i];

                std::copy_n(offspring.genes_ptr(i), dim, combined.genes_ptr(N + i));
                std::copy_n(offspring.objectives_ptr(i), n_obj, combined.objectives_ptr(N + i));
                combined.flags[N + i] = offspring.flags[i];
            }
            combined.pop_size = 2 * N;

            // Select best N using NSGA-II environmental selection
            auto selected = replacer.replace(combined, N);

            // Copy selected back to population
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < dim; ++j) {
                    pop.gene(i, j) = combined.gene(selected[i], j);
                }
                for (int o = 0; o < n_obj; ++o) {
                    pop.objective(i, o) = combined.objective(selected[i], o);
                }
                pop.flags[i] = combined.flags[selected[i]];
            }
        }
    }
};

} // namespace ea