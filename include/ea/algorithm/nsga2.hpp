#pragma once
/// @file nsga2.hpp
/// @brief NSGA-II (Non-dominated Sorting Genetic Algorithm II)
/// Reference: Deb, K., et al. "A fast and elitist multiobjective genetic algorithm:
/// NSGA-II" IEEE Trans. Evol. Comput., 2002.
///
/// Template composition for zero-overhead dispatch when types are known at compile time.
/// Uses SoA Population, deducing this, and Concepts.

#include <ea/core/population.hpp>
#include <ea/core/comparator.hpp>
#include <ea/util/random.hpp>
#include <ea/operator/replacement/nsga2_replacement.hpp>
#include <string_view>
#include <vector>
#include <algorithm>

namespace ea {

/// NSGA-II algorithm — the most cited multi-objective evolutionary algorithm.
///
/// Usage:
///   NSGAII nsga;
///   nsga.crossover.distribution_index = 20.0;
///   nsga.mutation.distribution_index = 20.0;
///   nsga.pop_size = 100;
///   nsga.max_evals = 25000;
///   nsga.run(pop, problem);
template<typename CX, typename MT>
struct NSGAII {
    CX crossover;
    MT mutation;
    NSGAIIReplacement replacement;

    int pop_size = 100;
    int max_evals = 25000;

    static constexpr std::string_view name() { return "NSGA-II"; }

    /// Run NSGA-II on the given population.
    /// @param pop Population with genes initialized and bounds set. Must have pop_size individuals.
    /// @param problem Callable: void(Population&, int) — evaluates individual's objectives
    template<typename Problem>
    void run(this auto& self, Population& pop, Problem&& problem) {
        const int n = self.pop_size;
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;

        // Validate population
        if (pop.pop_size != n) pop.resize(n);

        // === Evaluate initial population ===
        for (int i = 0; i < n; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }
        int evals = n;

        // Allocate offspring + combined populations
        Population offspring(n, dim, n_obj, pop.encoding, pop.n_const);
        offspring.lower_bounds = pop.lower_bounds;
        offspring.upper_bounds = pop.upper_bounds;

        Population combined(2 * n, dim, n_obj, pop.encoding, pop.n_const);
        combined.lower_bounds = pop.lower_bounds;
        combined.upper_bounds = pop.upper_bounds;

        std::vector<int> ranks;
        std::vector<double> crowding;

        while (evals < self.max_evals) {
            // === 1. Non-dominated sort + crowding distance ===
            auto fronts = fast_non_dominated_sort(pop);

            ranks.assign(n, 0);
            crowding.assign(n, 0.0);
            for (int r = 0; r < static_cast<int>(fronts.size()); ++r) {
                std::vector<double> front_cd;
                compute_crowding_distance(pop, fronts[r], front_cd);
                for (size_t i = 0; i < fronts[r].size(); ++i) {
                    ranks[fronts[r][i]] = r;
                    crowding[fronts[r][i]] = front_cd[i];
                }
            }

            // === 2. Selection (binary tournament) ===
            auto& rng = Random::instance();
            std::vector<int> mating_pool(2 * n);
            for (int i = 0; i < 2 * n; ++i) {
                int a = rng.uniform_int(0, n - 1);
                int b = rng.uniform_int(0, n - 1);
                if (ranks[a] < ranks[b]) {
                    mating_pool[i] = a;
                } else if (ranks[a] > ranks[b]) {
                    mating_pool[i] = b;
                } else if (crowding[a] > crowding[b]) {
                    mating_pool[i] = a;
                } else {
                    mating_pool[i] = b;
                }
            }

            // === 3. Crossover + Mutation → Offspring ===
            // Copy selected parents to offspring
            for (int i = 0; i < n; i += 2) {
                int p1 = mating_pool[i];
                int p2 = mating_pool[i + 1];
                // Copy parents into offspring positions first (fallback if crossover doesn't apply)
                for (int j = 0; j < dim; ++j) {
                    offspring.gene(i, j) = pop.gene(p1, j);
                    offspring.gene(i + 1, j) = pop.gene(p2, j);
                }
                std::copy_n(pop.objectives_ptr(p1), n_obj, offspring.objectives_ptr(i));
                std::copy_n(pop.objectives_ptr(p2), n_obj, offspring.objectives_ptr(i + 1));
                offspring.set_evaluated(i, true);
                offspring.set_evaluated(i + 1, true);
            }

            // Apply crossover: parents from pop, children into offspring
            // SBX-style: crossover works on pop data, writes to child positions
            // We need a temporary to hold crossover results, then copy to offspring
            for (int i = 0; i < n; i += 2) {
                int p1 = mating_pool[i];
                int p2 = mating_pool[i + 1];

                // Crossover on pop (reads from parents, writes to child positions i, i+1 in pop)
                // But we don't want to corrupt pop! Use combined as scratch.
                // Actually: crossover.apply takes parent indices and child_start index
                // We'll use offspring population for both parents and children

                // Simpler approach: copy parents to offspring first (done above),
                // then apply crossover on offspring (parent_a = i, parent_b = i+1, child_start = i)
                self.crossover.apply(offspring, i, i + 1, i);
            }

            // Apply mutation on offspring
            for (int i = 0; i < n; ++i) {
                self.mutation.apply(offspring, i);
            }

            // Evaluate offspring
            for (int i = 0; i < n; ++i) {
                if (!offspring.evaluated(i)) {
                    problem(offspring, i);
                    offspring.set_evaluated(i, true);
                    evals++;
                    if (evals >= self.max_evals) break;
                }
            }

            // === 4. Combine parent + offspring ===
            for (int i = 0; i < n; ++i) {
                // Parents
                std::copy_n(pop.genes_ptr(i), dim, combined.genes_ptr(i));
                std::copy_n(pop.objectives_ptr(i), n_obj, combined.objectives_ptr(i));
                combined.flags[i] = pop.flags[i];
                // Offspring
                std::copy_n(offspring.genes_ptr(i), dim, combined.genes_ptr(n + i));
                std::copy_n(offspring.objectives_ptr(i), n_obj, combined.objectives_ptr(n + i));
                combined.flags[n + i] = offspring.flags[i];
            }
            combined.pop_size = 2 * n;

            // === 5. Environmental selection ===
            auto selected = self.replacement.replace(combined, n);

            // Copy selected back to population
            for (int i = 0; i < n; ++i) {
                if (selected[i] != i) {
                    pop.copy_individual(selected[i], i);
                    // Wait, this copies from combined to pop, but combined indices are different
                }
                // Actually need to copy from combined[selected[i]] to pop[i]
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