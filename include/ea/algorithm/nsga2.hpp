#pragma once
/// @file nsga2.hpp
/// @brief NSGA-II (Non-dominated Sorting Genetic Algorithm II)
/// Reference: Deb, K., et al. "A fast and elitist multiobjective genetic algorithm: NSGA-II"
/// IEEE Trans. Evol. Comput., 2002.

#include <ea/core/population.hpp>
#include <ea/core/comparator.hpp>
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/operator/selection/binary_tournament.hpp>
#include <ea/util/random.hpp>
#include <string_view>
#include <vector>

namespace ea {

/// NSGA-II algorithm — template-based composition for zero-overhead dispatch.
template<typename CX = SBXCrossover,
         typename MT = PolynomialMutation,
         typename SEL = BinaryTournamentSelection>
struct NSGAII {
    CX crossover;
    MT mutation;
    SEL selection;

    int pop_size = 100;
    int max_evals = 25000;
    double crossover_prob = 0.9;
    double mutation_prob = -1.0; // auto: 1/dim

    static constexpr std::string_view name() { return "NSGA-II"; }

    /// Run NSGA-II on the given population.
    /// The population must be initialized (genes set, bounds set, dim/n_obj set).
    /// @param pop Population (modified in-place with final Pareto front)
    /// @param problem A callable that evaluates objectives: void(Population&, int)
    template<typename Problem>
    void run(this auto& self, Population& pop, Problem&& problem) {
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;

        // Allocate offspring population
        Population offspring(self.pop_size, dim, n_obj, pop.encoding, pop.n_const);
        offspring.lower_bounds = pop.lower_bounds;
        offspring.upper_bounds = pop.upper_bounds;

        // Rank and crowding distance storage
        std::vector<int> ranks(pop.pop_size);
        std::vector<double> crowding_dist(pop.pop_size);

        int evals = pop.pop_size; // Initial population evaluated

        // Evaluate initial population
        for (int i = 0; i < pop.pop_size; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }

        while (evals < self.max_evals) {
            // 1. Non-dominated sort + crowding distance
            auto fronts = fast_non_dominated_sort(pop);

            // Assign ranks
            for (int rank = 0; rank < static_cast<int>(fronts.size()); ++rank) {
                for (int idx : fronts[rank]) {
                    ranks[idx] = rank;
                }
            }

            // Compute crowding distance for each front
            for (auto& front : fronts) {
                compute_crowding_distance(pop, front, crowding_dist);
                // Map back to individual indices
                // crowding_dist is indexed by position in front, not individual
                // We need per-individual crowding distance for tournament
            }

            // Create per-individual crowding distance
            std::vector<double> ind_crowding(pop.pop_size, 0.0);
            for (auto& front : fronts) {
                std::vector<double> front_cd;
                compute_crowding_distance(pop, front, front_cd);
                for (size_t i = 0; i < front.size(); ++i) {
                    ind_crowding[front[i]] = front_cd[i];
                }
            }

            // 2. Selection — binary tournament
            std::vector<int> mating_pool;
            self.selection.select(pop, mating_pool, ranks, ind_crowding);

            // 3. Crossover + Mutation → offspring
            offspring.pop_size = self.pop_size;
            offspring.resize(self.pop_size);

            for (int i = 0; i < self.pop_size; i += 2) {
                int p1 = mating_pool[i];
                int p2 = mating_pool[i + 1];

                // Copy parents to offspring positions (in case crossover doesn't apply)
                pop.copy_individual(p1, i);
                pop.copy_individual(p2, i + 1);

                // We need to work on offspring, not pop
                // Copy parents into offspring first
                offspring.copy_individual(p1, i);
                offspring.copy_individual(p2, i + 1);

                // Crossover
                self.crossover.apply(pop, p1, p2, i);

                // ... Actually, let's use offspring directly
                // Copy parents into offspring
                for (int j = 0; j < dim; ++j) {
                    offspring.gene(i, j) = pop.gene(p1, j);
                    offspring.gene(i + 1, j) = pop.gene(p2, j);
                }

                // Crossover modifies offspring positions i and i+1
                // We need a combined population approach
            }

            // Simplified approach: work with combined parent+offspring population
            // ... (see below)
        }
    }
};

} // namespace ea