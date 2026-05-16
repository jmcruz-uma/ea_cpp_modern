#pragma once
/// @file nsga3.hpp
/// @brief NSGA-III (Non-dominated Sorting Genetic Algorithm III)
///
/// NSGA-III is designed for many-objective optimization problems (>3 objectives).
/// It replaces NSGA-II's crowding distance with reference-point-based niching
/// to maintain diversity in high-dimensional objective spaces.
///
/// Key differences from NSGA-II:
/// - Uses Das-Dennis reference points for diversity preservation
/// - Normalizes objectives using ideal point and intercepts
/// - Niching-based selection from critical front (instead of crowding distance)
/// - No crowding distance computation
///
/// Reference: K. Deb and H. Jain, "An Evolutionary Many-Objective Optimization
/// Algorithm Using Reference-Point-Based Nondominated Sorting Approach, Part I:
/// Solving Problems With Box Constraints," IEEE TEVC, vol. 18, no. 4, 2014.
///
/// jMetal reference classes:
///   - org.uma.jmetal.algorithm.multiobjective.nsgaiii.NSGAIII
///   - org.uma.jmetal.algorithm.multiobjective.nsgaiii.NSGAIIIBuilder
///   - org.uma.jmetal.algorithm.multiobjective.nsgaiii.util.EnvironmentalSelection
///   - org.uma.jmetal.algorithm.multiobjective.nsgaiii.util.ReferencePoint
///   - org.uma.jmetal.component.catalogue.ea.replacement.impl.NSGAIIIReplacement
///
/// Template composition for zero-overhead dispatch when types are known at compile time.
/// Uses SoA Population, deducing this, and Concepts.

#include <ea/core/population.hpp>
#include <ea/core/comparator.hpp>
#include <ea/util/random.hpp>
#include <ea/util/reference_point.hpp>
#include <ea/operator/replacement/nsga3_replacement.hpp>
#include <string_view>
#include <vector>
#include <algorithm>

namespace ea {

/// NSGA-III algorithm — many-objective evolutionary algorithm with reference-point niching.
///
/// Usage:
///   NSGAIII nsga3;
///   nsga3.crossover.distribution_index = 20.0;
///   nsga3.mutation.distribution_index = 20.0;
///   nsga3.divisions = 12;        // Das-Dennis divisions
///   nsga3.max_evals = 25000;
///   nsga3.run(pop, problem);
///
/// For problems with many objectives (>5), use two-layer reference points:
///   nsga3.divisions = 3;
///   nsga3.second_layer_divisions = 2;
template<typename CX, typename MT>
struct NSGAIII {
    CX crossover;
    MT mutation;
    NSGAIIIReplacement replacement;

    /// Number of divisions for Das-Dennis reference point generation (p in the paper)
    int divisions = 12;

    /// Divisions for inner layer (0 = disabled, single layer)
    /// Recommended for problems with >5 objectives.
    int second_layer_divisions = 0;

    /// Maximum number of function evaluations
    int max_evals = 25000;

    static constexpr std::string_view name() { return "NSGA-III"; }

    /// Run NSGA-III on the given population.
    /// @param pop Population with genes initialized and bounds set.
    /// @param problem Callable: void(Population&, int) — evaluates individual's objectives
    template<typename Problem>
    void run(this auto& self, Population& pop, Problem&& problem) {
        const int n_obj = pop.n_obj;
        const int dim = pop.dim;

        // --- Step 1: Generate reference points ---
        std::vector<ReferencePoint> reference_points;
        if (self.second_layer_divisions > 0) {
            generate_reference_points_two_layer(
                reference_points, n_obj, self.divisions, self.second_layer_divisions);
        } else {
            generate_reference_points_das_dennis(
                reference_points, n_obj, self.divisions);
        }

        // Set up replacement operator
        self.replacement.set_reference_points(std::move(reference_points));

        // Compute actual population size (round up to multiple of 4)
        int pop_size = compute_population_size(
            static_cast<int>(self.replacement.reference_points.size()));

        // Ensure population matches computed size
        if (pop.pop_size != pop_size) pop.resize(pop_size);

        // --- Step 2: Evaluate initial population ---
        for (int i = 0; i < pop_size; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }
        int evals = pop_size;

        // Allocate offspring + combined populations
        Population offspring(pop_size, dim, n_obj, pop.encoding, pop.n_const);
        offspring.lower_bounds = pop.lower_bounds;
        offspring.upper_bounds = pop.upper_bounds;

        Population combined(2 * pop_size, dim, n_obj, pop.encoding, pop.n_const);
        combined.lower_bounds = pop.lower_bounds;
        combined.upper_bounds = pop.upper_bounds;

        while (evals < self.max_evals) {
            // === 1. Non-dominated sort (no crowding distance for NSGA-III) ===
            auto fronts = fast_non_dominated_sort(pop);

            // NSGA-III doesn't use crowding distance — it uses reference points
            // Ranks are sufficient for selection
            std::vector<int> ranks(pop_size, 0);
            for (int r = 0; r < static_cast<int>(fronts.size()); ++r) {
                for (int idx : fronts[r]) {
                    ranks[idx] = r;
                }
            }

            // === 2. Selection (binary tournament) ===
            auto& rng = Random::instance();
            std::vector<int> mating_pool(2 * pop_size);
            for (int i = 0; i < 2 * pop_size; ++i) {
                int a = rng.uniform_int(0, pop_size - 1);
                int b = rng.uniform_int(0, pop_size - 1);
                if (ranks[a] < ranks[b]) {
                    mating_pool[i] = a;
                } else if (ranks[a] > ranks[b]) {
                    mating_pool[i] = b;
                } else {
                    // Same rank — random tiebreak (no crowding distance in NSGA-III)
                    mating_pool[i] = (rng.uniform() < 0.5) ? a : b;
                }
            }

            // === 3. Crossover + Mutation → Offspring ===
            for (int i = 0; i < pop_size; i += 2) {
                int p1 = mating_pool[i];
                int p2 = mating_pool[i + 1];

                // Copy parents into offspring positions first
                for (int j = 0; j < dim; ++j) {
                    offspring.gene(i, j) = pop.gene(p1, j);
                    offspring.gene(i + 1, j) = pop.gene(p2, j);
                }
                std::copy_n(pop.objectives_ptr(p1), n_obj, offspring.objectives_ptr(i));
                std::copy_n(pop.objectives_ptr(p2), n_obj, offspring.objectives_ptr(i + 1));
                offspring.set_evaluated(i, true);
                offspring.set_evaluated(i + 1, true);
            }

            // Apply crossover
            for (int i = 0; i < pop_size; i += 2) {
                self.crossover.apply(offspring, i, i + 1, i);
            }

            // Apply mutation
            for (int i = 0; i < pop_size; ++i) {
                self.mutation.apply(offspring, i);
            }

            // Evaluate offspring
            for (int i = 0; i < pop_size; ++i) {
                if (!offspring.evaluated(i)) {
                    problem(offspring, i);
                    offspring.set_evaluated(i, true);
                    evals++;
                    if (evals >= self.max_evals) break;
                }
            }

            // === 4. Combine parent + offspring ===
            for (int i = 0; i < pop_size; ++i) {
                // Parents
                std::copy_n(pop.genes_ptr(i), dim, combined.genes_ptr(i));
                std::copy_n(pop.objectives_ptr(i), n_obj, combined.objectives_ptr(i));
                combined.flags[i] = pop.flags[i];
                // Offspring
                std::copy_n(offspring.genes_ptr(i), dim, combined.genes_ptr(pop_size + i));
                std::copy_n(offspring.objectives_ptr(i), n_obj, combined.objectives_ptr(pop_size + i));
                combined.flags[pop_size + i] = offspring.flags[i];
            }
            combined.pop_size = 2 * pop_size;

            // === 5. Environmental selection (NSGA-III replacement) ===
            auto selected = self.replacement.replace(combined, pop_size);

            // Copy selected back to population
            for (int i = 0; i < pop_size; ++i) {
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
