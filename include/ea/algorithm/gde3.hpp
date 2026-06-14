#pragma once
/// @file gde3.hpp
/// @brief GDE3 (Generalized Differential Evolution 3)
/// Reference: Kukkonen, S. & Lampinen, J. "GDE3: A third step of differential
/// evolution", EMO 2005.
///
/// Differential evolution variant for multi-objective optimization.
/// Uses DE/rand/1/bin crossover with per-individual dominance test and
/// NSGA-II environmental selection.

#include <algorithm>
#include <ea/core/comparator.hpp>
#include <ea/core/concepts.hpp>
#include <ea/core/population.hpp>
#include <ea/operator/crossover/de.hpp>
#include <ea/operator/replacement/nsga2_replacement.hpp>
#include <ea/operator/selection/differential_evolution.hpp>
#include <ea/util/random.hpp>
#include <string_view>
#include <vector>

namespace ea {

/// GDE3 — Generalized Differential Evolution 3.
///
/// Each generation:
///   1. For each individual i, select r1, r2, r3 (random, distinct, ≠ i).
///   2. Create trial vector via DE/rand/1/bin:
///        trial[j] = r3[j] + F*(r1[j] − r2[j])  if rand<CR or j==jrand
///                   i[j]                          otherwise
///   3. Per-individual dominance test (parent vs trial):
///        parent dominates → keep only parent
///        trial dominates  → keep only trial
///        non-dominated    → keep both
///   4. NSGA-II environmental selection from the combined pool (N to 2N) → N.
struct GDE3 {
    DECrossover crossover; ///< DE operator — default: CR=0.5, F=0.5, Rand1Bin

    int pop_size  = 100;
    int max_evals = 25000;

    static constexpr std::string_view name() { return "GDE3"; }

    template <EvalFunctor F>
    void run(this auto& self, Population<>& pop, F&& problem) {
        const int dim   = pop.dim;
        const int n_obj = pop.n_obj;
        const int N     = self.pop_size;

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

        // Trial vector population (1 per individual per generation)
        Population<> offspring(N, dim, n_obj, pop.n_const);
        offspring.lower_bounds = pop.lower_bounds;
        offspring.upper_bounds = pop.upper_bounds;

        // Combined pool: at most 2N (may be smaller after dominance filtering)
        Population<> combined(2 * N, dim, n_obj, pop.n_const);
        combined.lower_bounds = pop.lower_bounds;
        combined.upper_bounds = pop.upper_bounds;

        DifferentialEvolutionSelection de_selector; // selects 3 distinct ≠ current
        NSGAIIReplacement replacer;
        auto& rng = Random::instance();
        std::vector<int> parents;

        while (evals < self.max_evals) {
            // ── Step 1: Generate trial vectors (DE/rand/1/bin) ─────────────────
            // Matches jMetal GDE3.selection() + GDE3.reproduction():
            //   selection: selectionOperator.setIndex(i); parents = execute(population)
            //   reproduction: currentSolution = population[i]; child = crossover.execute(parents)
            // child = copy of population[i] with some genes replaced by:
            //   mutant = parents[2] + F*(parents[0] − parents[1])
            for (int i = 0; i < N; ++i) {
                de_selector.select(pop, parents, i); // r1=parents[0], r2=parents[1], r3=parents[2]
                int r1 = parents[0], r2 = parents[1], r3 = parents[2];

                int jrand = rng.uniform_int(0, dim - 1);
                for (int j = 0; j < dim; ++j) {
                    double val;
                    if (rng.uniform() < self.crossover.cr || j == jrand) {
                        val = pop.gene(r3, j) +
                              self.crossover.f * (pop.gene(r1, j) - pop.gene(r2, j));
                        val = std::clamp(val, pop.lower_bounds[j], pop.upper_bounds[j]);
                    } else {
                        val = pop.gene(i, j);
                    }
                    offspring.gene(i, j) = val;
                }
                offspring.set_evaluated(i, false);
            }

            // ── Step 2: Evaluate trial vectors ─────────────────────────────────
            for (int i = 0; i < N && evals < self.max_evals; ++i) {
                problem(offspring, i);
                offspring.set_evaluated(i, true);
                ++evals;
            }

            // ── Step 3: Per-individual dominance test → build combined pool ────
            // Matches jMetal GDE3.replacement():
            //   for each i:
            //     result = dominanceComparator.compare(population[i], child[i])
            //     -1 → parent dominates → tmpList.add(parent)
            //      1 → child dominates  → tmpList.add(child)
            //      0 → non-dominated    → tmpList.add(child); tmpList.add(parent)
            int combined_size = 0;
            for (int i = 0; i < N; ++i) {
                // Inline cross-population Pareto dominance check:
                //   pop[i] (a) vs offspring[i] (b)
                bool a_leq_b = true, b_leq_a = true;
                for (int o = 0; o < n_obj; ++o) {
                    double fa = pop.objective(i, o);
                    double fb = offspring.objective(i, o);
                    if (fa > fb) a_leq_b = false;
                    if (fb > fa) b_leq_a = false;
                }
                bool a_eq_b = a_leq_b && b_leq_a;
                if (a_eq_b) {
                    // Identical objectives → non-dominated (add both)
                    combined.copy_from_other(offspring, i, combined_size++);
                    combined.copy_from_other(pop,      i, combined_size++);
                } else if (a_leq_b) {
                    // parent (a) dominates offspring (b)
                    combined.copy_from_other(pop, i, combined_size++);
                } else if (b_leq_a) {
                    // offspring (b) dominates parent (a)
                    combined.copy_from_other(offspring, i, combined_size++);
                } else {
                    // Non-dominated (each is better in at least one objective)
                    combined.copy_from_other(offspring, i, combined_size++);
                    combined.copy_from_other(pop,      i, combined_size++);
                }
            }
            combined.pop_size = combined_size;

            // ── Step 4: NSGA-II selection from combined (N..2N) → N ────────────
            // Matches jMetal: RankingAndCrowdingSelection.execute(tmpList) with N
            auto selected = replacer.replace(combined, N);

            for (int i = 0; i < N; ++i)
                pop.copy_from_other(combined, selected[i], i);
        }
    }
};

} // namespace ea
