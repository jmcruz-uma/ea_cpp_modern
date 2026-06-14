#pragma once
/// @file nsga3.hpp
/// @brief NSGA-III (Many-objective NSGA-II with reference points)
/// Reference: Deb, K. & Jain, H. "An Evolutionary Many-Objective Optimization
/// Algorithm Using Reference-Point-Based Nondominated Sorting Approach", IEEE TEVC 2014.
///
/// Bug history (corrected vs previous implementation):
///   1. Environmental selection used crowding distance (NSGA-II); must use
///      NSGAIIIReplacement — reference-point-based niching.
///   2. pop_size hardcoded to 100; must derive from reference point count:
///      compute_population_size(reference_points.size()).
///   3. reference_points vector generated but never passed to any selection.
///   4. Private generate_recursive duplicated generate_reference_points_das_dennis.
///   5. ideal_point maintained across generations but never used (dead code).

#include <algorithm>
#include <ea/core/comparator.hpp>
#include <ea/core/concepts.hpp>
#include <ea/core/population.hpp>
#include <ea/operator/replacement/nsga3_replacement.hpp>
#include <ea/util/random.hpp>
#include <ea/util/reference_point.hpp>
#include <iostream>
#include <numeric>
#include <string_view>
#include <vector>

namespace ea {

/// NSGA-III algorithm — many-objective optimization using reference points.
template <Crossover CX, Mutation MT> struct NSGAIII {
    CX crossover;
    MT mutation;

    /// Population size. 0 = auto-derive from reference point count (recommended).
    /// jMetal: setMaxPopulationSize(compute_population_size(referencePoints.size()))
    int pop_size = 0;
    int max_evals = 25000;
    int divisions = 12;       ///< Das-Dennis outer divisions
    int inner_divisions = 0;  ///< Second-layer divisions (0 = single layer)

    static constexpr std::string_view name() { return "NSGA-III"; }

    template <EvalFunctor F> void run(this auto& self, Population<>& pop, F&& problem) {
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;

        // Generate reference points (Das-Dennis, single or two-layer)
        std::vector<ReferencePoint> ref_pts;
        if (self.inner_divisions > 0)
            generate_reference_points_two_layer(ref_pts, n_obj, self.divisions,
                                                self.inner_divisions);
        else
            generate_reference_points_das_dennis(ref_pts, n_obj, self.divisions);

        // Derive population size from reference point count, matching jMetal:
        //   while (pop_size % 4 > 0) pop_size++;
        int n = self.pop_size;
        if (n <= 0)
            n = compute_population_size(static_cast<int>(ref_pts.size()));
        self.pop_size = n;

        if (pop.pop_size != n)
            pop.resize(n);

        // Evaluate initial population
        for (int i = 0; i < n; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }
        int evals = n;

        // Allocate offspring and combined populations
        Population<> offspring(n, dim, n_obj, pop.n_const);
        offspring.lower_bounds = pop.lower_bounds;
        offspring.upper_bounds = pop.upper_bounds;

        Population<> combined(2 * n, dim, n_obj, pop.n_const);
        combined.lower_bounds = pop.lower_bounds;
        combined.upper_bounds = pop.upper_bounds;

        // NSGAIIIReplacement is initialized once and reused every generation.
        // It clears per-generation state (member counts, potential members)
        // internally at the start of each replace() call.
        NSGAIIIReplacement replacement(ref_pts);

        auto& rng = Random::instance();

        while (evals < self.max_evals) {

            // ── Mating selection (binary tournament by NDS rank) ───────────
            auto fronts = fast_non_dominated_sort(pop);
            std::vector<int> ranks(n, 0);
            for (int r = 0; r < (int)fronts.size(); ++r)
                for (int idx : fronts[r])
                    ranks[idx] = r;

            std::vector<int> pool(2 * n);
            for (int i = 0; i < 2 * n; ++i) {
                int a = rng.uniform_int(0, n - 1);
                int b = rng.uniform_int(0, n - 1);
                pool[i] = (ranks[a] <= ranks[b]) ? a : b;
            }

            // ── Crossover + Mutation → Offspring ──────────────────────────
            for (int i = 0; i < n; i += 2) {
                int p1 = pool[i], p2 = pool[i + 1];
                std::copy_n(pop.genes_ptr(p1), dim, offspring.genes_ptr(i));
                std::copy_n(pop.genes_ptr(p2), dim, offspring.genes_ptr(i + 1));
                std::copy_n(pop.objectives_ptr(p1), n_obj, offspring.objectives_ptr(i));
                std::copy_n(pop.objectives_ptr(p2), n_obj, offspring.objectives_ptr(i + 1));
                offspring.set_evaluated(i, true);
                offspring.set_evaluated(i + 1, true);
            }

            for (int i = 0; i < n; i += 2)
                self.crossover.apply(offspring, i, i + 1, i);

            for (int i = 0; i < n; ++i)
                self.mutation.apply(offspring, i);

            // ── Evaluate offspring ─────────────────────────────────────────
            for (int i = 0; i < n && evals < self.max_evals; ++i) {
                problem(offspring, i);
                offspring.set_evaluated(i, true);
                ++evals;
            }

            // ── Combine parent + offspring ─────────────────────────────────
            for (int i = 0; i < n; ++i) {
                std::copy_n(pop.genes_ptr(i), dim, combined.genes_ptr(i));
                std::copy_n(pop.objectives_ptr(i), n_obj, combined.objectives_ptr(i));
                combined.flags[i] = pop.flags[i];
                std::copy_n(offspring.genes_ptr(i), dim, combined.genes_ptr(n + i));
                std::copy_n(offspring.objectives_ptr(i), n_obj, combined.objectives_ptr(n + i));
                combined.flags[n + i] = offspring.flags[i];
            }
            combined.pop_size = 2 * n;

            // ── Environmental selection (reference-point niching) ──────────
            auto sel = replacement.replace(combined, {}, n);

            for (int i = 0; i < n; ++i) {
                std::copy_n(combined.genes_ptr(sel[i]), dim, pop.genes_ptr(i));
                std::copy_n(combined.objectives_ptr(sel[i]), n_obj, pop.objectives_ptr(i));
                pop.flags[i] = combined.flags[sel[i]];
            }
        }
    }
};

} // namespace ea
