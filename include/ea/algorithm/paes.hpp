#pragma once
/// @file paes.hpp
/// @brief PAES (Pareto Archived Evolution Strategy)
/// Reference: Knowles, J.D. & Corne, D.W. "Approximating the nondominated front
/// using the Pareto Archived Evolution Strategy", Evolutionary Computation, 2000.
///
/// (1+1)-ES with adaptive grid archive. Mutates a single parent, compares the
/// mutant with the parent using dominance, and updates both the current solution
/// and the external archive based on grid density.

#include <algorithm>
#include <cmath>
#include <ea/core/comparator.hpp>
#include <ea/core/concepts.hpp>
#include <ea/core/population.hpp>
#include <ea/util/adaptive_grid.hpp>
#include <ea/util/random.hpp>
#include <string_view>
#include <vector>

namespace ea {

/// PAES — Pareto Archived Evolution Strategy.
/// (1+1)-ES with adaptive grid archive for diversity.
/// Template param: MT — mutation operator only (no crossover in PAES).
template <Mutation MT> struct PAES {
    MT mutation;

    int max_evals = 25000;
    int archive_size = 100; ///< External archive capacity
    int bisections = 5;     ///< Grid bisections per objective

    static constexpr std::string_view name() { return "PAES"; }

    /// Run PAES.
    template <EvalFunctor F> void run(this auto& self, Population<>& pop, F&& problem) {
        const int n_obj = pop.n_obj;
        const int dim = pop.dim;

        // PAES works with a single individual — we use pop index 0
        if (pop.pop_size < 1)
            pop.resize(1);

        // Initialize if needed
        auto& rng = Random::instance();
        if (!pop.evaluated(0)) {
            for (int j = 0; j < dim; ++j)
                pop.gene(0, j) = rng.uniform(pop.lower_bounds[j], pop.upper_bounds[j]);
            problem(pop, 0);
            pop.set_evaluated(0, true);
        }

        int evals = 1;

        // Create adaptive grid archive and add initial solution
        AdaptiveGridArchive archive(self.archive_size, self.bisections, n_obj);
        {
            std::vector<double> init_genes(dim), init_obj(n_obj);
            for (int j = 0; j < dim; ++j) init_genes[j] = pop.gene(0, j);
            for (int o = 0; o < n_obj; ++o) init_obj[o] = pop.objective(0, o);
            archive.add(init_genes, init_obj);
        }

        // Working populations: current solution and its mutant
        Population<> current(1, dim, n_obj, pop.n_const);
        current.lower_bounds = pop.lower_bounds;
        current.upper_bounds = pop.upper_bounds;
        for (int j = 0; j < dim; ++j) current.gene(0, j) = pop.gene(0, j);
        for (int o = 0; o < n_obj; ++o) current.objective(0, o) = pop.objective(0, o);
        current.flags[0] = pop.flags[0];

        Population<> mutant(1, dim, n_obj, pop.n_const);
        mutant.lower_bounds = pop.lower_bounds;
        mutant.upper_bounds = pop.upper_bounds;

        // Pre-allocate work vectors outside the loop to avoid per-iteration heap allocs
        std::vector<double> mut_genes(dim), mut_obj(n_obj), curr_obj(n_obj);

        while (evals < self.max_evals) {
            // === 1. Mutate current solution ===
            for (int j = 0; j < dim; ++j)
                mutant.gene(0, j) = current.gene(0, j);
            mutant.set_evaluated(0, false);
            self.mutation.apply(mutant, 0);

            // === 2. Evaluate mutant ===
            problem(mutant, 0);
            mutant.set_evaluated(0, true);
            evals++;

            // === 3. Capture objective/gene vectors for archive operations ===
            for (int j = 0; j < dim; ++j) mut_genes[j] = mutant.gene(0, j);
            for (int o = 0; o < n_obj; ++o) mut_obj[o] = mutant.objective(0, o);
            for (int o = 0; o < n_obj; ++o) curr_obj[o] = current.objective(0, o);

            // === 4. Replacement (jMetal PAES.replacement) ===
            auto dom = self.compare_dominance(mutant, 0, current, 0);

            if (dom == Dominance::Dominates) {
                // Mutant dominates current: accept mutant, unconditionally add to archive
                for (int j = 0; j < dim; ++j) current.gene(0, j) = mutant.gene(0, j);
                for (int o = 0; o < n_obj; ++o) current.objective(0, o) = mutant.objective(0, o);
                current.flags[0] = mutant.flags[0];
                archive.add(mut_genes, mut_obj);
            } else if (dom == Dominance::Equal) {
                // Non-dominated: try to add mutant to archive; if added, compare grid
                // density — accept mutant as new current only if it occupies a less
                // crowded cell than current (promotes diversity, jMetal archive.comparator).
                bool added = archive.add(mut_genes, mut_obj);
                if (added) {
                    // add() guarantees density_is_fresh() == true after return:
                    // grid_occupancy_ is accurate via O(1) incremental update
                    // (steady-state) or O(N) rebuild (bounds expansion only).
                    int curr_loc  = archive.find_grid_location_for(curr_obj);
                    int mut_loc   = archive.find_grid_location_for(mut_obj);
                    int curr_dens = archive.grid_occupancy_[curr_loc];
                    int mut_dens  = archive.grid_occupancy_[mut_loc];
                    if (curr_dens > mut_dens) {
                        for (int j = 0; j < dim; ++j) current.gene(0, j) = mutant.gene(0, j);
                        for (int o = 0; o < n_obj; ++o) current.objective(0, o) = mutant.objective(0, o);
                        current.flags[0] = mutant.flags[0];
                    }
                }
            }
            // Mutant dominated: keep current, do not add to archive
        }

        // === Copy archive back to population (result = archive, like jMetal) ===
        int result_size = archive.size();
        pop.resize(result_size);
        for (int i = 0; i < result_size; ++i) {
            const auto& ag = archive.genes(i);
            const auto& ao = archive.objectives(i);
            for (int j = 0; j < dim; ++j) pop.gene(i, j) = ag[j];
            for (int o = 0; o < n_obj; ++o) pop.objective(i, o) = ao[o];
            pop.set_evaluated(i, true);
        }
        pop.pop_size = result_size;
    }

private:
    /// Compare dominance between two individuals in potentially different populations.
    Dominance compare_dominance(this auto&, const Population<>& pop_a, int a, const Population<>& pop_b,
                                int b) {
        const int n_obj = pop_a.n_obj;
        bool a_dominates_b = false;
        bool b_dominates_a = false;

        for (int o = 0; o < n_obj; ++o) {
            double fa = pop_a.objective(a, o);
            double fb = pop_b.objective(b, o);
            if (fa < fb)
                a_dominates_b = true;
            else if (fb < fa)
                b_dominates_a = true;
            if (a_dominates_b && b_dominates_a)
                return Dominance::Equal;
        }

        if (a_dominates_b && !b_dominates_a)
            return Dominance::Dominates;
        if (b_dominates_a && !a_dominates_b)
            return Dominance::Dominated;
        return Dominance::Equal;
    }
};

} // namespace ea
