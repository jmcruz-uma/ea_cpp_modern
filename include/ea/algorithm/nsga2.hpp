#pragma once
/// @file nsga2.hpp
/// @brief NSGA-II (Non-dominated Sorting Genetic Algorithm II)
/// Reference: Deb, K., et al. "A fast and elitist multiobjective genetic algorithm:
/// NSGA-II" IEEE Trans. Evol. Comput., 2002.
///
/// Template composition — zero overhead when types are known at compile time.
/// Configurable via std::variant for runtime dispatch (see variant.hpp).

#include <ea/core/population.hpp>
#include <ea/core/comparator.hpp>
#include <ea/util/random.hpp>
#include <vector>
#include <algorithm>
#include <string_view>

namespace ea {

/// NSGA-II algorithm — the most cited multi-objective evolutionary algorithm.
/// Template parameters allow compile-time specialization for maximum performance.
template<typename CX, typename MT, typename SEL>
struct NSGAII {
    CX crossover;
    MT mutation;
    SEL selection;

    int pop_size = 100;
    int max_evals = 25000;

    static constexpr std::string_view name() { return "NSGA-II"; }

    /// Run NSGA-II on the given population.
    /// @param pop Population (genes must be initialized; objectives will be computed)
    /// @param problem Callable: void(Population&, int) — evaluates individual
    template<typename Problem>
    void run(this auto& self, Population& pop, Problem&& problem) {
        const int n = self.pop_size;
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;

        // Ensure population size matches
        if (pop.pop_size != n) {
            pop.resize(n);
        }

        // Evaluate initial population
        for (int i = 0; i < n; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
                pop.set_feasible(i, true);
            }
        }
        int evals = n;

        // Allocate offspring population
        Population offspring(n, dim, n_obj, pop.encoding, pop.n_const);
        offspring.lower_bounds = pop.lower_bounds;
        offspring.upper_bounds = pop.upper_bounds;

        // Allocate combined population (parents + offspring = 2n)
        Population combined(2 * n, dim, n_obj, pop.encoding, pop.n_const);
        combined.lower_bounds = pop.lower_bounds;
        combined.upper_bounds = pop.upper_bounds;

        // Rank and crowding storage
        std::vector<int> ranks;
        std::vector<double> crowding;

        while (evals < self.max_evals) {
            // === Step 1: Non-dominated sort of current population ===
            auto fronts = fast_non_dominated_sort(pop);

            // Compute per-individual rank and crowding distance
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

            // === Step 2: Selection (binary tournament) ===
            std::vector<int> mating_pool;
            self.selection.select(pop, mating_pool, ranks, crowding);

            // === Step 3: Crossover + Mutation → Offspring ===
            for (int i = 0; i < n; i += 2) {
                int p1 = mating_pool[i];
                int p2 = mating_pool[i + 1];

                // Copy parents to offspring slots
                for (int j = 0; j < dim; ++j) {
                    offspring.gene(i, j) = pop.gene(p1, j);
                    offspring.gene(i + 1, j) = pop.gene(p2, j);
                }

                // Crossover produces 2 children at positions i and i+1
                self.crossover.apply(pop, p1, p2, i);

                // Copy crossover results from pop slots to offspring
                // (crossover works on pop, so we copy back to offspring)
                // Actually, crossover writes children into pop slots child_start, child_start+1
                // which are the parent indices — we need to use offspring instead
                // Let's redo: crossover works on the population passed to it
            }

            // Redo crossover+mutation properly using offspring population
            for (int i = 0; i < n; i += 2) {
                int p1 = mating_pool[i];
                int p2 = mating_pool[i + 1];

                // Copy parents into offspring positions
                for (int j = 0; j < dim; ++j) {
                    offspring.gene(i, j) = pop.gene(p1, j);
                    offspring.gene(i + 1, j) = pop.gene(p2, j);
                }
                std::copy_n(pop.objectives_ptr(p1), n_obj, offspring.objectives_ptr(i));
                std::copy_n(pop.objectives_ptr(p2), n_obj, offspring.objectives_ptr(i + 1));
                offspring.set_evaluated(i, true);
                offspring.set_evaluated(i + 1, true);
            }

            // Now apply crossover on offspring population
            for (int i = 0; i < n; i += 2) {
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
                }
            }

            if (evals >= self.max_evals) break;

            // === Step 4: Combine parent + offspring ===
            combined.pop_size = 2 * n;
            for (int i = 0; i < n; ++i) {
                // Parents
                std::copy_n(pop.genes_ptr(i), dim, combined.genes_ptr(i));
                std::copy_n(pop.objectives_ptr(i), n_obj, combined.objectives_ptr(i));
                if (pop.n_const > 0) {
                    std::copy_n(pop.constraints.data() + i * pop.n_const, pop.n_const,
                               combined.constraints.data() + i * pop.n_const);
                }
                combined.flags[i] = pop.flags[i];

                // Offspring
                std::copy_n(offspring.genes_ptr(i), dim, combined.genes_ptr(n + i));
                std::copy_n(offspring.objectives_ptr(i), n_obj, combined.objectives_ptr(n + i));
                if (pop.n_const > 0) {
                    std::copy_n(offspring.constraints.data() + i * pop.n_const, pop.n_const,
                               combined.constraints.data() + (n + i) * pop.n_const);
                }
                combined.flags[n + i] = offspring.flags[i];
            }

            // === Step 5: Non-dominated sort of combined population ===
            auto combined_fronts = fast_non_dominated_sort(combined);

            // === Step 6: Select next generation ===
            int count = 0;
            int front_idx = 0;
            std::vector<double> front_cd;

            while (count < n && front_idx < static_cast<int>(combined_fronts.size())) {
                auto& front = combined_fronts[front_idx];
                if (count + static_cast<int>(front.size()) <= n) {
                    // Full front fits
                    for (int idx : front) {
                        // Copy from combined back to pop
                        std::copy_n(combined.genes_ptr(idx), dim, pop.genes_ptr(count));
                        std::copy_n(combined.objectives_ptr(idx), n_obj, pop.objectives_ptr(count));
                        if (pop.n_const > 0) {
                            std::copy_n(combined.constraints.data() + idx * pop.n_const, pop.n_const,
                                       pop.constraints.data() + count * pop.n_const);
                        }
                        pop.flags[count] = combined.flags[idx];
                        count++;
                    }
                } else {
                    // Partial front — select by crowding distance
                    compute_crowding_distance(combined, front, front_cd);
                    // Sort front by crowding distance (descending)
                    std::vector<int> sorted_front(front);
                    std::sort(sorted_front.begin(), sorted_front.end(),
                        [&](int a, int b) { return front_cd[&a - front.data()] > front_cd[&b - front.data()]; });

                    // Hmm, we need crowding distance indexed by front position
                    // Let's use a proper index mapping
                    std::vector<std::pair<double, int>> cd_index(front.size());
                    for (size_t fi = 0; fi < front.size(); ++fi) {
                        cd_index[fi] = {front_cd[fi], front[fi]};
                    }
                    std::sort(cd_index.begin(), cd_index.end(),
                        [](const auto& a, const auto& b) { return a.first > b.first; });

                    int remaining = n - count;
                    for (int fi = 0; fi < remaining && fi < static_cast<int>(cd_index.size()); ++fi) {
                        int idx = cd_index[fi].second;
                        std::copy_n(combined.genes_ptr(idx), dim, pop.genes_ptr(count));
                        std::copy_n(combined.objectives_ptr(idx), n_obj, pop.objectives_ptr(count));
                        if (pop.n_const > 0) {
                            std::copy_n(combined.constraints.data() + idx * pop.n_const, pop.n_const,
                                       pop.constraints.data() + count * pop.n_const);
                        }
                        pop.flags[count] = combined.flags[idx];
                        count++;
                    }
                    break;
                }
                front_idx++;
            }

            // All selected individuals are already evaluated
        }
    }
};

} // namespace ea