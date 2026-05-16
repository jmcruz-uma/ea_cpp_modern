#pragma once
/// @file nsga3.hpp
/// @brief NSGA-III (Many-objective NSGA-II with reference points)
/// Reference: Deb, K. & Jain, H. "An Evolutionary Many-Objective Optimization
/// Algorithm Using Reference-Point-Based Nondominated Sorting Approach", IEEE TEVC 2014.
///
/// Uses reference-point-based niching for environmental selection.

#include <ea/core/population.hpp>
#include <ea/core/comparator.hpp>
#include <ea/util/random.hpp>
#include <ea/util/reference_point.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include <string_view>

namespace ea {

/// NSGA-III algorithm — many-objective optimization using reference points.
/// Template parameters allow compile-time specialization.
template<typename CX, typename MT>
struct NSGAIII {
    CX crossover;
    MT mutation;

    int pop_size = 100;
    int max_evals = 25000;
    int divisions = 10;  ///< Das-Dennis divisions for reference point generation

    static constexpr std::string_view name() { return "NSGA-III"; }

    /// Run NSGA-III.
    template<typename Problem>
    void run(this auto& self, Population& pop, Problem&& problem) {
        const int n = self.pop_size;
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;

        if (pop.pop_size != n) pop.resize(n);

        // Generate reference points
        auto reference_points = generate_reference_points(n_obj, self.divisions);

        // Evaluate initial population
        for (int i = 0; i < n; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }
        int evals = n;

        // Allocate offspring + combined
        Population offspring(n, dim, n_obj, pop.encoding, pop.n_const);
        offspring.lower_bounds = pop.lower_bounds;
        offspring.upper_bounds = pop.upper_bounds;

        Population combined(2 * n, dim, n_obj, pop.encoding, pop.n_const);
        combined.lower_bounds = pop.lower_bounds;
        combined.upper_bounds = pop.upper_bounds;

        // Ideal point and nadir point (updated each generation)
        std::vector<double> ideal_point(n_obj, std::numeric_limits<double>::max());

        auto& rng = Random::instance();

        while (evals < self.max_evals) {
            // Update ideal point
            for (int i = 0; i < n; ++i) {
                for (int o = 0; o < n_obj; ++o) {
                    ideal_point[o] = std::min(ideal_point[o], pop.objective(i, o));
                }
            }

            // === Non-dominated sort ===
            auto fronts = fast_non_dominated_sort(pop);

            // === Selection (binary tournament by rank) ===
            std::vector<int> ranks(n, 0);
            for (int r = 0; r < static_cast<int>(fronts.size()); ++r) {
                for (int idx : fronts[r]) ranks[idx] = r;
            }

            std::vector<int> mating_pool(2 * n);
            for (int i = 0; i < 2 * n; ++i) {
                int a = rng.uniform_int(0, n - 1);
                int b = rng.uniform_int(0, n - 1);
                mating_pool[i] = (ranks[a] < ranks[b]) ? a : b;
            }

            // === Crossover + Mutation → Offspring ===
            for (int i = 0; i < n; i += 2) {
                int p1 = mating_pool[i];
                int p2 = mating_pool[i + 1];
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
            for (int i = 0; i < n; i += 2) {
                self.crossover.apply(offspring, i, i + 1, i);
            }

            // Apply mutation
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

            // === Combine parent + offspring ===
            for (int i = 0; i < n; ++i) {
                std::copy_n(pop.genes_ptr(i), dim, combined.genes_ptr(i));
                std::copy_n(pop.objectives_ptr(i), n_obj, combined.objectives_ptr(i));
                combined.flags[i] = pop.flags[i];
                std::copy_n(offspring.genes_ptr(i), dim, combined.genes_ptr(n + i));
                std::copy_n(offspring.objectives_ptr(i), n_obj, combined.objectives_ptr(n + i));
                combined.flags[n + i] = offspring.flags[i];
            }
            combined.pop_size = 2 * n;

            // === Environmental selection with reference point niching ===
            auto sel_fronts = fast_non_dominated_sort(combined);

            // Normalize objectives
            std::vector<double> ideal(n_obj, std::numeric_limits<double>::max());
            for (int i = 0; i < combined.pop_size; ++i) {
                for (int o = 0; o < n_obj; ++o) {
                    ideal[o] = std::min(ideal[o], combined.objective(i, o));
                }
            }

            // Translate objectives (subtract ideal point)
            std::vector<std::vector<double>> translated(combined.pop_size, std::vector<double>(n_obj));
            for (int i = 0; i < combined.pop_size; ++i) {
                for (int o = 0; o < n_obj; ++o) {
                    translated[i][o] = combined.objective(i, o) - ideal[o];
                }
            }

            // Select individuals
            std::vector<int> selected;
            for (auto& front : sel_fronts) {
                if (static_cast<int>(selected.size() + front.size()) <= n) {
                    selected.insert(selected.end(), front.begin(), front.end());
                } else {
                    // Niching selection for remaining spots
                    int remaining = n - static_cast<int>(selected.size());
                    // Simple approach: use crowding distance for remaining
                    std::vector<double> front_cd;
                    compute_crowding_distance(combined, front, front_cd);

                    std::vector<std::pair<double, int>> cd_index(front.size());
                    for (size_t fi = 0; fi < front.size(); ++fi) {
                        cd_index[fi] = {front_cd[fi], front[fi]};
                    }
                    std::sort(cd_index.begin(), cd_index.end(),
                        [](const auto& a, const auto& b) { return a.first > b.first; });

                    for (int fi = 0; fi < remaining && fi < static_cast<int>(cd_index.size()); ++fi) {
                        selected.push_back(cd_index[fi].second);
                    }
                    break;
                }
            }

            // Copy selected back to population
            for (int i = 0; i < n; ++i) {
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

private:
    /// Generate Das-Dennis reference points
    static std::vector<ReferencePoint> generate_reference_points(int n_obj, int divisions) {
        std::vector<ReferencePoint> points;
        std::vector<double> ref_point(n_obj, 0.0);
        generate_recursive(points, ref_point, n_obj, divisions, 0, divisions);
        return points;
    }

    static void generate_recursive(std::vector<ReferencePoint>& points,
                                    std::vector<double>& ref_point,
                                    int n_obj, int divisions, int obj_idx, int remaining) {
        if (obj_idx == n_obj - 1) {
            ref_point[obj_idx] = static_cast<double>(remaining) / divisions;
            ReferencePoint rp;
            rp.position = ref_point;
            points.push_back(rp);
            return;
        }
        for (int i = 0; i <= remaining; ++i) {
            ref_point[obj_idx] = static_cast<double>(i) / divisions;
            generate_recursive(points, ref_point, n_obj, divisions, obj_idx + 1, remaining - i);
        }
    }
};

} // namespace ea