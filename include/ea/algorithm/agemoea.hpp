#pragma once
/// @file agemoea.hpp
/// @brief AGE-MOEA (Adaptive Geometry Estimation Multi-objective Evolutionary Algorithm)
/// Reference: Panichella, A. "An adaptive evolutionary algorithm based on non-Euclidean
/// geometry for many-objective optimization", GECCO 2019.
///
/// Uses geometry-based environmental selection to estimate the shape of the Pareto front
/// and compute survival scores based on distance to a local geometry estimate.

#include <algorithm>
#include <cmath>
#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <limits>
#include <numeric>
#include <string_view>
#include <vector>

namespace ea {

/// AGE-MOEA — Adaptive Geometry Estimation based MOEA.
/// Template composition for crossover and mutation.
template <typename CX, typename MT> struct AGEMOEA {
    CX crossover;
    MT mutation;

    int pop_size = 100;
    int max_evals = 25000;

    static constexpr std::string_view name() { return "AGE-MOEA"; }

    /// Run AGE-MOEA.
    template <typename Problem> void run(this auto& self, Population<>& pop, Problem&& problem) {
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

        // Environmental selection on initial population (reduce if needed)
        self.environmental_selection(pop, n);

        // Allocate offspring population
        Population<> offspring(n, dim, n_obj, pop.n_const);
        offspring.lower_bounds = pop.lower_bounds;
        offspring.upper_bounds = pop.upper_bounds;

        Population<> combined(2 * n, dim, n_obj, pop.n_const);
        combined.lower_bounds = pop.lower_bounds;
        combined.upper_bounds = pop.upper_bounds;

        auto& rng = Random::instance();

        while (evals < self.max_evals) {
            // === 1. Binary tournament selection ===
            std::vector<int> mating_pool(2 * n);
            for (int i = 0; i < 2 * n; ++i) {
                int a = rng.uniform_int(0, n - 1);
                int b = rng.uniform_int(0, n - 1);
                // For selection, use simple random (could use rank)
                mating_pool[i] = a;
            }

            // === 2. Crossover + Mutation → Offspring ===
            for (int i = 0; i < n; i += 2) {
                int p1 = mating_pool[i];
                int p2 = mating_pool[i + 1];
                for (int j = 0; j < dim; ++j) {
                    offspring.gene(i, j) = pop.gene(p1, j);
                    offspring.gene(i + 1, j) = pop.gene(p2, j);
                }
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
                    if (evals >= self.max_evals)
                        break;
                }
            }
            if (evals >= self.max_evals)
                break;

            // === 3. Combine parent + offspring ===
            for (int i = 0; i < n; ++i) {
                std::copy_n(pop.genes_ptr(i), dim, combined.genes_ptr(i));
                std::copy_n(pop.objectives_ptr(i), n_obj, combined.objectives_ptr(i));
                combined.flags[i] = pop.flags[i];
                std::copy_n(offspring.genes_ptr(i), dim, combined.genes_ptr(n + i));
                std::copy_n(offspring.objectives_ptr(i), n_obj, combined.objectives_ptr(n + i));
                combined.flags[n + i] = offspring.flags[i];
            }
            combined.pop_size = 2 * n;

            // === 4. AGE-MOEA environmental selection ===
            self.environmental_selection(combined, n);

            // Copy selected back to population
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < dim; ++j) {
                    pop.gene(i, j) = combined.gene(i, j);
                }
                for (int o = 0; o < n_obj; ++o) {
                    pop.objective(i, o) = combined.objective(i, o);
                }
                pop.flags[i] = combined.flags[i];
            }
        }
    }

private:
    /// AGE-MOEA environmental selection with geometry estimation.
    /// Sorts by non-dominated rank, then uses survival score based on
    /// geometry estimation of the Pareto front shape.
    void environmental_selection(this auto& self, Population<>& pop, int target_size) {
        const int n_obj = pop.n_obj;
        const int n = pop.pop_size;

        if (n <= target_size)
            return;

        // Step 1: Fast non-dominated sort
        auto fronts = fast_non_dominated_sort(pop);

        // Step 2: Select full fronts
        std::vector<int> selected;
        std::vector<int> critical_front; // Last front that doesn't fit

        for (auto& front : fronts) {
            if (static_cast<int>(selected.size() + front.size()) <= target_size) {
                selected.insert(selected.end(), front.begin(), front.end());
            } else {
                critical_front = front;
                break;
            }
        }

        if (critical_front.empty()) {
            // All fit — reorder population
            self.reorder_population(pop, selected);
            pop.pop_size = target_size;
            return;
        }

        // Step 3: Compute geometry-based survival scores for critical front
        // First, compute the estimated shape (p) of the Pareto front
        double p = self.estimate_geometry(pop, selected, critical_front, n_obj);

        // Compute survival scores for critical front individuals
        std::vector<double> scores(critical_front.size());
        for (size_t i = 0; i < critical_front.size(); ++i) {
            scores[i] = self.survival_score(pop, critical_front[i], selected, n_obj, p);
        }

        // Sort critical front by score (higher = better, keep individuals with highest scores)
        std::vector<std::pair<double, int>> score_index;
        score_index.reserve(critical_front.size());
        for (size_t i = 0; i < critical_front.size(); ++i) {
            score_index.push_back({scores[i], critical_front[i]});
        }
        std::sort(score_index.begin(), score_index.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });

        int remaining = target_size - static_cast<int>(selected.size());
        for (int i = 0; i < remaining && i < static_cast<int>(score_index.size()); ++i) {
            selected.push_back(score_index[i].second);
        }

        self.reorder_population(pop, selected);
        pop.pop_size = target_size;
    }

    /// Estimate the geometry (shape) parameter p of the Pareto front.
    /// p = 1: linear, p = 2: convex, p > 2: more convex, p < 1: concave.
    double estimate_geometry(this auto& self, Population<>& pop, const std::vector<int>& selected,
                             const std::vector<int>& critical_front, int n_obj) {
        // Combine selected + critical front for estimation
        std::vector<int> all_individuals = selected;
        all_individuals.insert(all_individuals.end(), critical_front.begin(), critical_front.end());

        if (all_individuals.size() < 2)
            return 1.0;

        // Compute ideal and nadir points
        std::vector<double> ideal(n_obj, std::numeric_limits<double>::max());
        std::vector<double> nadir(n_obj, -std::numeric_limits<double>::max());

        for (int idx : all_individuals) {
            for (int o = 0; o < n_obj; ++o) {
                ideal[o] = std::min(ideal[o], pop.objective(idx, o));
                nadir[o] = std::max(nadir[o], pop.objective(idx, o));
            }
        }

        // Compute average distance to the ideal point
        double total_dist = 0.0;
        for (int idx : all_individuals) {
            double dist = 0.0;
            for (int o = 0; o < n_obj; ++o) {
                double range = nadir[o] - ideal[o];
                if (range < 1e-12)
                    range = 1e-12;
                double normalized = (pop.objective(idx, o) - ideal[o]) / range;
                dist += normalized * normalized;
            }
            total_dist += std::sqrt(dist);
        }
        double avg_dist = total_dist / all_individuals.size();

        // Estimate p from average distance: closer to ideal = more convex (higher p)
        // Simple heuristic: p = 1 / avg_dist (closer = larger p)
        if (avg_dist < 1e-12)
            return 1.0;
        double p = 1.0 / avg_dist;
        return std::clamp(p, 0.1, 10.0);
    }

    /// Compute survival score for an individual based on geometry estimation.
    /// Uses the Minkowski distance with estimated p to measure contribution.
    double survival_score(this auto& self, Population<>& pop, int idx,
                          const std::vector<int>& selected, int n_obj, double p) {
        if (selected.empty()) {
            return std::numeric_limits<double>::infinity();
        }

        // Compute ideal point from selected
        std::vector<double> ideal(n_obj, std::numeric_limits<double>::max());
        for (int s : selected) {
            for (int o = 0; o < n_obj; ++o) {
                ideal[o] = std::min(ideal[o], pop.objective(s, o));
            }
        }

        // Compute Minkowski distance from ideal to individual
        double dist = 0.0;
        for (int o = 0; o < n_obj; ++o) {
            double d = pop.objective(idx, o) - ideal[o];
            if (d < 0)
                d = 0; // Shouldn't happen for ideal
            dist += std::pow(d, p);
        }
        dist = std::pow(dist, 1.0 / p);

        // Higher distance = less crowded = better survival score
        // But also prefer solutions that extend the front
        // Simple: just return distance (larger = better)
        if (dist < 1e-12) {
            // Compute crowding-like distance to nearest neighbor
            double min_neighbor_dist = std::numeric_limits<double>::infinity();
            for (int s : selected) {
                double nd = 0.0;
                for (int o = 0; o < n_obj; ++o) {
                    double d = pop.objective(idx, o) - pop.objective(s, o);
                    nd += d * d;
                }
                min_neighbor_dist = std::min(min_neighbor_dist, std::sqrt(nd));
            }
            return min_neighbor_dist;
        }

        return dist;
    }

    /// Reorder population so selected individuals are in positions [0, selected.size()).
    void reorder_population(this auto&, Population<>& pop, const std::vector<int>& selected) {
        Population<> temp(pop.pop_size, pop.dim, pop.n_obj, pop.n_const);
        temp.lower_bounds = pop.lower_bounds;
        temp.upper_bounds = pop.upper_bounds;

        // Copy selected to temp
        for (size_t i = 0; i < selected.size(); ++i) {
            int src = selected[i];
            for (int j = 0; j < pop.dim; ++j) {
                temp.gene(static_cast<int>(i), j) = pop.gene(src, j);
            }
            for (int o = 0; o < pop.n_obj; ++o) {
                temp.objective(static_cast<int>(i), o) = pop.objective(src, o);
            }
            temp.flags[i] = pop.flags[src];
        }

        // Copy back to pop
        for (size_t i = 0; i < selected.size(); ++i) {
            for (int j = 0; j < pop.dim; ++j) {
                pop.gene(static_cast<int>(i), j) = temp.gene(static_cast<int>(i), j);
            }
            for (int o = 0; o < pop.n_obj; ++o) {
                pop.objective(static_cast<int>(i), o) = temp.objective(static_cast<int>(i), o);
            }
            pop.flags[i] = temp.flags[i];
        }
    }
};

} // namespace ea
