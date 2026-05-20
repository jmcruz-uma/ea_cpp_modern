#pragma once
/// @file smsemoa.hpp
/// @brief SMS-EMOA (S-Metric Selection Evolutionary Multi-Objective Algorithm)
/// Reference: Emmerich, M., Beume, N., & Naujoks, B. "An EMO Algorithm Using the
/// Hypervolume Measure as Selection Indicator", EMO 2005.
///
/// Uses hypervolume contribution for environmental selection.

#include <algorithm>
#include <cmath>
#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <ea/indicator/hypervolume.hpp>
#include <ea/util/random.hpp>
#include <limits>
#include <vector>

namespace ea {

/// SMS-EMOA — Hypervolume-based selection MOEA.
/// Each generation: create offspring, merge with population, remove individual
/// with smallest hypervolume contribution.
template <typename CX, typename MT> struct SMSEMOA {
    CX crossover;
    MT mutation;

    int pop_size = 100;
    int max_evals = 25000;

    static constexpr std::string_view name() { return "SMS-EMOA"; }

    /// Run SMS-EMOA.
    template <typename Problem> void run(this auto& self, Population& pop, Problem&& problem) {
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;

        if (pop.pop_size != self.pop_size)
            pop.resize(self.pop_size);

        // Evaluate initial population
        for (int i = 0; i < self.pop_size; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }
        int evals = self.pop_size;

        // Reference point for hypervolume (worst case + margin)
        std::vector<double> ref_point(n_obj);
        compute_reference_point(pop, ref_point, 0.1);

        // Combined population (pop + offspring)
        Population combined(self.pop_size + 1, dim, n_obj, pop.encoding, pop.n_const);
        combined.lower_bounds = pop.lower_bounds;
        combined.upper_bounds = pop.upper_bounds;

        auto& rng = Random::instance();

        while (evals < self.max_evals) {
            // Selection — binary tournament
            int p1 = rng.uniform_int(0, self.pop_size - 1);
            int p2 = rng.uniform_int(0, self.pop_size - 1);

            // Create offspring by crossover + mutation
            Population offspring(2, dim, n_obj, pop.encoding, pop.n_const);
            offspring.lower_bounds = pop.lower_bounds;
            offspring.upper_bounds = pop.upper_bounds;

            // Copy parents to offspring
            for (int j = 0; j < dim; ++j) {
                offspring.gene(0, j) = pop.gene(p1, j);
                offspring.gene(1, j) = pop.gene(p2, j);
            }

            // Crossover
            self.crossover.apply(pop, p1, p2, 0);

            // Copy crossover result to a temp individual
            Population child(1, dim, n_obj, pop.encoding, pop.n_const);
            child.lower_bounds = pop.lower_bounds;
            child.upper_bounds = pop.upper_bounds;
            for (int j = 0; j < dim; ++j) {
                child.gene(0, j) = pop.gene(p1, j); // Use first child
            }

            // Mutation
            self.mutation.apply(child, 0);

            // Evaluate offspring
            problem(child, 0);
            child.set_evaluated(0, true);
            evals++;
            if (evals >= self.max_evals)
                break;

            // Update reference point if needed
            for (int o = 0; o < n_obj; ++o) {
                if (child.objective(0, o) > ref_point[o] - 0.1) {
                    ref_point[o] = child.objective(0, o) + 0.1;
                    // Update ref for all existing too
                    for (int i = 0; i < self.pop_size; ++i) {
                        ref_point[o] = std::max(ref_point[o], pop.objective(i, o) + 0.1);
                    }
                }
            }

            // Combine population + offspring
            combined.pop_size = self.pop_size + 1;
            for (int i = 0; i < self.pop_size; ++i) {
                std::copy_n(pop.genes_ptr(i), dim, combined.genes_ptr(i));
                std::copy_n(pop.objectives_ptr(i), n_obj, combined.objectives_ptr(i));
                combined.flags[i] = pop.flags[i];
            }
            std::copy_n(child.genes_ptr(0), dim, combined.genes_ptr(self.pop_size));
            std::copy_n(child.objectives_ptr(0), n_obj, combined.objectives_ptr(self.pop_size));
            combined.flags[self.pop_size] = child.flags[0];

            // Remove individual with smallest hypervolume contribution
            int worst = find_worst_by_hv_contribution(combined, ref_point);

            // Replace worst with offspring (or keep if worst is the offspring)
            if (worst != self.pop_size) {
                // Remove worst by replacing with offspring position
                if (worst < self.pop_size) {
                    // Remove worst: shift everything after worst one position left
                    for (int i = worst; i < self.pop_size - 1; ++i) {
                        std::copy_n(combined.genes_ptr(i + 1), dim, combined.genes_ptr(i));
                        std::copy_n(combined.objectives_ptr(i + 1), n_obj,
                                    combined.objectives_ptr(i));
                        combined.flags[i] = combined.flags[i + 1];
                    }
                    // Add offspring at the end
                    std::copy_n(child.genes_ptr(0), dim, combined.genes_ptr(self.pop_size - 1));
                    std::copy_n(child.objectives_ptr(0), n_obj,
                                combined.objectives_ptr(self.pop_size - 1));
                    combined.flags[self.pop_size - 1] = child.flags[0];
                }
            }
            // If worst is the offspring, we keep the current population as-is

            // Copy back from combined (first pop_size individuals)
            for (int i = 0; i < self.pop_size; ++i) {
                std::copy_n(combined.genes_ptr(i), dim, pop.genes_ptr(i));
                std::copy_n(combined.objectives_ptr(i), n_obj, pop.objectives_ptr(i));
                pop.flags[i] = combined.flags[i];
            }
        }
    }

private:
    static void compute_reference_point(const Population& pop, std::vector<double>& ref_point,
                                        double margin) {
        const int n_obj = pop.n_obj;
        for (int o = 0; o < n_obj; ++o) {
            double max_obj = -std::numeric_limits<double>::max();
            for (int i = 0; i < pop.pop_size; ++i) {
                max_obj = std::max(max_obj, pop.objective(i, o));
            }
            ref_point[o] = max_obj + margin;
        }
    }

    static int find_worst_by_hv_contribution(Population& combined,
                                             const std::vector<double>& ref_point) {
        const int n = combined.pop_size;
        const int n_obj = combined.n_obj;

        // First do non-dominated sort; worst from worst front
        auto fronts = fast_non_dominated_sort(combined);

        // If last front has only one individual, remove it
        if (!fronts.empty() && fronts.back().size() == 1) {
            return fronts.back()[0];
        }

        // Otherwise, from the worst front, remove the one with smallest HV contribution
        auto& last_front = fronts.back();
        if (last_front.size() == 1) {
            return last_front[0];
        }

        // Compute HV contribution for each individual in last front
        // HV contribution of i = HV(all) - HV(all \ {i})
        // Simplified: for 2D, compute contribution directly
        if (n_obj == 2) {
            // Sort by first objective
            std::vector<int> sorted_front(last_front);
            std::sort(sorted_front.begin(), sorted_front.end(), [&](int a, int b) {
                return combined.objective(a, 0) < combined.objective(b, 0);
            });

            // Compute contribution: area between this point and neighbors
            double min_contribution = std::numeric_limits<double>::max();
            int worst_idx = sorted_front[0];

            for (size_t i = 0; i < sorted_front.size(); ++i) {
                double contribution;
                if (sorted_front.size() == 1) {
                    contribution = std::numeric_limits<double>::max();
                } else if (i == 0) {
                    // Leftmost point: contribution = (f1_right - f1_this) * (ref[1] - f2_this)
                    double f1_right = combined.objective(sorted_front[i + 1], 0);
                    contribution = (f1_right - combined.objective(sorted_front[i], 0)) *
                                   (ref_point[1] - combined.objective(sorted_front[i], 1));
                } else if (i == sorted_front.size() - 1) {
                    // Rightmost point: contribution = (ref[0] - f1_this) * (f2_left - f2_this)
                    double f2_left = combined.objective(sorted_front[i - 1], 1);
                    contribution = (ref_point[0] - combined.objective(sorted_front[i], 0)) *
                                   (f2_left - combined.objective(sorted_front[i], 1));
                } else {
                    // Middle point
                    double f1_right = combined.objective(sorted_front[i + 1], 0);
                    double f2_left = combined.objective(sorted_front[i - 1], 1);
                    contribution = (f1_right - combined.objective(sorted_front[i], 0)) *
                                   (f2_left - combined.objective(sorted_front[i], 1));
                }

                if (contribution < min_contribution) {
                    min_contribution = contribution;
                    worst_idx = sorted_front[i];
                }
            }

            return worst_idx;
        }

        // Higher dimensions: use general HV computation
        // Compute full HV, then subtract HV without each individual
        // This is expensive — simplified: just remove from worst front
        return last_front[0];
    }
};

} // namespace ea