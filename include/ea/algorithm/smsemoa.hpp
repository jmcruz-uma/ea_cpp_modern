#pragma once
/// @file smsemoa.hpp
/// @brief SMS-EMOA (S-Metric Selection Evolutionary Multi-Objective Algorithm)
/// Reference: Emmerich, M., Beume, N., & Naujoks, B. "An EMO Algorithm Using the
/// Hypervolume Measure as Selection Indicator", EMO 2005.
///
/// Matches jMetal 7.4 SMSEMOA.java + SMSEMOABuilder.java:
///   - Selection : RandomSelection (2 random parents, with replacement)
///   - Crossover : SBX applied to offspring copy, uses only offspring[0]
///   - Mutation  : PolynomialMutation on offspring[0]
///   - Replacement: NDS + HV contribution on worst front; offset = 100.0

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
/// Each generation: create one offspring, merge with population (N+1),
/// run NDS, remove individual with smallest HV contribution in worst front.
template <typename CX, typename MT> struct SMSEMOA {
    CX crossover;
    MT mutation;

    int pop_size  = 100;
    int max_evals = 25000;

    static constexpr std::string_view name() { return "SMS-EMOA"; }

    template <typename Problem> void run(this auto& self, Population<>& pop, Problem&& problem) {
        const int dim   = pop.dim;
        const int n_obj = pop.n_obj;

        if (pop.pop_size != self.pop_size)
            pop.resize(self.pop_size);

        for (int i = 0; i < self.pop_size; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }
        int evals = self.pop_size;

        // Combined buffer (pop + 1 offspring); reused every generation
        Population<> combined(self.pop_size + 1, dim, n_obj, pop.n_const);
        combined.lower_bounds = pop.lower_bounds;
        combined.upper_bounds = pop.upper_bounds;

        std::vector<double> ref_point(n_obj, 0.0);
        auto& rng = Random::instance();

        while (evals < self.max_evals) {
            // ── Selection: 2 random parents (jMetal RandomSelection) ──────────
            int p1 = rng.uniform_int(0, self.pop_size - 1);
            int p2 = rng.uniform_int(0, self.pop_size - 1);

            // ── Offspring: copy parents, crossover, mutation ──────────────────
            // SBX requires parent genes already in child_start / child_start+1
            Population<> offspring(2, dim, n_obj, pop.n_const);
            offspring.lower_bounds = pop.lower_bounds;
            offspring.upper_bounds = pop.upper_bounds;
            for (int j = 0; j < dim; ++j) {
                offspring.gene(0, j) = pop.gene(p1, j);
                offspring.gene(1, j) = pop.gene(p2, j);
            }
            self.crossover.apply(offspring, 0, 1, 0);  // reads offspring[0,1], writes offspring[0,1]
            self.mutation.apply(offspring, 0);

            problem(offspring, 0);
            offspring.set_evaluated(0, true);
            evals++;

            // ── Build combined (N+1) ──────────────────────────────────────────
            combined.pop_size = self.pop_size + 1;
            for (int i = 0; i < self.pop_size; ++i) {
                std::copy_n(pop.genes_ptr(i),       dim,   combined.genes_ptr(i));
                std::copy_n(pop.objectives_ptr(i),  n_obj, combined.objectives_ptr(i));
                combined.flags[i] = pop.flags[i];
            }
            std::copy_n(offspring.genes_ptr(0),      dim,   combined.genes_ptr(self.pop_size));
            std::copy_n(offspring.objectives_ptr(0), n_obj, combined.objectives_ptr(self.pop_size));
            combined.flags[self.pop_size] = offspring.flags[0];

            // ── Reference point: max over joint population + offset (jMetal DEFAULT_OFFSET=100) ─
            for (int o = 0; o < n_obj; ++o) {
                double max_o = -std::numeric_limits<double>::max();
                for (int i = 0; i <= self.pop_size; ++i)
                    max_o = std::max(max_o, combined.objective(i, o));
                ref_point[o] = max_o + 100.0;
            }

            // ── Remove individual with smallest HV contribution in worst front ─
            int worst = find_worst_by_hv_contribution(combined, ref_point);

            // ── Copy combined without worst back to pop ────────────────────────
            int dst = 0;
            for (int i = 0; i <= self.pop_size; ++i) {
                if (i == worst) continue;
                std::copy_n(combined.genes_ptr(i),      dim,   pop.genes_ptr(dst));
                std::copy_n(combined.objectives_ptr(i), n_obj, pop.objectives_ptr(dst));
                pop.flags[dst] = combined.flags[i];
                ++dst;
            }
        }
    }

private:
    static int find_worst_by_hv_contribution(Population<>& combined,
                                             const std::vector<double>& ref_point) {
        const int n_obj = combined.n_obj;
        auto fronts = fast_non_dominated_sort(combined);

        auto& last_front = fronts.back();

        // Single individual in worst front — remove it directly
        if (last_front.size() == 1)
            return last_front[0];

        if (n_obj == 2)
            return least_contributor_2d(combined, last_front, ref_point);

        return least_contributor_wfg(combined, last_front, ref_point);
    }

    // O(n) exact 2D exclusive HV contribution for a non-dominated front.
    // Formula: contribution_i = width_i × height_i where
    //   width_i  = f1_{i+1} - f1_i  (or ref[0] - f1_i for rightmost)
    //   height_i = f2_{i-1} - f2_i  (or ref[1] - f2_i for leftmost)
    // sorted by f1 ascending.
    static int least_contributor_2d(const Population<>& pop,
                                    const std::vector<int>& front,
                                    const std::vector<double>& ref_point) {
        std::vector<int> sorted(front);
        std::sort(sorted.begin(), sorted.end(),
                  [&](int a, int b) { return pop.objective(a, 0) < pop.objective(b, 0); });

        const int n = static_cast<int>(sorted.size());
        double min_contrib = std::numeric_limits<double>::max();
        int worst = sorted[0];

        for (int i = 0; i < n; ++i) {
            double f1_this = pop.objective(sorted[i], 0);
            double f2_this = pop.objective(sorted[i], 1);
            double contrib;
            if (i == 0) {
                double f1_next = pop.objective(sorted[i + 1], 0);
                contrib = (f1_next - f1_this) * (ref_point[1] - f2_this);
            } else if (i == n - 1) {
                double f2_prev = pop.objective(sorted[i - 1], 1);
                contrib = (ref_point[0] - f1_this) * (f2_prev - f2_this);
            } else {
                double f1_next = pop.objective(sorted[i + 1], 0);
                double f2_prev = pop.objective(sorted[i - 1], 1);
                contrib = (f1_next - f1_this) * (f2_prev - f2_this);
            }
            if (contrib < min_contrib) {
                min_contrib = contrib;
                worst = sorted[i];
            }
        }
        return worst;
    }

    // O(n × WFG) HV contribution for n_obj ≥ 3.
    // contribution_i = HV(front) - HV(front \ {i}); remove the minimum.
    static int least_contributor_wfg(const Population<>& pop,
                                     const std::vector<int>& front,
                                     const std::vector<double>& ref_point) {
        const int n = static_cast<int>(front.size());
        double total_hv = hypervolume(pop, front, ref_point);

        double min_contrib = std::numeric_limits<double>::max();
        int worst = front[0];

        std::vector<int> reduced(n - 1);
        for (int k = 0; k < n; ++k) {
            int dst = 0;
            for (int j = 0; j < n; ++j)
                if (j != k) reduced[dst++] = front[j];
            double hv_without = reduced.empty() ? 0.0 : hypervolume(pop, reduced, ref_point);
            double contrib = total_hv - hv_without;
            if (contrib < min_contrib) {
                min_contrib = contrib;
                worst = front[k];
            }
        }
        return worst;
    }
};

} // namespace ea
