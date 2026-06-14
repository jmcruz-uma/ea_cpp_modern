#pragma once
/// @file smsemoa.hpp
/// @brief SMS-EMOA (S-Metric Selection Evolutionary Multi-Objective Algorithm)
/// Reference: Emmerich, M., Beume, N., & Naujoks, B. "An EMO Algorithm Using the
/// Hypervolume Measure as Selection Indicator", EMO 2005.
///
/// Matches jMetal 7.4 SMSEMOAReplacement (jmetal-component):
///   - Selection   : RandomSelection (2 random parents)
///   - Crossover   : SBX on offspring copy, uses only offspring[0]
///   - Mutation    : PolynomialMutation on offspring[0]
///   - Replacement : incremental NDS + HV contribution on worst front
///   - Ref point   : max(joint) × 1.1  (REFERENCE_POINT_OFFSET_FACTOR)
///
/// Performance: incremental front insertion avoids O(N²) full NDS per step,
/// matching jMetal's O(N) amortized per-generation complexity.

#include <algorithm>
#include <cmath>
#include <ea/core/comparator.hpp>
#include <ea/core/concepts.hpp>
#include <ea/core/population.hpp>
#include <ea/indicator/hypervolume.hpp>
#include <ea/util/random.hpp>
#include <limits>
#include <vector>

namespace ea {

template <Crossover CX, Mutation MT> struct SMSEMOA {
    CX crossover;
    MT mutation;

    int pop_size  = 100;
    int max_evals = 25000;

    static constexpr std::string_view name() { return "SMS-EMOA"; }

    template <EvalFunctor F> void run(this auto& self, Population<>& pop, F&& problem) {
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

        // Working buffer: combined[0..N-1] mirrors current population,
        // combined[N] is the offspring slot. Kept persistent to avoid
        // O(N×dim) copies per generation.
        Population<> combined(self.pop_size + 1, dim, n_obj, pop.n_const);
        combined.lower_bounds = pop.lower_bounds;
        combined.upper_bounds = pop.upper_bounds;

        for (int i = 0; i < self.pop_size; ++i) {
            std::copy_n(pop.genes_ptr(i),      dim,   combined.genes_ptr(i));
            std::copy_n(pop.objectives_ptr(i), n_obj, combined.objectives_ptr(i));
            combined.flags[i] = pop.flags[i];
        }

        // Initial NDS on combined[0..N-1] — only once.
        combined.pop_size = self.pop_size;
        auto fronts = fast_non_dominated_sort(combined);
        combined.pop_size = self.pop_size + 1;

        // Pre-allocated offspring scratch buffer.
        Population<> offspring_buf(2, dim, n_obj, pop.n_const);
        offspring_buf.lower_bounds = combined.lower_bounds;
        offspring_buf.upper_bounds = combined.upper_bounds;

        std::vector<double> ref_point(n_obj, 0.0);
        auto& rng = Random::instance();

        while (evals < self.max_evals) {
            // ── Selection: 2 random parents ──────────────────────────────────
            int p1 = rng.uniform_int(0, self.pop_size - 1);
            int p2 = rng.uniform_int(0, self.pop_size - 1);

            // ── Offspring: copy parents, crossover, mutation ──────────────────
            for (int j = 0; j < dim; ++j) {
                offspring_buf.gene(0, j) = combined.gene(p1, j);
                offspring_buf.gene(1, j) = combined.gene(p2, j);
            }
            self.crossover.apply(offspring_buf, 0, 1, 0);
            self.mutation.apply(offspring_buf, 0);

            // Write offspring to the extra slot combined[pop_size]
            std::copy_n(offspring_buf.genes_ptr(0),      dim,   combined.genes_ptr(self.pop_size));
            std::copy_n(offspring_buf.objectives_ptr(0), n_obj, combined.objectives_ptr(self.pop_size));
            combined.flags[self.pop_size] = offspring_buf.flags[0];

            problem(combined, self.pop_size);
            combined.set_evaluated(self.pop_size, true);
            evals++;

            // ── Reference point: max over N+1 joint individuals × 1.1 ─────────
            for (int o = 0; o < n_obj; ++o) {
                double max_o = -std::numeric_limits<double>::max();
                for (int i = 0; i <= self.pop_size; ++i)
                    max_o = std::max(max_o, combined.objective(i, o));
                ref_point[o] = max_o * 1.1;
            }

            // ── Incremental NDS: insert combined[pop_size] into cached fronts ─
            // O(N) amortized — matches SMSEMOAReplacement.insertOffspringIntoFronts
            insert_incremental(combined, fronts, self.pop_size);

            // ── Find worst in last front ──────────────────────────────────────
            auto& last_front = fronts.back();
            int worst;
            if (last_front.size() == 1) {
                worst = last_front[0];
            } else if (n_obj == 2) {
                worst = least_contributor_2d(combined, last_front, ref_point);
            } else {
                worst = least_contributor_wfg(combined, last_front, ref_point);
            }

            // ── Remove worst, keep offspring: update fronts + swap data ───────
            // Step 1: remove `worst` from fronts
            remove_from_fronts(fronts, worst);
            // Step 2: if offspring (pop_size) is not the worst, move it to slot `worst`
            if (worst != self.pop_size) {
                rename_in_fronts(fronts, self.pop_size, worst);
                for (int j = 0; j < dim; ++j)
                    std::swap(combined.gene(worst, j), combined.gene(self.pop_size, j));
                for (int o = 0; o < n_obj; ++o)
                    std::swap(combined.objective(worst, o), combined.objective(self.pop_size, o));
                std::swap(combined.flags[worst], combined.flags[self.pop_size]);
            }
            while (!fronts.empty() && fronts.back().empty())
                fronts.pop_back();
            // combined[0..pop_size-1] is now the new population; fronts are up-to-date
        }

        // Copy combined[0..N-1] back to pop
        for (int i = 0; i < self.pop_size; ++i) {
            std::copy_n(combined.genes_ptr(i),      dim,   pop.genes_ptr(i));
            std::copy_n(combined.objectives_ptr(i), n_obj, pop.objectives_ptr(i));
            pop.flags[i] = combined.flags[i];
        }
    }

private:
    // ── Incremental front insertion ───────────────────────────────────────────
    // Inserts new_idx into the first front it is not dominated by.
    // Displaced members cascade to subsequent fronts (same logic as
    // SMSEMOAReplacement.insertOffspringIntoFronts + propagateDisplacedSolutions).
    static void insert_incremental(Population<>& pop,
                                   std::vector<std::vector<int>>& fronts,
                                   int new_idx) {
        for (int fi = 0; fi < (int)fronts.size(); ++fi) {
            // Is new_idx dominated by anyone in this front?
            bool dominated = false;
            for (int idx : fronts[fi]) {
                if (dominates(pop, idx, new_idx)) { dominated = true; break; }
            }
            if (dominated) continue;

            // Partition: survivors stay in fi, dominated by new_idx are displaced
            auto& front = fronts[fi];
            int split = (int)(std::partition(front.begin(), front.end(), [&](int idx) {
                return !dominates(pop, new_idx, idx);
            }) - front.begin());

            std::vector<int> displaced(front.begin() + split, front.end());
            front.resize(split);
            front.push_back(new_idx);

            if (!displaced.empty())
                propagate_displaced(pop, fronts, fi + 1, std::move(displaced));
            return;
        }
        fronts.push_back({new_idx});
    }

    // Cascades displaced solutions from front start_fi onward.
    // Matches SMSEMOAReplacement.propagateDisplacedSolutions.
    static void propagate_displaced(Population<>& pop,
                                    std::vector<std::vector<int>>& fronts,
                                    int start_fi,
                                    std::vector<int> pending) {
        int fi = start_fi;
        while (!pending.empty()) {
            if (fi >= (int)fronts.size()) {
                fronts.push_back(std::move(pending));
                return;
            }
            auto& front = fronts[fi];
            // Remove from front anything dominated by a pending solution
            int split = (int)(std::partition(front.begin(), front.end(), [&](int idx) {
                for (int p : pending) if (dominates(pop, p, idx)) return false;
                return true;
            }) - front.begin());
            std::vector<int> next_pending(front.begin() + split, front.end());
            front.resize(split);
            // Prepend pending to front (they are mutually non-dominated and
            // non-dominated w.r.t. surviving front members)
            front.insert(front.end(), pending.begin(), pending.end());
            pending = std::move(next_pending);
            ++fi;
        }
    }

    // a strictly dominates b (minimization)
    static bool dominates(const Population<>& pop, int a, int b) {
        bool any_better = false;
        for (int o = 0; o < pop.n_obj; ++o) {
            if (pop.objective(a, o) > pop.objective(b, o)) return false;
            if (pop.objective(a, o) < pop.objective(b, o)) any_better = true;
        }
        return any_better;
    }

    static void remove_from_fronts(std::vector<std::vector<int>>& fronts, int idx) {
        for (auto& front : fronts) {
            auto it = std::find(front.begin(), front.end(), idx);
            if (it != front.end()) { front.erase(it); return; }
        }
    }

    static void rename_in_fronts(std::vector<std::vector<int>>& fronts,
                                  int old_idx, int new_idx) {
        for (auto& front : fronts)
            for (int& i : front)
                if (i == old_idx) { i = new_idx; return; }
    }

    // ── HV contribution ───────────────────────────────────────────────────────

    // O(n log n) exact 2D exclusive HV contribution for a non-dominated front.
    static int least_contributor_2d(const Population<>& pop,
                                    const std::vector<int>& front,
                                    const std::vector<double>& ref_point) {
        std::vector<int> sorted(front);
        std::sort(sorted.begin(), sorted.end(),
                  [&](int a, int b) { return pop.objective(a, 0) < pop.objective(b, 0); });

        const int n = (int)sorted.size();
        double min_contrib = std::numeric_limits<double>::max();
        int worst = sorted[0];

        for (int i = 0; i < n; ++i) {
            double f1 = pop.objective(sorted[i], 0);
            double f2 = pop.objective(sorted[i], 1);
            double contrib;
            if (i == 0) {
                contrib = (pop.objective(sorted[1], 0) - f1) * (ref_point[1] - f2);
            } else if (i == n - 1) {
                contrib = (ref_point[0] - f1) * (pop.objective(sorted[i-1], 1) - f2);
            } else {
                contrib = (pop.objective(sorted[i+1], 0) - f1)
                        * (pop.objective(sorted[i-1], 1) - f2);
            }
            if (contrib < min_contrib) { min_contrib = contrib; worst = sorted[i]; }
        }
        return worst;
    }

    // O(n × WFG) HV contribution for n_obj ≥ 3.
    static int least_contributor_wfg(const Population<>& pop,
                                     const std::vector<int>& front,
                                     const std::vector<double>& ref_point) {
        const int n = (int)front.size();
        double total_hv = hypervolume(pop, front, ref_point);
        double min_contrib = std::numeric_limits<double>::max();
        int worst = front[0];

        std::vector<int> reduced(n - 1);
        for (int k = 0; k < n; ++k) {
            int dst = 0;
            for (int j = 0; j < n; ++j)
                if (j != k) reduced[dst++] = front[j];
            double contrib = total_hv - (reduced.empty() ? 0.0 : hypervolume(pop, reduced, ref_point));
            if (contrib < min_contrib) { min_contrib = contrib; worst = front[k]; }
        }
        return worst;
    }
};

} // namespace ea
