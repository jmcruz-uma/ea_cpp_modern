#pragma once
/// @file agemoea.hpp
/// @brief AGE-MOEA (Adaptive Geometry Estimation Multi-objective Evolutionary Algorithm)
/// Reference: Panichella, A. "An adaptive evolutionary algorithm based on non-Euclidean
/// geometry for many-objective optimization", GECCO 2019.
///
/// Faithfully mirrors jMetal 7.4 AGEMOEAEnvironmentalSelection:
///   - Survival score for front[0]: extreme points (+inf), then greedy diversity
///     using pairwise Minkowski distances normalised by proximity to utopia.
///   - Convergence score for fronts rank > 0: Minkowski distance to utopia / intercepts.
///   - Binary tournament: rank first, survival score second (higher = better).

#include <algorithm>
#include <cassert>
#include <cmath>
#include <ea/core/comparator.hpp>
#include <ea/core/concepts.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <limits>
#include <numeric>
#include <string_view>
#include <vector>

namespace ea {

template <Crossover CX, Mutation MT> struct AGEMOEA {
    CX crossover;
    MT mutation;

    int pop_size  = 100;
    int max_evals = 25000;

    static constexpr std::string_view name() { return "AGE-MOEA"; }

    template <EvalFunctor F> void run(this auto& self, Population<>& pop, F&& problem) {
        const int n     = self.pop_size;
        const int dim   = pop.dim;
        const int n_obj = pop.n_obj;

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

        // Assign initial survival scores (used by first generation's tournament)
        self.scores_.assign(n, 0.0);
        self.ranks_.assign(n, 0);
        self.environmental_selection(pop, n); // scores/ranks computed, no pruning

        Population<> offspring(n, dim, n_obj, pop.n_const);
        offspring.lower_bounds = pop.lower_bounds;
        offspring.upper_bounds = pop.upper_bounds;

        Population<> combined(2 * n, dim, n_obj, pop.n_const);
        combined.lower_bounds = pop.lower_bounds;
        combined.upper_bounds = pop.upper_bounds;

        auto& rng = Random::instance();

        while (evals < self.max_evals) {
            // 1. Binary tournament selection (rank first, survival score second, higher = better)
            std::vector<int> mating_pool(2 * n);
            for (int i = 0; i < 2 * n; ++i) {
                int a = rng.uniform_int(0, n - 1);
                int b = rng.uniform_int(0, n - 1);
                if (self.ranks_[a] < self.ranks_[b]) {
                    mating_pool[i] = a;
                } else if (self.ranks_[a] > self.ranks_[b]) {
                    mating_pool[i] = b;
                } else if (self.scores_[a] >= self.scores_[b]) {
                    mating_pool[i] = a;
                } else {
                    mating_pool[i] = b;
                }
            }

            // 2. Crossover + Mutation → Offspring
            for (int i = 0; i < n; i += 2) {
                int p1 = mating_pool[i], p2 = mating_pool[i + 1];
                for (int j = 0; j < dim; ++j) {
                    offspring.gene(i,     j) = pop.gene(p1, j);
                    offspring.gene(i + 1, j) = pop.gene(p2, j);
                }
                offspring.set_evaluated(i,     false);
                offspring.set_evaluated(i + 1, false);
            }
            for (int i = 0; i < n; i += 2)
                self.crossover.apply(offspring, i, i + 1, i);
            for (int i = 0; i < n; ++i)
                self.mutation.apply(offspring, i);

            // Evaluate offspring
            for (int i = 0; i < n; ++i) {
                problem(offspring, i);
                offspring.set_evaluated(i, true);
                ++evals;
                if (evals >= self.max_evals) break;
            }
            if (evals >= self.max_evals) break;

            // 3. Combine parent + offspring → combined[2n]
            for (int i = 0; i < n; ++i) {
                std::copy_n(pop.genes_ptr(i),      dim,   combined.genes_ptr(i));
                std::copy_n(pop.objectives_ptr(i),  n_obj, combined.objectives_ptr(i));
                combined.flags[i] = pop.flags[i];
                std::copy_n(offspring.genes_ptr(i),     dim,   combined.genes_ptr(n + i));
                std::copy_n(offspring.objectives_ptr(i), n_obj, combined.objectives_ptr(n + i));
                combined.flags[n + i] = offspring.flags[i];
            }
            combined.pop_size = 2 * n;

            // 4. AGE-MOEA environmental selection → combined shrinks to n
            self.environmental_selection(combined, n);

            // Copy selected back to population
            for (int i = 0; i < n; ++i) {
                std::copy_n(combined.genes_ptr(i),      dim,   pop.genes_ptr(i));
                std::copy_n(combined.objectives_ptr(i),  n_obj, pop.objectives_ptr(i));
                pop.flags[i] = combined.flags[i];
            }
        }
    }

private:
    // State carried across iterations (computed in environmental_selection, used in tournament)
    double              P_         = 1.0; // geometry exponent estimated from front[0]
    std::vector<double> intercepts_;      // hyperplane intercepts from front[0]
    std::vector<double> scores_;          // survival/convergence score per individual [0..pop_size)
    std::vector<int>    ranks_;           // NDS rank per individual [0..pop_size)

    // =========================================================================
    // Geometry utilities
    // =========================================================================

    /// Perpendicular distance from point P to line defined by two points A and B.
    static double point2line_dist(const double* P, const std::vector<double>& A,
                                  const std::vector<double>& B, int m) {
        std::vector<double> pa(m), ba(m);
        for (int i = 0; i < m; ++i) { pa[i] = P[i] - A[i]; ba[i] = B[i] - A[i]; }
        double t_num = 0.0, t_den = 0.0;
        for (int i = 0; i < m; ++i) { t_num += pa[i] * ba[i]; t_den += ba[i] * ba[i]; }
        double t = (t_den > 1e-15) ? t_num / t_den : 0.0;
        double d = 0.0;
        for (int i = 0; i < m; ++i) { double x = pa[i] - t * ba[i]; d += x * x; }
        return std::sqrt(d);
    }

    /// Minkowski distance between two arrays of length m.
    static double minkowski_dist(const double* a, const double* b, double p, int m) {
        double sum = 0.0;
        for (int i = 0; i < m; ++i) sum += std::pow(std::abs(a[i] - b[i]), p);
        return std::pow(sum, 1.0 / p);
    }

    /// For each objective axis i, find the solution in `front` closest to that axis
    /// (minimum perpendicular distance). Returns pop-level indices.
    std::vector<int> find_extremes(Population<>& pop, const std::vector<int>& front,
                                   int n_obj) const {
        std::vector<double> A(n_obj, 0.0);
        std::vector<int> extremes;
        extremes.reserve(n_obj);
        for (int o = 0; o < n_obj; ++o) {
            std::vector<double> B(n_obj, 0.0);
            B[o] = 1.0;
            double best   = std::numeric_limits<double>::infinity();
            int    best_i = front[0];
            for (int idx : front) {
                double d = point2line_dist(pop.objectives_ptr(idx), A, B, n_obj);
                if (d < best) { best = d; best_i = idx; }
            }
            extremes.push_back(best_i);
        }
        return extremes;
    }

    /// Gaussian elimination solving Ax = b; returns x.
    static std::vector<double> gaussian_elimination(std::vector<std::vector<double>> A,
                                                    std::vector<double>              b) {
        int N = static_cast<int>(A.size());
        for (int i = 0; i < N; ++i) A[i].push_back(b[i]); // augmented
        for (int base = 0; base < N - 1; ++base) {
            for (int target = base + 1; target < N; ++target) {
                if (std::abs(A[base][base]) < 1e-15) continue;
                double ratio = A[target][base] / A[base][base];
                for (int term = 0; term <= N; ++term)
                    A[target][term] -= A[base][term] * ratio;
            }
        }
        std::vector<double> x(N, 0.0);
        for (int i = N - 1; i >= 0; --i) {
            for (int known = i + 1; known < N; ++known)
                A[i][N] -= A[i][known] * x[known];
            x[i] = (std::abs(A[i][i]) > 1e-15) ? A[i][N] / A[i][i] : 0.0;
        }
        return x;
    }

    /// Construct hyperplane through extreme points → intercepts.
    /// Handles duplicate extreme points (same individual selected for two objectives).
    std::vector<double> construct_hyperplane(Population<>& pop,
                                             const std::vector<int>& extremes, int n_obj) {
        // Detect duplicate extreme points
        bool duplicate = false;
        for (int i = 0; !duplicate && i < n_obj; ++i)
            for (int j = i + 1; !duplicate && j < n_obj; ++j)
                if (extremes[i] == extremes[j]) duplicate = true;

        if (duplicate) {
            std::vector<double> intercepts(n_obj);
            for (int f = 0; f < n_obj; ++f)
                intercepts[f] = pop.objective(extremes[f], f);
            return intercepts;
        }

        std::vector<std::vector<double>> A(n_obj, std::vector<double>(n_obj));
        std::vector<double>              bb(n_obj, 1.0);
        for (int i = 0; i < n_obj; ++i)
            for (int j = 0; j < n_obj; ++j)
                A[i][j] = pop.objective(extremes[i], j);

        auto x = gaussian_elimination(A, bb);
        std::vector<double> intercepts(n_obj);
        for (int f = 0; f < n_obj; ++f)
            intercepts[f] = (std::abs(x[f]) > 1e-15) ? 1.0 / x[f] : 1e15;
        return intercepts;
    }

    /// Estimate geometry exponent P from the normalised front.
    /// Mirrors jMetal computeGeometry(): central point = arg min distance to diagonal,
    /// then p = log(M) / log(1/average).
    double compute_geometry(const std::vector<std::vector<double>>& norm,
                            const std::vector<int>& extreme_local, int n_obj) {
        int m = static_cast<int>(norm.size());
        std::vector<double> utopia(n_obj, 0.0), nadir(n_obj, 1.0);

        // Distance to diagonal (utopia→nadir line) for non-extreme points
        std::vector<double> d(m, std::numeric_limits<double>::infinity());
        for (int i = 0; i < m; ++i) {
            bool is_extreme = (std::find(extreme_local.begin(), extreme_local.end(), i)
                               != extreme_local.end());
            if (!is_extreme)
                d[i] = point2line_dist(norm[i].data(), utopia, nadir, n_obj);
        }

        // Central point = min distance to diagonal
        int central = 0;
        double min_d = std::numeric_limits<double>::infinity();
        for (int i = 0; i < m; ++i)
            if (d[i] < min_d) { min_d = d[i]; central = i; }

        double sum = 0.0;
        for (int o = 0; o < n_obj; ++o) sum += norm[central][o];
        double average = sum / n_obj;

        if (average <= 0.0 || average >= 1.0) return 1.0;
        double p = std::log(static_cast<double>(n_obj)) / std::log(1.0 / average);
        if (std::isnan(p) || std::isinf(p) || p < 0.1) p = 1.0;
        return p;
    }

    // =========================================================================
    // Score computation (jMetal AGEMOEAEnvironmentalSelection)
    // =========================================================================

    /// Compute survival scores for front[0].
    /// Also sets P_ and intercepts_ for reuse by assign_convergence_scores().
    /// out_scores[i] corresponds to front[i] (local index).
    void compute_survival_scores(Population<>& pop, const std::vector<int>& front,
                                 int n_obj, std::vector<double>& out_scores) {
        P_ = 1.0;
        intercepts_.clear();
        int m = static_cast<int>(front.size());
        out_scores.assign(m, 0.0);

        if (m == 0) return;

        // Ideal point from front[0]
        std::vector<double> ideal(n_obj, std::numeric_limits<double>::max());
        for (int idx : front)
            for (int o = 0; o < n_obj; ++o)
                ideal[o] = std::min(ideal[o], pop.objective(idx, o));

        // Extreme points (pop-level indices)
        auto extremes_pop = find_extremes(pop, front, n_obj);

        // Map extreme pop-indices → front-local indices; mark +inf and add to selected
        std::vector<int> selected_local; // local indices of selected (extreme) points
        for (int ep : extremes_pop) {
            for (int li = 0; li < m; ++li) {
                if (front[li] == ep) {
                    bool already = (std::find(selected_local.begin(), selected_local.end(), li)
                                    != selected_local.end());
                    if (!already) selected_local.push_back(li);
                    out_scores[li] = std::numeric_limits<double>::infinity();
                    break;
                }
            }
        }

        // Small front: intercepts = max per objective; no full geometry computation
        if (m <= n_obj) {
            intercepts_.assign(n_obj, 0.0);
            for (int idx : front)
                for (int o = 0; o < n_obj; ++o)
                    intercepts_[o] = std::max(intercepts_[o], pop.objective(idx, o));
            return;
        }

        // Construct hyperplane → intercepts_
        intercepts_ = construct_hyperplane(pop, extremes_pop, n_obj);

        // Normalize front: (obj - ideal) / (intercept - ideal)
        std::vector<std::vector<double>> norm(m, std::vector<double>(n_obj));
        for (int i = 0; i < m; ++i) {
            for (int o = 0; o < n_obj; ++o) {
                double range = intercepts_[o] - ideal[o];
                if (std::abs(range) < 1e-12) range = 1e-12;
                norm[i][o] = (pop.objective(front[i], o) - ideal[o]) / range;
            }
        }

        // Geometry P_
        std::vector<int> extreme_local;
        for (int ep : extremes_pop) {
            for (int li = 0; li < m; ++li)
                if (front[li] == ep) { extreme_local.push_back(li); break; }
        }
        P_ = compute_geometry(norm, extreme_local, n_obj);

        // Proximity score nn[i] = minkowski(norm[i], utopia, P_)
        std::vector<double> utopia(n_obj, 0.0);
        std::vector<double> nn(m);
        for (int i = 0; i < m; ++i)
            nn[i] = minkowski_dist(norm[i].data(), utopia.data(), P_, n_obj);

        // Pairwise Minkowski distances (upper triangle, symmetric before normalization)
        std::vector<std::vector<double>> distances(m, std::vector<double>(m, 0.0));
        for (int i = 0; i < m - 1; ++i)
            for (int j = i + 1; j < m; ++j)
                distances[i][j] = distances[j][i] =
                    minkowski_dist(norm[i].data(), norm[j].data(), P_, n_obj);

        // Normalize distances[i][j] by nn[i] (matches jMetal: distances[i][j] /= nn[i])
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < m; ++j)
                if (nn[i] > 1e-15) distances[i][j] /= nn[i];

        // Greedy diversity score assignment: remaining → selected
        std::vector<int> remaining;
        remaining.reserve(m - selected_local.size());
        for (int i = 0; i < m; ++i)
            if (std::find(selected_local.begin(), selected_local.end(), i) == selected_local.end())
                remaining.push_back(i);

        while (!remaining.empty()) {
            double best_score = -1.0;
            int    best_li    = remaining[0];
            for (int i : remaining) {
                double min1 = std::numeric_limits<double>::infinity();
                double min2 = std::numeric_limits<double>::infinity();
                for (int j : selected_local) {
                    double dij = distances[i][j];
                    if (dij < min1) { min2 = min1; min1 = dij; }
                    else if (dij < min2) { min2 = dij; }
                }
                double score = min1 + min2;
                if (score >= best_score) { best_score = score; best_li = i; }
            }
            out_scores[best_li] = best_score;
            remaining.erase(std::find(remaining.begin(), remaining.end(), best_li));
            selected_local.push_back(best_li);
        }
    }

    /// Assign convergence (proximity) scores to a non-first-rank front.
    /// Uses stored P_ and intercepts_ from the last compute_survival_scores() call.
    /// Mirrors jMetal assignConvergenceScore(): normalise by intercepts only, then
    /// compute Minkowski distance to utopia [0,...,0].
    void assign_convergence_scores(Population<>& pop, const std::vector<int>& front,
                                   int n_obj, std::vector<double>& out_scores) {
        int m = static_cast<int>(front.size());
        out_scores.resize(m);
        std::vector<double> utopia(n_obj, 0.0);
        for (int i = 0; i < m; ++i) {
            std::vector<double> norm(n_obj);
            for (int o = 0; o < n_obj; ++o) {
                double denom = (!intercepts_.empty() && intercepts_[o] > 1e-15)
                                   ? intercepts_[o] : 1.0;
                norm[o] = pop.objective(front[i], o) / denom;
            }
            out_scores[i] = minkowski_dist(norm.data(), utopia.data(), P_, n_obj);
        }
    }

    // =========================================================================
    // Environmental selection
    // =========================================================================

    /// AGE-MOEA environmental selection: reduce pop to target_size.
    /// Updates scores_ and ranks_ (used by tournament in next generation).
    void environmental_selection(Population<>& pop, int target_size) {
        const int n_obj = pop.n_obj;
        const int n     = pop.pop_size;

        // Non-dominated sort
        auto fronts = fast_non_dominated_sort(pop);

        // Rank lookup indexed by pop position
        std::vector<int> temp_ranks(n, 0);
        for (int r = 0; r < static_cast<int>(fronts.size()); ++r)
            for (int idx : fronts[r])
                temp_ranks[idx] = r;

        // Collect whole fronts until we exceed target_size
        std::vector<int> selected;
        selected.reserve(target_size);
        int critical_rank = -1;

        for (int r = 0; r < static_cast<int>(fronts.size()); ++r) {
            const auto& fr = fronts[r];
            if (static_cast<int>(selected.size() + fr.size()) <= target_size) {
                selected.insert(selected.end(), fr.begin(), fr.end());
            } else {
                critical_rank = r;
                break;
            }
        }

        // All-individual score table (indexed by pop position, filled as needed)
        std::vector<double> all_scores(n, 0.0);

        // Always compute survival scores for front[0] (sets P_ and intercepts_)
        {
            std::vector<double> s0;
            compute_survival_scores(pop, fronts[0], n_obj, s0);
            for (int i = 0; i < static_cast<int>(fronts[0].size()); ++i)
                all_scores[fronts[0][i]] = s0[i];
        }

        if (critical_rank == -1) {
            // All individuals fit — assign convergence scores to fronts > 0
            for (int r = 1; r < static_cast<int>(fronts.size()); ++r) {
                std::vector<double> cs;
                assign_convergence_scores(pop, fronts[r], n_obj, cs);
                for (int i = 0; i < static_cast<int>(fronts[r].size()); ++i)
                    all_scores[fronts[r][i]] = cs[i];
            }
            // selected already contains all individuals; reorder and update
            reorder_population(pop, selected, all_scores, temp_ranks);
            pop.pop_size = target_size;
            return;
        }

        if (critical_rank == 0) {
            // front[0] is the critical front: sort by survival score (higher first)
            const auto& fr = fronts[0];
            std::vector<int> order(fr.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](int a, int b) {
                return all_scores[fr[a]] > all_scores[fr[b]];
            });
            selected.clear();
            int take = std::min(target_size, static_cast<int>(fr.size()));
            for (int k = 0; k < take; ++k) selected.push_back(fr[order[k]]);
        } else {
            // Critical front has rank > 0: assign convergence scores, sort, take remainder
            const auto& cf = fronts[critical_rank];
            std::vector<double> cs;
            assign_convergence_scores(pop, cf, n_obj, cs);
            for (int i = 0; i < static_cast<int>(cf.size()); ++i)
                all_scores[cf[i]] = cs[i];

            std::vector<int> order(cf.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](int a, int b) {
                return cs[a] > cs[b];
            });
            int slots = target_size - static_cast<int>(selected.size());
            for (int k = 0; k < slots && k < static_cast<int>(cf.size()); ++k)
                selected.push_back(cf[order[k]]);
        }

        reorder_population(pop, selected, all_scores, temp_ranks);
        pop.pop_size = target_size;
    }

    /// Reorder population so that selected[i] moves to position i.
    /// Also updates scores_ and ranks_ arrays.
    void reorder_population(Population<>& pop, const std::vector<int>& selected,
                            const std::vector<double>& all_scores,
                            const std::vector<int>&    all_ranks) {
        int m = static_cast<int>(selected.size());
        Population<> tmp(m, pop.dim, pop.n_obj, pop.n_const);
        tmp.lower_bounds = pop.lower_bounds;
        tmp.upper_bounds = pop.upper_bounds;

        for (int i = 0; i < m; ++i) {
            int src = selected[i];
            std::copy_n(pop.genes_ptr(src),      pop.dim,   tmp.genes_ptr(i));
            std::copy_n(pop.objectives_ptr(src),  pop.n_obj, tmp.objectives_ptr(i));
            tmp.flags[i] = pop.flags[src];
        }
        for (int i = 0; i < m; ++i) {
            std::copy_n(tmp.genes_ptr(i),      pop.dim,   pop.genes_ptr(i));
            std::copy_n(tmp.objectives_ptr(i),  pop.n_obj, pop.objectives_ptr(i));
            pop.flags[i] = tmp.flags[i];
        }

        scores_.resize(m);
        ranks_.resize(m);
        for (int i = 0; i < m; ++i) {
            scores_[i] = all_scores[selected[i]];
            ranks_[i]  = all_ranks[selected[i]];
        }
    }
};

} // namespace ea
