#pragma once
/// @file mu_lambda_es.hpp
/// @brief (µ,λ) Evolution Strategy for single-objective real-valued optimisation.
///
/// Replicates the algorithmic behaviour of Bio4Res/ea (ea_cpp_original):
///   - k-tournament selection (by objective[0], minimisation)
///   - BLX-style crossover: arity=2, one child per pair (prob-gated)
///   - Gaussian mutation: one gene, multiplicative perturbation (prob-gated)
///   - Comma replacement: best min(λ,µ) offspring; if λ<µ fill with best parents
///
/// Parameters come from JSON via MuLambdaConfig (see mu_lambda_runner.cpp).
///
/// Usage with StatsCollector:
///   ea::MuLambdaES<BLXAlphaBetaCrossover, GaussianMutation> es;
///   es.pop_size = 100; es.offspring_size = 99; es.max_evals = 10'000'000;
///   es.tournament_k = 2;
///   es.crossover.crossover_probability = 0.9;
///   es.crossover.alpha = 0.1; es.crossover.beta = 0.1;
///   es.mutation.mutation_rate = 0.6321; es.mutation.step_size = 1.0;
///   StatsCollector stats;
///   stats.new_run(seed);
///   es.run(pop, problem, stats);
///   stats.close_run();

#include <algorithm>
#include <cmath>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <numeric>
#include <string_view>
#include <vector>

namespace ea {

/// (µ,λ) Evolution Strategy — single-objective, real-valued.
///
/// @tparam CX  Crossover operator (arity=2, one child per pair).
/// @tparam MT  Mutation operator (per-individual).
template <typename CX, typename MT>
struct MuLambdaES {
    CX crossover;
    MT mutation;

    int pop_size      = 100;         ///< µ — parent population size
    int offspring_size = 99;         ///< λ — offspring per generation
    int max_evals     = 10'000'000;  ///< budget in function evaluations
    int tournament_k  = 2;           ///< tournament size for parent selection

    static constexpr std::string_view name() { return "MuLambdaES"; }

    /// Run without stats collection.
    template <typename Problem>
    void run(this auto& self, Population<>& pop, Problem&& problem) {
        struct NoStats {
            void take(const Population<>&, uint64_t) {}
        } no_stats;
        self.run(pop, problem, no_stats);
    }

    /// Run with stats collection.
    /// @param pop     Population<> initialised and with bounds set. Size == pop_size.
    /// @param problem Callable void(Population<>&, int idx) — evaluates individual.
    /// @param stats   Object with void take(const Population<>&, uint64_t evals) method.
    template <typename Problem, typename Stats>
    void run(this auto& self, Population<>& pop, Problem&& problem, Stats& stats) {
        const int mu  = self.pop_size;
        const int lam = self.offspring_size;
        const int dim = pop.dim;

        if (pop.pop_size != mu)
            pop.resize(mu);

        // === Evaluate initial population ===
        for (int i = 0; i < mu; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }
        uint64_t evals = static_cast<uint64_t>(mu);
        stats.take(pop, evals);

        // --- Pre-allocate all buffers once outside the main loop ---
        // Eliminates ~303K malloc/free calls at 10M evals (next + off_order + par_order).
        Population<> offspring(lam, dim, pop.n_obj, pop.n_const);
        offspring.lower_bounds = pop.lower_bounds;
        offspring.upper_bounds = pop.upper_bounds;

        Population<> next(mu, dim, pop.n_obj, pop.n_const);
        next.lower_bounds = pop.lower_bounds;
        next.upper_bounds = pop.upper_bounds;

        std::vector<int> off_order(lam);
        std::vector<int> par_order(lam < mu ? mu : 0);

        const int n_from_offspring = std::min(lam, mu);

        auto& rng = Random::instance();

        while (evals < static_cast<uint64_t>(self.max_evals)) {
            // === 1. Generate offspring ===
            for (int c = 0; c < lam; ++c) {
                int p1 = tournament_select(pop, self.tournament_k, rng);
                int p2 = tournament_select(pop, self.tournament_k, rng);
                apply_crossover_one_child(self.crossover, pop, p1, p2, offspring, c);
                self.mutation.apply(offspring, c);
                offspring.set_evaluated(c, false);
            }

            // === 2. Evaluate offspring ===
            for (int c = 0; c < lam && evals < static_cast<uint64_t>(self.max_evals); ++c) {
                if (!offspring.evaluated(c)) {
                    problem(offspring, c);
                    offspring.set_evaluated(c, true);
                    ++evals;
                }
            }

            // === 3. Comma replacement (reuses pre-allocated buffers) ===
            std::iota(off_order.begin(), off_order.end(), 0);
            std::sort(off_order.begin(), off_order.end(), [&](int a, int b) {
                return offspring.objective(a, 0) < offspring.objective(b, 0);
            });

            if (lam < mu) {
                std::iota(par_order.begin(), par_order.end(), 0);
                std::sort(par_order.begin(), par_order.end(), [&](int a, int b) {
                    return pop.objective(a, 0) < pop.objective(b, 0);
                });
            }

            for (int i = 0; i < n_from_offspring; ++i)
                next.copy_from_other(offspring, off_order[i], i);

            for (int i = n_from_offspring; i < mu; ++i)
                next.copy_from_other(pop, par_order[i - n_from_offspring], i);

            // Swap buffers: pop gets the new generation, next keeps a valid
            // (reusable) buffer for the following iteration — no heap allocation.
            std::swap(pop.genes,      next.genes);
            std::swap(pop.objectives, next.objectives);
            std::swap(pop.flags,      next.flags);

            stats.take(pop, evals);
        }
    }

private:
    /// k-tournament selection on objective[0] (minimisation).
    static int tournament_select(const Population<>& pop, int k, Random& rng) {
        int best = rng.uniform_int(0, pop.pop_size - 1);
        for (int i = 1; i < k; ++i) {
            int cand = rng.uniform_int(0, pop.pop_size - 1);
            if (pop.objective(cand, 0) < pop.objective(best, 0))
                best = cand;
        }
        return best;
    }

    /// Apply crossover producing exactly ONE child at offspring[child_idx].
    ///
    /// ONE coin flip per pair (not per gene) — mirrors both ea_cpp_original
    /// VariationOperator::apply() and jMetal BLXAlphaCrossover.doCrossover() line 98.
    /// If the flip fails, child = verbatim copy of parent p1 (all genes).
    /// If it fires, BLX-alpha is applied to every gene.
    template <typename CX2>
    static void apply_crossover_one_child(CX2& cx, const Population<>& parents,
                                          int p1, int p2,
                                          Population<>& offspring, int child_idx) {
        auto& rng = Random::instance();
        const int dim = parents.dim;

        // One gate per pair — copy parent p1 if crossover does not fire
        if (!rng.coin_flip(cx.crossover_probability)) {
            for (int j = 0; j < dim; ++j)
                offspring.gene(child_idx, j) = parents.gene(p1, j);
            offspring.set_evaluated(child_idx, false);
            return;
        }

        for (int j = 0; j < dim; ++j) {
            double a = parents.gene(p1, j);
            double b = parents.gene(p2, j);

            double lo = parents.lower_bounds[j];
            double hi = parents.upper_bounds[j];
            double mn = std::min(a, b);
            double mx = std::max(a, b);
            double d  = mx - mn;

            // BLX-alpha exploration (symmetric: alpha == beta)
            double c_min = mn - cx.alpha * d;
            double c_max = mx + cx.beta  * d;
            double val   = c_min + rng.uniform() * (c_max - c_min);
            offspring.gene(child_idx, j) = std::clamp(val, lo, hi);
        }
        offspring.set_evaluated(child_idx, false);
    }
};

} // namespace ea
