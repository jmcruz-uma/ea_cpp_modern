#pragma once
/// @file moead_de.hpp
/// @brief MOEA/D-DE — MOEA/D with Differential Evolution operators.
/// @par Reference Li & Zhang, "Multiobjective Optimization Problems With Complicated
///            Pareto Sets, MOEA/D and NSGA-II", IEEE Trans. Evol. Comput., 2009.
///
/// Key differences from base MOEA/D:
/// - Uses DE/rand/1/bin crossover instead of SBX
/// - Uses polynomial mutation on the DE trial vector
/// - Maintains the decomposition + neighborhood structure
/// - Generally better for problems with complicated Pareto sets

#include <algorithm>
#include <cmath>
#include <ea/core/concepts.hpp>
#include <ea/core/population.hpp>
#include <ea/operator/crossover/de.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/util/aggregation.hpp>
#include <ea/util/random.hpp>
#include <functional>
#include <limits>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>

namespace ea {

/// MOEA/D-DE algorithm — decomposition-based MOEA with DE operators.
///
/// Usage:
///   MOEAD_DE moead_de;
///   moead_de.pop_size = 100;
///   moead_de.max_evals = 25000;
///   moead_de.aggregation_type = AggregationType::Tchebycheff;
///   moead_de.run(pop, problem);
///
/// The algorithm maintains:
/// - Weight vectors and neighborhoods (same as MOEA/D)
/// - Ideal/nadir points for normalization
/// - DE crossover: target + 2 difference vectors from neighborhood
/// - Polynomial mutation on DE trial vector
struct MOEAD_DE {
    DECrossover crossover;      ///< DE crossover (default: rand/1/bin, CR=1.0, F=0.5)
    PolynomialMutation mutation; ///< Polynomial mutation on trial vector

    int pop_size = 100;             ///< Number of subproblems (N)
    int max_evals = 25000;          ///< Maximum function evaluations
    int neighbor_size = -1;         ///< Neighborhood size (T), -1 = pop_size/10
    int max_replaced_solutions = 2; ///< Max replacements per offspring (nr)
    double neighborhood_prob = 0.9; ///< Probability of using neighborhood (delta)
    AggregationType aggregation_type = AggregationType::Tchebycheff;
    bool normalize = false; ///< Whether to normalize

    static constexpr std::string_view name() { return "MOEA/D-DE"; }

    // ---- Internal state ----
    std::vector<std::vector<double>> lambda_;
    std::vector<std::vector<int>> neighborhood_;
    std::vector<double> ideal_point_;
    std::vector<double> nadir_point_;
    int evals_ = 0;

    enum class NeighborType : uint8_t { Neighbor, WholePopulation };

    /// Run MOEA/D-DE on the given population.
    /// @param pop Population<> with genes initialized and bounds set.
    /// @param problem Callable: void(Population<>&, int) — evaluates individual
    template <EvalFunctor F> void run(this auto& self, Population<>& pop, F&& problem) {
        const int n = self.pop_size;
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;

        if (pop.pop_size != n)
            pop.resize(n);

        int T = self.neighbor_size < 0 ? n / 10 : self.neighbor_size;
        if (T < 2) T = 2;
        if (T > n - 1) T = n - 1;

        // Allocate state
        self.lambda_.assign(n, std::vector<double>(n_obj, 0.0));
        self.neighborhood_.assign(n, std::vector<int>(T, 0));
        self.ideal_point_.assign(n_obj, std::numeric_limits<double>::max());
        self.nadir_point_.assign(n_obj, std::numeric_limits<double>::lowest());
        self.evals_ = 0;

        self.initialize_uniform_weights(n, n_obj);
        self.initialize_neighborhood(n, T);

        // Evaluate initial population
        for (int i = 0; i < n; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
            self.update_ideal_point(pop, i);
            self.update_nadir_point(pop, i);
        }
        self.evals_ = n;

        // Scratch population for trial vector — allocated once, reused every iteration
        Population<> trial_pop(1, dim, n_obj, pop.n_const);
        trial_pop.lower_bounds = pop.lower_bounds;
        trial_pop.upper_bounds = pop.upper_bounds;

        // Main loop
        while (self.evals_ < self.max_evals) {
            auto permutation = Random::instance().permutation(n);

            for (int i = 0; i < n && self.evals_ < self.max_evals; ++i) {
                int sub_problem_id = permutation[i];
                NeighborType ntype = self.choose_neighbor_type();

                // === DE/rand/1: select 2 distinct neighbors as difference vector sources ===
                // Matches jMetal: matingSelection(subProblemId, 2, neighbourType)
                auto parents = self.mating_selection(sub_problem_id, 2, ntype, n, T);
                // parents[0] = r1, parents[1] = r2
                // trial = current + F*(r1 - r2), BIN crossover with current as base

                // Create trial vector via DE crossover + clamp to bounds
                // jMetal: DifferentialEvolutionCrossover.repairVariableValues() clips after DE
                auto& rng = Random::instance();
                int j_rand = rng.uniform_int(0, dim - 1);
                for (int j = 0; j < dim; ++j) {
                    double lb = trial_pop.lower_bounds[j];
                    double ub = trial_pop.upper_bounds[j];
                    double val;
                    if (rng.coin_flip(self.crossover.cr) || j == j_rand) {
                        val = pop.gene(sub_problem_id, j) +
                              self.crossover.f *
                                  (pop.gene(parents[0], j) - pop.gene(parents[1], j));
                    } else {
                        val = pop.gene(sub_problem_id, j);
                    }
                    trial_pop.gene(0, j) = std::clamp(val, lb, ub);
                }

                self.mutation.apply(trial_pop, 0);

                // Evaluate offspring
                problem(trial_pop, 0);
                trial_pop.set_evaluated(0, true);
                self.evals_++;

                self.update_ideal_point(trial_pop, 0);
                self.update_nadir_point(trial_pop, 0);

                self.update_neighborhood(pop, trial_pop, sub_problem_id, ntype, n, T);
            }
        }
    }

private:
    void initialize_uniform_weights(this auto& self, int n, int n_obj) {
        if (n_obj == 2) {
            for (int i = 0; i < n; ++i) {
                double a = 1.0 * i / (n - 1);
                self.lambda_[i][0] = a;
                self.lambda_[i][1] = 1.0 - a;
            }
        } else {
            auto& rng = Random::instance();
            for (int i = 0; i < n; ++i) {
                double sum = 0.0;
                for (int o = 0; o < n_obj; ++o) {
                    self.lambda_[i][o] = rng.uniform();
                    sum += self.lambda_[i][o];
                }
                for (int o = 0; o < n_obj; ++o) {
                    self.lambda_[i][o] /= sum;
                }
            }
        }
    }

    void initialize_neighborhood(this auto& self, int n, int T) {
        std::vector<double> dists(n);
        std::vector<int> indices(n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double sum = 0.0;
                for (size_t k = 0; k < self.lambda_[i].size(); ++k) {
                    double d = self.lambda_[i][k] - self.lambda_[j][k];
                    sum += d * d;
                }
                dists[j] = std::sqrt(sum);
                indices[j] = j;
            }
            std::partial_sort(indices.begin(), indices.begin() + T, indices.end(),
                              [&dists](int a, int b) { return dists[a] < dists[b]; });
            for (int t = 0; t < T; ++t) {
                self.neighborhood_[i][t] = indices[t];
            }
        }
    }

    NeighborType choose_neighbor_type(this auto& self) {
        return Random::instance().coin_flip(self.neighborhood_prob) ? NeighborType::Neighbor
                                                                    : NeighborType::WholePopulation;
    }

    std::vector<int> mating_selection(this auto& self, int subproblem_id,
                                      int number, NeighborType ntype, int n, int T) {
        std::vector<int> selected;
        selected.reserve(number);
        auto& rng = Random::instance();
        while (static_cast<int>(selected.size()) < number) {
            int candidate;
            if (ntype == NeighborType::Neighbor) {
                int idx = rng.uniform_int(0, T - 1);
                candidate = self.neighborhood_[subproblem_id][idx];
            } else {
                candidate = rng.uniform_int(0, n - 1);
            }
            bool dup = false;
            for (int s : selected) {
                if (s == candidate) { dup = true; break; }
            }
            if (!dup) {
                selected.push_back(candidate);
            }
        }
        return selected;
    }

    double fitness(this auto& self, const Population<>& pop, int idx, int w_idx, int n_obj) {
        AggregationFunction aggr;
        aggr.type = self.aggregation_type;
        aggr.normalize = self.normalize;
        return aggr.compute(pop.objectives_ptr(idx), self.lambda_[w_idx].data(),
                            self.ideal_point_.data(), self.nadir_point_.data(), n_obj);
    }

    void update_ideal_point(this auto& self, const Population<>& pop, int idx) {
        for (int o = 0; o < pop.n_obj; ++o) {
            double val = pop.objective(idx, o);
            if (val < self.ideal_point_[o]) self.ideal_point_[o] = val;
        }
    }

    void update_nadir_point(this auto& self, const Population<>& pop, int idx) {
        for (int o = 0; o < pop.n_obj; ++o) {
            double val = pop.objective(idx, o);
            if (val > self.nadir_point_[o]) self.nadir_point_[o] = val;
        }
    }

    void update_neighborhood(this auto& self, Population<>& pop, const Population<>& offspring,
                             int sub_problem_id, NeighborType ntype, int n, int T) {
        int size = (ntype == NeighborType::Neighbor) ? T : n;
        int n_obj = pop.n_obj;

        std::vector<int> perm(size);
        std::iota(perm.begin(), perm.end(), 0);
        Random::instance().shuffle(std::span<int>(perm.data(), perm.size()));

        int replaced = 0;
        for (int i = 0; i < size; ++i) {
            int k;
            if (ntype == NeighborType::Neighbor) {
                k = self.neighborhood_[sub_problem_id][perm[i]];
            } else {
                k = perm[i];
            }

            double f1 = self.fitness(pop, k, k, n_obj);
            double f2 = self.fitness(offspring, 0, k, n_obj);

            if (f2 < f1) {
                for (int j = 0; j < pop.dim; ++j) {
                    pop.gene(k, j) = offspring.gene(0, j);
                }
                for (int o = 0; o < n_obj; ++o) {
                    pop.objective(k, o) = offspring.objective(0, o);
                }
                pop.set_evaluated(k, true);
                replaced++;
            }
            if (replaced >= self.max_replaced_solutions) break;
        }
    }
};

} // namespace ea