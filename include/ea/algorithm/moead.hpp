#pragma once
/// @file moead.hpp
/// @brief MOEA/D (Multi-Objective Evolutionary Algorithm based on Decomposition)
/// Reference: Zhang & Li, "MOEA/D: A Multiobjective Evolutionary Algorithm Based
/// on Decomposition", IEEE Trans. Evol. Comput., 2007.
///
/// This is the classic MOEA/D with Tchebycheff or Weighted Sum aggregation.
/// Uses the SoA Population, deducing this, and Concepts.
///
/// jMetal reference classes:
/// - org.uma.jmetal.algorithm.multiobjective.moead.AbstractMOEAD
/// - org.uma.jmetal.algorithm.multiobjective.moead.MOEAD
/// - org.uma.jmetal.component.algorithm.multiobjective.MOEADBuilder
///
/// Key differences from jMetal:
/// - SoA instead of AoS (list of Solution objects)
/// - Deducing this instead of virtual methods
/// - std::vector instead of ArrayList
/// - Batch-friendly evaluation

#include <ea/core/population.hpp>
#include <ea/core/concepts.hpp>
#include <ea/util/random.hpp>
#include <ea/util/aggregation.hpp>
#include <string_view>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>
#include <string>

namespace ea {

/// MOEA/D algorithm — decompose a multi-objective problem into a set of
/// single-objective subproblems and solve them cooperatively.
///
/// Usage:
///   MOEAD<SBXCrossover, PolynomialMutation> moead;
///   moead.crossover.distribution_index = 20.0;
///   moead.mutation.distribution_index = 20.0;
///   moead.pop_size = 100;
///   moead.max_evals = 25000;
///   moead.aggregation_type = AggregationType::Tchebycheff;
///   moead.run(pop, problem);
///
/// @tparam CX Crossover operator type (must satisfy Crossover concept)
/// @tparam MT Mutation operator type (must satisfy Mutation concept)
///
/// The algorithm maintains:
/// - A set of weight vectors (lambda), one per subproblem
/// - A neighborhood for each subproblem (nearest weight vectors)
/// - An ideal point (best value per objective)
///
/// Each iteration:
/// 1. Random permutation of subproblem indices
/// 2. For each subproblem:
///    a. Select parents from neighborhood or population
///    b. Apply crossover + mutation
///    c. Evaluate offspring
///    d. Update ideal point
///    e. Replace neighbors that are improved by the offspring
template<typename CX, typename MT>
struct MOEAD {
    CX crossover;
    MT mutation;

    int pop_size = 100;              ///< Number of subproblems (N)
    int max_evals = 25000;           ///< Maximum number of function evaluations
    int neighbor_size = -1;          ///< Neighborhood size (T), -1 = pop_size/10
    int max_replaced_solutions = 2;  ///< Max replacements per offspring (nr)
    double neighborhood_prob = 0.9;  ///< Probability of using neighborhood (delta)
    AggregationType aggregation_type = AggregationType::Tchebycheff;
    bool normalize = false;          ///< Whether to normalize using ideal/nadir
    std::string weight_directory;    ///< Directory with weight vector files (empty = generate)

    static constexpr std::string_view name() { return "MOEA/D"; }

    // ---- Internal state (allocated during run) ----
    std::vector<std::vector<double>> lambda_;       ///< Weight vectors [pop_size][n_obj]
    std::vector<std::vector<int>> neighborhood_;    ///< Neighborhood indices [pop_size][T]
    std::vector<double> ideal_point_;                ///< Ideal point [n_obj]
    std::vector<double> nadir_point_;                ///< Nadir point [n_obj]
    int evals_ = 0;

    enum class NeighborType : uint8_t { Neighbor, Population };

    /// Run MOEA/D on the given population.
    /// @param pop Population with genes initialized and bounds set.
    /// @param problem Callable: void(Population&, int) — evaluates individual's objectives
    template<typename Problem>
    void run(this auto& self, Population& pop, Problem&& problem) {
        const int n = self.pop_size;
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;

        if (pop.pop_size != n) pop.resize(n);

        // Resolve neighbor size
        int T = self.neighbor_size < 0 ? n / 10 : self.neighbor_size;
        if (T < 2) T = 2;
        if (T > n - 1) T = n - 1;

        // Allocate state
        self.lambda_.assign(n, std::vector<double>(n_obj, 0.0));
        self.neighborhood_.assign(n, std::vector<int>(T, 0));
        self.ideal_point_.assign(n_obj, std::numeric_limits<double>::max());
        self.nadir_point_.assign(n_obj, std::numeric_limits<double>::lowest());
        self.evals_ = 0;

        // === Initialize weight vectors ===
        self.initialize_uniform_weights(n, n_obj);

        // === Initialize neighborhoods ===
        self.initialize_neighborhood(n, T);

        // === Evaluate initial population ===
        for (int i = 0; i < n; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
            self.update_ideal_point(pop, i);
            self.update_nadir_point(pop, i);
        }
        self.evals_ = n;

        // Scratch population for offspring (single offspring per subproblem per iteration)
        // We need at least 2 slots for arity-2 crossover operators (SBX writes 2 children)
        // And up to what DE crossover needs. Use a scratch population with enough room.
        Population scratch(std::max(4, CX::arity() + 1), dim, n_obj, pop.encoding, pop.n_const);
        scratch.lower_bounds = pop.lower_bounds;
        scratch.upper_bounds = pop.upper_bounds;

        // Main loop
        while (self.evals_ < self.max_evals) {
            // Random permutation of subproblem indices
            auto permutation = Random::instance().permutation(n);

            for (int i = 0; i < n && self.evals_ < self.max_evals; ++i) {
                int sub_problem_id = permutation[i];

                // Choose neighbor type
                NeighborType neighbor_type = self.choose_neighbor_type();

                // Parent selection
                auto parents = self.parent_selection(sub_problem_id, neighbor_type, n, T);

                // Copy parents to scratch population positions
                for (size_t p = 0; p < parents.size(); ++p) {
                    for (int j = 0; j < dim; ++j) {
                        scratch.gene(static_cast<int>(p), j) = pop.gene(parents[p], j);
                    }
                    std::copy_n(pop.objectives_ptr(parents[p]), n_obj,
                                scratch.objectives_ptr(static_cast<int>(p)));
                    scratch.set_evaluated(static_cast<int>(p), true);
                }

                // Apply crossover on scratch (writes child at child_start)
                // For SBX: arity()=2, writes to child_start and child_start+1
                // We only need child_start=0 (single offspring for MOEA/D)
                if constexpr (CX::arity() == 2) {
                    // Two-parent crossover: apply to first two parents in scratch
                    self.crossover.apply(scratch, 0, 1, 0);
                } else if constexpr (CX::arity() == 4) {
                    // DE-style: needs 4 parents (target + 3 for difference vectors)
                    // But we may only have 2-3 parents selected.
                    // For MOEA/D with DE, typically use: target + 2 random from neighborhood
                    // Simplified: use parents[0]=target, parents[1]=r1, parents[2]=r2
                    // Need a 4th parent, duplicate one
                    if (parents.size() >= 3) {
                        // Copy parent 3 to scratch position 3 if needed
                        if (parents.size() == 3) {
                            for (int j = 0; j < dim; ++j) {
                                scratch.gene(3, j) = scratch.gene(2, j);
                            }
                        }
                        const int parent_indices[4] = {0, 1, 2, 3};
                        self.crossover.apply(scratch, parent_indices, 4, 0, 0);
                    } else {
                        // Not enough parents for DE — copy parent 0 as-is
                        // (child already copied at position 0)
                    }
                }

                // Apply mutation on scratch position 0 (the primary child)
                self.mutation.apply(scratch, 0);

                // Evaluate offspring
                if (!scratch.evaluated(0)) {
                    problem(scratch, 0);
                    scratch.set_evaluated(0, true);
                }
                self.evals_++;

                // Update ideal and nadir points
                self.update_ideal_point(scratch, 0);
                self.update_nadir_point(scratch, 0);

                // Update neighborhood
                self.update_neighborhood(pop, scratch, sub_problem_id, neighbor_type, n, T);
            }
        }
    }

private:
    // ============================================================
    // Weight vector initialization
    // ============================================================

    /// Initialize uniform weight vectors.
    /// For 2-objective problems, generate evenly distributed weights.
    /// For >2 objectives, try to load from file; if unavailable, generate simplex lattice.
    void initialize_uniform_weights(this auto& self, int n, int n_obj) {
        if (n_obj == 2) {
            // Evenly distributed on the unit interval
            for (int i = 0; i < n; ++i) {
                double a = 1.0 * i / (n - 1);
                self.lambda_[i][0] = a;
                self.lambda_[i][1] = 1.0 - a;
            }
        } else if (n_obj == 3) {
            // Simplex lattice with H divisions where C(H+2,2) ≈ n
            int H = static_cast<int>(std::sqrt(2.0 * n)) - 1;
            if (H < 1) H = 1;
            int count = 0;
            for (int i = 0; i <= H && count < n; ++i) {
                for (int j = 0; j <= H - i && count < n; ++j) {
                    int k = H - i - j;
                    self.lambda_[count][0] = i / static_cast<double>(H);
                    self.lambda_[count][1] = j / static_cast<double>(H);
                    self.lambda_[count][2] = k / static_cast<double>(H);
                    ++count;
                }
            }
            // Fill remaining with random weights if count < n
            auto& rng = Random::instance();
            while (count < n) {
                double sum = 0.0;
                for (int o = 0; o < n_obj; ++o) {
                    self.lambda_[count][o] = rng.uniform();
                    sum += self.lambda_[count][o];
                }
                for (int o = 0; o < n_obj; ++o) {
                    self.lambda_[count][o] /= sum;
                }
                ++count;
            }
        } else {
            // Try to load from file first
            if (!self.weight_directory.empty()) {
                std::string dataFileName = "W" + std::to_string(n_obj) + "D_" +
                                           std::to_string(n) + ".dat";
                std::string path = self.weight_directory + "/" + dataFileName;
                std::ifstream file(path);
                if (file.is_open()) {
                    int i = 0;
                    std::string line;
                    while (std::getline(file, line) && i < n) {
                        std::istringstream iss(line);
                        int j = 0;
                        double value;
                        while (iss >> value && j < n_obj) {
                            self.lambda_[i][j] = value;
                            ++j;
                        }
                        ++i;
                    }
                    if (i == n) return; // Successfully loaded
                }
            }
            // Fallback: generate random normalized weight vectors
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

    // ============================================================
    // Neighborhood initialization
    // ============================================================

    /// Initialize neighborhoods based on Euclidean distance between weight vectors.
    /// For each subproblem i, find the T nearest weight vectors.
    void initialize_neighborhood(this auto& self, int n, int T) {
        std::vector<double> dists(n);
        std::vector<int> indices(n);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                dists[j] = self.dist_vector(
                    self.lambda_[i].data(), self.lambda_[j].data(),
                    static_cast<int>(self.lambda_[i].size()));
                indices[j] = j;
            }
            // Sort indices by distance (ascending)
            std::partial_sort(indices.begin(), indices.begin() + T, indices.end(),
                [&dists](int a, int b) { return dists[a] < dists[b]; });
            for (int t = 0; t < T; ++t) {
                self.neighborhood_[i][t] = indices[t];
            }
        }
    }

    /// Euclidean distance between two weight vectors.
    static double dist_vector(const double* a, const double* b, int n_obj) {
        double sum = 0.0;
        for (int i = 0; i < n_obj; ++i) {
            double d = a[i] - b[i];
            sum += d * d;
        }
        return std::sqrt(sum);
    }

    // ============================================================
    // Neighbor type selection
    // ============================================================

    /// Choose whether to use neighborhood or whole population for mating.
    NeighborType choose_neighbor_type(this auto& self) {
        return Random::instance().coin_flip(self.neighborhood_prob)
                   ? NeighborType::Neighbor
                   : NeighborType::Population;
    }

    // ============================================================
    // Parent selection
    // ============================================================

    /// Select parents for crossover.
    /// Returns vector of parent indices.
    std::vector<int> parent_selection(this auto& self, int sub_problem_id,
                                      NeighborType neighbor_type, int n, int T) {
        auto mating_pool = self.mating_selection(sub_problem_id, 2, neighbor_type, n, T);
        // Add current subproblem as third "parent" (for DE-style crossover)
        mating_pool.push_back(sub_problem_id);
        return mating_pool;
    }

    /// Mating selection: select distinct individuals from neighborhood or population.
    std::vector<int> mating_selection(this auto& self, int subproblem_id,
                                      int number_to_select, NeighborType neighbor_type,
                                      int n, int T) {
        std::vector<int> selected;
        selected.reserve(number_to_select);
        auto& rng = Random::instance();

        int pool_size = (neighbor_type == NeighborType::Neighbor) ? T : n;

        while (static_cast<int>(selected.size()) < number_to_select) {
            int candidate;
            if (neighbor_type == NeighborType::Neighbor) {
                int idx = rng.uniform_int(0, T - 1);
                candidate = self.neighborhood_[subproblem_id][idx];
            } else {
                candidate = rng.uniform_int(0, n - 1);
            }

            bool duplicate = false;
            for (int s : selected) {
                if (s == candidate) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                selected.push_back(candidate);
            }
        }

        return selected;
    }

    // ============================================================
    // Fitness computation
    // ============================================================

    /// Compute scalar fitness for individual i in population pop
    /// using weight vector w_idx.
    double fitness(this auto& self, const Population& pop, int individual_idx,
                   int w_idx, int n_obj) {
        AggregationFunction aggr;
        aggr.type = self.aggregation_type;
        aggr.normalize = self.normalize;

        return aggr.compute(
            pop.objectives_ptr(individual_idx),
            self.lambda_[w_idx].data(),
            self.ideal_point_.data(),
            self.nadir_point_.data(),
            n_obj);
    }

    // ============================================================
    // Ideal / Nadir point updates
    // ============================================================

    /// Update ideal point with individual's objectives.
    void update_ideal_point(this auto& self, const Population& pop, int idx) {
        for (int o = 0; o < pop.n_obj; ++o) {
            double val = pop.objective(idx, o);
            if (val < self.ideal_point_[o]) {
                self.ideal_point_[o] = val;
            }
        }
    }

    /// Update nadir point with individual's objectives.
    void update_nadir_point(this auto& self, const Population& pop, int idx) {
        for (int o = 0; o < pop.n_obj; ++o) {
            double val = pop.objective(idx, o);
            if (val > self.nadir_point_[o]) {
                self.nadir_point_[o] = val;
            }
        }
    }

    // ============================================================
    // Neighborhood update
    // ============================================================

    /// Update neighborhood with the new offspring.
    /// Replace at most max_replaced_solutions neighbors that are improved.
    void update_neighborhood(this auto& self, Population& pop,
                             const Population& offspring,
                             int sub_problem_id, NeighborType neighbor_type,
                             int n, int T) {
        int size = (neighbor_type == NeighborType::Neighbor) ? T : n;
        int n_obj = pop.n_obj;

        std::vector<int> perm(size);
        std::iota(perm.begin(), perm.end(), 0);
        Random::instance().shuffle(std::span<int>(perm.data(), perm.size()));

        int replaced = 0;
        for (int i = 0; i < size; ++i) {
            int k;
            if (neighbor_type == NeighborType::Neighbor) {
                k = self.neighborhood_[sub_problem_id][perm[i]];
            } else {
                k = perm[i];
            }

            double f1 = self.fitness(pop, k, k, n_obj);
            double f2 = self.fitness(offspring, 0, k, n_obj);

            if (f2 < f1) {
                // Copy offspring to position k in population
                for (int j = 0; j < pop.dim; ++j) {
                    pop.gene(k, j) = offspring.gene(0, j);
                }
                for (int o = 0; o < n_obj; ++o) {
                    pop.objective(k, o) = offspring.objective(0, o);
                }
                pop.set_evaluated(k, true);
                replaced++;
            }

            if (replaced >= self.max_replaced_solutions) {
                break;
            }
        }
    }
};

} // namespace ea
