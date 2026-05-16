#pragma once
/// @file spea2.hpp
/// @brief SPEA2 (Strength Pareto Evolutionary Algorithm 2)
/// Reference: Zitzler, E., Laumanns, M., & Thiele, L. "SPEA2: Improving the
/// Strength Pareto Evolutionary Algorithm", TIK Report 103, 2001.
///
/// Uses Strength Raw Fitness density estimator for environmental selection.

#include <ea/core/population.hpp>
#include <ea/core/comparator.hpp>
#include <ea/util/random.hpp>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

namespace ea {

/// SPEA2 algorithm — Strength Pareto Evolutionary Algorithm 2.
/// Maintains an external archive of non-dominated solutions.
/// Uses strength-based fitness and k-th nearest neighbor density.
template<typename CX, typename MT>
struct SPEA2 {
    CX crossover;
    MT mutation;

    int pop_size = 100;
    int archive_size = 100;  // = pop_size for standard SPEA2
    int max_evals = 25000;
    int k_nearest = 1;        // k for density estimation (default: 1)

    static constexpr std::string_view name() { return "SPEA2"; }

    /// Run SPEA2 on the given population.
    template<typename Problem>
    void run(this auto& self, Population& pop, Problem&& problem) {
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;
        const int N = self.pop_size;

        // === Initialize population ===
        if (pop.pop_size != N) pop.resize(N);
        for (int i = 0; i < N; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }

        // Archive population
        Population archive(N, dim, n_obj, pop.encoding, pop.n_const);
        archive.lower_bounds = pop.lower_bounds;
        archive.upper_bounds = pop.upper_bounds;
        int archive_count = 0;

        // Union population (population + archive)
        Population union_pop(2 * N, dim, n_obj, pop.encoding, pop.n_const);
        union_pop.lower_bounds = pop.lower_bounds;
        union_pop.upper_bounds = pop.upper_bounds;

        // Fitness values (raw fitness + density)
        std::vector<double> fitness(2 * N);

        int evals = N;

        while (evals < self.max_evals) {
            // === Combine population and archive ===
            int union_size = N + archive_count;
            for (int i = 0; i < N; ++i) {
                std::copy_n(pop.genes_ptr(i), dim, union_pop.genes_ptr(i));
                std::copy_n(pop.objectives_ptr(i), n_obj, union_pop.objectives_ptr(i));
                union_pop.flags[i] = pop.flags[i];
            }
            for (int i = 0; i < archive_count; ++i) {
                std::copy_n(archive.genes_ptr(i), dim, union_pop.genes_ptr(N + i));
                std::copy_n(archive.objectives_ptr(i), n_obj, union_pop.objectives_ptr(N + i));
                union_pop.flags[N + i] = archive.flags[i];
            }
            union_pop.pop_size = union_size;

            // === Compute fitness ===
            compute_spea2_fitness(union_pop, fitness);

            // === Environmental selection → new archive ===
            archive_count = spea2_environmental_selection(union_pop, fitness, archive, N);

            // === Selection + Crossover + Mutation → new population ===
            Population offspring(N, dim, n_obj, pop.encoding, pop.n_const);
            offspring.lower_bounds = pop.lower_bounds;
            offspring.upper_bounds = pop.upper_bounds;

            auto& rng = Random::instance();
            for (int i = 0; i < N; i += 2) {
                // Binary tournament selection from archive
                int p1 = spea2_tournament(archive, archive_count, fitness, N, rng);
                int p2 = spea2_tournament(archive, archive_count, fitness, N, rng);

                // Copy parents to offspring
                for (int j = 0; j < dim; ++j) {
                    offspring.gene(i, j) = archive.gene(p1, j);
                    offspring.gene(i + 1, j) = archive.gene(p2, j);
                }
                std::copy_n(archive.objectives_ptr(p1), n_obj, offspring.objectives_ptr(i));
                std::copy_n(archive.objectives_ptr(p2), n_obj, offspring.objectives_ptr(i + 1));

                // Apply crossover
                self.crossover.apply(offspring, i, i + 1, i);

                // Apply mutation
                self.mutation.apply(offspring, i);
                self.mutation.apply(offspring, i + 1);
            }

            // Evaluate offspring
            for (int i = 0; i < N; ++i) {
                if (!offspring.evaluated(i)) {
                    problem(offspring, i);
                    offspring.set_evaluated(i, true);
                    evals++;
                    if (evals >= self.max_evals) break;
                }
            }

            // Copy offspring to population
            for (int i = 0; i < N; ++i) {
                std::copy_n(offspring.genes_ptr(i), dim, pop.genes_ptr(i));
                std::copy_n(offspring.objectives_ptr(i), n_obj, pop.objectives_ptr(i));
                pop.flags[i] = offspring.flags[i];
            }
        }

        // Final: copy archive back to pop
        for (int i = 0; i < std::min(archive_count, N); ++i) {
            std::copy_n(archive.genes_ptr(i), dim, pop.genes_ptr(i));
            std::copy_n(archive.objectives_ptr(i), n_obj, pop.objectives_ptr(i));
            pop.flags[i] = archive.flags[i];
        }
    }

private:
    /// Compute SPEA2 fitness: S(i) = strength (number of solutions dominated by i)
    /// R(i) = sum of strengths of solutions dominating i
    /// D(i) = 1 / (sigma_k + 2) where sigma_k is distance to k-th nearest neighbor
    /// Fitness(i) = R(i) + D(i)  (lower is better; < 1 means non-dominated)
    static void compute_spea2_fitness(const Population& pop, std::vector<double>& fitness) {
        const int n = pop.pop_size;
        fitness.resize(n);

        // Step 1: Compute strength S(i)
        std::vector<int> strength(n, 0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                if (compare_dominance(pop, i, j) == Dominance::Dominates) {
                    strength[i]++;
                }
            }
        }

        // Step 2: Compute raw fitness R(i) = sum of S(j) for all j that dominate i
        for (int i = 0; i < n; ++i) {
            double raw = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                if (compare_dominance(pop, j, i) == Dominance::Dominates) {
                    raw += strength[j];
                }
            }
            fitness[i] = raw;
        }

        // Step 3: Compute density D(i) = 1 / (sigma_k + 2)
        // where sigma_k is the distance to the k-th nearest neighbor
        int k = std::min(n - 1, 1); // k = 1 by default

        for (int i = 0; i < n; ++i) {
            // Compute distances to all other solutions
            std::vector<double> distances;
            distances.reserve(n - 1);
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                double dist = 0.0;
                for (int o = 0; o < pop.n_obj; ++o) {
                    double diff = pop.objective(i, o) - pop.objective(j, o);
                    dist += diff * diff;
                }
                distances.push_back(std::sqrt(dist));
            }

            // Sort distances
            std::sort(distances.begin(), distances.end());

            // k-th nearest neighbor distance
            double sigma_k = distances[k];
            fitness[i] += 1.0 / (sigma_k + 2.0);
        }
    }

    /// SPEA2 environmental selection: select best N solutions from union
    static int spea2_environmental_selection(const Population& union_pop,
                                              const std::vector<double>& fitness,
                                              Population& archive, int archive_size) {
        const int n = union_pop.pop_size;
        const int dim = union_pop.dim;
        const int n_obj = union_pop.n_obj;

        // Collect non-dominated solutions (fitness < 1)
        std::vector<int> nondominated;
        std::vector<int> dominated;
        for (int i = 0; i < n; ++i) {
            if (fitness[i] < 1.0) {
                nondominated.push_back(i);
            } else {
                dominated.push_back(i);
            }
        }

        std::vector<int> selected;

        if (static_cast<int>(nondominated.size()) < archive_size) {
            // Not enough non-dominated — fill with best dominated
            selected = nondominated;
            std::sort(dominated.begin(), dominated.end(), [&](int a, int b) {
                return fitness[a] < fitness[b];
            });
            int remaining = archive_size - static_cast<int>(selected.size());
            for (int i = 0; i < remaining && i < static_cast<int>(dominated.size()); ++i) {
                selected.push_back(dominated[i]);
            }
        } else if (static_cast<int>(nondominated.size()) == archive_size) {
            selected = nondominated;
        } else {
            // Too many non-dominated — truncate using density
            // Remove solutions with smallest k-th nearest neighbor distance
            selected = nondominated;

            // Compute pairwise distances for non-dominated set
            while (static_cast<int>(selected.size()) > archive_size) {
                // Find solution with minimum distance to k-th nearest neighbor
                int worst_idx = -1;
                double min_dist = std::numeric_limits<double>::max();

                for (size_t si = 0; si < selected.size(); ++si) {
                    int i = selected[si];
                    std::vector<double> distances;
                    for (size_t sj = 0; sj < selected.size(); ++sj) {
                        if (si == sj) continue;
                        int j = selected[sj];
                        double dist = 0.0;
                        for (int o = 0; o < n_obj; ++o) {
                            double diff = union_pop.objective(i, o) - union_pop.objective(j, o);
                            dist += diff * diff;
                        }
                        distances.push_back(std::sqrt(dist));
                    }
                    std::sort(distances.begin(), distances.end());
                    if (!distances.empty() && distances[0] < min_dist) {
                        min_dist = distances[0];
                        worst_idx = static_cast<int>(si);
                    }
                }

                if (worst_idx >= 0) {
                    selected.erase(selected.begin() + worst_idx);
                } else {
                    break;
                }
            }
        }

        // Copy selected to archive
        int count = std::min(static_cast<int>(selected.size()), archive_size);
        for (int i = 0; i < count; ++i) {
            for (int j = 0; j < dim; ++j) {
                archive.gene(i, j) = union_pop.gene(selected[i], j);
            }
            for (int o = 0; o < n_obj; ++o) {
                archive.objective(i, o) = union_pop.objective(selected[i], o);
            }
            archive.flags[i] = union_pop.flags[selected[i]];
        }

        return count;
    }

    /// Binary tournament selection based on SPEA2 fitness
    static int spea2_tournament(const Population& archive, int archive_count,
                                 const std::vector<double>& fitness, int offset,
                                 Random& rng) {
        int a = rng.uniform_int(0, archive_count - 1);
        int b = rng.uniform_int(0, archive_count - 1);
        // Lower fitness is better
        return (fitness[a + offset] < fitness[b + offset]) ? a : b;
    }
};

} // namespace ea