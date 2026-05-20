#pragma once
/// @file spea2_replacement.hpp
/// @brief SPEA2 Replacement — archive truncation based on strength + density.
///
/// Reference: Zitzler, E., Laumanns, M., & Thiele, L. "SPEA2: Improving the
/// Strength Pareto Evolutionary Algorithm", TIK Report 103, 2001.
///
/// SPEA2 environmental selection from a combined population of parents + offspring:
/// 1. Compute SPEA2 fitness (strength + k-NN density) for all
/// 2. Select all non-dominated individuals (fitness < 1.0)
/// 3. If too few, fill with best dominated
/// 4. If too many, truncate by removing individuals with smallest k-th NN distance

#include <algorithm>
#include <cmath>
#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <limits>
#include <vector>

namespace ea {

/// SPEA2 environmental selection (replacement operator).
/// Selects archive_size individuals from combined population using
/// strength-based fitness and k-th nearest neighbor density.
struct SPEA2Replacement {
    int archive_size = 100; ///< Target archive size
    int k_nearest = 1;      ///< k for density estimation (default: 1)

    /// Perform SPEA2 environmental selection.
    /// @param combined      Combined population (parents + offspring)
    /// @param offspring_indices  Indices of offspring in combined (unused, for API compat)
    /// @param target_size   Target population size (typically = archive_size)
    /// @return Indices of selected individuals
    std::vector<int> replace(this SPEA2Replacement& self, Population& combined,
                               const std::vector<int>& offspring_indices, int target_size) {
        (void)offspring_indices;
        const int n = combined.pop_size;
        const int n_obj = combined.n_obj;

        // Compute SPEA2 fitness
        std::vector<double> fitness(n);
        compute_spea2_fitness(combined, fitness);

        // Separate non-dominated (fitness < 1.0) and dominated
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

        if (static_cast<int>(nondominated.size()) < target_size) {
            // Not enough non-dominated — fill with best dominated by fitness
            selected = nondominated;
            std::sort(dominated.begin(), dominated.end(),
                      [&fitness](int a, int b) { return fitness[a] < fitness[b]; });
            int remaining = target_size - static_cast<int>(selected.size());
            for (int i = 0; i < remaining && i < static_cast<int>(dominated.size()); ++i) {
                selected.push_back(dominated[i]);
            }
        } else if (static_cast<int>(nondominated.size()) == target_size) {
            selected = nondominated;
        } else {
            // Too many non-dominated — truncate using k-th NN distance
            selected = nondominated;

            while (static_cast<int>(selected.size()) > target_size) {
                // Find individual with minimum k-th nearest neighbor distance
                int worst_idx = -1;
                double min_dist = std::numeric_limits<double>::max();

                for (size_t si = 0; si < selected.size(); ++si) {
                    int i = selected[si];
                    std::vector<double> distances;
                    distances.reserve(selected.size() - 1);

                    for (size_t sj = 0; sj < selected.size(); ++sj) {
                        if (si == sj)
                            continue;
                        int j = selected[sj];
                        double dist = 0.0;
                        for (int o = 0; o < n_obj; ++o) {
                            double diff = combined.objective(i, o) - combined.objective(j, o);
                            dist += diff * diff;
                        }
                        distances.push_back(std::sqrt(dist));
                    }

                    std::sort(distances.begin(), distances.end());
                    double kth_dist = distances[std::min(self.k_nearest,
                                                          static_cast<int>(distances.size()) - 1)];

                    if (kth_dist < min_dist) {
                        min_dist = kth_dist;
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

        return selected;
    }

private:
    /// Compute SPEA2 fitness: R(i) + D(i)
    /// R(i) = sum of strengths of dominators
    /// D(i) = 1 / (sigma_k + 2) where sigma_k = k-th NN distance
    static void compute_spea2_fitness(const Population& pop, std::vector<double>& fitness) {
        const int n = pop.pop_size;
        fitness.resize(n);

        // Step 1: strength S(i) = number of individuals dominated by i
        std::vector<int> strength(n, 0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j)
                    continue;
                if (compare_dominance(pop, i, j) == Dominance::Dominates) {
                    strength[i]++;
                }
            }
        }

        // Step 2: raw fitness R(i) = sum of S(j) for all j that dominate i
        for (int i = 0; i < n; ++i) {
            double raw = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i == j)
                    continue;
                if (compare_dominance(pop, j, i) == Dominance::Dominates) {
                    raw += strength[j];
                }
            }
            fitness[i] = raw;
        }

        // Step 3: density D(i) = 1 / (sigma_k + 2)
        int k = std::min(n - 1, 1);
        for (int i = 0; i < n; ++i) {
            std::vector<double> distances;
            distances.reserve(n - 1);
            for (int j = 0; j < n; ++j) {
                if (i == j)
                    continue;
                double dist = 0.0;
                for (int o = 0; o < pop.n_obj; ++o) {
                    double diff = pop.objective(i, o) - pop.objective(j, o);
                    dist += diff * diff;
                }
                distances.push_back(std::sqrt(dist));
            }
            std::sort(distances.begin(), distances.end());
            double sigma_k = distances[k];
            fitness[i] += 1.0 / (sigma_k + 2.0);
        }
    }
};

} // namespace ea
