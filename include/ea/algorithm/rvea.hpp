#pragma once
/// @file rvea.hpp
/// @brief RVEA (Reference Vector-guided Evolutionary Algorithm)
///
/// Reference: R. Cheng et al. "A Reference Vector Guided Evolutionary Algorithm
/// for Many-Objective Optimization" IEEE Trans. Evol. Comput., 2016.
///
/// Key components:
/// - Reference vectors uniformly distributed on unit simplex
/// - Angle-penalized distance (APD) for selection
/// - Adaptive reference vector adjustment

#include <ea/core/population.hpp>
#include <ea/core/comparator.hpp>
#include <ea/util/random.hpp>
#include <ea/util/reference_point.hpp>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <unordered_map>

namespace ea {

/// RVEA — Reference Vector-guided Evolutionary Algorithm.
/// Designed for many-objective optimization (3+ objectives).
///
/// Usage:
///   RVEA<ea::SBXCrossover, ea::PolynomialMutation> rvea;
///   rvea.pop_size = 100;
///   rvea.max_evals = 25000;
///   rvea.run(pop, problem);
template<typename CX, typename MT>
struct RVEA {
    CX crossover;
    MT mutation;

    int pop_size = 100;
    int max_evals = 25000;
    int num_reference_points = 100;  ///< Number of reference vectors
    double alpha = 2.0;              ///< APD penalty parameter
    int freq_adapt = 10;             ///< Frequency of reference vector adaptation (generations)

    static constexpr std::string_view name() { return "RVEA"; }

    void run(this auto& self, Population& pop, auto&& problem) {
        const int n = self.pop_size;
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;

        if (pop.pop_size != n) pop.resize(n);

        // === Evaluate initial population ===
        for (int i = 0; i < n; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }
        int evals = n;

        // === Generate reference vectors ===
        std::vector<ReferencePoint> ref_points;
        if (n_obj == 2) {
            // For 2 objectives, generate uniformly on the line
            for (int i = 0; i < self.num_reference_points; ++i) {
                ReferencePoint rp(2);
                double w0 = static_cast<double>(i) / (self.num_reference_points - 1);
                rp.position[0] = w0;
                rp.position[1] = 1.0 - w0;
                ref_points.push_back(rp);
            }
        } else {
            // For 3+ objectives, use Das-Dennis
            int divisions = static_cast<int>(std::cbrt(self.num_reference_points)) + 1;
            generate_reference_points_das_dennis(ref_points, n_obj, divisions);
        }

        // Normalize reference vectors
        for (auto& rp : ref_points) {
            double norm = 0.0;
            for (double v : rp.position) norm += v * v;
            norm = std::sqrt(norm);
            if (norm > 1e-12) {
                for (double& v : rp.position) v /= norm;
            }
        }

        // === Evolution loop ===
        auto& rng = Random::instance();
        int generation = 0;

        while (evals < self.max_evals) {
            generation++;

            // 1. Generate offspring via crossover + mutation
            Population offspring(n, dim, n_obj);
            offspring.lower_bounds = pop.lower_bounds;
            offspring.upper_bounds = pop.upper_bounds;

            for (int i = 0; i < n; i += 2) {
                int p1 = rng.uniform_int(0, n - 1);
                int p2 = rng.uniform_int(0, n - 1);

                // Copy parents
                for (int j = 0; j < dim; ++j) {
                    offspring.gene(i, j) = pop.gene(p1, j);
                    if (i + 1 < n) offspring.gene(i + 1, j) = pop.gene(p2, j);
                }

                // Crossover
                self.crossover.apply(offspring, i, (i + 1 < n) ? i + 1 : i, i);

                // Mutate
                self.mutation.apply(offspring, i);
                if (i + 1 < n) self.mutation.apply(offspring, i + 1);
            }

            // Evaluate offspring
            for (int i = 0; i < n; ++i) {
                if (!offspring.evaluated(i)) {
                    problem(offspring, i);
                    offspring.set_evaluated(i, true);
                    evals++;
                    if (evals >= self.max_evals) break;
                }
            }

            // 2. Combine parent + offspring
            Population combined(2 * n, dim, n_obj);
            combined.lower_bounds = pop.lower_bounds;
            combined.upper_bounds = pop.upper_bounds;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < dim; ++j) {
                    combined.gene(i, j) = pop.gene(i, j);
                    combined.gene(n + i, j) = offspring.gene(i, j);
                }
                for (int o = 0; o < n_obj; ++o) {
                    combined.objective(i, o) = pop.objective(i, o);
                    combined.objective(n + i, o) = offspring.objective(i, o);
                }
            }
            combined.pop_size = 2 * n;

            // 3. Find ideal and nadir points
            std::vector<double> ideal(n_obj, std::numeric_limits<double>::infinity());
            std::vector<double> nadir(n_obj, std::numeric_limits<double>::lowest());
            for (int i = 0; i < 2 * n; ++i) {
                for (int o = 0; o < n_obj; ++o) {
                    double val = combined.objective(i, o);
                    if (val < ideal[o]) ideal[o] = val;
                    if (val > nadir[o]) nadir[o] = val;
                }
            }

            // 4. Normalize objectives
            for (int i = 0; i < 2 * n; ++i) {
                for (int o = 0; o < n_obj; ++o) {
                    double range = nadir[o] - ideal[o];
                    if (range > 1e-12) {
                        combined.objective(i, o) = (combined.objective(i, o) - ideal[o]) / range;
                    } else {
                        combined.objective(i, o) = 0.0;
                    }
                }
            }

            // 5. Assign each individual to nearest reference vector
            std::vector<int> assigned_ref(2 * n, -1);
            std::vector<double> min_angles(2 * n, std::numeric_limits<double>::infinity());

            for (int i = 0; i < 2 * n; ++i) {
                double best_angle = std::numeric_limits<double>::infinity();
                int best_ref = 0;

                for (int r = 0; r < static_cast<int>(ref_points.size()); ++r) {
                    double dot = 0.0;
                    double norm_f = 0.0;
                    for (int o = 0; o < n_obj; ++o) {
                        dot += combined.objective(i, o) * ref_points[r].position[o];
                        norm_f += combined.objective(i, o) * combined.objective(i, o);
                    }
                    norm_f = std::sqrt(norm_f);
                    double norm_ref = 1.0; // Already normalized

                    double cos_angle = 0.0;
                    if (norm_f > 1e-12) {
                        cos_angle = dot / (norm_f * norm_ref);
                    }
                    // Clamp to [-1, 1] for acos
                    cos_angle = std::clamp(cos_angle, -1.0, 1.0);
                    double angle = std::acos(cos_angle);

                    if (angle < best_angle) {
                        best_angle = angle;
                        best_ref = r;
                    }
                }

                assigned_ref[i] = best_ref;
                min_angles[i] = best_angle;
            }

            // 6. Compute APD (Angle-Penalized Distance) for selection
            std::vector<double> apd_values(2 * n);
            double t = static_cast<double>(evals) / self.max_evals;

            for (int i = 0; i < 2 * n; ++i) {
                double norm_f = 0.0;
                for (int o = 0; o < n_obj; ++o) {
                    norm_f += combined.objective(i, o) * combined.objective(i, o);
                }
                norm_f = std::sqrt(norm_f);

                // APD = norm_f * (1 + self.alpha * t * min_angles[i])
                apd_values[i] = norm_f * (1.0 + self.alpha * t * min_angles[i]);
            }

            // 7. Environmental selection: pick best n by APD
            std::vector<int> selected_indices;
            selected_indices.reserve(n);

            // Count how many per reference vector
            std::unordered_map<int, std::vector<std::pair<double, int>>> ref_members;
            for (int i = 0; i < 2 * n; ++i) {
                ref_members[assigned_ref[i]].push_back({apd_values[i], i});
            }

            // Select best from each reference vector
            std::vector<bool> selected(2 * n, false);
            for (auto& [ref_idx, members] : ref_members) {
                std::sort(members.begin(), members.end());
                for (const auto& [apd, idx] : members) {
                    if (static_cast<int>(selected_indices.size()) < n) {
                        selected_indices.push_back(idx);
                        selected[idx] = true;
                    }
                }
            }

            // Fill remaining slots if any (shouldn't happen with enough ref vectors)
            if (static_cast<int>(selected_indices.size()) < n) {
                for (int i = 0; i < 2 * n && static_cast<int>(selected_indices.size()) < n; ++i) {
                    if (!selected[i]) {
                        selected_indices.push_back(i);
                    }
                }
            }

            // Copy selected back to population
            for (int i = 0; i < n; ++i) {
                int src = selected_indices[i];
                for (int j = 0; j < dim; ++j) {
                    pop.gene(i, j) = combined.gene(src, j);
                }
                for (int o = 0; o < n_obj; ++o) {
                    pop.objective(i, o) = combined.objective(src, o);
                }
                pop.set_evaluated(i, true);
            }

            // 8. Adaptive reference vector adjustment (every freq_adapt generations)
            if (generation % self.freq_adapt == 0) {
                // Scale reference vectors based on objective ranges
                for (auto& rp : ref_points) {
                    for (int o = 0; o < n_obj; ++o) {
                        rp.position[o] *= (nadir[o] - ideal[o]);
                    }
                    // Renormalize
                    double norm = 0.0;
                    for (double v : rp.position) norm += v * v;
                    norm = std::sqrt(norm);
                    if (norm > 1e-12) {
                        for (double& v : rp.position) v /= norm;
                    }
                }
            }
        }
    }
};

} // namespace ea
