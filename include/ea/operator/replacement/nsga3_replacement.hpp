#pragma once
/// @file nsga3_replacement.hpp
/// @brief NSGA-III Environmental Selection — reference-point-based niching.
///
/// Implements Algorithm 4 from the NSGA-III paper:
/// 1. Normalize objectives using ideal point and intercepts
/// 2. Associate individuals with reference points
/// 3. Niching-based selection from critical front
///
/// Reference: Deb & Jain, IEEE TEVC 2014.

#include <algorithm>
#include <cmath>
#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <ea/util/reference_point.hpp>
#include <limits>
#include <vector>

namespace ea {

/// NSGA-III environmental selection (replacement operator).
/// Uses reference-point-based niching for many-objective selection.
struct NSGAIIIReplacement {
    std::vector<ReferencePoint> reference_points;
    int n_obj = 0;
    int pop_size = 100;

    NSGAIIIReplacement() = default;
    explicit NSGAIIIReplacement(std::vector<ReferencePoint> refs)
        : reference_points(std::move(refs)) {}

    /// Perform NSGA-III environmental selection on a combined population.
    std::vector<int> replace(this NSGAIIIReplacement& self, Population<>& combined, int target_size);

private:
    /// Compute ideal point from front
    static void compute_ideal_point(Population<>& pop, const std::vector<int>& front,
                                    std::vector<double>& ideal, int n_obj) {
        ideal.assign(n_obj, std::numeric_limits<double>::max());
        for (int idx : front) {
            for (int o = 0; o < n_obj; ++o) {
                ideal[o] = std::min(ideal[o], pop.objective(idx, o));
            }
        }
    }

    /// Compute ASF (Achievement Scalarizing Function) for extreme point finding
    static double asf(const std::vector<double>& objectives, int axis, int n_obj) {
        double max_val = -std::numeric_limits<double>::max();
        for (int i = 0; i < n_obj; ++i) {
            double weight = (i == axis) ? 1.0 : 1e-6;
            double val = objectives[i] / weight;
            max_val = std::max(max_val, val);
        }
        return max_val;
    }

    /// Compute perpendicular distance from a point to a reference direction
    static double perpendicular_distance(const double* point, const double* direction, int n_obj) {
        // Project point onto direction
        double dot = 0.0, dir_norm_sq = 0.0;
        for (int i = 0; i < n_obj; ++i) {
            dot += point[i] * direction[i];
            dir_norm_sq += direction[i] * direction[i];
        }

        // Perpendicular component
        double dist_sq = 0.0;
        for (int i = 0; i < n_obj; ++i) {
            double proj = (dir_norm_sq > 1e-14) ? (dot / dir_norm_sq) * direction[i] : 0.0;
            double diff = point[i] - proj;
            dist_sq += diff * diff;
        }
        return std::sqrt(dist_sq);
    }

    /// Find which reference point is closest to a normalized point
    static int find_closest_ref(const double* normalized, const std::vector<ReferencePoint>& refs,
                                int n_obj) {
        int closest = 0;
        double min_dist = std::numeric_limits<double>::max();
        for (size_t r = 0; r < refs.size(); ++r) {
            double dist = 0.0;
            for (int o = 0; o < n_obj; ++o) {
                double diff = normalized[o] - refs[r].position[o];
                dist += diff * diff;
            }
            if (dist < min_dist) {
                min_dist = dist;
                closest = static_cast<int>(r);
            }
        }
        return closest;
    }
};

inline std::vector<int> NSGAIIIReplacement::replace(this NSGAIIIReplacement& self,
                                                    Population<>& combined, int target_size) {
    const int n_obj = combined.n_obj;

    // Non-dominated sort
    auto fronts = fast_non_dominated_sort(combined);

    // Step 1: Fill with complete fronts until one doesn't fit
    std::vector<int> selected;
    int last_front_idx = -1;

    for (int f = 0; f < static_cast<int>(fronts.size()); ++f) {
        if (static_cast<int>(selected.size() + fronts[f].size()) <= target_size) {
            selected.insert(selected.end(), fronts[f].begin(), fronts[f].end());
        } else {
            last_front_idx = f;
            break;
        }
    }

    // If we already have enough, return
    if (static_cast<int>(selected.size()) >= target_size || last_front_idx < 0) {
        selected.resize(target_size);
        return selected;
    }

    // Step 2: Niching selection for the last front
    auto& last_front = fronts[last_front_idx];
    int remaining = target_size - static_cast<int>(selected.size());

    // Compute ideal point from all selected + last front
    std::vector<double> ideal(n_obj, std::numeric_limits<double>::max());
    for (int idx : selected) {
        for (int o = 0; o < n_obj; ++o) {
            ideal[o] = std::min(ideal[o], combined.objective(idx, o));
        }
    }
    for (int idx : last_front) {
        for (int o = 0; o < n_obj; ++o) {
            ideal[o] = std::min(ideal[o], combined.objective(idx, o));
        }
    }

    // Normalize objectives (subtract ideal point)
    // Then associate each individual with nearest reference point
    std::vector<int> association(last_front.size());
    std::vector<int> niche_counts(self.reference_points.size(), 0);

    // Count niche members from already-selected individuals
    for (int idx : selected) {
        std::vector<double> normalized(n_obj);
        for (int o = 0; o < n_obj; ++o) {
            normalized[o] = combined.objective(idx, o) - ideal[o];
        }
        int ref = find_closest_ref(normalized.data(), self.reference_points, n_obj);
        niche_counts[ref]++;
    }

    // Associate last front individuals
    for (size_t i = 0; i < last_front.size(); ++i) {
        std::vector<double> normalized(n_obj);
        for (int o = 0; o < n_obj; ++o) {
            normalized[o] = combined.objective(last_front[i], o) - ideal[o];
        }
        association[i] = find_closest_ref(normalized.data(), self.reference_points, n_obj);
    }

    // Niching: pick from last front, preferring reference points with fewer members
    std::vector<int> selected_from_front;
    std::vector<bool> chosen(last_front.size(), false);

    for (int pick = 0; pick < remaining; ++pick) {
        // Find reference point with minimum niche count
        int min_count = std::numeric_limits<int>::max();
        std::vector<int> min_refs;

        for (size_t r = 0; r < self.reference_points.size(); ++r) {
            if (niche_counts[r] < min_count) {
                min_count = niche_counts[r];
                min_refs.clear();
                min_refs.push_back(static_cast<int>(r));
            } else if (niche_counts[r] == min_count) {
                min_refs.push_back(static_cast<int>(r));
            }
        }

        // Pick a random ref from minimum-niche refs
        auto& rng = Random::instance();
        int chosen_ref = min_refs[rng.uniform_int(0, static_cast<int>(min_refs.size()) - 1)];

        // Find closest individual from last front associated with this ref
        int best_idx = -1;
        double best_dist = std::numeric_limits<double>::max();
        for (size_t i = 0; i < last_front.size(); ++i) {
            if (chosen[i])
                continue;
            if (association[i] != chosen_ref)
                continue;

            std::vector<double> normalized(n_obj);
            for (int o = 0; o < n_obj; ++o) {
                normalized[o] = combined.objective(last_front[i], o) - ideal[o];
            }
            double dist = perpendicular_distance(
                normalized.data(), self.reference_points[chosen_ref].position.data(), n_obj);

            if (dist < best_dist) {
                best_dist = dist;
                best_idx = static_cast<int>(i);
            }
        }

        if (best_idx >= 0) {
            selected_from_front.push_back(last_front[best_idx]);
            chosen[best_idx] = true;
            niche_counts[chosen_ref]++;
        } else {
            // No individual for this ref — pick any remaining
            for (size_t i = 0; i < last_front.size(); ++i) {
                if (!chosen[i]) {
                    selected_from_front.push_back(last_front[i]);
                    chosen[i] = true;
                    niche_counts[association[i]]++;
                    break;
                }
            }
        }
    }

    selected.insert(selected.end(), selected_from_front.begin(), selected_from_front.end());
    return selected;
}

} // namespace ea