#pragma once
/// @file nsga3_replacement.hpp
/// @brief NSGA-III Environmental Selection — reference-point-based niching.
///
/// Implements Algorithm 4 from the NSGA-III paper:
/// 1. Normalize objectives using ideal point and intercepts
/// 2. Associate individuals with reference points
/// 3. Niching-based selection from critical front
///
/// Reference: K. Deb and H. Jain, "An Evolutionary Many-Objective Optimization
/// Algorithm Using Reference-Point-Based Nondominated Sorting Approach, Part I:
/// Solving Problems With Box Constraints," IEEE TEVC, vol. 18, no. 4, 2014.
///
/// jMetal reference classes:
///   - org.uma.jmetal.algorithm.multiobjective.nsgaiii.util.EnvironmentalSelection
///   - org.uma.jmetal.component.catalogue.ea.replacement.impl.NSGAIIIReplacement

#include <ea/core/population.hpp>
#include <ea/core/comparator.hpp>
#include <ea/util/reference_point.hpp>
#include <ea/util/random.hpp>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <random>

namespace ea {

/// NSGA-III environmental selection (replacement operator).
/// Takes a combined population of 2*pop_size individuals,
/// performs fast non-dominated sort, then uses reference-point-based
/// niching to select the best pop_size individuals.
///
/// jMetal equivalent: EnvironmentalSelection + NSGAIIIReplacement
struct NSGAIIIReplacement {
    /// Reference points for niching
    std::vector<ReferencePoint> reference_points;

    /// Number of objectives
    int n_obj = 0;

    /// Target population size
    int pop_size = 100;

    // --- Temporary storage (reused across calls) ---
    std::vector<double> ideal_point_;
    std::vector<double> intercepts_;
    std::vector<double> translated_;        // Flat [n * n_obj] for translated objectives
    std::vector<double> normalized_;        // Flat [n * n_obj] for normalized objectives
    std::vector<std::vector<double>> extreme_points_;  // For hyperplane construction
    std::vector<int> niche_counts_;        // Per-reference-point member counts

    NSGAIIIReplacement() = default;

    explicit NSGAIIIReplacement(std::vector<ReferencePoint> ref_points)
        : reference_points(std::move(ref_points)) {}

    /// Set reference points and update internal state
    void set_reference_points(this auto& self, std::vector<ReferencePoint> ref_points) {
        self.reference_points = std::move(ref_points);
    }

    /// Initialize dimensions and allocate temporary storage
    void initialize(this auto& self, int num_objectives, int target_pop_size) {
        self.n_obj = num_objectives;
        self.pop_size = target_pop_size;
        self.ideal_point_.resize(num_objectives);
        self.intercepts_.resize(num_objectives);
        self.extreme_points_.resize(num_objectives, std::vector<double>(num_objectives));
    }

    /// Perform NSGA-III environmental selection on combined population.
    /// @param combined Combined population of size 2*pop_size
    /// @param target_pop_size Target population size
    /// @return Vector of selected individual indices
    std::vector<int> replace(this auto& self, Population& combined, int target_pop_size) {
        // Step 1: Fast non-dominated sort
        auto fronts = fast_non_dominated_sort(combined);

        std::vector<int> selected;
        selected.reserve(target_pop_size);

        int last_complete_front = 0;

        // Add complete fronts
        for (; last_complete_front < static_cast<int>(fronts.size()); ++last_complete_front) {
            if (static_cast<int>(selected.size() + fronts[last_complete_front].size()) <= target_pop_size) {
                selected.insert(selected.end(), fronts[last_complete_front].begin(), fronts[last_complete_front].end());
            } else {
                break;
            }
        }

        if (static_cast<int>(selected.size()) == target_pop_size) {
            return selected;
        }

        // We need to select from the critical front using reference points
        int remaining = target_pop_size - static_cast<int>(selected.size());

        // Build list of all fronts (up to and including critical front)
        std::vector<std::vector<int>> all_fronts;
        for (int f = 0; f < last_complete_front; ++f) {
            all_fronts.push_back(fronts[f]);
        }
        if (last_complete_front < static_cast<int>(fronts.size())) {
            all_fronts.push_back(fronts[last_complete_front]);
        }

        self.initialize(combined.n_obj, target_pop_size);

        // Step 2: Normalize objectives (Algorithm 2 from paper)
        self.normalize_objectives(combined, all_fronts);

        // Step 3: Associate and niching (Algorithm 3 + 4 from paper)
        auto niching_selected = self.niching_selection(combined, all_fronts, remaining);
        selected.insert(selected.end(), niching_selected.begin(), niching_selected.end());

        return selected;
    }

    /// Compact selected individuals into the first pop_size positions.
    void compact(this auto& self, Population& combined, int pop_size) {
        auto selected = self.replace(combined, pop_size);

        for (int i = 0; i < pop_size; ++i) {
            if (selected[i] != i) {
                combined.copy_individual(selected[i], i);
            }
        }
    }

private:
    // ============================================================
    // Objective Normalization (Algorithm 2)
    // ============================================================

    /// Compute ideal point from the first front
    void compute_ideal_point(this auto& self, Population& pop, const std::vector<int>& first_front) {
        for (int f = 0; f < self.n_obj; ++f) {
            self.ideal_point_[f] = std::numeric_limits<double>::max();
        }

        for (int idx : first_front) {
            for (int f = 0; f < self.n_obj; ++f) {
                self.ideal_point_[f] = std::min(self.ideal_point_[f], pop.objective(idx, f));
            }
        }
    }

    /// Translate objectives by subtracting ideal point.
    /// Stores translated objectives for all individuals in all fronts.
    void translate_objectives(this auto& self, Population& pop,
                              const std::vector<std::vector<int>>& fronts) {
        int total_individuals = 0;
        for (const auto& front : fronts) {
            total_individuals += static_cast<int>(front.size());
        }
        self.translated_.resize(total_individuals * self.n_obj);

        int offset = 0;
        for (const auto& front : fronts) {
            for (int idx : front) {
                for (int f = 0; f < self.n_obj; ++f) {
                    self.translated_[offset * self.n_obj + f] = pop.objective(idx, f) - self.ideal_point_[f];
                }
                ++offset;
            }
        }
    }

    /// Achievement Scalarizing Function (ASF) for finding extreme points
    double asf(const std::vector<double>& translated, int target_idx) const {
        double max_ratio = std::numeric_limits<double>::lowest();
        for (int i = 0; i < self.n_obj; ++i) {
            double weight = (i == target_idx) ? 1.0 : 1e-6;
            max_ratio = std::max(max_ratio, translated[i] / weight);
        }
        return max_ratio;
    }

    /// Find extreme points (one per objective) from the first front
    /// Uses ASF to find individuals that maximize each objective axis
    void find_extreme_points(this auto& self, Population& pop,
                           const std::vector<int>& first_front) {
        // Get translated objectives for first front individuals
        std::vector<std::vector<double>> translated_front;
        translated_front.reserve(first_front.size());

        for (int idx : first_front) {
            std::vector<double> translated(self.n_obj);
            for (int f = 0; f < self.n_obj; ++f) {
                translated[f] = pop.objective(idx, f) - self.ideal_point_[f];
            }
            translated_front.push_back(std::move(translated));
        }

        for (int f = 0; f < self.n_obj; ++f) {
            double min_asf = std::numeric_limits<double>::max();
            int extreme_idx = first_front[0];

            for (size_t i = 0; i < translated_front.size(); ++i) {
                double val = asf(translated_front[i], f);
                if (val < min_asf) {
                    min_asf = val;
                    extreme_idx = first_front[i];
                }
            }

            // Store extreme point coordinates (translated)
            for (int j = 0; j < self.n_obj; ++j) {
                self.extreme_points_[f][j] = pop.objective(extreme_idx, j) - self.ideal_point_[j];
            }
        }
    }

    /// Construct hyperplane from extreme points to find intercepts.
    /// Uses Gaussian elimination to solve Ax = 1 where A is extreme points matrix.
    /// Falls back to extreme point values if hyperplane is degenerate.
    void construct_hyperplane(this auto& self) {
        // Check for duplicate extreme points (degenerate case)
        bool degenerate = false;
        for (int i = 0; i < self.n_obj && !degenerate; ++i) {
            for (int j = i + 1; j < self.n_obj; ++j) {
                bool same = true;
                for (int k = 0; k < self.n_obj; ++k) {
                    if (std::abs(self.extreme_points_[i][k] - self.extreme_points_[j][k]) > 1e-10) {
                        same = false;
                        break;
                    }
                }
                if (same) {
                    degenerate = true;
                    break;
                }
            }
        }

        if (degenerate) {
            // Fallback: use extreme point values as intercepts
            for (int f = 0; f < self.n_obj; ++f) {
                self.intercepts_[f] = self.extreme_points_[f][f];
                if (self.intercepts_[f] < 1e-10) self.intercepts_[f] = 1e-10;
            }
            return;
        }

        // Gaussian elimination to solve A * x = b where b = [1, 1, ..., 1]
        // A is the matrix of extreme points (translated)
        std::vector<std::vector<double>> a = self.extreme_points_;
        std::vector<double> b(self.n_obj, 1.0);

        // Forward elimination (convert to upper triangular)
        for (int base = 0; base < self.n_obj - 1; ++base) {
            for (int target = base + 1; target < self.n_obj; ++target) {
                double ratio = a[target][base] / a[base][base];
                for (int term = base; term < self.n_obj; ++term) {
                    a[target][term] -= a[base][term] * ratio;
                }
                b[target] -= b[base] * ratio;
            }
        }

        // Back substitution
        std::vector<double> x(self.n_obj, 0.0);
        for (int i = self.n_obj - 1; i >= 0; --i) {
            double known_sum = 0.0;
            for (int j = i + 1; j < self.n_obj; ++j) {
                known_sum += a[i][j] * x[j];
            }
            x[i] = (b[i] - known_sum) / a[i][i];
        }

        // Compute intercepts as 1 / x[i]
        for (int f = 0; f < self.n_obj; ++f) {
            if (std::abs(x[f]) > 1e-10) {
                self.intercepts_[f] = 1.0 / x[f];
            } else {
                self.intercepts_[f] = self.extreme_points_[f][f];
            }
            if (self.intercepts_[f] < 1e-10) self.intercepts_[f] = 1e-10;
        }
    }

    /// Normalize all objectives by subtracting ideal point and dividing by intercepts.
    void normalize_objectives(this auto& self, Population& pop,
                              const std::vector<std::vector<int>>& fronts) {
        // Step 1: Compute ideal point
        self.compute_ideal_point(pop, fronts[0]);

        // Step 2: Translate objectives
        self.translate_objectives(pop, fronts);

        // Step 3: Find extreme points and construct hyperplane
        self.find_extreme_points(pop, fronts[0]);
        self.construct_hyperplane();

        // Step 4: Normalize
        int total_individuals = 0;
        for (const auto& front : fronts) {
            total_individuals += static_cast<int>(front.size());
        }
        self.normalized_.resize(total_individuals * self.n_obj);

        int offset = 0;
        for (const auto& front : fronts) {
            for (int idx : front) {
                for (int f = 0; f < self.n_obj; ++f) {
                    double translated = pop.objective(idx, f) - self.ideal_point_[f];
                    double range = self.intercepts_[f] - self.ideal_point_[f];
                    if (std::abs(range) > 1e-10) {
                        self.normalized_[offset * self.n_obj + f] = translated / range;
                    } else {
                        self.normalized_[offset * self.n_obj + f] = translated / 1e-10;
                    }
                }
                ++offset;
            }
        }
    }

    // ============================================================
    // Reference Point Association and Niching
    // ============================================================

    /// Compute perpendicular distance from a point to a reference direction.
    /// This is the distance from the point to the line through the origin
    /// in the direction of the reference point.
    double perpendicular_distance(const std::vector<double>& direction,
                                   const double* point) const {
        double numerator = 0.0;
        double denominator = 0.0;

        for (int i = 0; i < self.n_obj; ++i) {
            numerator += direction[i] * point[i];
            denominator += direction[i] * direction[i];
        }

        if (denominator < 1e-10) {
            return std::numeric_limits<double>::max();
        }

        double k = numerator / denominator;

        double dist = 0.0;
        for (int i = 0; i < self.n_obj; ++i) {
            double diff = k * direction[i] - point[i];
            dist += diff * diff;
        }

        return std::sqrt(dist);
    }

    /// Find the closest reference point for a normalized objective vector.
    int find_closest_reference_point(const double* normalized) const {
        int closest = 0;
        double min_dist = std::numeric_limits<double>::max();

        for (size_t r = 0; r < self.reference_points.size(); ++r) {
            double dist = self.perpendicular_distance(self.reference_points[r].position, normalized);
            if (dist < min_dist) {
                min_dist = dist;
                closest = static_cast<int>(r);
            }
        }

        return closest;
    }

    /// Associate individuals with reference points and build niche counts.
    /// @return Mapping: reference_point_idx -> list of (individual_idx, distance) pairs for critical front
    void associate_and_build_niche_counts(this auto& self, Population& pop,
                                         const std::vector<std::vector<int>>& fronts) {
        // Clear reference point state
        for (auto& rp : self.reference_points) {
            rp.clear();
        }

        self.niche_counts_.assign(self.reference_points.size(), 0);

        int offset = 0;
        for (size_t f = 0; f < fronts.size(); ++f) {
            for (int idx : fronts[f]) {
                int closest_rp = self.find_closest_reference_point(
                    &self.normalized_[offset * self.n_obj]);

                if (f + 1 < fronts.size()) {
                    // Individual on a complete front: just increment niche count
                    self.niche_counts_[closest_rp]++;
                } else {
                    // Individual on critical front: add as potential member
                    double dist = self.perpendicular_distance(
                        self.reference_points[closest_rp].position,
                        &self.normalized_[offset * self.n_obj]);
                    self.reference_points[closest_rp].add_potential_member(idx, dist);
                }
                ++offset;
            }
        }

        // Copy niche counts to reference points and sort potential members
        for (size_t r = 0; r < self.reference_points.size(); ++r) {
            self.reference_points[r].member_count = self.niche_counts_[r];
            self.reference_points[r].sort_potential_members();
        }
    }

    /// Niching-based selection from the critical front (Algorithm 4).
    std::vector<int> niching_selection(this auto& self, Population& pop,
                                        const std::vector<std::vector<int>>& fronts,
                                        int remaining) {
        // Associate all individuals with reference points
        self.associate_and_build_niche_counts(pop, fronts);

        std::vector<int> selected;
        selected.reserve(remaining);

        auto& rng = Random::instance().engine();

        // Group reference points by niche count
        // Map: niche_count -> list of reference point indices
        std::map<int, std::vector<int>> niche_groups;
        for (size_t i = 0; i < self.reference_points.size(); ++i) {
            if (self.reference_points[i].has_potential_members()) {
                niche_groups[self.niche_counts_[i]].push_back(static_cast<int>(i));
            }
        }

        while (static_cast<int>(selected.size()) < remaining && !niche_groups.empty()) {
            // Find reference points with minimum niche count
            auto it = niche_groups.begin();
            auto& min_niche_refs = it->second;

            // Randomly select one reference point with minimum niche count
            std::uniform_int_distribution<size_t> dist(0, min_niche_refs.size() - 1);
            size_t rp_idx_in_group = dist(rng);
            int ref_point_idx = min_niche_refs[rp_idx_in_group];

            auto& rp = self.reference_points[ref_point_idx];

            if (rp.has_potential_members()) {
                int chosen_idx;
                if (self.niche_counts_[ref_point_idx] == 0) {
                    // Select closest member if no members yet
                    chosen_idx = rp.remove_closest_member();
                } else {
                    // Select random member otherwise
                    chosen_idx = rp.remove_random_member(rng);
                }

                if (chosen_idx >= 0) {
                    selected.push_back(chosen_idx);

                    // Update niche count
                    int old_count = self.niche_counts_[ref_point_idx];
                    self.niche_counts_[ref_point_idx]++;
                    int new_count = self.niche_counts_[ref_point_idx];

                    // Move reference point to new group
                    // Remove from old group
                    min_niche_refs.erase(min_niche_refs.begin() + rp_idx_in_group);
                    if (min_niche_refs.empty()) {
                        niche_groups.erase(it);
                    }

                    // Add to new group if still has potential members
                    if (rp.has_potential_members()) {
                        niche_groups[new_count].push_back(ref_point_idx);
                    }
                }
            } else {
                // No more potential members, remove from group
                min_niche_refs.erase(min_niche_refs.begin() + rp_idx_in_group);
                if (min_niche_refs.empty()) {
                    niche_groups.erase(it);
                }
            }
        }

        return selected;
    }
};

} // namespace ea
