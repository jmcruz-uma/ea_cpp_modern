#pragma once
/// @file reference_point.hpp
/// @brief Reference point generation and management for NSGA-III.
/// Uses the Das-Dennis method for generating well-spread reference directions.
///
/// Reference: K. Deb and H. Jain, "An Evolutionary Many-Objective Optimization
/// Algorithm Using Reference-Point-Based Nondominated Sorting Approach, Part I:
/// Solving Problems With Box Constraints," IEEE TEVC, vol. 18, no. 4, 2014.
///
/// jMetal reference:
///   - org.uma.jmetal.algorithm.multiobjective.nsgaiii.util.ReferencePoint
///   - org.uma.jmetal.component.catalogue.ea.replacement.impl.nsgaiii.ReferencePoint

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

namespace ea {

/// A reference point (direction) in objective space for NSGA-III.
/// Tracks associated individuals and niche counts for environmental selection.
///
/// jMetal equivalent: ReferencePoint.java (both old and component versions)
struct ReferencePoint {
    /// Position of the reference point in objective space (normalized, sums to 1)
    std::vector<double> position;

    /// Number of currently associated selected individuals
    int member_count = 0;

    /// Indices of potential members (individuals on the critical front) with distances
    struct MemberInfo {
        int individual_index = -1; ///< Index in combined population
        double distance = 0.0;     ///< Perpendicular distance to this reference point

        bool operator<(const MemberInfo& other) const {
            return distance < other.distance; // For sorting by ascending distance
        }
    };
    std::vector<MemberInfo> potential_members;

    ReferencePoint() = default;

    explicit ReferencePoint(int num_objectives)
        : position(num_objectives, 0.0) {}

    /// Copy constructor (resets member tracking)
    ReferencePoint(const ReferencePoint& other)
        : position(other.position)
        , member_count(0) {}

    /// Add a selected member to this reference point
    void add_member(this auto& self) { self.member_count++; }

    /// Clear all member associations (reset niche count and potential members)
    void clear(this auto& self) {
        self.member_count = 0;
        self.potential_members.clear();
    }

    /// Add a potential member from the critical front
    void add_potential_member(this auto& self, int idx, double distance) {
        self.potential_members.push_back({idx, distance});
    }

    /// Sort potential members by distance (ascending — closest first)
    void sort_potential_members(this auto& self) {
        std::sort(self.potential_members.begin(), self.potential_members.end());
    }

    /// Check if there are potential members available
    bool has_potential_members(this auto& self) { return !self.potential_members.empty(); }

    /// Remove and return the closest potential member
    int remove_closest_member(this auto& self) {
        if (self.potential_members.empty())
            return -1;
        int idx = self.potential_members.front().individual_index;
        self.potential_members.erase(self.potential_members.begin());
        return idx;
    }

    /// Remove a random potential member, return its individual index
    template <typename Rng> int remove_random_member(this auto& self, Rng& rng) {
        if (self.potential_members.empty())
            return -1;
        std::uniform_int_distribution<size_t> dist(0, self.potential_members.size() - 1);
        size_t i = dist(rng);
        int idx = self.potential_members[i].individual_index;
        self.potential_members.erase(self.potential_members.begin() + i);
        return idx;
    }
};

/// Generate reference points using the Das-Dennis method.
/// Creates uniformly distributed reference points on the normalized hyperplane
/// (sum of coordinates = 1) using simplex-lattice design.
///
/// @param[out] reference_points  Output vector of generated reference points
/// @param[in]  num_objectives   Number of objectives (dimension of the simplex)
/// @param[in]  num_divisions    Number of divisions per dimension (p in the paper)
///
/// jMetal equivalent: ReferencePoint.generateReferencePoints()
/// Recursive implementation following the original NSGA-III paper.
inline void generate_reference_points_das_dennis(std::vector<ReferencePoint>& reference_points,
                                                 int num_objectives, int num_divisions) {

    ReferencePoint point(num_objectives);

    // Recursive lambda to enumerate all combinations
    auto generate_recursive = [&](this auto& self, std::vector<ReferencePoint>& out,
                                  ReferencePoint& ref, int objectives_left, int divisions_left,
                                  int element) -> void {
        if (element == num_objectives - 1) {
            ref.position[element] = static_cast<double>(divisions_left) / num_divisions;
            out.emplace_back(ref);
        } else {
            for (int i = 0; i <= divisions_left; ++i) {
                ref.position[element] = static_cast<double>(i) / num_divisions;
                self(out, ref, objectives_left, divisions_left - i, element + 1);
            }
        }
    };

    generate_recursive(reference_points, point, num_objectives, num_divisions, 0);
}

/// Generate two-layer reference points (outer + inner layer).
/// Recommended for problems with many objectives (>5) as per NSGA-III paper.
/// Inner layer is shrunk towards the center (0.5, ..., 0.5) and renormalized.
///
/// @param[out] reference_points     Output vector
/// @param[in]  num_objectives      Number of objectives
/// @param[in]  outer_divisions     Divisions for outer layer
/// @param[in]  inner_divisions     Divisions for inner layer (0 = single layer)
///
/// jMetal equivalent: NSGAIII constructor with secondLayerDivisions
inline void generate_reference_points_two_layer(std::vector<ReferencePoint>& reference_points,
                                                int num_objectives, int outer_divisions,
                                                int inner_divisions) {

    // Generate outer layer
    generate_reference_points_das_dennis(reference_points, num_objectives, outer_divisions);

    if (inner_divisions > 0) {
        // Generate inner layer
        std::vector<ReferencePoint> inner_points;
        generate_reference_points_das_dennis(inner_points, num_objectives, inner_divisions);

        const double shrink_factor = 0.5;

        // Shrink inner points towards center and renormalize
        for (auto& rp : inner_points) {
            double sum = 0.0;
            for (double& coord : rp.position) {
                coord = shrink_factor + coord * shrink_factor;
                sum += coord;
            }
            // Renormalize to sum to 1
            for (double& coord : rp.position) {
                coord /= sum;
            }
            reference_points.push_back(std::move(rp));
        }
    }
}

/// Compute the population size needed for NSGA-III.
/// Rounds up to the nearest multiple of 4 (for crossover pairing).
///
/// @param num_reference_points Number of reference points
/// @return Population<> size (>= num_reference_points, divisible by 4)
inline int compute_population_size(int num_reference_points) {
    int pop_size = num_reference_points;
    while (pop_size % 4 > 0) {
        ++pop_size;
    }
    return pop_size;
}

} // namespace ea
