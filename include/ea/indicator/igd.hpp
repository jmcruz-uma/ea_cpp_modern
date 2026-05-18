#pragma once
/// @file igd.hpp
/// @brief Inverted Generational Distance (IGD) indicator.
///
/// Reference: jMetal
///   jmetal-core/src/main/java/org/uma/jmetal/qualityindicator/impl/InvertedGenerationalDistance.java
///
/// IGD measures the average distance from each point in the reference front
/// to its closest point in the approximated front. Lower is better.

#include <cmath>
#include <ea/core/population.hpp>
#include <limits>
#include <vector>

namespace ea {

/// Compute Euclidean distance between two objective vectors.
inline double euclidean_distance(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

/// Find the minimum distance from a point to any point in a front.
inline double distance_to_closest(const std::vector<double>& point,
                                  const std::vector<std::vector<double>>& front) {
    double min_dist = std::numeric_limits<double>::infinity();
    for (const auto& f : front) {
        double dist = euclidean_distance(point, f);
        if (dist < min_dist)
            min_dist = dist;
    }
    return min_dist;
}

/// Inverted Generational Distance indicator.
///
/// @param pop            Population containing the approximated front
/// @param indices        Indices of individuals in the approximated front
/// @param reference_front Reference Pareto front
/// @param p              Power parameter (default 2.0 for Euclidean)
/// @return IGD value (lower is better)
inline double igd(const Population& pop, const std::vector<int>& indices,
                  const std::vector<std::vector<double>>& reference_front, double p = 2.0) {
    if (reference_front.empty() || indices.empty())
        return 0.0;

    // Extract front points from population
    std::vector<std::vector<double>> front;
    front.reserve(indices.size());
    for (int idx : indices) {
        std::vector<double> point(pop.n_obj);
        for (int o = 0; o < pop.n_obj; ++o)
            point[o] = pop.objective(idx, o);
        front.push_back(point);
    }

    double sum = 0.0;
    for (const auto& ref_point : reference_front) {
        double min_dist = distance_to_closest(ref_point, front);
        sum += std::pow(min_dist, p);
    }

    sum = std::pow(sum, 1.0 / p);
    return sum / static_cast<double>(reference_front.size());
}

/// IGD from two vector fronts.
inline double igd(const std::vector<std::vector<double>>& front,
                  const std::vector<std::vector<double>>& reference_front, double p = 2.0) {
    if (reference_front.empty() || front.empty())
        return 0.0;

    double sum = 0.0;
    for (const auto& ref_point : reference_front) {
        double min_dist = distance_to_closest(ref_point, front);
        sum += std::pow(min_dist, p);
    }

    sum = std::pow(sum, 1.0 / p);
    return sum / static_cast<double>(reference_front.size());
}

struct IGDIndicator {
    double p = 2.0;

    double compute(this IGDIndicator& self, const Population& pop, const std::vector<int>& indices,
                   const std::vector<std::vector<double>>& reference_front) {
        return igd(pop, indices, reference_front, self.p);
    }

    double compute(this IGDIndicator& self, const std::vector<std::vector<double>>& front,
                   const std::vector<std::vector<double>>& reference_front) {
        return igd(front, reference_front, self.p);
    }
};

} // namespace ea
