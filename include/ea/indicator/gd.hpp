#pragma once
/// @file gd.hpp
/// @brief Generational Distance (GD) indicator.
///
/// Reference: jMetal
///   jmetal-core/src/main/java/org/uma/jmetal/qualityindicator/impl/GenerationalDistance.java
///
/// GD measures the average distance from each point in the approximated front
/// to its closest point in the reference front. Lower is better.

#include <cmath>
#include <ea/core/population.hpp>
#include <limits>
#include <vector>

namespace ea {

/// Generational Distance indicator.
///
/// @param pop            Population<> containing the approximated front
/// @param indices        Indices of individuals in the approximated front
/// @param reference_front Reference Pareto front
/// @param p              Power parameter (default 2.0 for Euclidean)
/// @return GD value (lower is better)
inline double gd(const Population<>& pop, const std::vector<int>& indices,
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
    for (const auto& point : front) {
        double min_dist = std::numeric_limits<double>::infinity();
        for (const auto& ref_point : reference_front) {
            double dist = 0.0;
            for (size_t i = 0; i < point.size(); ++i) {
                double diff = point[i] - ref_point[i];
                dist += diff * diff;
            }
            dist = std::sqrt(dist);
            if (dist < min_dist)
                min_dist = dist;
        }
        sum += std::pow(min_dist, p);
    }

    sum = std::pow(sum, 1.0 / p);
    return sum / static_cast<double>(front.size());
}

/// GD from two vector fronts.
inline double gd(const std::vector<std::vector<double>>& front,
                 const std::vector<std::vector<double>>& reference_front, double p = 2.0) {
    if (reference_front.empty() || front.empty())
        return 0.0;

    double sum = 0.0;
    for (const auto& point : front) {
        double min_dist = std::numeric_limits<double>::infinity();
        for (const auto& ref_point : reference_front) {
            double dist = 0.0;
            for (size_t i = 0; i < point.size(); ++i) {
                double diff = point[i] - ref_point[i];
                dist += diff * diff;
            }
            dist = std::sqrt(dist);
            if (dist < min_dist)
                min_dist = dist;
        }
        sum += std::pow(min_dist, p);
    }

    sum = std::pow(sum, 1.0 / p);
    return sum / static_cast<double>(front.size());
}

struct GDIndicator {
    double p = 2.0;

    double compute(this GDIndicator& self, const Population<>& pop, const std::vector<int>& indices,
                   const std::vector<std::vector<double>>& reference_front) {
        return gd(pop, indices, reference_front, self.p);
    }

    double compute(this GDIndicator& self, const std::vector<std::vector<double>>& front,
                   const std::vector<std::vector<double>>& reference_front) {
        return gd(front, reference_front, self.p);
    }
};

} // namespace ea
