#pragma once
/// @file spread.hpp
/// @brief Spread (Δ) indicator — diversity metric for Pareto front approximation.
///
/// Reference: Deb et al. "A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II"
///   IEEE Transactions on Evolutionary Computation, 2002.
///
///   jmetal-core/src/main/java/org/uma/jmetal/qualityindicator/impl/Spread.java
///
/// Spread measures the diversity of the approximated front by computing the
/// sum of distances between consecutive solutions, comparing them against the
/// average distance. A perfect distribution gives Spread = 0. Lower is better.
///
/// Formula:
///   Spread = (df + dl + Σ|di - d̄|) / (df + dl + (N-1)*d̄)
///
/// where:
///   df = distance from first extreme of approximated front to first extreme of reference front
///   dl = distance from last extreme of approximated front to last extreme of reference front
///   di = distance between consecutive solutions in the approximated front (sorted)
///   d̄ = average of all di
///   N = number of solutions in the approximated front

#include <ea/core/population.hpp>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

namespace ea {
namespace detail {

/// Compute Euclidean distance between two objective vectors.
inline double euclidean_dist(const std::vector<double>& a,
                              const std::vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

/// Sort a front by the first objective dimension (ascending).
inline std::vector<std::vector<double>> sort_by_first_objective(
    const std::vector<std::vector<double>>& front) {
    auto sorted = front;
    std::sort(sorted.begin(), sorted.end(),
              [](const auto& a, const auto& b) { return a[0] < b[0]; });
    return sorted;
}

/// Find extreme points (minimum and maximum of first objective) in a front.
/// Returns {min_first_obj_point, max_first_obj_point}.
inline std::vector<std::vector<double>> find_extremes(
    const std::vector<std::vector<double>>& front) {
    if (front.empty()) return {};
    
    auto min_it = std::min_element(front.begin(), front.end(),
        [](const auto& a, const auto& b) { return a[0] < b[0]; });
    auto max_it = std::max_element(front.begin(), front.end(),
        [](const auto& a, const auto& b) { return a[0] < b[0]; });
    
    return {*min_it, *max_it};
}

/// Core Spread computation from two sorted fronts.
inline double spread_impl(const std::vector<std::vector<double>>& sorted_front,
                           const std::vector<std::vector<double>>& sorted_ref) {
    if (sorted_front.size() < 2) return 0.0;

    // Find extreme points
    auto front_extremes = find_extremes(sorted_front);
    auto ref_extremes = find_extremes(sorted_ref);

    // df = distance from first extreme of front to first extreme of reference
    double df = euclidean_dist(front_extremes[0], ref_extremes[0]);
    // dl = distance from last extreme of front to last extreme of reference
    double dl = euclidean_dist(front_extremes[1], ref_extremes[1]);

    // Compute distances between consecutive solutions in the approximated front
    std::vector<double> distances;
    distances.reserve(sorted_front.size() - 1);
    for (size_t i = 1; i < sorted_front.size(); ++i) {
        distances.push_back(euclidean_dist(sorted_front[i], sorted_front[i - 1]));
    }

    // Compute average distance d̄
    double d_bar = 0.0;
    for (double d : distances) d_bar += d;
    d_bar /= static_cast<double>(distances.size());

    if (d_bar == 0.0) return 0.0;

    // Compute Σ|di - d̄|
    double sum_deviation = 0.0;
    for (double d : distances) {
        sum_deviation += std::abs(d - d_bar);
    }

    // Spread = (df + dl + Σ|di - d̄|) / (df + dl + (N-1)*d̄)
    double numerator = df + dl + sum_deviation;
    double denominator = df + dl + (sorted_front.size() - 1) * d_bar;

    return numerator / denominator;
}

} // namespace detail

/// Spread (Δ) indicator from two vector fronts.
inline double spread(const std::vector<std::vector<double>>& front,
                     const std::vector<std::vector<double>>& reference_front) {
    if (reference_front.empty() || front.size() < 2) return 0.0;

    auto sorted_front = detail::sort_by_first_objective(front);
    auto sorted_ref = detail::sort_by_first_objective(reference_front);

    return detail::spread_impl(sorted_front, sorted_ref);
}

/// Spread (Δ) indicator.
///
/// @param pop             Population containing the approximated front
/// @param indices         Indices of individuals in the approximated front
/// @param reference_front Reference Pareto front
/// @return Spread value (lower is better, 0 = perfect distribution)
inline double spread(const Population& pop,
                     const std::vector<int>& indices,
                     const std::vector<std::vector<double>>& reference_front) {
    if (reference_front.empty() || indices.size() < 2) return 0.0;

    // Extract front points from population
    std::vector<std::vector<double>> front;
    front.reserve(indices.size());
    for (int idx : indices) {
        std::vector<double> point(pop.n_obj);
        for (int o = 0; o < pop.n_obj; ++o) point[o] = pop.objective(idx, o);
        front.push_back(point);
    }

    return spread(front, reference_front);
}

struct SpreadIndicator {
    double compute(this SpreadIndicator& self,
                   const Population& pop,
                   const std::vector<int>& indices,
                   const std::vector<std::vector<double>>& reference_front) {
        return spread(pop, indices, reference_front);
    }

    double compute(this SpreadIndicator& self,
                   const std::vector<std::vector<double>>& front,
                   const std::vector<std::vector<double>>& reference_front) {
        return spread(front, reference_front);
    }
};

} // namespace ea
