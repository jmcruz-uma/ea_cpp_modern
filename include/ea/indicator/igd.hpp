#pragma once
/// @file igd.hpp
/// @brief Inverted Generational Distance (IGD) quality indicator.
/// Measures both convergence and diversity of an approximation set.

#include <vector>
#include <cmath>
#include <algorithm>
#include <ea/core/population.hpp>

namespace ea {

/// Compute Inverted Generational Distance between an approximation set and a reference front.
/// IGD = (1/|P*|) * sum_{p in P*} min_{a in A} dist(p, a)
/// Lower is better. Requires a well-distributed reference front P*.
/// @param approximation The approximation set (population with evaluated objectives)
/// @param reference_front Vector of reference points (each a vector of objectives)
/// @return IGD value
inline double compute_igd(const Population& approximation,
                          const std::vector<std::vector<double>>& reference_front) {
    if (reference_front.empty() || approximation.pop_size == 0) return 0.0;

    double sum = 0.0;
    int n_obj = approximation.n_obj;

    for (const auto& ref_point : reference_front) {
        double min_dist = std::numeric_limits<double>::max();
        for (int i = 0; i < approximation.pop_size; ++i) {
            double dist = 0.0;
            for (int o = 0; o < n_obj; ++o) {
                double diff = approximation.objective(i, o) - ref_point[o];
                dist += diff * diff;
            }
            min_dist = std::min(min_dist, std::sqrt(dist));
        }
        sum += min_dist;
    }

    return sum / static_cast<double>(reference_front.size());
}

/// Compute Generational Distance between an approximation set and a reference front.
/// GD = (1/|A|) * sqrt(sum_{a in A} min_{p in P*} dist(a, p)^2)
/// @param approximation The approximation set
/// @param reference_front Vector of reference points
/// @return GD value
inline double compute_gd(const Population& approximation,
                          const std::vector<std::vector<double>>& reference_front) {
    if (reference_front.empty() || approximation.pop_size == 0) return 0.0;

    int n_obj = approximation.n_obj;
    double sum = 0.0;

    for (int i = 0; i < approximation.pop_size; ++i) {
        double min_dist = std::numeric_limits<double>::max();
        for (const auto& ref_point : reference_front) {
            double dist = 0.0;
            for (int o = 0; o < n_obj; ++o) {
                double diff = approximation.objective(i, o) - ref_point[o];
                dist += diff * diff;
            }
            min_dist = std::min(min_dist, dist);
        }
        sum += min_dist;
    }

    return std::sqrt(sum) / static_cast<double>(approximation.pop_size);
}

/// Compute Spread (Δ) indicator for 2-objective problems.
/// Measures how well-distributed the approximation set is along the Pareto front.
/// Δ = (d_f + d_l + sum |d_i - d_avg|) / (d_f + d_l + (n-1)*d_avg)
/// @param approximation The approximation set (must be 2 objectives)
/// @param reference_front Vector of reference points (for extreme points)
/// @return Spread value (0 = perfect distribution)
inline double compute_spread(const Population& approximation,
                              const std::vector<std::vector<double>>& reference_front) {
    if (approximation.n_obj != 2 || approximation.pop_size < 2) return 0.0;

    int n = approximation.pop_size;

    // Find extreme points in reference front
    double f1_min_ref = std::numeric_limits<double>::max();
    double f2_min_ref = std::numeric_limits<double>::max();
    for (const auto& pt : reference_front) {
        f1_min_ref = std::min(f1_min_ref, pt[0]);
        f2_min_ref = std::min(f2_min_ref, pt[1]);
    }

    // Sort approximation by first objective
    std::vector<int> sorted(n);
    std::iota(sorted.begin(), sorted.end(), 0);
    std::sort(sorted.begin(), sorted.end(), [&](int a, int b) {
        return approximation.objective(a, 0) < approximation.objective(b, 0);
    });

    // Distance to extreme points
    double df = std::sqrt(std::pow(approximation.objective(sorted[0], 0) - f1_min_ref, 2) +
                           std::pow(approximation.objective(sorted[0], 1) - 0, 2));
    double dl = std::sqrt(std::pow(approximation.objective(sorted[n-1], 0) - 0, 2) +
                           std::pow(approximation.objective(sorted[n-1], 1) - f2_min_ref, 2));

    // Consecutive distances
    std::vector<double> distances(n - 1);
    double sum_dist = 0.0;
    for (int i = 0; i < n - 1; ++i) {
        double d = 0.0;
        for (int o = 0; o < 2; ++o) {
            double diff = approximation.objective(sorted[i+1], o) - approximation.objective(sorted[i], o);
            d += diff * diff;
        }
        distances[i] = std::sqrt(d);
        sum_dist += distances[i];
    }

    double d_avg = sum_dist / (n - 1);
    if (d_avg < 1e-14) return 0.0;

    double sum_dev = 0.0;
    for (double d : distances) {
        sum_dev += std::abs(d - d_avg);
    }

    return (df + dl + sum_dev) / (df + dl + (n - 1) * d_avg);
}

} // namespace ea