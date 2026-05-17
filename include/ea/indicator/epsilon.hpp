#pragma once
/// @file epsilon.hpp
/// @brief Epsilon Indicator (Additive and Multiplicative)
///
/// Reference: Zitzler et al. "Performance Assessment of Multiobjective Optimizers:
/// An Analysis and Review" IEEE Trans. Evol. Comput., 2003.
///
/// The additive epsilon indicator measures the minimum factor by which the
/// approximated front must be shifted to weakly dominate the reference front.
/// Lower is better (0 = perfect).

#include <ea/core/population.hpp>
#include <vector>
#include <cmath>
#include <limits>

namespace ea {

/// Additive Epsilon Indicator (I_epsilon+).
/// Measures the minimum value to add to each objective of the approximated front
/// so that it weakly dominates the reference front.
/// @param pop Population with approximated front
/// @param indices Indices of individuals in the approximated front
/// @param reference_front Reference Pareto front
/// @return Additive epsilon value (lower is better, 0 = perfect)
inline double epsilon_additive(const Population& pop,
                                const std::vector<int>& indices,
                                const std::vector<std::vector<double>>& reference_front) {
    if (reference_front.empty() || indices.empty()) return 0.0;

    double epsilon = 0.0;
    int n_obj = pop.n_obj;

    for (const auto& ref_point : reference_front) {
        double min_epsilon_for_point = std::numeric_limits<double>::infinity();

        for (int idx : indices) {
            double point_epsilon = std::numeric_limits<double>::lowest();
            for (int o = 0; o < n_obj; ++o) {
                double val = pop.objective(idx, o);
                double needed = ref_point[o] - val; // How much we need to add
                if (needed > point_epsilon) point_epsilon = needed;
            }
            if (point_epsilon < min_epsilon_for_point) {
                min_epsilon_for_point = point_epsilon;
            }
        }

        if (min_epsilon_for_point > epsilon) {
            epsilon = min_epsilon_for_point;
        }
    }

    return std::max(0.0, epsilon);
}

/// Multiplicative Epsilon Indicator (I_epsilon*).
/// Measures the minimum factor by which to multiply each objective.
/// @return Multiplicative epsilon value (lower is better, 1 = perfect)
inline double epsilon_multiplicative(const Population& pop,
                                      const std::vector<int>& indices,
                                      const std::vector<std::vector<double>>& reference_front) {
    if (reference_front.empty() || indices.empty()) return 1.0;

    double epsilon = 0.0;
    int n_obj = pop.n_obj;

    for (const auto& ref_point : reference_front) {
        double min_epsilon_for_point = std::numeric_limits<double>::infinity();

        for (int idx : indices) {
            double point_epsilon = std::numeric_limits<double>::lowest();
            for (int o = 0; o < n_obj; ++o) {
                double val = pop.objective(idx, o);
                if (val <= 0 || ref_point[o] <= 0) {
                    // Fallback to additive for non-positive values
                    double needed = ref_point[o] - val;
                    if (needed > point_epsilon) point_epsilon = needed;
                } else {
                    double ratio = ref_point[o] / val;
                    if (ratio > point_epsilon) point_epsilon = ratio;
                }
            }
            if (point_epsilon < min_epsilon_for_point) {
                min_epsilon_for_point = point_epsilon;
            }
        }

        if (min_epsilon_for_point > epsilon) {
            epsilon = min_epsilon_for_point;
        }
    }

    return std::max(1.0, epsilon);
}

/// Epsilon Indicator struct for convenience.
struct EpsilonIndicator {
    bool additive = true; ///< true = additive, false = multiplicative

    double compute(this EpsilonIndicator& self,
                   const Population& pop,
                   const std::vector<int>& indices,
                   const std::vector<std::vector<double>>& reference_front) {
        if (self.additive) {
            return epsilon_additive(pop, indices, reference_front);
        } else {
            return epsilon_multiplicative(pop, indices, reference_front);
        }
    }
};

} // namespace ea
