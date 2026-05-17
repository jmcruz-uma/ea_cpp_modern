#pragma once
/// @file r2.hpp
/// @brief R2 Indicator (R2-metric)
///
/// Reference: J. Knowles and D. Corne, "Properties of an Adaptive Archiving Algorithm for
/// Storing Nondominated Vectors" IEEE Trans. Evol. Comput., 2003.
///
/// The R2 indicator measures the utility of an approximated front with respect to
/// a set of reference vectors (weights). Lower is better.
/// Uses Tchebycheff utility function.

#include <ea/core/population.hpp>
#include <vector>
#include <cmath>
#include <limits>

namespace ea {

/// R2 Indicator.
/// Computes the average utility difference between the reference front and the approximated front.
/// @param pop Population with approximated front
/// @param indices Indices of individuals in the approximated front
/// @param reference_front Reference Pareto front
/// @param weights Vector of weight vectors (one per direction)
/// @param ideal Ideal point (minimum of reference front)
/// @param nadir Nadir point (maximum of reference front)
/// @return R2 value (lower is better)
inline double r2(const Population& pop,
                 const std::vector<int>& indices,
                 const std::vector<std::vector<double>>& reference_front,
                 const std::vector<std::vector<double>>& weights,
                 const std::vector<double>& ideal,
                 const std::vector<double>& nadir) {
    if (reference_front.empty() || indices.empty() || weights.empty()) return 0.0;

    int n_obj = pop.n_obj;
    if (ideal.size() != static_cast<size_t>(n_obj) || nadir.size() != static_cast<size_t>(n_obj)) {
        return 0.0; // Mismatched dimensions
    }

    double r2_value = 0.0;

    for (const auto& weight : weights) {
        // Best utility in reference front for this weight
        double best_ref = std::numeric_limits<double>::infinity();
        for (const auto& ref_point : reference_front) {
            double max_val = std::numeric_limits<double>::lowest();
            for (int o = 0; o < n_obj; ++o) {
                double range = nadir[o] - ideal[o];
                if (range < 1e-12) range = 1.0;
                double normalized = (ref_point[o] - ideal[o]) / range;
                double utility = weight[o] * std::abs(normalized);
                if (utility > max_val) max_val = utility;
            }
            if (max_val < best_ref) best_ref = max_val;
        }

        // Best utility in approximated front for this weight
        double best_approx = std::numeric_limits<double>::infinity();
        for (int idx : indices) {
            double max_val = std::numeric_limits<double>::lowest();
            for (int o = 0; o < n_obj; ++o) {
                double range = nadir[o] - ideal[o];
                if (range < 1e-12) range = 1.0;
                double normalized = (pop.objective(idx, o) - ideal[o]) / range;
                double utility = weight[o] * std::abs(normalized);
                if (utility > max_val) max_val = utility;
            }
            if (max_val < best_approx) best_approx = max_val;
        }

        r2_value += std::abs(best_approx - best_ref);
    }

    return r2_value / static_cast<double>(weights.size());
}

/// R2 Indicator with auto-computed ideal and nadir.
inline double r2(const Population& pop,
                 const std::vector<int>& indices,
                 const std::vector<std::vector<double>>& reference_front,
                 const std::vector<std::vector<double>>& weights) {
    if (reference_front.empty()) return 0.0;

    int n_obj = static_cast<int>(reference_front[0].size());
    std::vector<double> ideal(n_obj, std::numeric_limits<double>::infinity());
    std::vector<double> nadir(n_obj, std::numeric_limits<double>::lowest());

    for (const auto& point : reference_front) {
        for (int o = 0; o < n_obj; ++o) {
            if (point[o] < ideal[o]) ideal[o] = point[o];
            if (point[o] > nadir[o]) nadir[o] = point[o];
        }
    }

    return r2(pop, indices, reference_front, weights, ideal, nadir);
}

/// R2 Indicator struct for convenience.
struct R2Indicator {
    std::vector<std::vector<double>> weights; ///< Weight vectors
    std::vector<double> ideal;                     ///< Ideal point (optional)
    std::vector<double> nadir;                     ///< Nadir point (optional)

    double compute(this R2Indicator& self,
                   const Population& pop,
                   const std::vector<int>& indices,
                   const std::vector<std::vector<double>>& reference_front) {
        if (self.weights.empty()) {
            // Generate default weights (uniform)
            int n_obj = pop.n_obj;
            int n_weights = 100;
            self.weights = generate_uniform_weights(n_obj, n_weights);
        }

        if (self.ideal.empty() || self.nadir.empty()) {
            return r2(pop, indices, reference_front, self.weights);
        } else {
            return r2(pop, indices, reference_front, self.weights, self.ideal, self.nadir);
        }
    }

private:
    static std::vector<std::vector<double>> generate_uniform_weights(int n_obj, int n_weights) {
        std::vector<std::vector<double>> weights;
        if (n_obj == 2) {
            for (int i = 0; i < n_weights; ++i) {
                double w0 = static_cast<double>(i) / (n_weights - 1);
                weights.push_back({w0, 1.0 - w0});
            }
        } else {
            // Simple random weights for higher dimensions
            for (int i = 0; i < n_weights; ++i) {
                std::vector<double> w(n_obj, 0.0);
                double sum = 0.0;
                for (int o = 0; o < n_obj; ++o) {
                    w[o] = static_cast<double>(rand()) / RAND_MAX;
                    sum += w[o];
                }
                for (int o = 0; o < n_obj; ++o) w[o] /= sum;
                weights.push_back(w);
            }
        }
        return weights;
    }
};

} // namespace ea
