#pragma once
/// @file aggregation.hpp
/// @brief Aggregation functions for decomposition-based MOEAs (MOEA/D).
/// Reference: Zhang & Li, "MOEA/D: A Multiobjective Evolutionary Algorithm
/// Based on Decomposition", IEEE TEVC 2007.
///
/// jMetal reference:
/// - org.uma.jmetal.util.aggregationfunction.impl.WeightedSum
/// - org.uma.jmetal.util.aggregationfunction.impl.Tschebyscheff
/// - org.uma.jmetal.util.aggregationfunction.impl.PenaltyBoundaryIntersection

#include <cmath>
#include <algorithm>
#include <string_view>
#include <stdexcept>

namespace ea {

/// Type of aggregation function used in MOEA/D.
/// PBI is also known as Penalty-based Boundary Intersection.
enum class AggregationType : uint8_t {
    WeightedSum,  ///< Weighted sum aggregation
    Tchebycheff,  ///< Tchebycheff (minimax) aggregation
    PBI,          ///< Penalty Boundary Intersection (not implemented yet)
};

/// Parse aggregation type from string.
constexpr AggregationType aggregation_type_from_string(std::string_view s) {
    if (s == "WEIGHTED_SUM" || s == "WS" || s == "AGG") {
        return AggregationType::WeightedSum;
    } else if (s == "TCHEBYCHEFF" || s == "TCHE" || s == "TCH") {
        return AggregationType::Tchebycheff;
    } else if (s == "PBI") {
        return AggregationType::PBI;
    }
    throw std::invalid_argument("Invalid aggregation type string");
}

// ============================================================
// Weighted Sum Aggregation
// ============================================================

/// Weighted sum aggregation: sum_i(weight[i] * (objective[i] - ideal[i])).
///
/// Reference: org.uma.jmetal.util.aggregationfunction.impl.WeightedSum
///
/// Optionally normalizes using nadir point:
/// value = (obj[i] - ideal[i]) / (nadir[i] - ideal[i] + epsilon)
struct WeightedSum {
    bool normalize = false;  ///< Whether to normalize using ideal/nadir points
    double epsilon = 1e-6;   ///< Small value to avoid division by zero

    /// Compute weighted sum aggregation.
    /// @param obj Objective values (size n_obj)
    /// @param weight Weight vector (size n_obj), sum of weights = 1.0
    /// @param ideal Ideal point (size n_obj)
    /// @param nadir Nadir point (size n_obj), only used if normalize = true
    /// @param n_obj Number of objectives
    double compute(this const WeightedSum& self,
                   const double* obj, const double* weight,
                   const double* ideal, const double* nadir, int n_obj) {
        double sum = 0.0;
        for (int n = 0; n < n_obj; ++n) {
            double value;
            if (self.normalize) {
                double range = nadir[n] - ideal[n];
                value = (obj[n] - ideal[n]) / (range + self.epsilon);
            } else {
                value = obj[n] - ideal[n];
            }
            sum += weight[n] * value;
        }
        return sum;
    }
};

// ============================================================
// Tchebycheff Aggregation
// ============================================================

/// Tchebycheff (minimax) aggregation:
/// max_i(|objective[i] - ideal[i]| * weight[i]).
///
/// Reference: org.uma.jmetal.util.aggregationfunction.impl.Tschebyscheff
///
/// If weight[i] == 0, the term becomes 0.0001 * diff (jMetal convention).
/// Optionally normalizes using nadir point.
struct Tchebycheff {
    bool normalize = false;  ///< Whether to normalize using ideal/nadir points
    double epsilon = 1e-6;     ///< Small value to avoid division by zero

    /// Compute Tchebycheff aggregation.
    /// @param obj Objective values (size n_obj)
    /// @param weight Weight vector (size n_obj), sum of weights = 1.0
    /// @param ideal Ideal point (size n_obj)
    /// @param nadir Nadir point (size n_obj), only used if normalize = true
    /// @param n_obj Number of objectives
    double compute(this const Tchebycheff& self,
                   const double* obj, const double* weight,
                   const double* ideal, const double* nadir, int n_obj) {
        double max_fun = -1.0e+30;

        for (int n = 0; n < n_obj; ++n) {
            double diff;
            if (self.normalize) {
                double range = nadir[n] - ideal[n];
                diff = std::abs((obj[n] - ideal[n]) / (range + self.epsilon));
            } else {
                diff = std::abs(obj[n] - ideal[n]);
            }

            double feval;
            if (weight[n] == 0.0) {
                feval = 0.0001 * diff;
            } else {
                feval = diff * weight[n];
            }
            if (feval > max_fun) {
                max_fun = feval;
            }
        }

        return max_fun;
    }
};

// ============================================================
// Penalty Boundary Intersection (PBI)
// ============================================================

/// Penalty Boundary Intersection aggregation.
///
/// Reference: org.uma.jmetal.util.aggregationfunction.impl.PenaltyBoundaryIntersection
///
/// fitness = d1 + theta * d2 where:
///   d1 = |sum_i((obj[i] - ideal[i]) * weight[i]) / ||weight||
///   d2 = ||(obj - ideal) - d1 * (weight / ||weight||)||
struct PenaltyBoundaryIntersection {
    double theta = 5.0;      ///< Penalty parameter controlling diversity
    bool normalize = false;  ///< Whether to normalize using ideal/nadir points
    double epsilon = 1e-6;   ///< Small value to avoid division by zero

    /// Compute PBI aggregation.
    /// @param obj Objective values (size n_obj)
    /// @param weight Weight vector (size n_obj), sum of weights = 1.0
    /// @param ideal Ideal point (size n_obj)
    /// @param nadir Nadir point (size n_obj), only used if normalize = true
    /// @param n_obj Number of objectives
    double compute(this const PenaltyBoundaryIntersection& self,
                   const double* obj, const double* weight,
                   const double* ideal, const double* nadir, int n_obj) {
        double d1 = 0.0;
        double d2 = 0.0;
        double nl = 0.0;

        for (int i = 0; i < n_obj; ++i) {
            double diff;
            if (self.normalize) {
                double range = nadir[i] - ideal[i];
                diff = (obj[i] - ideal[i]) / (range + self.epsilon);
            } else {
                diff = obj[i] - ideal[i];
            }
            d1 += diff * weight[i];
            nl += weight[i] * weight[i];
        }
        nl = std::sqrt(nl);
        d1 = std::abs(d1) / nl;

        for (int i = 0; i < n_obj; ++i) {
            double diff;
            if (self.normalize) {
                double range = nadir[i] - ideal[i];
                diff = (obj[i] - ideal[i]) / (range + self.epsilon);
            } else {
                diff = obj[i] - ideal[i];
            }
            double t = diff - d1 * (weight[i] / nl);
            d2 += t * t;
        }
        d2 = std::sqrt(d2);

        return d1 + self.theta * d2;
    }
};

// ============================================================
// Generic dispatch
// ============================================================

/// Generic aggregation function that dispatches to the concrete implementation.
/// Useful for runtime configuration via variant.
struct AggregationFunction {
    AggregationType type = AggregationType::Tchebycheff;
    bool normalize = false;
    double theta = 5.0;    ///< Only used for PBI
    double epsilon = 1e-6;

    /// Compute aggregation value for the given objective vector.
    double compute(this const AggregationFunction& self,
                   const double* obj, const double* weight,
                   const double* ideal, const double* nadir, int n_obj) {
        switch (self.type) {
            case AggregationType::WeightedSum: {
                WeightedSum aggr{self.normalize, self.epsilon};
                return aggr.compute(obj, weight, ideal, nadir, n_obj);
            }
            case AggregationType::Tchebycheff: {
                Tchebycheff aggr{self.normalize, self.epsilon};
                return aggr.compute(obj, weight, ideal, nadir, n_obj);
            }
            case AggregationType::PBI: {
                PenaltyBoundaryIntersection aggr{self.theta, self.normalize, self.epsilon};
                return aggr.compute(obj, weight, ideal, nadir, n_obj);
            }
        }
        // unreachable
        return std::numeric_limits<double>::quiet_NaN();
    }
};

} // namespace ea
