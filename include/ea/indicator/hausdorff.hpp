#pragma once
/// @file hausdorff.hpp
/// @brief Average Hausdorff Distance
///
/// Reference: Schutze et al. "Averaging on the Broyden Class of Quasi-Newton Methods"
///
/// The Average Hausdorff Distance combines Generational Distance (GD) and
/// Inverted Generational Distance (IGD) to measure both convergence and diversity.
/// It is defined as: max(GD, IGD) or alternatively (GD + IGD) / 2.
/// Lower is better.

#include <algorithm>
#include <ea/core/population.hpp>
#include <ea/indicator/gd.hpp>
#include <ea/indicator/igd.hpp>
#include <vector>

namespace ea {

/// Average Hausdorff Distance.
/// Combines GD and IGD to measure both convergence and diversity.
/// @param pop Population with approximated front
/// @param indices Indices of individuals in the approximated front
/// @param reference_front Reference Pareto front
/// @param p Power parameter for GD/IGD (default 2.0)
/// @return Hausdorff distance (lower is better)
inline double average_hausdorff(const Population& pop, const std::vector<int>& indices,
                                const std::vector<std::vector<double>>& reference_front,
                                double p = 2.0) {
    double gd_val = gd(pop, indices, reference_front, p);
    double igd_val = igd(pop, indices, reference_front);

    // Average Hausdorff: max(GD, IGD)
    // Alternative: (GD + IGD) / 2
    return std::max(gd_val, igd_val);
}

/// Average Hausdorff Distance from vector fronts.
inline double average_hausdorff(const std::vector<std::vector<double>>& front,
                                const std::vector<std::vector<double>>& reference_front,
                                double p = 2.0) {
    double gd_val = gd(front, reference_front, p);
    double igd_val = igd(front, reference_front);
    return std::max(gd_val, igd_val);
}

/// Average Hausdorff Distance struct for convenience.
struct AverageHausdorffIndicator {
    double p = 2.0;
    bool use_max = true; ///< true = max(GD, IGD), false = (GD + IGD)/2

    double compute(this AverageHausdorffIndicator& self, const Population& pop,
                   const std::vector<int>& indices,
                   const std::vector<std::vector<double>>& reference_front) {
        double gd_val = gd(pop, indices, reference_front, self.p);
        double igd_val = igd(pop, indices, reference_front);

        if (self.use_max) {
            return std::max(gd_val, igd_val);
        } else {
            return (gd_val + igd_val) / 2.0;
        }
    }
};

} // namespace ea
