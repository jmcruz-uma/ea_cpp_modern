#pragma once
/// @file comparator.hpp
/// @brief Dominance and crowding distance comparators for multi-objective optimization.

#include <ea/core/population.hpp>
#include <cmath>
#include <algorithm>
#include <limits>

namespace ea {

/// Pareto dominance comparison result.
enum class Dominance : int8_t {
    Dominates = -1,  ///< a dominates b
    Equal = 0,       ///< a and b are non-dominated
    Dominated = 1    ///< b dominates a
};

/// Compare two individuals for Pareto dominance.
/// @return Dominance::Dominates if a dominates b,
///         Dominance::Dominated if b dominates a,
///         Dominance::Equal otherwise
inline Dominance compare_dominance(const Population& pop, int a, int b) {
    bool a_dominates_b = false;
    bool b_dominates_a = false;

    for (int o = 0; o < pop.n_obj; ++o) {
        double fa = pop.objective(a, o);
        double fb = pop.objective(b, o);
        if (fa < fb) a_dominates_b = true;
        else if (fb < fa) b_dominates_a = true;

        // Early exit: if both dominate in different objectives, they're non-dominated
        if (a_dominates_b && b_dominates_a) return Dominance::Equal;
    }

    if (a_dominates_b && !b_dominates_a) return Dominance::Dominates;
    if (b_dominates_a && !a_dominates_b) return Dominance::Dominated;
    return Dominance::Equal;
}

/// Compute crowding distances for a set of individuals (for NSGA-II).
/// Matches jMetal's CrowdingDistanceDensityEstimator exactly:
/// - Normalized by (max_obj - min_obj) per objective
/// - Boundary solutions get infinite distance
/// - Distances are assigned to the correct individuals (not sorted positions)
/// @param pop Population
/// @param indices Indices of individuals in the front
/// @param crowding_dist Output crowding distances indexed by individual position in `indices`
inline void compute_crowding_distance(const Population& pop,
                                       const std::vector<int>& indices,
                                       std::vector<double>& crowding_dist) {
    const int n = static_cast<int>(indices.size());
    crowding_dist.assign(n, 0.0);

    if (n <= 2) {
        for (int i = 0; i < n; ++i)
            crowding_dist[i] = std::numeric_limits<double>::infinity();
        return;
    }

    for (int o = 0; o < pop.n_obj; ++o) {
        // Sort indices by objective o
        std::vector<int> sorted(indices);
        std::sort(sorted.begin(), sorted.end(),
            [&pop, o](int a, int b) {
                return pop.objective(a, o) < pop.objective(b, o);
            });

        double f_min = pop.objective(sorted[0], o);
        double f_max = pop.objective(sorted[n - 1], o);

        if (f_min == f_max) continue;

        double range = f_max - f_min;

        // Boundary solutions: find their position in `indices` and set to infinity
        // sorted[0] is the individual with min objective, sorted[n-1] with max
        int idx_min = sorted[0];
        int idx_max = sorted[n - 1];

        // Map individual index back to position in `indices` array
        for (int pos = 0; pos < n; ++pos) {
            if (indices[pos] == idx_min || indices[pos] == idx_max) {
                crowding_dist[pos] = std::numeric_limits<double>::infinity();
            }
        }

        // Intermediate solutions: add normalized distance
        for (int i = 1; i < n - 1; ++i) {
            double f_prev = pop.objective(sorted[i - 1], o);
            double f_next = pop.objective(sorted[i + 1], o);
            double distance = (f_next - f_prev) / range;

            // Find position of sorted[i] in `indices`
            int individual = sorted[i];
            for (int pos = 0; pos < n; ++pos) {
                if (indices[pos] == individual) {
                    crowding_dist[pos] += distance;
                    break;
                }
            }
        }
    }
}

/// Fast non-dominated sort. Returns fronts as vectors of individual indices.
/// @param pop Population (must have objectives evaluated)
/// @return Vector of fronts, front[0] is the best (non-dominated) front
inline std::vector<std::vector<int>> fast_non_dominated_sort(Population& pop) {
    const int n = pop.pop_size;
    std::vector<int> domination_count(n, 0);
    std::vector<std::vector<int>> dominated_by(n);
    std::vector<std::vector<int>> fronts;
    fronts.push_back({});

    for (int p = 0; p < n; ++p) {
        for (int q = 0; q < n; ++q) {
            if (p == q) continue;
            auto dom = compare_dominance(pop, p, q);
            if (dom == Dominance::Dominates) {
                dominated_by[p].push_back(q);
            } else if (dom == Dominance::Dominated) {
                domination_count[p]++;
            }
        }
        if (domination_count[p] == 0) {
            fronts[0].push_back(p);
        }
    }

    int k = 0;
    while (!fronts[k].empty()) {
        std::vector<int> next_front;
        for (int p : fronts[k]) {
            for (int q : dominated_by[p]) {
                domination_count[q]--;
                if (domination_count[q] == 0) {
                    next_front.push_back(q);
                }
            }
        }
        k++;
        fronts.push_back(std::move(next_front));
    }

    // Remove last empty front
    if (!fronts.back().empty() == false) {
        fronts.pop_back();
    }

    return fronts;
}

/// Compute overall constraint violation for an individual.
/// Returns the sum of constraint violations (negative values = feasible).
inline double constraint_violation(const Population& pop, int idx) {
    if (pop.n_const == 0) return 0.0;
    double violation = 0.0;
    for (int c = 0; c < pop.n_const; ++c) {
        double val = pop.constraint(idx, c);
        if (val < 0.0) violation += val; // Only count violations
    }
    return violation;
}

} // namespace ea