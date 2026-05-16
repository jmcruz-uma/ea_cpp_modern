#pragma once
/// @file dtlz.hpp
/// @brief DTLZ benchmark problem suite for multi-objective optimization.
///
/// Reference: jMetal jmetal-problem/src/main/java/org/uma/jmetal/problem/multiobjective/dtlz/
///
/// The DTLZ (Deb-Thiele-Laumanns-Zitzler) test suite is scalable in both
/// number of objectives and number of variables. All problems share a common
/// structure with position variables (first M-1) and distance variables (remaining k).
/// Bounds: [0, 1] for all variables.

#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <span>
#include <cmath>
#include <vector>
#include <numeric>

namespace ea {

// ============================================================
// DTLZ1 — linear Pareto-optimal front
// ============================================================

/// DTLZ1: linear Pareto-optimal front.
/// Default: 7 variables, 3 objectives.
/// k = dim - num_obj + 1 distance variables.
struct DTLZ1 {
    int dim_ = 7;  ///< Number of decision variables
    int n_obj_ = 3; ///< Number of objectives
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    DTLZ1() : DTLZ1(7, 3) {}

    DTLZ1(int dim, int n_obj) : dim_(dim), n_obj_(n_obj) {
        if (dim_ < n_obj_) dim_ = n_obj_; // k >= 1 required
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    int num_objectives() const { return n_obj_; }
    int dimension() const { return dim_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        int k = dim_ - n_obj_ + 1;
        double g = 0.0;
        for (int i = dim_ - k; i < dim_; ++i) {
            double xi = pop.gene(idx, i);
            g += (xi - 0.5) * (xi - 0.5) - std::cos(20.0 * M_PI * (xi - 0.5));
        }
        g = 100.0 * (k + g);

        // f_i = (1 + g) * 0.5 * product(x_j for j=0..M-(i+1)-1)
        // For i > 0: f_i *= (1 - x_{M-(i+1)})
        for (int i = 0; i < n_obj_; ++i) {
            double fi = (1.0 + g) * 0.5;
            for (int j = 0; j < n_obj_ - (i + 1); ++j) {
                fi *= pop.gene(idx, j);
            }
            if (i != 0) {
                fi *= (1.0 - pop.gene(idx, n_obj_ - (i + 1)));
            }
            pop.objective(idx, i) = fi;
        }
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// DTLZ2 — spherical Pareto-optimal front
// ============================================================

/// DTLZ2: spherical Pareto-optimal front.
/// Default: 12 variables, 3 objectives.
struct DTLZ2 {
    int dim_ = 12;
    int n_obj_ = 3;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    DTLZ2() : DTLZ2(12, 3) {}

    DTLZ2(int dim, int n_obj) : dim_(dim), n_obj_(n_obj) {
        if (dim_ < n_obj_) dim_ = n_obj_;
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    int num_objectives() const { return n_obj_; }
    int dimension() const { return dim_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        int k = dim_ - n_obj_ + 1;
        double g = 0.0;
        for (int i = dim_ - k; i < dim_; ++i) {
            double xi = pop.gene(idx, i);
            g += (xi - 0.5) * (xi - 0.5);
        }

        for (int i = 0; i < n_obj_; ++i) {
            double fi = 1.0 + g;
            for (int j = 0; j < n_obj_ - (i + 1); ++j) {
                fi *= std::cos(pop.gene(idx, j) * 0.5 * M_PI);
            }
            if (i != 0) {
                int aux = n_obj_ - (i + 1);
                fi *= std::sin(pop.gene(idx, aux) * 0.5 * M_PI);
            }
            pop.objective(idx, i) = fi;
        }
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

} // namespace ea
