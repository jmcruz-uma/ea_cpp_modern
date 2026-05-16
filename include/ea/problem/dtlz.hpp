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

// ============================================================
// DTLZ3 — many local Pareto-optimal fronts (multi-modal)
// ============================================================

/// DTLZ3: many local Pareto-optimal fronts (multi-modal).
/// Same structure as DTLZ2 but g uses a multi-modal function.
/// Default: 12 variables, 3 objectives.
struct DTLZ3 {
    int dim_ = 12;
    int n_obj_ = 3;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    DTLZ3() : DTLZ3(12, 3) {}

    DTLZ3(int dim, int n_obj) : dim_(dim), n_obj_(n_obj) {
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
            g += (xi - 0.5) * (xi - 0.5) - std::cos(20.0 * M_PI * (xi - 0.5));
        }
        g = 100.0 * (k + g);

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

// ============================================================
// DTLZ4 — biased density towards some objective planes
// ============================================================

/// DTLZ4: biased density towards some objective planes.
/// Same structure as DTLZ2 but x_j is raised to power alpha=100.
/// Default: 12 variables, 3 objectives.
struct DTLZ4 {
    int dim_ = 12;
    int n_obj_ = 3;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    DTLZ4() : DTLZ4(12, 3) {}

    DTLZ4(int dim, int n_obj) : dim_(dim), n_obj_(n_obj) {
        if (dim_ < n_obj_) dim_ = n_obj_;
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    int num_objectives() const { return n_obj_; }
    int dimension() const { return dim_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        constexpr double alpha = 100.0;
        int k = dim_ - n_obj_ + 1;
        double g = 0.0;
        for (int i = dim_ - k; i < dim_; ++i) {
            double xi = pop.gene(idx, i);
            g += (xi - 0.5) * (xi - 0.5);
        }

        for (int i = 0; i < n_obj_; ++i) {
            double fi = 1.0 + g;
            for (int j = 0; j < n_obj_ - (i + 1); ++j) {
                fi *= std::cos(std::pow(pop.gene(idx, j), alpha) * 0.5 * M_PI);
            }
            if (i != 0) {
                int aux = n_obj_ - (i + 1);
                fi *= std::sin(std::pow(pop.gene(idx, aux), alpha) * 0.5 * M_PI);
            }
            pop.objective(idx, i) = fi;
        }
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// DTLZ5 — Pareto front is a curve (reduced dimensionality)
// ============================================================

/// DTLZ5: Pareto front is a curve (reduced dimensionality).
/// Theta mapping reduces the effective dimensionality of the Pareto front.
/// Default: 12 variables, 3 objectives.
struct DTLZ5 {
    int dim_ = 12;
    int n_obj_ = 3;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    DTLZ5() : DTLZ5(12, 3) {}

    DTLZ5(int dim, int n_obj) : dim_(dim), n_obj_(n_obj) {
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

        double t = M_PI / (4.0 * (1.0 + g));

        std::vector<double> theta(n_obj_ - 1);
        theta[0] = pop.gene(idx, 0) * 0.5 * M_PI;
        for (int i = 1; i < n_obj_ - 1; ++i) {
            theta[i] = t * (1.0 + 2.0 * g * pop.gene(idx, i));
        }

        for (int i = 0; i < n_obj_; ++i) {
            double fi = 1.0 + g;
            for (int j = 0; j < n_obj_ - (i + 1); ++j) {
                fi *= std::cos(theta[j]);
            }
            if (i != 0) {
                int aux = n_obj_ - (i + 1);
                fi *= std::sin(theta[aux]);
            }
            pop.objective(idx, i) = fi;
        }
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// DTLZ6 — degenerate version of DTLZ5
// ============================================================

/// DTLZ6: degenerate version of DTLZ5.
/// Uses g = sum(x_i^0.1) instead of the quadratic g from DTLZ5.
/// Default: 12 variables, 3 objectives.
struct DTLZ6 {
    int dim_ = 12;
    int n_obj_ = 3;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    DTLZ6() : DTLZ6(12, 3) {}

    DTLZ6(int dim, int n_obj) : dim_(dim), n_obj_(n_obj) {
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
            g += std::pow(pop.gene(idx, i), 0.1);
        }

        double t = M_PI / (4.0 * (1.0 + g));

        std::vector<double> theta(n_obj_ - 1);
        theta[0] = pop.gene(idx, 0) * 0.5 * M_PI;
        for (int i = 1; i < n_obj_ - 1; ++i) {
            theta[i] = t * (1.0 + 2.0 * g * pop.gene(idx, i));
        }

        for (int i = 0; i < n_obj_; ++i) {
            double fi = 1.0 + g;
            for (int j = 0; j < n_obj_ - (i + 1); ++j) {
                fi *= std::cos(theta[j]);
            }
            if (i != 0) {
                int aux = n_obj_ - (i + 1);
                fi *= std::sin(theta[aux]);
            }
            pop.objective(idx, i) = fi;
        }
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// DTLZ7 — disconnected Pareto front
// ============================================================

/// DTLZ7: disconnected Pareto front.
/// Default: 22 variables, 3 objectives.
/// f_i = x_i for i = 0..M-2
/// g = 1 + (9/k) * sum(x_i for i=M-1..n-1)
/// h = M - sum( f_i/(1+g) * (1 + sin(3*pi*f_i)) for i=0..M-2 )
/// f_{M-1} = (1+g) * h
struct DTLZ7 {
    int dim_ = 22;
    int n_obj_ = 3;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    DTLZ7() : DTLZ7(22, 3) {}

    DTLZ7(int dim, int n_obj) : dim_(dim), n_obj_(n_obj) {
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

        // Copy first M-1 objectives directly from position variables
        for (int i = 0; i < n_obj_ - 1; ++i) {
            pop.objective(idx, i) = pop.gene(idx, i);
        }

        double g = 0.0;
        for (int i = dim_ - k; i < dim_; ++i) {
            g += pop.gene(idx, i);
        }
        g = 1.0 + (9.0 * g) / k;

        double h = 0.0;
        for (int i = 0; i < n_obj_ - 1; ++i) {
            double fi = pop.objective(idx, i);
            h += (fi / (1.0 + g)) * (1.0 + std::sin(3.0 * M_PI * fi));
        }
        h = static_cast<double>(n_obj_) - h;

        pop.objective(idx, n_obj_ - 1) = (1.0 + g) * h;
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

} // namespace ea
