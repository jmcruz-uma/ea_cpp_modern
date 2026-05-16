#pragma once
/// @file zdt.hpp
/// @brief ZDT benchmark problem suite for multi-objective optimization.
///
/// Reference: jMetal jmetal-problem/src/main/java/org/uma/jmetal/problem/multiobjective/zdt/
///
/// The ZDT (Zitzler-Deb-Thiele) test suite contains problems with two objectives
/// and a known Pareto-optimal front. All problems share a common structure:
///   f1(x) = x1
///   g(x)  = function of remaining variables
///   h(f1, g) = shape function defining the front
///   f2(x) = h(f1, g) * g(x)
///
/// Bounds: [0, 1] for all variables.

#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <span>
#include <cmath>
#include <vector>
#include <numeric>

namespace ea {

// ============================================================
// ZDT1 — convex Pareto-optimal front
// ============================================================

/// ZDT1: convex Pareto-optimal front.
/// g(x) = 1 + 9/(n-1) * sum(x_i for i=2..n)
/// h(f1, g) = 1 - sqrt(f1/g)
struct ZDT1 {
    int dim_ = 30; ///< Number of decision variables (default: 30)
    std::vector<double> lower_bounds_; ///< [dim] lower bounds
    std::vector<double> upper_bounds_; ///< [dim] upper bounds

    ZDT1() : ZDT1(30) {}

    explicit ZDT1(int dim) : dim_(dim) {
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        double g = eval_g(pop, idx);
        double f1 = pop.gene(idx, 0);
        double h = eval_h(f1, g);
        pop.objective(idx, 0) = f1;
        pop.objective(idx, 1) = h * g;
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }

protected:
    double eval_g(const Population& pop, int idx) const {
        double g = 0.0;
        for (int i = 1; i < dim_; ++i) {
            g += pop.gene(idx, i);
        }
        return 1.0 + 9.0 / (dim_ - 1) * g;
    }

    static double eval_h(double f1, double g) {
        return 1.0 - std::sqrt(f1 / g);
    }
};

// ============================================================
// ZDT2 — non-convex Pareto-optimal front
// ============================================================

/// ZDT2: non-convex Pareto-optimal front.
/// Same g(x) as ZDT1, but h(f1, g) = 1 - (f1/g)^2.
struct ZDT2 {
    int dim_ = 30;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    ZDT2() : ZDT2(30) {}

    explicit ZDT2(int dim) : dim_(dim) {
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        double g = eval_g(pop, idx);
        double f1 = pop.gene(idx, 0);
        double h = eval_h(f1, g);
        pop.objective(idx, 0) = f1;
        pop.objective(idx, 1) = h * g;
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }

protected:
    double eval_g(const Population& pop, int idx) const {
        double g = 0.0;
        for (int i = 1; i < dim_; ++i) {
            g += pop.gene(idx, i);
        }
        return 1.0 + 9.0 / (dim_ - 1) * g;
    }

    static double eval_h(double f1, double g) {
        return 1.0 - std::pow(f1 / g, 2.0);
    }
};

// ============================================================
// ZDT3 — disconnected Pareto-optimal front
// ============================================================

/// ZDT3: disconnected Pareto-optimal front.
/// Same g(x) as ZDT1, but h(f1, g) = 1 - sqrt(f1/g) - (f1/g)*sin(10*pi*f1).
struct ZDT3 {
    int dim_ = 30;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    ZDT3() : ZDT3(30) {}

    explicit ZDT3(int dim) : dim_(dim) {
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        double g = eval_g(pop, idx);
        double f1 = pop.gene(idx, 0);
        double h = eval_h(f1, g);
        pop.objective(idx, 0) = f1;
        pop.objective(idx, 1) = h * g;
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }

protected:
    double eval_g(const Population& pop, int idx) const {
        double g = 0.0;
        for (int i = 1; i < dim_; ++i) {
            g += pop.gene(idx, i);
        }
        return 1.0 + 9.0 / (dim_ - 1) * g;
    }

    static double eval_h(double f1, double g) {
        double x = f1 / g;
        return 1.0 - std::sqrt(x) - x * std::sin(10.0 * M_PI * f1);
    }
};

} // namespace ea
