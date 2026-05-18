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

#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <numeric>
#include <span>
#include <vector>

namespace ea {

// ============================================================
// ZDT1 — convex Pareto-optimal front
// ============================================================

/// ZDT1: convex Pareto-optimal front.
/// g(x) = 1 + 9/(n-1) * sum(x_i for i=2..n)
/// h(f1, g) = 1 - sqrt(f1/g)
struct ZDT1 {
    int dim_ = 30;                     ///< Number of decision variables (default: 30)
    std::vector<double> lower_bounds_; ///< [dim] lower bounds
    std::vector<double> upper_bounds_; ///< [dim] upper bounds

    ZDT1()
        : ZDT1(30) {}

    explicit ZDT1(int dim)
        : dim_(dim) {
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
        for (int i = start; i < start + count; ++i)
            evaluate(pop, i);
    }

protected:
    double eval_g(const Population& pop, int idx) const {
        double g = 0.0;
        for (int i = 1; i < dim_; ++i) {
            g += pop.gene(idx, i);
        }
        return 1.0 + 9.0 / (dim_ - 1) * g;
    }

    static double eval_h(double f1, double g) { return 1.0 - std::sqrt(f1 / g); }
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

    ZDT2()
        : ZDT2(30) {}

    explicit ZDT2(int dim)
        : dim_(dim) {
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
        for (int i = start; i < start + count; ++i)
            evaluate(pop, i);
    }

protected:
    double eval_g(const Population& pop, int idx) const {
        double g = 0.0;
        for (int i = 1; i < dim_; ++i) {
            g += pop.gene(idx, i);
        }
        return 1.0 + 9.0 / (dim_ - 1) * g;
    }

    static double eval_h(double f1, double g) { return 1.0 - std::pow(f1 / g, 2.0); }
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

    ZDT3()
        : ZDT3(30) {}

    explicit ZDT3(int dim)
        : dim_(dim) {
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
        for (int i = start; i < start + count; ++i)
            evaluate(pop, i);
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

// ============================================================
// ZDT4 — multi-modal ZDT variant with many local fronts
// ============================================================

/// ZDT4: multi-modal ZDT variant with many local fronts.
/// 10 variables: first in [0,1], rest in [-5,5].
/// g(x) = 1 + 10*(n-1) + sum(x_i^2 - 10*cos(4*pi*x_i) for i=2..n)
/// h(f1, g) = 1 - sqrt(f1/g)
struct ZDT4 {
    int dim_ = 10;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    ZDT4()
        : ZDT4(10) {}

    explicit ZDT4(int dim)
        : dim_(dim) {
        lower_bounds_.reserve(dim_);
        upper_bounds_.reserve(dim_);
        lower_bounds_.push_back(0.0);
        upper_bounds_.push_back(1.0);
        for (int i = 1; i < dim_; ++i) {
            lower_bounds_.push_back(-5.0);
            upper_bounds_.push_back(5.0);
        }
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
        for (int i = start; i < start + count; ++i)
            evaluate(pop, i);
    }

protected:
    double eval_g(const Population& pop, int idx) const {
        double g = 0.0;
        for (int i = 1; i < dim_; ++i) {
            double xi = pop.gene(idx, i);
            g += xi * xi - 10.0 * std::cos(4.0 * M_PI * xi);
        }
        double constant = 1.0 + 10.0 * (dim_ - 1);
        return g + constant;
    }

    static double eval_h(double f1, double g) { return 1.0 - std::sqrt(f1 / g); }
};

// ============================================================
// ZDT5 — deceptive ZDT variant (binary encoding mapped to double 0/1)
// ============================================================

/// ZDT5: deceptive ZDT variant.
/// Mapped to double encoding where values are 0 or 1.
/// 11 variables: first has 30 bits (mapped to one double), rest have 5 bits each.
/// f1 = 1 + u(x1)
/// g(x) = sum(v(u(xi))) for i=2..n
/// h(f1, g) = 1 / f1
/// where u = bit count (for us: just the double value, 0..1, scaled by bits)
/// v(y) = 2+y if y<5, else 1
///
/// Since we use real encoding, each "bit group" is a single double in [0,1].
/// u(x) = 30 * x for the first variable, 5 * x for the rest.
struct ZDT5 {
    int dim_ = 11;
    std::vector<int> bits_per_variable_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    ZDT5()
        : ZDT5(11) {}

    explicit ZDT5(int dim)
        : dim_(dim) {
        bits_per_variable_.reserve(dim_);
        bits_per_variable_.push_back(30);
        for (int i = 1; i < dim_; ++i) {
            bits_per_variable_.push_back(5);
        }
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        double u1 = u(pop.gene(idx, 0), 0);
        double f1 = 1.0 + u1;
        double g = eval_g(pop, idx);
        double h = eval_h(f1, g);
        pop.objective(idx, 0) = f1;
        pop.objective(idx, 1) = h * g;
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i)
            evaluate(pop, i);
    }

protected:
    /// u(x) = bits * x  (maps [0,1] double to expected bit count)
    double u(double x, int var_idx) const { return bits_per_variable_[var_idx] * x; }

    double eval_v(double y) const { return (y < 5.0) ? (2.0 + y) : 1.0; }

    double eval_g(const Population& pop, int idx) const {
        double g = 0.0;
        for (int i = 1; i < dim_; ++i) {
            double ui = u(pop.gene(idx, i), i);
            g += eval_v(ui);
        }
        return g;
    }

    static double eval_h(double f1, double /*g*/) { return 1.0 / f1; }
};

// ============================================================
// ZDT6 — non-uniform ZDT variant with non-convex Pareto front
// ============================================================

/// ZDT6: non-uniform ZDT variant with non-convex Pareto front.
/// f1(x1) = 1 - exp(-4*x1) * sin^6(6*pi*x1)
/// g(x) = 1 + 9 * [ (sum(x_i for i=2..n) / (n-1)) ^ 0.25 ]
/// h(f1, g) = 1 - (f1/g)^2
struct ZDT6 {
    int dim_ = 10;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    ZDT6()
        : ZDT6(10) {}

    explicit ZDT6(int dim)
        : dim_(dim) {
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        double x1 = pop.gene(idx, 0);
        double f1 = 1.0 - std::exp(-4.0 * x1) * std::pow(std::sin(6.0 * M_PI * x1), 6);
        double g = eval_g(pop, idx);
        double h = eval_h(f1, g);
        pop.objective(idx, 0) = f1;
        pop.objective(idx, 1) = h * g;
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i)
            evaluate(pop, i);
    }

protected:
    double eval_g(const Population& pop, int idx) const {
        double g = 0.0;
        for (int i = 1; i < dim_; ++i) {
            g += pop.gene(idx, i);
        }
        g = g / (dim_ - 1);
        g = std::pow(g, 0.25);
        g = 9.0 * g;
        g = 1.0 + g;
        return g;
    }

    static double eval_h(double f1, double g) { return 1.0 - std::pow(f1 / g, 2.0); }
};

} // namespace ea
