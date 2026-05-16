#pragma once
/// @file wfg.hpp
/// @brief WFG benchmark problem suite for multi-objective optimization.
///
/// Reference: jMetal jmetal-problem/src/main/java/org/uma/jmetal/problem/multiobjective/wfg/
///
/// The WFG (Walking Fish Group) test suite is scalable in both number of
/// objectives and number of variables. Each problem uses:
/// - k position-related parameters (default: 2*(M-1))
/// - l distance-related parameters (default: 20)
/// - Variables are bounded to [0, 2*(i+1)] for variable i.
///
/// The evaluation pipeline:
/// 1. Normalize variables to [0, 1]
/// 2. Apply problem-specific transformations t1, t2, ...
/// 3. Calculate x vector from transformed y
/// 4. Apply shape functions to get objectives.

#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <span>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cassert>

namespace ea {

// ============================================================
// WFG transformation primitives
// ============================================================

namespace wfg_transform {

constexpr double EPSILON = 1.0e-10;

/// Clamp a value to [0, 1] with epsilon tolerance.
inline double correct_to_01(double a) {
    constexpr double min = 0.0;
    constexpr double max = 1.0;
    constexpr double min_epsilon = min - EPSILON;
    constexpr double max_epsilon = max + EPSILON;

    if ((a <= min && a >= min_epsilon) || (a >= min && a <= min_epsilon)) {
        return min;
    }
    if ((a >= max && a <= max_epsilon) || (a <= max && a >= max_epsilon)) {
        return max;
    }
    return a;
}

/// bPoly: power transformation. alpha > 0.
inline double b_poly(double y, double alpha) {
    return correct_to_01(std::pow(y, alpha));
}

/// bFlat: flat region transformation.
inline double b_flat(double y, double A, double B, double C) {
    double tmp1 = std::min(0.0, std::floor(y - B)) * A * (B - y) / B;
    double tmp2 = std::min(0.0, std::floor(C - y)) * (1.0 - A) * (y - C) / (1.0 - C);
    return correct_to_01(A + tmp1 - tmp2);
}

/// sLinear: linear transformation.
inline double s_linear(double y, double A) {
    return correct_to_01(std::abs(y - A) / std::abs(std::floor(A - y) + A));
}

/// sDecept: deceptive transformation.
inline double s_decept(double y, double A, double B, double C) {
    double tmp1 = std::floor(y - A + B) * (1.0 - C + (A - B) / B) / (A - B);
    double tmp2 = std::floor(A + B - y) * (1.0 - C + (1.0 - A - B) / B) / (1.0 - A - B);
    double tmp = std::abs(y - A) - B;
    return correct_to_01(1.0 + tmp * (tmp1 + tmp2 + 1.0 / B));
}

/// sMulti: multi-modal transformation.
inline double s_multi(double y, int A, int B, double C) {
    double denom = 2.0 * (std::floor(C - y) + C);
    double t = (0.5 - std::abs(y - C) / denom);
    double tmp1 = (4.0 * A + 2.0) * M_PI * t;
    double tmp2 = 4.0 * B * std::pow(std::abs(y - C) / denom, 2.0);
    return correct_to_01((1.0 + std::cos(tmp1) + tmp2) / (B + 2.0));
}

/// rSum: weighted sum reduction.
inline double r_sum(const std::vector<double>& y, const std::vector<double>& w) {
    assert(y.size() == w.size());
    double tmp1 = 0.0, tmp2 = 0.0;
    for (size_t i = 0; i < y.size(); ++i) {
        tmp1 += y[i] * w[i];
        tmp2 += w[i];
    }
    return correct_to_01(tmp1 / tmp2);
}

/// rNonsep: non-separable reduction.
inline double r_nonsep(const std::vector<double>& y, int A) {
    double tmp = std::ceil(A / 2.0);
    double denominator = y.size() * tmp * (1.0 + 2.0 * A - 2.0 * tmp) / A;
    double numerator = 0.0;
    for (size_t j = 0; j < y.size(); ++j) {
        numerator += y[j];
        for (int k = 0; k <= A - 2; ++k) {
            numerator += std::abs(y[j] - y[(j + k + 1) % y.size()]);
        }
    }
    return correct_to_01(numerator / denominator);
}

/// bParam: parameter-dependent transformation.
inline double b_param(double y, double u, double A, double B, double C) {
    double v = A - (1.0 - 2.0 * u) * std::abs(std::floor(0.5 - u) + A);
    double exp = B + (C - B) * v;
    return correct_to_01(std::pow(y, exp));
}

} // namespace wfg_transform

// ============================================================
// WFG shape functions
// ============================================================

namespace wfg_shape {

/// Linear shape function. m is 1-based index.
inline double linear(const std::vector<double>& x, int m) {
    double result = 1.0;
    int M = static_cast<int>(x.size());
    for (int i = 1; i <= M - m; ++i) {
        result *= x[i - 1];
    }
    if (m != 1) {
        result *= (1.0 - x[M - m]);
    }
    return result;
}

/// Convex shape function. m is 1-based index.
inline double convex(const std::vector<double>& x, int m) {
    double result = 1.0;
    int M = static_cast<int>(x.size());
    for (int i = 1; i <= M - m; ++i) {
        result *= (1.0 - std::cos(x[i - 1] * M_PI * 0.5));
    }
    if (m != 1) {
        result *= (1.0 - std::sin(x[M - m] * M_PI * 0.5));
    }
    return result;
}

/// Concave shape function. m is 1-based index.
inline double concave(const std::vector<double>& x, int m) {
    double result = 1.0;
    int M = static_cast<int>(x.size());
    for (int i = 1; i <= M - m; ++i) {
        result *= std::sin(x[i - 1] * M_PI * 0.5);
    }
    if (m != 1) {
        result *= std::cos(x[M - m] * M_PI * 0.5);
    }
    return result;
}

/// Mixed shape function (used by WFG1 last objective).
inline double mixed(const std::vector<double>& x, int A, double alpha) {
    double tmp = std::cos(2.0 * A * M_PI * x[0] + M_PI * 0.5) / (2.0 * A * M_PI);
    return std::pow(1.0 - x[0] - tmp, alpha);
}

/// Disc shape function (used by WFG2 last objective).
inline double disc(const std::vector<double>& x, int A, double alpha, double beta) {
    double tmp = std::cos(A * std::pow(x[0], beta) * M_PI);
    return 1.0 - std::pow(x[0], alpha) * tmp * tmp;
}

} // namespace wfg_shape

// ============================================================
// WFG base helper functions
// ============================================================

namespace wfg_detail {

/// Normalize variables. Bound for variable i is 2*(i+1).
inline std::vector<double> normalise(const std::vector<double>& z) {
    std::vector<double> result(z.size());
    for (size_t i = 0; i < z.size(); ++i) {
        double bound = 2.0 * (i + 1);
        result[i] = z[i] / bound;
        result[i] = wfg_transform::correct_to_01(result[i]);
    }
    return result;
}

/// Calculate x vector from transformed t vector.
/// a[i] values are 0 or 1. m is number of objectives.
inline std::vector<double> calculate_x(const std::vector<double>& t,
                                          const std::vector<int>& a,
                                          int m) {
    std::vector<double> x(m);
    for (int i = 0; i < m - 1; ++i) {
        x[i] = std::max(t[m - 1], static_cast<double>(a[i])) * (t[i] - 0.5) + 0.5;
    }
    x[m - 1] = t[m - 1];
    return x;
}

/// Default s vector: s[i] = 2*(i+1) for 1-based indexing → s[0] = 2, s[1] = 4, ...
inline std::vector<int> default_s(int m) {
    std::vector<int> s(m);
    for (int i = 0; i < m; ++i) {
        s[i] = 2 * (i + 1);
    }
    return s;
}

} // namespace wfg_detail

// ============================================================
// WFG1 — convex, mixed
// ============================================================

struct WFG1 {
    int k_;  ///< position variables
    int l_;  ///< distance variables
    int m_;  ///< objectives
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;
    std::vector<int> s_;
    std::vector<int> a_;
    int d_ = 1;

    WFG1(int m = 3, int k = -1, int l = 20)
        : l_(l), m_(m) {
        k_ = (k < 0) ? 2 * (m_ - 1) : k;
        int dim = k_ + l_;
        lower_bounds_.resize(dim);
        upper_bounds_.resize(dim);
        for (int i = 0; i < dim; ++i) {
            lower_bounds_[i] = 0.0;
            upper_bounds_[i] = 2.0 * (i + 1);
        }
        s_ = wfg_detail::default_s(m_);
        a_.assign(m_ - 1, 1);
    }

    int num_objectives() const { return m_; }
    int dimension() const { return k_ + l_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        int n = k_ + l_;
        std::vector<double> z(n);
        for (int i = 0; i < n; ++i) {
            z[i] = pop.gene(idx, i);
        }

        std::vector<double> y = wfg_detail::normalise(z);

        // t1: sLinear on distance variables
        for (int i = k_; i < n; ++i) {
            y[i] = wfg_transform::s_linear(y[i], 0.35);
        }

        // t2: bFlat on distance variables
        for (int i = k_; i < n; ++i) {
            y[i] = wfg_transform::b_flat(y[i], 0.8, 0.75, 0.85);
        }

        // t3: bPoly on all variables
        for (int i = 0; i < n; ++i) {
            y[i] = wfg_transform::b_poly(y[i], 0.02);
        }

        // t4: rSum reduction to M variables
        std::vector<double> t(m_);
        std::vector<double> w(n);
        for (int i = 0; i < n; ++i) {
            w[i] = 2.0 * (i + 1);
        }
        for (int i = 1; i <= m_ - 1; ++i) {
            int head = (i - 1) * k_ / (m_ - 1);
            int tail = i * k_ / (m_ - 1) - 1;
            std::vector<double> sub_y(y.begin() + head, y.begin() + tail + 1);
            std::vector<double> sub_w(w.begin() + head, w.begin() + tail + 1);
            t[i - 1] = wfg_transform::r_sum(sub_y, sub_w);
        }
        {
            std::vector<double> sub_y(y.begin() + k_, y.end());
            std::vector<double> sub_w(w.begin() + k_, w.end());
            t[m_ - 1] = wfg_transform::r_sum(sub_y, sub_w);
        }

        std::vector<double> x = wfg_detail::calculate_x(t, a_, m_);

        // Shape functions
        for (int i = 1; i <= m_ - 1; ++i) {
            pop.objective(idx, i - 1) = d_ * x[m_ - 1] + s_[i - 1] * wfg_shape::convex(x, i);
        }
        pop.objective(idx, m_ - 1) = d_ * x[m_ - 1] + s_[m_ - 1] * wfg_shape::mixed(x, 5, 1.0);
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// WFG2 — convex, disconnected
// ============================================================

struct WFG2 {
    int k_;
    int l_;
    int m_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;
    std::vector<int> s_;
    std::vector<int> a_;
    int d_ = 1;

    WFG2(int m = 3, int k = -1, int l = 20)
        : l_(l), m_(m) {
        k_ = (k < 0) ? 2 * (m_ - 1) : k;
        int dim = k_ + l_;
        lower_bounds_.resize(dim);
        upper_bounds_.resize(dim);
        for (int i = 0; i < dim; ++i) {
            lower_bounds_[i] = 0.0;
            upper_bounds_[i] = 2.0 * (i + 1);
        }
        s_ = wfg_detail::default_s(m_);
        a_.assign(m_ - 1, 1);
    }

    int num_objectives() const { return m_; }
    int dimension() const { return k_ + l_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        int n = k_ + l_;
        std::vector<double> z(n);
        for (int i = 0; i < n; ++i) {
            z[i] = pop.gene(idx, i);
        }

        std::vector<double> y = wfg_detail::normalise(z);

        // t1: sLinear on distance variables
        for (int i = k_; i < n; ++i) {
            y[i] = wfg_transform::s_linear(y[i], 0.35);
        }

        // t2: rNonsep on pairs of distance variables
        int l_prime = n - k_;
        std::vector<double> y2(n);
        for (int i = 0; i < k_; ++i) {
            y2[i] = y[i];
        }
        for (int i = 1; i <= l_prime / 2; ++i) {
            int head = k_ + 2 * (i - 1);
            int tail = k_ + 2 * (i - 1) + 1;
            std::vector<double> sub_y{y[head], y[tail]};
            y2[k_ + i - 1] = wfg_transform::r_nonsep(sub_y, 2);
        }
        y2.resize(k_ + l_prime / 2);

        // t3: rSum reduction to M variables
        int new_n = static_cast<int>(y2.size());
        std::vector<double> t(m_);
        std::vector<double> w(new_n, 1.0);
        for (int i = 1; i <= m_ - 1; ++i) {
            int head = (i - 1) * k_ / (m_ - 1);
            int tail = i * k_ / (m_ - 1) - 1;
            std::vector<double> sub_y(y2.begin() + head, y2.begin() + tail + 1);
            std::vector<double> sub_w(w.begin() + head, w.begin() + tail + 1);
            t[i - 1] = wfg_transform::r_sum(sub_y, sub_w);
        }
        {
            std::vector<double> sub_y(y2.begin() + k_, y2.end());
            std::vector<double> sub_w(w.begin() + k_, w.end());
            t[m_ - 1] = wfg_transform::r_sum(sub_y, sub_w);
        }

        std::vector<double> x = wfg_detail::calculate_x(t, a_, m_);

        for (int i = 1; i <= m_ - 1; ++i) {
            pop.objective(idx, i - 1) = d_ * x[m_ - 1] + s_[i - 1] * wfg_shape::convex(x, i);
        }
        pop.objective(idx, m_ - 1) = d_ * x[m_ - 1] + s_[m_ - 1] * wfg_shape::disc(x, 5, 1.0, 1.0);
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// WFG3 — linear, degenerate
// ============================================================

struct WFG3 {
    int k_;
    int l_;
    int m_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;
    std::vector<int> s_;
    std::vector<int> a_;
    int d_ = 1;

    WFG3(int m = 3, int k = -1, int l = 20)
        : l_(l), m_(m) {
        k_ = (k < 0) ? 2 * (m_ - 1) : k;
        int dim = k_ + l_;
        lower_bounds_.resize(dim);
        upper_bounds_.resize(dim);
        for (int i = 0; i < dim; ++i) {
            lower_bounds_[i] = 0.0;
            upper_bounds_[i] = 2.0 * (i + 1);
        }
        s_ = wfg_detail::default_s(m_);
        a_.assign(m_ - 1, 0);
        if (!a_.empty()) a_[0] = 1;
    }

    int num_objectives() const { return m_; }
    int dimension() const { return k_ + l_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        int n = k_ + l_;
        std::vector<double> z(n);
        for (int i = 0; i < n; ++i) {
            z[i] = pop.gene(idx, i);
        }

        std::vector<double> y = wfg_detail::normalise(z);

        // t1: sLinear on distance variables
        for (int i = k_; i < n; ++i) {
            y[i] = wfg_transform::s_linear(y[i], 0.35);
        }

        // t2: rNonsep on pairs of distance variables
        int l_prime = n - k_;
        std::vector<double> y2(n);
        for (int i = 0; i < k_; ++i) {
            y2[i] = y[i];
        }
        for (int i = 1; i <= l_prime / 2; ++i) {
            int head = k_ + 2 * (i - 1);
            int tail = head + 1;
            std::vector<double> sub_y{y[head], y[tail]};
            y2[k_ + i - 1] = wfg_transform::r_nonsep(sub_y, 2);
        }
        y2.resize(k_ + l_prime / 2);

        // t3: rSum reduction to M variables
        int new_n = static_cast<int>(y2.size());
        std::vector<double> t(m_);
        std::vector<double> w(new_n, 1.0);
        for (int i = 1; i <= m_ - 1; ++i) {
            int head = (i - 1) * k_ / (m_ - 1);
            int tail = i * k_ / (m_ - 1) - 1;
            std::vector<double> sub_y(y2.begin() + head, y2.begin() + tail + 1);
            std::vector<double> sub_w(w.begin() + head, w.begin() + tail + 1);
            t[i - 1] = wfg_transform::r_sum(sub_y, sub_w);
        }
        {
            std::vector<double> sub_y(y2.begin() + k_, y2.end());
            std::vector<double> sub_w(w.begin() + k_, w.end());
            t[m_ - 1] = wfg_transform::r_sum(sub_y, sub_w);
        }

        std::vector<double> x = wfg_detail::calculate_x(t, a_, m_);

        for (int i = 1; i <= m_; ++i) {
            pop.objective(idx, i - 1) = d_ * x[m_ - 1] + s_[i - 1] * wfg_shape::linear(x, i);
        }
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// WFG4 — multi-modal, separable
// ============================================================

struct WFG4 {
    int k_;
    int l_;
    int m_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;
    std::vector<int> s_;
    std::vector<int> a_;
    int d_ = 1;

    WFG4(int m = 3, int k = -1, int l = 20)
        : l_(l), m_(m) {
        k_ = (k < 0) ? 2 * (m_ - 1) : k;
        int dim = k_ + l_;
        lower_bounds_.resize(dim);
        upper_bounds_.resize(dim);
        for (int i = 0; i < dim; ++i) {
            lower_bounds_[i] = 0.0;
            upper_bounds_[i] = 2.0 * (i + 1);
        }
        s_ = wfg_detail::default_s(m_);
        a_.assign(m_ - 1, 1);
    }

    int num_objectives() const { return m_; }
    int dimension() const { return k_ + l_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        int n = k_ + l_;
        std::vector<double> z(n);
        for (int i = 0; i < n; ++i) {
            z[i] = pop.gene(idx, i);
        }

        std::vector<double> y = wfg_detail::normalise(z);

        // t1: sMulti on all variables
        for (int i = 0; i < n; ++i) {
            y[i] = wfg_transform::s_multi(y[i], 30, 10, 0.35);
        }

        // t2: rSum reduction to M variables
        std::vector<double> t(m_);
        std::vector<double> w(n, 1.0);
        for (int i = 1; i <= m_ - 1; ++i) {
            int head = (i - 1) * k_ / (m_ - 1);
            int tail = i * k_ / (m_ - 1) - 1;
            std::vector<double> sub_y(y.begin() + head, y.begin() + tail + 1);
            std::vector<double> sub_w(w.begin() + head, w.begin() + tail + 1);
            t[i - 1] = wfg_transform::r_sum(sub_y, sub_w);
        }
        {
            std::vector<double> sub_y(y.begin() + k_, y.end());
            std::vector<double> sub_w(w.begin() + k_, w.end());
            t[m_ - 1] = wfg_transform::r_sum(sub_y, sub_w);
        }

        std::vector<double> x = wfg_detail::calculate_x(t, a_, m_);

        for (int i = 1; i <= m_; ++i) {
            pop.objective(idx, i - 1) = d_ * x[m_ - 1] + s_[i - 1] * wfg_shape::concave(x, i);
        }
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// WFG5 — deceptive, separable
// ============================================================

struct WFG5 {
    int k_;
    int l_;
    int m_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;
    std::vector<int> s_;
    std::vector<int> a_;
    int d_ = 1;

    WFG5(int m = 3, int k = -1, int l = 20)
        : l_(l), m_(m) {
        k_ = (k < 0) ? 2 * (m_ - 1) : k;
        int dim = k_ + l_;
        lower_bounds_.resize(dim);
        upper_bounds_.resize(dim);
        for (int i = 0; i < dim; ++i) {
            lower_bounds_[i] = 0.0;
            upper_bounds_[i] = 2.0 * (i + 1);
        }
        s_ = wfg_detail::default_s(m_);
        a_.assign(m_ - 1, 1);
    }

    int num_objectives() const { return m_; }
    int dimension() const { return k_ + l_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        int n = k_ + l_;
        std::vector<double> z(n);
        for (int i = 0; i < n; ++i) {
            z[i] = pop.gene(idx, i);
        }

        std::vector<double> y = wfg_detail::normalise(z);

        // t1: sDecept on all variables
        for (int i = 0; i < n; ++i) {
            y[i] = wfg_transform::s_decept(y[i], 0.35, 0.001, 0.05);
        }

        // t2: rSum reduction to M variables
        std::vector<double> t(m_);
        std::vector<double> w(n, 1.0);
        for (int i = 1; i <= m_ - 1; ++i) {
            int head = (i - 1) * k_ / (m_ - 1);
            int tail = i * k_ / (m_ - 1) - 1;
            std::vector<double> sub_y(y.begin() + head, y.begin() + tail + 1);
            std::vector<double> sub_w(w.begin() + head, w.begin() + tail + 1);
            t[i - 1] = wfg_transform::r_sum(sub_y, sub_w);
        }
        {
            std::vector<double> sub_y(y.begin() + k_, y.end());
            std::vector<double> sub_w(w.begin() + k_, w.end());
            t[m_ - 1] = wfg_transform::r_sum(sub_y, sub_w);
        }

        std::vector<double> x = wfg_detail::calculate_x(t, a_, m_);

        for (int i = 1; i <= m_; ++i) {
            pop.objective(idx, i - 1) = d_ * x[m_ - 1] + s_[i - 1] * wfg_shape::concave(x, i);
        }
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// WFG6 — non-separable
// ============================================================

struct WFG6 {
    int k_;
    int l_;
    int m_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;
    std::vector<int> s_;
    std::vector<int> a_;
    int d_ = 1;

    WFG6(int m = 3, int k = -1, int l = 20)
        : l_(l), m_(m) {
        k_ = (k < 0) ? 2 * (m_ - 1) : k;
        int dim = k_ + l_;
        lower_bounds_.resize(dim);
        upper_bounds_.resize(dim);
        for (int i = 0; i < dim; ++i) {
            lower_bounds_[i] = 0.0;
            upper_bounds_[i] = 2.0 * (i + 1);
        }
        s_ = wfg_detail::default_s(m_);
        a_.assign(m_ - 1, 1);
    }

    int num_objectives() const { return m_; }
    int dimension() const { return k_ + l_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        int n = k_ + l_;
        std::vector<double> z(n);
        for (int i = 0; i < n; ++i) {
            z[i] = pop.gene(idx, i);
        }

        std::vector<double> y = wfg_detail::normalise(z);

        // t1: sLinear on distance variables
        for (int i = k_; i < n; ++i) {
            y[i] = wfg_transform::s_linear(y[i], 0.35);
        }

        // t2: rNonsep reduction to M variables
        std::vector<double> t(m_);
        for (int i = 1; i <= m_ - 1; ++i) {
            int head = (i - 1) * k_ / (m_ - 1);
            int tail = i * k_ / (m_ - 1) - 1;
            int group_size = tail - head + 1;
            std::vector<double> sub_y(y.begin() + head, y.begin() + tail + 1);
            t[i - 1] = wfg_transform::r_nonsep(sub_y, group_size);
        }
        {
            std::vector<double> sub_y(y.begin() + k_, y.end());
            t[m_ - 1] = wfg_transform::r_nonsep(sub_y, n - k_);
        }

        std::vector<double> x = wfg_detail::calculate_x(t, a_, m_);

        for (int i = 1; i <= m_; ++i) {
            pop.objective(idx, i - 1) = d_ * x[m_ - 1] + s_[i - 1] * wfg_shape::concave(x, i);
        }
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// WFG7 — biased parameter, non-separable
// ============================================================

struct WFG7 {
    int k_;
    int l_;
    int m_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;
    std::vector<int> s_;
    std::vector<int> a_;
    int d_ = 1;

    WFG7(int m = 3, int k = -1, int l = 20)
        : l_(l), m_(m) {
        k_ = (k < 0) ? 2 * (m_ - 1) : k;
        int dim = k_ + l_;
        lower_bounds_.resize(dim);
        upper_bounds_.resize(dim);
        for (int i = 0; i < dim; ++i) {
            lower_bounds_[i] = 0.0;
            upper_bounds_[i] = 2.0 * (i + 1);
        }
        s_ = wfg_detail::default_s(m_);
        a_.assign(m_ - 1, 1);
    }

    int num_objectives() const { return m_; }
    int dimension() const { return k_ + l_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        int n = k_ + l_;
        std::vector<double> z(n);
        for (int i = 0; i < n; ++i) {
            z[i] = pop.gene(idx, i);
        }

        std::vector<double> y = wfg_detail::normalise(z);

        // t1: bParam on position variables
        std::vector<double> w(n, 1.0);
        for (int i = 0; i < k_; ++i) {
            std::vector<double> sub_y(y.begin() + i + 1, y.end());
            std::vector<double> sub_w(w.begin() + i + 1, w.end());
            double aux = wfg_transform::r_sum(sub_y, sub_w);
            y[i] = wfg_transform::b_param(y[i], aux, 0.98 / 49.98, 0.02, 50);
        }

        // t2: sLinear on distance variables
        for (int i = k_; i < n; ++i) {
            y[i] = wfg_transform::s_linear(y[i], 0.35);
        }

        // t3: rSum reduction to M variables
        std::vector<double> t(m_);
        for (int i = 1; i <= m_ - 1; ++i) {
            int head = (i - 1) * k_ / (m_ - 1);
            int tail = i * k_ / (m_ - 1) - 1;
            std::vector<double> sub_y(y.begin() + head, y.begin() + tail + 1);
            std::vector<double> sub_w(w.begin() + head, w.begin() + tail + 1);
            t[i - 1] = wfg_transform::r_sum(sub_y, sub_w);
        }
        {
            std::vector<double> sub_y(y.begin() + k_, y.end());
            std::vector<double> sub_w(w.begin() + k_, w.end());
            t[m_ - 1] = wfg_transform::r_sum(sub_y, sub_w);
        }

        std::vector<double> x = wfg_detail::calculate_x(t, a_, m_);

        for (int i = 1; i <= m_; ++i) {
            pop.objective(idx, i - 1) = d_ * x[m_ - 1] + s_[i - 1] * wfg_shape::concave(x, i);
        }
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// WFG8 — biased parameter, non-separable, asymmetric
// ============================================================

struct WFG8 {
    int k_;
    int l_;
    int m_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;
    std::vector<int> s_;
    std::vector<int> a_;
    int d_ = 1;

    WFG8(int m = 3, int k = -1, int l = 20)
        : l_(l), m_(m) {
        k_ = (k < 0) ? 2 * (m_ - 1) : k;
        int dim = k_ + l_;
        lower_bounds_.resize(dim);
        upper_bounds_.resize(dim);
        for (int i = 0; i < dim; ++i) {
            lower_bounds_[i] = 0.0;
            upper_bounds_[i] = 2.0 * (i + 1);
        }
        s_ = wfg_detail::default_s(m_);
        a_.assign(m_ - 1, 1);
    }

    int num_objectives() const { return m_; }
    int dimension() const { return k_ + l_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        int n = k_ + l_;
        std::vector<double> z(n);
        for (int i = 0; i < n; ++i) {
            z[i] = pop.gene(idx, i);
        }

        std::vector<double> y = wfg_detail::normalise(z);

        // t1: bParam on distance variables
        std::vector<double> w(n, 1.0);
        for (int i = k_; i < n; ++i) {
            std::vector<double> sub_y(y.begin(), y.begin() + i);
            std::vector<double> sub_w(w.begin(), w.begin() + i);
            double aux = wfg_transform::r_sum(sub_y, sub_w);
            y[i] = wfg_transform::b_param(y[i], aux, 0.98 / 49.98, 0.02, 50);
        }

        // t2: sLinear on distance variables
        for (int i = k_; i < n; ++i) {
            y[i] = wfg_transform::s_linear(y[i], 0.35);
        }

        // t3: rSum reduction to M variables
        std::vector<double> t(m_);
        for (int i = 1; i <= m_ - 1; ++i) {
            int head = (i - 1) * k_ / (m_ - 1);
            int tail = i * k_ / (m_ - 1) - 1;
            std::vector<double> sub_y(y.begin() + head, y.begin() + tail + 1);
            std::vector<double> sub_w(w.begin() + head, w.begin() + tail + 1);
            t[i - 1] = wfg_transform::r_sum(sub_y, sub_w);
        }
        {
            std::vector<double> sub_y(y.begin() + k_, y.end());
            std::vector<double> sub_w(w.begin() + k_, w.end());
            t[m_ - 1] = wfg_transform::r_sum(sub_y, sub_w);
        }

        std::vector<double> x = wfg_detail::calculate_x(t, a_, m_);

        for (int i = 1; i <= m_; ++i) {
            pop.objective(idx, i - 1) = d_ * x[m_ - 1] + s_[i - 1] * wfg_shape::concave(x, i);
        }
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// WFG9 — deceptive, biased parameter, non-separable
// ============================================================

struct WFG9 {
    int k_;
    int l_;
    int m_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;
    std::vector<int> s_;
    std::vector<int> a_;
    int d_ = 1;

    WFG9(int m = 3, int k = -1, int l = 20)
        : l_(l), m_(m) {
        k_ = (k < 0) ? 2 * (m_ - 1) : k;
        int dim = k_ + l_;
        lower_bounds_.resize(dim);
        upper_bounds_.resize(dim);
        for (int i = 0; i < dim; ++i) {
            lower_bounds_[i] = 0.0;
            upper_bounds_[i] = 2.0 * (i + 1);
        }
        s_ = wfg_detail::default_s(m_);
        a_.assign(m_ - 1, 1);
    }

    int num_objectives() const { return m_; }
    int dimension() const { return k_ + l_; }

    std::span<const double> lower_bounds() const { return std::span(lower_bounds_); }
    std::span<const double> upper_bounds() const { return std::span(upper_bounds_); }

    void evaluate(Population& pop, int idx) const {
        int n = k_ + l_;
        std::vector<double> z(n);
        for (int i = 0; i < n; ++i) {
            z[i] = pop.gene(idx, i);
        }

        std::vector<double> y = wfg_detail::normalise(z);

        // t1: bParam on all but last variable
        std::vector<double> w(n, 1.0);
        for (int i = 0; i < n - 1; ++i) {
            std::vector<double> sub_y(y.begin() + i + 1, y.end());
            std::vector<double> sub_w(w.begin() + i + 1, w.end());
            double aux = wfg_transform::r_sum(sub_y, sub_w);
            y[i] = wfg_transform::b_param(y[i], aux, 0.98 / 49.98, 0.02, 50);
        }

        // t2: sDecept on position, sMulti on distance
        for (int i = 0; i < k_; ++i) {
            y[i] = wfg_transform::s_decept(y[i], 0.35, 0.001, 0.05);
        }
        for (int i = k_; i < n; ++i) {
            y[i] = wfg_transform::s_multi(y[i], 30, 95, 0.35);
        }

        // t3: rNonsep reduction to M variables
        std::vector<double> t(m_);
        for (int i = 1; i <= m_ - 1; ++i) {
            int head = (i - 1) * k_ / (m_ - 1);
            int tail = i * k_ / (m_ - 1) - 1;
            int group_size = tail - head + 1;
            std::vector<double> sub_y(y.begin() + head, y.begin() + tail + 1);
            t[i - 1] = wfg_transform::r_nonsep(sub_y, group_size);
        }
        {
            std::vector<double> sub_y(y.begin() + k_, y.end());
            t[m_ - 1] = wfg_transform::r_nonsep(sub_y, n - k_);
        }

        std::vector<double> x = wfg_detail::calculate_x(t, a_, m_);

        for (int i = 1; i <= m_; ++i) {
            pop.objective(idx, i - 1) = d_ * x[m_ - 1] + s_[i - 1] * wfg_shape::concave(x, i);
        }
    }

    void evaluate_batch(Population& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

} // namespace ea
