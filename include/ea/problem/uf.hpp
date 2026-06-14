#pragma once
/// @file uf.hpp
/// @brief UF benchmark problem suite (CEC 2009) for multi-objective optimization.
///
/// Reference: jMetal jmetal-problem/src/main/java/org/uma/jmetal/problem/multiobjective/uf/
///
/// UF1–UF7: 2 objectives, 30 variables (default).
/// UF8–UF10: 3 objectives, 30 variables (default).
///
/// Bounds: x[0] in [0,1]; remaining in [-1,1], [-2,2], or [0,1] depending on problem.
///
/// Audited against jMetal 7.4 — all formulas verified correct.

#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <span>
#include <vector>

namespace ea {

// ============================================================
// UF1 — 2 objectives
// ============================================================

struct UF1 {
    int dim_ = 30;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    UF1() : UF1(30) {}
    explicit UF1(int dim) : dim_(dim) {
        lower_bounds_.reserve(dim_);
        upper_bounds_.reserve(dim_);
        lower_bounds_.push_back(0.0);
        upper_bounds_.push_back(1.0);
        for (int i = 1; i < dim_; ++i) {
            lower_bounds_.push_back(-1.0);
            upper_bounds_.push_back(1.0);
        }
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const {
        const double x0 = pop.gene(idx, 0);
        double sum1 = 0.0, sum2 = 0.0;
        int count1 = 0, count2 = 0;
        for (int j = 2; j <= dim_; ++j) {
            const double yj = pop.gene(idx, j - 1) - std::sin(6.0 * M_PI * x0 + j * M_PI / dim_);
            const double yj2 = yj * yj;
            if (j % 2 == 0) { sum2 += yj2; ++count2; }
            else             { sum1 += yj2; ++count1; }
        }
        pop.objective(idx, 0) = x0 + 2.0 * sum1 / count1;
        pop.objective(idx, 1) = 1.0 - std::sqrt(x0) + 2.0 * sum2 / count2;
    }

    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// UF2 — 2 objectives
// ============================================================

struct UF2 {
    int dim_ = 30;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    UF2() : UF2(30) {}
    explicit UF2(int dim) : dim_(dim) {
        lower_bounds_.reserve(dim_);
        upper_bounds_.reserve(dim_);
        lower_bounds_.push_back(0.0);
        upper_bounds_.push_back(1.0);
        for (int i = 1; i < dim_; ++i) {
            lower_bounds_.push_back(-1.0);
            upper_bounds_.push_back(1.0);
        }
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const {
        const double x0 = pop.gene(idx, 0);
        double sum1 = 0.0, sum2 = 0.0;
        int count1 = 0, count2 = 0;
        for (int j = 2; j <= dim_; ++j) {
            const double term = 0.3 * x0 * x0 * std::cos(24.0 * M_PI * x0 + 4.0 * j * M_PI / dim_)
                                + 0.6 * x0;
            const double theta = 6.0 * M_PI * x0 + j * M_PI / dim_;
            if (j % 2 == 0) {
                const double yj = pop.gene(idx, j - 1) - term * std::sin(theta);
                sum2 += yj * yj; ++count2;
            } else {
                const double yj = pop.gene(idx, j - 1) - term * std::cos(theta);
                sum1 += yj * yj; ++count1;
            }
        }
        pop.objective(idx, 0) = x0 + 2.0 * sum1 / count1;
        pop.objective(idx, 1) = 1.0 - std::sqrt(x0) + 2.0 * sum2 / count2;
    }

    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// UF3 — 2 objectives, all variables in [0,1]
// ============================================================

struct UF3 {
    int dim_ = 30;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    UF3() : UF3(30) {}
    explicit UF3(int dim) : dim_(dim) {
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const {
        const double x0 = pop.gene(idx, 0);
        double sum1 = 0.0, sum2 = 0.0;
        double prod1 = 1.0, prod2 = 1.0;
        int count1 = 0, count2 = 0;
        for (int j = 2; j <= dim_; ++j) {
            const double yj = pop.gene(idx, j - 1)
                              - std::pow(x0, 0.5 * (1.0 + 3.0 * (j - 2.0) / (dim_ - 2.0)));
            const double pj = std::cos(20.0 * yj * M_PI / std::sqrt(static_cast<double>(j)));
            if (j % 2 == 0) { sum2 += yj * yj; prod2 *= pj; ++count2; }
            else             { sum1 += yj * yj; prod1 *= pj; ++count1; }
        }
        pop.objective(idx, 0) = x0 + 2.0 * (4.0 * sum1 - 2.0 * prod1 + 2.0) / count1;
        pop.objective(idx, 1) = 1.0 - std::sqrt(x0) + 2.0 * (4.0 * sum2 - 2.0 * prod2 + 2.0) / count2;
    }

    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// UF4 — 2 objectives, x[i>0] in [-2,2]
// ============================================================

struct UF4 {
    int dim_ = 30;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    UF4() : UF4(30) {}
    explicit UF4(int dim) : dim_(dim) {
        lower_bounds_.reserve(dim_);
        upper_bounds_.reserve(dim_);
        lower_bounds_.push_back(0.0);
        upper_bounds_.push_back(1.0);
        for (int i = 1; i < dim_; ++i) {
            lower_bounds_.push_back(-2.0);
            upper_bounds_.push_back(2.0);
        }
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const {
        const double x0 = pop.gene(idx, 0);
        double sum1 = 0.0, sum2 = 0.0;
        int count1 = 0, count2 = 0;
        for (int j = 2; j <= dim_; ++j) {
            const double yj = pop.gene(idx, j - 1) - std::sin(6.0 * M_PI * x0 + j * M_PI / dim_);
            const double hj = std::abs(yj) / (1.0 + std::exp(2.0 * std::abs(yj)));
            if (j % 2 == 0) { sum2 += hj; ++count2; }
            else             { sum1 += hj; ++count1; }
        }
        pop.objective(idx, 0) = x0 + 2.0 * sum1 / count1;
        pop.objective(idx, 1) = 1.0 - x0 * x0 + 2.0 * sum2 / count2;
    }

    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// UF5 — 2 objectives; N=10, epsilon=0.1
// ============================================================

struct UF5 {
    int    dim_     = 30;
    int    n_       = 10;
    double epsilon_ = 0.1;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    UF5() : UF5(30, 10, 0.1) {}
    UF5(int dim, int n, double epsilon) : dim_(dim), n_(n), epsilon_(epsilon) {
        lower_bounds_.reserve(dim_);
        upper_bounds_.reserve(dim_);
        lower_bounds_.push_back(0.0);
        upper_bounds_.push_back(1.0);
        for (int i = 1; i < dim_; ++i) {
            lower_bounds_.push_back(-1.0);
            upper_bounds_.push_back(1.0);
        }
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const {
        const double x0 = pop.gene(idx, 0);
        double sum1 = 0.0, sum2 = 0.0;
        int count1 = 0, count2 = 0;
        for (int j = 2; j <= dim_; ++j) {
            const double yj = pop.gene(idx, j - 1) - std::sin(6.0 * M_PI * x0 + j * M_PI / dim_);
            const double hj = 2.0 * yj * yj - std::cos(4.0 * M_PI * yj) + 1.0;
            if (j % 2 == 0) { sum2 += hj; ++count2; }
            else             { sum1 += hj; ++count1; }
        }
        const double hterm = (0.5 / n_ + epsilon_) * std::abs(std::sin(2.0 * n_ * M_PI * x0));
        pop.objective(idx, 0) = x0 + hterm + 2.0 * sum1 / count1;
        pop.objective(idx, 1) = 1.0 - x0 + hterm + 2.0 * sum2 / count2;
    }

    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// UF6 — 2 objectives; N=2, epsilon=0.1
// ============================================================

struct UF6 {
    int    dim_     = 30;
    int    n_       = 2;
    double epsilon_ = 0.1;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    UF6() : UF6(30, 2, 0.1) {}
    UF6(int dim, int n, double epsilon) : dim_(dim), n_(n), epsilon_(epsilon) {
        lower_bounds_.reserve(dim_);
        upper_bounds_.reserve(dim_);
        lower_bounds_.push_back(0.0);
        upper_bounds_.push_back(1.0);
        for (int i = 1; i < dim_; ++i) {
            lower_bounds_.push_back(-1.0);
            upper_bounds_.push_back(1.0);
        }
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const {
        const double x0 = pop.gene(idx, 0);
        double sum1 = 0.0, sum2 = 0.0;
        double prod1 = 1.0, prod2 = 1.0;
        int count1 = 0, count2 = 0;
        for (int j = 2; j <= dim_; ++j) {
            const double yj = pop.gene(idx, j - 1) - std::sin(6.0 * M_PI * x0 + j * M_PI / dim_);
            const double pj = std::cos(20.0 * yj * M_PI / std::sqrt(static_cast<double>(j)));
            if (j % 2 == 0) { sum2 += yj * yj; prod2 *= pj; ++count2; }
            else             { sum1 += yj * yj; prod1 *= pj; ++count1; }
        }
        double hj = 2.0 * (0.5 / n_ + epsilon_) * std::sin(2.0 * n_ * M_PI * x0);
        if (hj < 0.0) hj = 0.0;
        pop.objective(idx, 0) = x0 + hj + 2.0 * (4.0 * sum1 - 2.0 * prod1 + 2.0) / count1;
        pop.objective(idx, 1) = 1.0 - x0 + hj + 2.0 * (4.0 * sum2 - 2.0 * prod2 + 2.0) / count2;
    }

    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// UF7 — 2 objectives; PS function: y0 = x0^0.2
// ============================================================

struct UF7 {
    int dim_ = 30;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    UF7() : UF7(30) {}
    explicit UF7(int dim) : dim_(dim) {
        lower_bounds_.reserve(dim_);
        upper_bounds_.reserve(dim_);
        lower_bounds_.push_back(0.0);
        upper_bounds_.push_back(1.0);
        for (int i = 1; i < dim_; ++i) {
            lower_bounds_.push_back(-1.0);
            upper_bounds_.push_back(1.0);
        }
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const {
        const double x0 = pop.gene(idx, 0);
        double sum1 = 0.0, sum2 = 0.0;
        int count1 = 0, count2 = 0;
        for (int j = 2; j <= dim_; ++j) {
            const double yj = pop.gene(idx, j - 1) - std::sin(6.0 * M_PI * x0 + j * M_PI / dim_);
            if (j % 2 == 0) { sum2 += yj * yj; ++count2; }
            else             { sum1 += yj * yj; ++count1; }
        }
        const double y0 = std::pow(x0, 0.2);
        pop.objective(idx, 0) = y0 + 2.0 * sum1 / count1;
        pop.objective(idx, 1) = 1.0 - y0 + 2.0 * sum2 / count2;
    }

    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// UF8 — 3 objectives; x[0],x[1] in [0,1], x[i>1] in [-2,2]
// ============================================================

struct UF8 {
    int dim_ = 30;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    UF8() : UF8(30) {}
    explicit UF8(int dim) : dim_(dim) {
        lower_bounds_.reserve(dim_);
        upper_bounds_.reserve(dim_);
        lower_bounds_.push_back(0.0); upper_bounds_.push_back(1.0);
        lower_bounds_.push_back(0.0); upper_bounds_.push_back(1.0);
        for (int i = 2; i < dim_; ++i) {
            lower_bounds_.push_back(-2.0);
            upper_bounds_.push_back(2.0);
        }
    }

    static constexpr int num_objectives() { return 3; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const {
        const double x0 = pop.gene(idx, 0);
        const double x1 = pop.gene(idx, 1);
        double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
        int count1 = 0, count2 = 0, count3 = 0;
        for (int j = 3; j <= dim_; ++j) {
            const double yj = pop.gene(idx, j - 1)
                              - 2.0 * x1 * std::sin(2.0 * M_PI * x0 + j * M_PI / dim_);
            const double yj2 = yj * yj;
            if      (j % 3 == 1) { sum1 += yj2; ++count1; }
            else if (j % 3 == 2) { sum2 += yj2; ++count2; }
            else                 { sum3 += yj2; ++count3; }
        }
        pop.objective(idx, 0) = std::cos(0.5 * M_PI * x0) * std::cos(0.5 * M_PI * x1) + 2.0 * sum1 / count1;
        pop.objective(idx, 1) = std::cos(0.5 * M_PI * x0) * std::sin(0.5 * M_PI * x1) + 2.0 * sum2 / count2;
        pop.objective(idx, 2) = std::sin(0.5 * M_PI * x0)                              + 2.0 * sum3 / count3;
    }

    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// UF9 — 3 objectives; epsilon=0.1
// ============================================================

struct UF9 {
    int    dim_     = 30;
    double epsilon_ = 0.1;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    UF9() : UF9(30, 0.1) {}
    UF9(int dim, double epsilon) : dim_(dim), epsilon_(epsilon) {
        lower_bounds_.reserve(dim_);
        upper_bounds_.reserve(dim_);
        lower_bounds_.push_back(0.0); upper_bounds_.push_back(1.0);
        lower_bounds_.push_back(0.0); upper_bounds_.push_back(1.0);
        for (int i = 2; i < dim_; ++i) {
            lower_bounds_.push_back(-2.0);
            upper_bounds_.push_back(2.0);
        }
    }

    static constexpr int num_objectives() { return 3; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const {
        const double x0 = pop.gene(idx, 0);
        const double x1 = pop.gene(idx, 1);
        double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
        int count1 = 0, count2 = 0, count3 = 0;
        for (int j = 3; j <= dim_; ++j) {
            const double yj = pop.gene(idx, j - 1)
                              - 2.0 * x1 * std::sin(2.0 * M_PI * x0 + j * M_PI / dim_);
            const double yj2 = yj * yj;
            if      (j % 3 == 1) { sum1 += yj2; ++count1; }
            else if (j % 3 == 2) { sum2 += yj2; ++count2; }
            else                 { sum3 += yj2; ++count3; }
        }
        double t = (1.0 + epsilon_) * (1.0 - 4.0 * (2.0 * x0 - 1.0) * (2.0 * x0 - 1.0));
        if (t < 0.0) t = 0.0;
        pop.objective(idx, 0) = 0.5 * (t + 2.0 * x0) * x1          + 2.0 * sum1 / count1;
        pop.objective(idx, 1) = 0.5 * (t - 2.0 * x0 + 2.0) * x1    + 2.0 * sum2 / count2;
        pop.objective(idx, 2) = 1.0 - x1                             + 2.0 * sum3 / count3;
    }

    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

// ============================================================
// UF10 — 3 objectives; Rastrigin-type hj
// ============================================================

struct UF10 {
    int dim_ = 30;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    UF10() : UF10(30) {}
    explicit UF10(int dim) : dim_(dim) {
        lower_bounds_.reserve(dim_);
        upper_bounds_.reserve(dim_);
        lower_bounds_.push_back(0.0); upper_bounds_.push_back(1.0);
        lower_bounds_.push_back(0.0); upper_bounds_.push_back(1.0);
        for (int i = 2; i < dim_; ++i) {
            lower_bounds_.push_back(-2.0);
            upper_bounds_.push_back(2.0);
        }
    }

    static constexpr int num_objectives() { return 3; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const {
        const double x0 = pop.gene(idx, 0);
        const double x1 = pop.gene(idx, 1);
        double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
        int count1 = 0, count2 = 0, count3 = 0;
        for (int j = 3; j <= dim_; ++j) {
            const double yj = pop.gene(idx, j - 1)
                              - 2.0 * x1 * std::sin(2.0 * M_PI * x0 + j * M_PI / dim_);
            const double hj = 4.0 * yj * yj - std::cos(8.0 * M_PI * yj) + 1.0;
            if      (j % 3 == 1) { sum1 += hj; ++count1; }
            else if (j % 3 == 2) { sum2 += hj; ++count2; }
            else                 { sum3 += hj; ++count3; }
        }
        pop.objective(idx, 0) = std::cos(0.5 * M_PI * x0) * std::cos(0.5 * M_PI * x1) + 2.0 * sum1 / count1;
        pop.objective(idx, 1) = std::cos(0.5 * M_PI * x0) * std::sin(0.5 * M_PI * x1) + 2.0 * sum2 / count2;
        pop.objective(idx, 2) = std::sin(0.5 * M_PI * x0)                              + 2.0 * sum3 / count3;
    }

    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

} // namespace ea
