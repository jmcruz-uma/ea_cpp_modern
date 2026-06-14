#pragma once
/// @file lz09.hpp
/// @brief LZ09 benchmark problem suite for multi-objective optimization.
///
/// Reference: H. Li and Q. Zhang, "Multiobjective optimization problem with
/// complicated Pareto sets, MOEA/D and NSGA-II", IEEE Trans. Evol. Comput.,
/// 12(2):284-302, April 2009.
///
/// jMetal reference: jmetal-problem/.../problem/multiobjective/lz09/
///
/// LZ09F1-5, F7-F9: 2 objectives. LZ09F6: 3 objectives.
/// Bounds: all variables in [0,1].
///
/// Audited against jMetal 7.4 — all formulas and configurations verified correct.

#include <cmath>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <span>
#include <vector>

namespace ea {

// ============================================================
// LZ09Base — shared alpha/beta/psfunc helpers
// ============================================================

struct LZ09Base {
    int nvar  = 30;
    int nobj  = 2;
    int ptype = 21; ///< PF shape type
    int dtype = 1;  ///< Distance type
    int ltype = 21; ///< PS shape type

    LZ09Base() = default;
    LZ09Base(int nvar_, int nobj_, int ptype_, int dtype_, int ltype_)
        : nvar(nvar_), nobj(nobj_), ptype(ptype_), dtype(dtype_), ltype(ltype_) {}

    void alphaFunction(double alpha[], const std::vector<double>& x, int dim) const {
        if (dim == 2) {
            if (ptype == 21) {
                alpha[0] = x[0];
                alpha[1] = 1.0 - std::sqrt(x[0]);
            } else if (ptype == 22) {
                alpha[0] = x[0];
                alpha[1] = 1.0 - x[0] * x[0];
            } else if (ptype == 23) {
                alpha[0] = x[0];
                alpha[1] = 1.0 - std::sqrt(alpha[0]) - alpha[0] * std::sin(10.0 * alpha[0] * alpha[0] * M_PI);
            } else if (ptype == 24) {
                alpha[0] = x[0];
                alpha[1] = 1.0 - x[0] - 0.05 * std::sin(4.0 * M_PI * x[0]);
            }
        } else {
            if (ptype == 31) {
                alpha[0] = std::cos(x[0] * M_PI / 2.0) * std::cos(x[1] * M_PI / 2.0);
                alpha[1] = std::cos(x[0] * M_PI / 2.0) * std::sin(x[1] * M_PI / 2.0);
                alpha[2] = std::sin(x[0] * M_PI / 2.0);
            } else if (ptype == 32) {
                alpha[0] = 1.0 - std::cos(x[0] * M_PI / 2.0) * std::cos(x[1] * M_PI / 2.0);
                alpha[1] = 1.0 - std::cos(x[0] * M_PI / 2.0) * std::sin(x[1] * M_PI / 2.0);
                alpha[2] = 1.0 - std::sin(x[0] * M_PI / 2.0);
            } else if (ptype == 33) {
                alpha[0] = x[0];
                alpha[1] = x[1];
                alpha[2] = 3.0 - (std::sin(3.0 * M_PI * x[0]) + std::sin(3.0 * M_PI * x[1]))
                           - 2.0 * (x[0] + x[1]);
            } else if (ptype == 34) {
                alpha[0] = x[0] * x[1];
                alpha[1] = x[0] * (1.0 - x[1]);
                alpha[2] = 1.0 - x[0];
            }
        }
    }

    double betaFunction(const std::vector<double>& x) const {
        const int dim = static_cast<int>(x.size());
        if (dim == 0) return 0.0;
        if (dtype == 1) {
            double s = 0.0;
            for (int i = 0; i < dim; ++i) s += x[i] * x[i];
            return 2.0 * s / dim;
        }
        if (dtype == 2) {
            double s = 0.0;
            for (int i = 0; i < dim; ++i) s += std::sqrt(static_cast<double>(i + 1)) * x[i] * x[i];
            return 2.0 * s / dim;
        }
        if (dtype == 3) {
            double s = 0.0;
            for (int i = 0; i < dim; ++i) {
                const double xx = 2.0 * x[i];
                s += xx * xx - std::cos(4.0 * M_PI * xx) + 1.0;
            }
            return 2.0 * s / dim;
        }
        if (dtype == 4) {
            double s = 0.0, p = 1.0;
            for (int i = 0; i < dim; ++i) {
                const double xx = 2.0 * x[i];
                s += xx * xx;
                p *= std::cos(10.0 * M_PI * xx / std::sqrt(static_cast<double>(i + 1)));
            }
            return 2.0 * (s - 2.0 * p + 2.0) / dim;
        }
        return 0.0;
    }

    double psfunc2(double x, double t1, int dim, int css) const {
        ++dim; // matches Java original
        const double xy = 2.0 * (x - 0.5);
        if (ltype == 21) {
            return xy - std::pow(t1, 0.5 * (nvar + 3.0 * dim - 8.0) / (nvar - 2.0));
        }
        if (ltype == 22) {
            return xy - std::sin(6.0 * M_PI * t1 + dim * M_PI / nvar);
        }
        if (ltype == 23) {
            const double theta = 6.0 * M_PI * t1 + dim * M_PI / nvar;
            const double ra = 0.8 * t1;
            return (css == 1) ? xy - ra * std::cos(theta) : xy - ra * std::sin(theta);
        }
        if (ltype == 24) {
            const double theta = 6.0 * M_PI * t1 + dim * M_PI / nvar;
            const double ra = 0.8 * t1;
            return (css == 1) ? xy - ra * std::cos(theta / 3.0) : xy - ra * std::sin(theta);
        }
        if (ltype == 25) {
            const double rho   = 0.8;
            const double phi   = M_PI * t1;
            const double theta = 6.0 * M_PI * t1 + dim * M_PI / nvar;
            if (css == 1) return xy - rho * std::sin(phi) * std::sin(theta);
            if (css == 2) return xy - rho * std::sin(phi) * std::cos(theta);
            return xy - rho * std::cos(phi);
        }
        if (ltype == 26) {
            const double theta = 6.0 * M_PI * t1 + dim * M_PI / nvar;
            const double ra = 0.3 * t1 * (t1 * std::cos(4.0 * theta) + 2.0);
            return (css == 1) ? xy - ra * std::cos(theta) : xy - ra * std::sin(theta);
        }
        return 0.0;
    }

    double psfunc3(double x, double t1, double t2, int dim) const {
        ++dim; // matches Java original
        const double xy = 4.0 * (x - 0.5);
        if (ltype == 31) {
            const double rate = static_cast<double>(dim) / nvar;
            return xy - 4.0 * (t1 * t1 * rate + t2 * (1.0 - rate)) + 2.0;
        }
        if (ltype == 32) {
            return xy - 2.0 * t2 * std::sin(2.0 * M_PI * t1 + dim * M_PI / nvar);
        }
        return 0.0;
    }

    void objective(const std::vector<double>& xVar, std::vector<double>& yObj) const {
        if (nobj == 2) {
            std::vector<double> aa, bb;
            if (ltype == 21 || ltype == 22 || ltype == 23 || ltype == 24 || ltype == 26) {
                for (int n = 1; n < nvar; ++n) {
                    if (n % 2 == 0) aa.push_back(psfunc2(xVar[n], xVar[0], n, 1));
                    else            bb.push_back(psfunc2(xVar[n], xVar[0], n, 2));
                }
            } else if (ltype == 25) {
                for (int n = 1; n < nvar; ++n) {
                    if      (n % 3 == 0) aa.push_back(psfunc2(xVar[n], xVar[0], n, 1));
                    else if (n % 3 == 1) bb.push_back(psfunc2(xVar[n], xVar[0], n, 2));
                    else {
                        const double c = psfunc2(xVar[n], xVar[0], n, 3);
                        if (n % 2 == 0) aa.push_back(c);
                        else            bb.push_back(c);
                    }
                }
            }
            const double g = betaFunction(aa);
            const double h = betaFunction(bb);
            double alpha[2] = {};
            alphaFunction(alpha, xVar, 2);
            yObj[0] = alpha[0] + h;
            yObj[1] = alpha[1] + g;
        }

        if (nobj == 3) {
            std::vector<double> aa, bb, cc;
            for (int n = 2; n < nvar; ++n) {
                const double a = psfunc3(xVar[n], xVar[0], xVar[1], n);
                if      (n % 3 == 0) aa.push_back(a);
                else if (n % 3 == 1) bb.push_back(a);
                else                 cc.push_back(a);
            }
            const double g = betaFunction(aa);
            const double h = betaFunction(bb);
            const double e = betaFunction(cc);
            double alpha[3] = {};
            alphaFunction(alpha, xVar, 3);
            yObj[0] = alpha[0] + h;
            yObj[1] = alpha[1] + g;
            yObj[2] = alpha[2] + e;
        }
    }
};

// ============================================================
// Concrete LZ09 problems
// ============================================================

// Helper: evaluate via LZ09Base, read nvar variables, write nobj objectives.
// Used identically by every concrete struct below.
#define LZ09_EVALUATE_BODY(NOBJ) \
    std::vector<double> x(dim_); \
    std::vector<double> y(NOBJ, 0.0); \
    for (int i = 0; i < dim_; ++i) x[i] = pop.gene(idx, i); \
    base_.objective(x, y); \
    for (int o = 0; o < NOBJ; ++o) pop.objective(idx, o) = y[o];

struct LZ09F1 {
    int dim_ = 10;
    LZ09Base base_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    LZ09F1() : LZ09F1(10) {}
    explicit LZ09F1(int dim) : dim_(dim), base_(dim, 2, 21, 1, 21) {
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const { LZ09_EVALUATE_BODY(2) }
    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

struct LZ09F2 {
    int dim_ = 30;
    LZ09Base base_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    LZ09F2() : LZ09F2(30) {}
    explicit LZ09F2(int dim) : dim_(dim), base_(dim, 2, 21, 1, 22) {
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const { LZ09_EVALUATE_BODY(2) }
    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

struct LZ09F3 {
    int dim_ = 30;
    LZ09Base base_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    LZ09F3() : LZ09F3(30) {}
    explicit LZ09F3(int dim) : dim_(dim), base_(dim, 2, 21, 1, 23) {
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const { LZ09_EVALUATE_BODY(2) }
    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

struct LZ09F4 {
    int dim_ = 30;
    LZ09Base base_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    LZ09F4() : LZ09F4(30) {}
    explicit LZ09F4(int dim) : dim_(dim), base_(dim, 2, 21, 1, 24) {
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const { LZ09_EVALUATE_BODY(2) }
    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

struct LZ09F5 {
    int dim_ = 30;
    LZ09Base base_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    LZ09F5() : LZ09F5(30) {}
    explicit LZ09F5(int dim) : dim_(dim), base_(dim, 2, 21, 1, 26) {
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const { LZ09_EVALUATE_BODY(2) }
    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

struct LZ09F6 {
    int dim_ = 10;
    LZ09Base base_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    LZ09F6() : LZ09F6(10) {}
    explicit LZ09F6(int dim) : dim_(dim), base_(dim, 3, 31, 1, 32) {
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 3; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const { LZ09_EVALUATE_BODY(3) }
    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

struct LZ09F7 {
    int dim_ = 10;
    LZ09Base base_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    LZ09F7() : LZ09F7(10) {}
    explicit LZ09F7(int dim) : dim_(dim), base_(dim, 2, 21, 3, 21) {
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const { LZ09_EVALUATE_BODY(2) }
    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

struct LZ09F8 {
    int dim_ = 10;
    LZ09Base base_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    LZ09F8() : LZ09F8(10) {}
    explicit LZ09F8(int dim) : dim_(dim), base_(dim, 2, 21, 4, 21) {
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const { LZ09_EVALUATE_BODY(2) }
    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

struct LZ09F9 {
    int dim_ = 30;
    LZ09Base base_;
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    LZ09F9() : LZ09F9(30) {}
    explicit LZ09F9(int dim) : dim_(dim), base_(dim, 2, 22, 1, 22) {
        lower_bounds_.assign(dim_, 0.0);
        upper_bounds_.assign(dim_, 1.0);
    }

    static constexpr int num_objectives() { return 2; }
    int dimension() const { return dim_; }
    std::span<const double> lower_bounds() const { return lower_bounds_; }
    std::span<const double> upper_bounds() const { return upper_bounds_; }

    void evaluate(Population<>& pop, int idx) const { LZ09_EVALUATE_BODY(2) }
    void evaluate_batch(Population<>& pop, int start, int count) const {
        for (int i = start; i < start + count; ++i) evaluate(pop, i);
    }
};

#undef LZ09_EVALUATE_BODY

} // namespace ea
