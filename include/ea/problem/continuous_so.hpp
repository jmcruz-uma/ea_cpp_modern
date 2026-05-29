#pragma once
/// @file continuous_so.hpp
/// @brief Single-objective continuous benchmark functions.
///
/// Matches the functions used in:
///   Cotta & Martínez-Cruz, "Evaluating the Impact of Hysteretic Phenomena and
///   Implementation Choices on Energy Consumption in Evolutionary Algorithms" (2025)
///
/// All functions are MINIMISATION. Domain is [-range, range]^d unless noted.
/// Implementations match Bio4Res/ea (ea_cpp_original) exactly, using SoA access.

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>

#include <ea/core/population.hpp>

namespace ea {

// ------------------------------------------------------------------
// Sphere  f(x) = Σ x_i²
// ------------------------------------------------------------------
struct Sphere {
    int dim;
    double lb, ub;

    Sphere(int d, double range) : dim(d), lb(-range), ub(range) {}

    int  dimension()    const { return dim; }
    int  num_objectives() const { return 1; }
    double lower_bound() const { return lb; }
    double upper_bound() const { return ub; }

    void evaluate(Population<>& pop, int idx) const {
        double s = 0.0;
        for (int j = 0; j < dim; ++j) {
            double v = pop.gene(idx, j);
            s += v * v;
        }
        pop.objective(idx, 0) = s;
        pop.set_evaluated(idx, true);
    }
};

// ------------------------------------------------------------------
// Rastrigin  f(x) = 10d + Σ (x_i² - 10 cos(2π x_i))
// ------------------------------------------------------------------
struct Rastrigin {
    int dim;
    double lb, ub;
    static constexpr double A = 10.0;

    Rastrigin(int d, double range) : dim(d), lb(-range), ub(range) {}

    int  dimension()    const { return dim; }
    int  num_objectives() const { return 1; }
    double lower_bound() const { return lb; }
    double upper_bound() const { return ub; }

    void evaluate(Population<>& pop, int idx) const {
        double s = A * dim;
        for (int j = 0; j < dim; ++j) {
            double v = pop.gene(idx, j);
            s += v * v - A * std::cos(2.0 * std::numbers::pi * v);
        }
        pop.objective(idx, 0) = s;
        pop.set_evaluated(idx, true);
    }
};

// ------------------------------------------------------------------
// Ackley  f(x) = -a·exp(-b·sqrt(Σx²/d)) - exp(Σcos(c·x_i)/d) + a + e
// ------------------------------------------------------------------
struct Ackley {
    int dim;
    double lb, ub;
    static constexpr double A = 20.0;
    static constexpr double B = 0.2;
    static constexpr double C = 2.0 * std::numbers::pi;

    Ackley(int d, double range) : dim(d), lb(-range), ub(range) {}

    int  dimension()    const { return dim; }
    int  num_objectives() const { return 1; }
    double lower_bound() const { return lb; }
    double upper_bound() const { return ub; }

    void evaluate(Population<>& pop, int idx) const {
        double sum1 = 0.0, sum2 = 0.0;
        for (int j = 0; j < dim; ++j) {
            double v = pop.gene(idx, j);
            sum1 += v * v;
            sum2 += std::cos(C * v);
        }
        pop.objective(idx, 0) =
            -A * std::exp(-B * std::sqrt(sum1 / dim))
            - std::exp(sum2 / dim)
            + A + std::numbers::e;
        pop.set_evaluated(idx, true);
    }
};

// ------------------------------------------------------------------
// Rosenbrock  f(x) = Σ [100(x_{i+1}-x_i²)² + (1-x_i)²]
// ------------------------------------------------------------------
struct Rosenbrock {
    int dim;
    double lb, ub;

    Rosenbrock(int d, double range) : dim(d), lb(-range), ub(range) {}

    int  dimension()    const { return dim; }
    int  num_objectives() const { return 1; }
    double lower_bound() const { return lb; }
    double upper_bound() const { return ub; }

    void evaluate(Population<>& pop, int idx) const {
        double s = 0.0;
        for (int j = 0; j < dim - 1; ++j) {
            double x  = pop.gene(idx, j);
            double xn = pop.gene(idx, j + 1);
            double t  = xn - x * x;
            double u  = 1.0 - x;
            s += 100.0 * t * t + u * u;
        }
        pop.objective(idx, 0) = s;
        pop.set_evaluated(idx, true);
    }
};

// ------------------------------------------------------------------
// Factory: build problem by name (matches ea_cpp_original problem strings)
// ------------------------------------------------------------------
struct SOProblemConfig {
    int    dim;
    double range;
};

} // namespace ea
