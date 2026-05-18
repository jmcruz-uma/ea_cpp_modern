#pragma once
/// @file random.hpp
/// @brief Random number generation utilities for EAs.
/// Thread-safe, cache-friendly, with SIMD-friendly batch generation.

#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

namespace ea {

/// Thread-local random number generator with common EA distributions.
class Random {
public:
    /// Get thread-local instance
    static Random& instance() {
        thread_local Random rng;
        return rng;
    }

    /// Set seed for reproducibility
    void seed(uint64_t s) { gen_.seed(s); }

    // --- Uniform distributions ---

    double uniform(double lo = 0.0, double hi = 1.0) {
        return std::uniform_real_distribution<double>(lo, hi)(gen_);
    }

    int uniform_int(int lo, int hi) { return std::uniform_int_distribution<int>(lo, hi)(gen_); }

    /// Fill a range with uniform random values — SIMD-friendly bulk generation
    void uniform_batch(double* out, int n, double lo = 0.0, double hi = 1.0) {
        std::uniform_real_distribution<double> dist(lo, hi);
        for (int i = 0; i < n; ++i) {
            out[i] = dist(gen_);
        }
    }

    /// Fill a range with uniform random values within per-element bounds
    void uniform_batch(double* out, int n, const double* lo, const double* hi) {
        for (int i = 0; i < n; ++i) {
            out[i] = std::uniform_real_distribution<double>(lo[i], hi[i])(gen_);
        }
    }

    // --- Normal distribution ---

    double normal(double mean = 0.0, double stddev = 1.0) {
        return std::normal_distribution<double>(mean, stddev)(gen_);
    }

    // --- Cauchy distribution (for Cauchy mutation) ---

    double cauchy(double location = 0.0, double scale = 1.0) {
        return std::cauchy_distribution<double>(location, scale)(gen_);
    }

    // --- Levy distribution (for Levy flight mutation) ---

    double levy(double scale = 1.0) {
        // Levy distribution: generate using Mantel's method
        double u = uniform(-0.5 * M_PI, 0.5 * M_PI);
        double v = uniform(0.0, 1.0);
        return scale * std::sin(u) / std::pow(std::cos(u), 1.0 / 1.5) *
               std::pow(std::cos((1.0 - 1.5) * u / (1.0 / 1.5)) / v, (1.5 - 1.0) / 1.0);
    }

    // --- Permutation generation ---

    /// Generate a random permutation of [0, n)
    std::vector<int> permutation(int n) {
        std::vector<int> perm(n);
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), gen_);
        return perm;
    }

    /// Shuffle a span in-place
    void shuffle(std::span<int> s) { std::shuffle(s.begin(), s.end(), gen_); }

    // --- Coin flip ---

    bool coin_flip(double p = 0.5) { return uniform() < p; }

    std::mt19937_64& engine() { return gen_; }
    Random() {
        std::random_device rd;
        gen_.seed(rd());
    }

    std::mt19937_64 gen_; // Mersenne Twister 64-bit — fast, high quality
};

} // namespace ea