#pragma once
/// @file population.hpp
/// @brief Structure of Arrays (SoA) population storage — cache-friendly, SIMD-ready.

#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <span>
#include <ea/core/encoding.hpp>

namespace ea {

/// Bit flags for individual state.
enum class IndFlags : uint8_t {
    None      = 0,
    Evaluated = 1 << 0,
    Feasible  = 1 << 1,
    Dominated = 1 << 2,
};

constexpr IndFlags operator|(IndFlags a, IndFlags b) {
    return static_cast<IndFlags>(static_cast<uint8_t>(a) | static_cast<uint8_t>(b));
}
constexpr IndFlags operator&(IndFlags a, IndFlags b) {
    return static_cast<IndFlags>(static_cast<uint8_t>(a) & static_cast<uint8_t>(b));
}
constexpr bool has_flag(IndFlags f, IndFlags flag) {
    return (f & flag) == flag;
}

/// Structure of Arrays population — all genes, objectives, and constraints
/// stored contiguously for maximum cache locality and SIMD vectorization.
struct Population {
    // --- SoA data ---
    std::vector<double> genes;          ///< [pop_size * dim] contiguos
    std::vector<double> objectives;    ///< [pop_size * n_obj] contiguos
    std::vector<double> constraints;   ///< [pop_size * n_const] contiguos
    std::vector<IndFlags> flags;       ///< [pop_size] state flags
    std::vector<double> lower_bounds;  ///< [dim] bounds
    std::vector<double> upper_bounds;  ///< [dim] bounds

    // --- Dimensions ---
    int pop_size = 0;   ///< Number of individuals
    int dim = 0;        ///< Number of decision variables
    int n_obj = 0;      ///< Number of objectives
    int n_const = 0;    ///< Number of constraints
    Encoding encoding = Encoding::Real;

    // --- Constructors ---
    Population() = default;

    Population(int pop_size, int dim, int n_obj,
               Encoding enc = Encoding::Real,
               int n_const = 0,
               double lb = 0.0, double ub = 1.0)
        : genes(pop_size * dim, 0.0)
        , objectives(pop_size * n_obj, 0.0)
        , constraints(pop_size * n_const, 0.0)
        , flags(pop_size, IndFlags::None)
        , lower_bounds(dim, lb)
        , upper_bounds(dim, ub)
        , pop_size(pop_size)
        , dim(dim)
        , n_obj(n_obj)
        , n_const(n_const)
        , encoding(enc)
    {}

    // --- Gene access ---
    double gene(int i, int j) const { return genes[i * dim + j]; }
    double& gene(int i, int j) { return genes[i * dim + j]; }

    /// Raw pointer to individual's genes — SIMD-friendly
    const double* genes_ptr(int i) const { return genes.data() + i * dim; }
    double* genes_ptr(int i) { return genes.data() + i * dim; }

    /// Mutable span over individual's genes
    std::span<double> genes_span(int i) { return {genes.data() + i * dim, static_cast<size_t>(dim)}; }
    std::span<const double> genes_span(int i) const { return {genes.data() + i * dim, static_cast<size_t>(dim)}; }

    // --- Objective access ---
    double objective(int i, int o) const { return objectives[i * n_obj + o]; }
    double& objective(int i, int o) { return objectives[i * n_obj + o]; }

    const double* objectives_ptr(int i) const { return objectives.data() + i * n_obj; }
    double* objectives_ptr(int i) { return objectives.data() + i * n_obj; }

    std::span<double> objectives_span(int i) { return {objectives.data() + i * n_obj, static_cast<size_t>(n_obj)}; }
    std::span<const double> objectives_span(int i) const { return {objectives.data() + i * n_obj, static_cast<size_t>(n_obj)}; }

    // --- Constraint access ---
    double constraint(int i, int c) const { return constraints[i * n_const + c]; }
    double& constraint(int i, int c) { return constraints[i * n_const + c]; }

    // --- Flags ---
    bool evaluated(int i) const { return has_flag(flags[i], IndFlags::Evaluated); }
    void set_evaluated(int i, bool v) {
        if (v) flags[i] = flags[i] | IndFlags::Evaluated;
        else   flags[i] = static_cast<IndFlags>(static_cast<uint8_t>(flags[i]) & ~static_cast<uint8_t>(IndFlags::Evaluated));
    }

    bool feasible(int i) const { return has_flag(flags[i], IndFlags::Feasible); }
    void set_feasible(int i, bool v) {
        if (v) flags[i] = flags[i] | IndFlags::Feasible;
        else   flags[i] = static_cast<IndFlags>(static_cast<uint8_t>(flags[i]) & ~static_cast<uint8_t>(IndFlags::Feasible));
    }

    // --- Utility ---
    /// Resize population (reallocate all arrays)
    void resize(int new_pop_size) {
        pop_size = new_pop_size;
        genes.resize(pop_size * dim, 0.0);
        objectives.resize(pop_size * n_obj, 0.0);
        constraints.resize(pop_size * n_const, 0.0);
        flags.resize(pop_size, IndFlags::None);
    }

    /// Copy individual src to dst
    void copy_individual(int src, int dst) {
        std::copy_n(genes_ptr(src), dim, genes_ptr(dst));
        std::copy_n(objectives_ptr(src), n_obj, objectives_ptr(dst));
        if (n_const > 0) {
            std::copy_n(constraints.data() + src * n_const, n_const,
                        constraints.data() + dst * n_const);
        }
        flags[dst] = flags[src];
    }

    /// Swap individuals i and j
    void swap_individuals(int i, int j) {
        if (i == j) return;
        // Swap genes
        std::swap_ranges(genes_ptr(i), genes_ptr(i) + dim, genes_ptr(j));
        // Swap objectives
        std::swap_ranges(objectives_ptr(i), objectives_ptr(i) + n_obj, objectives_ptr(j));
        // Swap constraints
        if (n_const > 0) {
            std::swap_ranges(constraints.data() + i * n_const,
                            constraints.data() + i * n_const + n_const,
                            constraints.data() + j * n_const);
        }
        std::swap(flags[i], flags[j]);
    }

    /// Invalidate all evaluations (mark as not evaluated)
    void invalidate_all() {
        std::fill(flags.begin(), flags.end(), IndFlags::None);
    }

    /// Total memory used in bytes
    [[nodiscard]] size_t memory_bytes() const {
        return genes.capacity() * sizeof(double)
             + objectives.capacity() * sizeof(double)
             + constraints.capacity() * sizeof(double)
             + flags.capacity() * sizeof(IndFlags)
             + lower_bounds.capacity() * sizeof(double)
             + upper_bounds.capacity() * sizeof(double);
    }
};

} // namespace ea