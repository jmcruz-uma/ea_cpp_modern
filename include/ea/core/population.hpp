#pragma once
/// @file population.hpp
/// @brief Structure of Arrays (SoA) population storage — cache-friendly, SIMD-ready.
///
/// Population<E> is parameterised on Encoding E (default: Encoding::Real).
/// The gene storage type is gene_t<E>:
///   Real → double (8 B), Binary → uint8_t (1 B), Integer/Permutation → int32_t (4 B).
///
/// For Encoding::Real this template instantiates to code identical to a plain
/// std::vector<double> implementation — zero overhead from the abstraction.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <ea/core/encoding.hpp>
#include <span>
#include <vector>

namespace ea {

/// Bit flags for individual state.
enum class IndFlags : uint8_t {
    None = 0,
    Evaluated = 1 << 0,
    Feasible = 1 << 1,
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
///
/// @tparam E  Encoding that determines the gene value type via gene_t<E>.
template <Encoding E = Encoding::Real>
struct Population {
    // --- Gene type ---
    using gene_type = gene_t<E>;
    static constexpr Encoding encoding_value = E;

    // --- SoA data ---
    std::vector<gene_type> genes;       ///< [pop_size × dim] contiguous
    std::vector<double>    objectives;  ///< [pop_size × n_obj] contiguous
    std::vector<double>    constraints; ///< [pop_size × n_const] contiguous
    std::vector<IndFlags>  flags;       ///< [pop_size] state flags
    std::vector<double>    lower_bounds; ///< [dim] per-variable lower bound
    std::vector<double>    upper_bounds; ///< [dim] per-variable upper bound

    // --- Dimensions ---
    int pop_size = 0; ///< Number of individuals
    int dim = 0;      ///< Number of decision variables
    int n_obj = 0;    ///< Number of objectives
    int n_const = 0;  ///< Number of constraints
    Encoding encoding = E; ///< Runtime encoding tag (always equals template param E)

    // --- Constructors ---
    Population() = default;

    Population(int pop_size, int dim, int n_obj, int n_const = 0,
               double lb = 0.0, double ub = 1.0)
        : genes(pop_size * dim, gene_type{})
        , objectives(pop_size * n_obj, 0.0)
        , constraints(pop_size * n_const, 0.0)
        , flags(pop_size, IndFlags::None)
        , lower_bounds(dim, lb)
        , upper_bounds(dim, ub)
        , pop_size(pop_size)
        , dim(dim)
        , n_obj(n_obj)
        , n_const(n_const)
        , encoding(E) {}

    // --- Gene access ---
    gene_type  gene(int i, int j) const { return genes[i * dim + j]; }
    gene_type& gene(int i, int j)       { return genes[i * dim + j]; }

    /// Raw pointer to individual's genes — SIMD-friendly
    const gene_type* genes_ptr(int i) const { return genes.data() + i * dim; }
    gene_type*       genes_ptr(int i)       { return genes.data() + i * dim; }

    /// Span over individual's genes
    std::span<gene_type>       genes_span(int i)       { return {genes.data() + i * dim, static_cast<size_t>(dim)}; }
    std::span<const gene_type> genes_span(int i) const { return {genes.data() + i * dim, static_cast<size_t>(dim)}; }

    // --- Objective access ---
    double  objective(int i, int o) const { return objectives[i * n_obj + o]; }
    double& objective(int i, int o)       { return objectives[i * n_obj + o]; }

    const double* objectives_ptr(int i) const { return objectives.data() + i * n_obj; }
    double*       objectives_ptr(int i)       { return objectives.data() + i * n_obj; }

    std::span<double>       objectives_span(int i)       { return {objectives.data() + i * n_obj, static_cast<size_t>(n_obj)}; }
    std::span<const double> objectives_span(int i) const { return {objectives.data() + i * n_obj, static_cast<size_t>(n_obj)}; }

    // --- Constraint access ---
    double  constraint(int i, int c) const { return constraints[i * n_const + c]; }
    double& constraint(int i, int c)       { return constraints[i * n_const + c]; }

    // --- Flags ---
    bool evaluated(int i) const { return has_flag(flags[i], IndFlags::Evaluated); }
    void set_evaluated(int i, bool v) {
        if (v)
            flags[i] = flags[i] | IndFlags::Evaluated;
        else
            flags[i] = static_cast<IndFlags>(static_cast<uint8_t>(flags[i]) &
                                             ~static_cast<uint8_t>(IndFlags::Evaluated));
    }

    bool feasible(int i) const { return has_flag(flags[i], IndFlags::Feasible); }
    void set_feasible(int i, bool v) {
        if (v)
            flags[i] = flags[i] | IndFlags::Feasible;
        else
            flags[i] = static_cast<IndFlags>(static_cast<uint8_t>(flags[i]) &
                                             ~static_cast<uint8_t>(IndFlags::Feasible));
    }

    // --- Utility ---
    /// Resize population (reallocate all arrays)
    void resize(int new_pop_size) {
        pop_size = new_pop_size;
        genes.resize(pop_size * dim, gene_type{});
        objectives.resize(pop_size * n_obj, 0.0);
        constraints.resize(pop_size * n_const, 0.0);
        flags.resize(pop_size, IndFlags::None);
    }

    /// Copy individual src_idx from another population into this population's dst_idx.
    void copy_from_other(const Population<E>& other, int src_idx, int dst_idx) {
        std::copy_n(other.genes_ptr(src_idx), dim, genes_ptr(dst_idx));
        std::copy_n(other.objectives_ptr(src_idx), n_obj, objectives_ptr(dst_idx));
        if (n_const > 0) {
            std::copy_n(other.constraints.data() + src_idx * n_const, n_const,
                        constraints.data() + dst_idx * n_const);
        }
        flags[dst_idx] = other.flags[src_idx];
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
        if (i == j)
            return;
        std::swap_ranges(genes_ptr(i), genes_ptr(i) + dim, genes_ptr(j));
        std::swap_ranges(objectives_ptr(i), objectives_ptr(i) + n_obj, objectives_ptr(j));
        if (n_const > 0) {
            std::swap_ranges(constraints.data() + i * n_const,
                             constraints.data() + i * n_const + n_const,
                             constraints.data() + j * n_const);
        }
        std::swap(flags[i], flags[j]);
    }

    /// Invalidate all evaluations (mark as not evaluated)
    void invalidate_all() { std::fill(flags.begin(), flags.end(), IndFlags::None); }

    /// Total memory used in bytes
    [[nodiscard]] size_t memory_bytes() const {
        return genes.capacity() * sizeof(gene_type)
             + objectives.capacity() * sizeof(double)
             + constraints.capacity() * sizeof(double)
             + flags.capacity() * sizeof(IndFlags)
             + lower_bounds.capacity() * sizeof(double)
             + upper_bounds.capacity() * sizeof(double);
    }
};

} // namespace ea
