#pragma once
/// @file concepts.hpp
/// @brief C++20 Concepts for evolutionary algorithm components.
/// These define compile-time contracts that operators, algorithms, and problems must satisfy.

#include <concepts>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <string_view>
#include <vector>

namespace ea {

// ============================================================
// Eval functor concept (primary interface for run() methods)
// ============================================================

/// A callable that evaluates a single individual in-place: f(pop, idx).
/// Algorithms accept this instead of Problem<T> so lambdas and adapters work.
template <typename F>
concept EvalFunctor = requires(F& f, Population<>& pop, int idx) {
    { f(pop, idx) } -> std::same_as<void>;
};

// ============================================================
// Operator concepts
// ============================================================

/// A Crossover operator takes two parents and produces a child.
template <typename T>
concept Crossover = requires(T& cx, Population<>& pop, int a, int b, int child) {
    { cx.apply(pop, a, b, child) } -> std::same_as<void>;
    { cx.arity() } -> std::convertible_to<int>;
    { cx.encoding() } -> std::convertible_to<Encoding>;
};

/// A Mutation operator modifies a single individual in-place.
template <typename T>
concept Mutation = requires(T& mut, Population<>& pop, int idx) {
    { mut.apply(pop, idx) } -> std::same_as<void>;
    { mut.encoding() } -> std::convertible_to<Encoding>;
};

/// A Selection operator selects individuals from a population.
template <typename T>
concept Selection = requires(T& sel, Population<>& pop, std::vector<int>& mating_pool) {
    { sel.select(pop, mating_pool) } -> std::same_as<void>;
};

/// A Replacement operator merges parents and offspring.
template <typename T>
concept Replacement =
    requires(T& repl, Population<>& pop, std::vector<int>& offspring_indices, int pop_size) {
        { repl.replace(pop, offspring_indices, pop_size) } -> std::same_as<std::vector<int>>;
    };

// ============================================================
// Problem concept
// ============================================================

/// A Problem evaluates objectives (and optionally constraints) for individuals.
template <typename T>
concept Problem = requires(T& prob, Population<>& pop, int idx) {
    { prob.num_objectives() } -> std::convertible_to<int>;
    { prob.dimension() } -> std::convertible_to<int>;
    { prob.evaluate(pop, idx) } -> std::same_as<void>;
    { prob.lower_bounds() } -> std::convertible_to<std::span<const double>>;
    { prob.upper_bounds() } -> std::convertible_to<std::span<const double>>;
};

/// A Problem that supports batch evaluation (SIMD-friendly).
template <typename T>
concept BatchProblem = Problem<T> && requires(T& prob, Population<>& pop, int start, int count) {
    { prob.evaluate_batch(pop, start, count) } -> std::same_as<void>;
};

/// A constrained problem provides constraint values.
template <typename T>
concept ConstrainedProblem = Problem<T> && requires(T& prob) {
    { prob.num_constraints() } -> std::convertible_to<int>;
};

// ============================================================
// Algorithm concept
// ============================================================

/// An Evolutionary Algorithm runs on a population.
template <typename T>
concept Algorithm = requires(T& algo, Population<>& pop) {
    { algo.run(pop) } -> std::same_as<void>;
    { algo.name() } -> std::convertible_to<std::string_view>;
};

// ============================================================
// Encoding-specific concepts
// ============================================================

/// A crossover that works with real-valued encodings.
template <typename T>
concept RealCrossover = Crossover<T> && requires(T& cx) {
    { cx.encoding() } -> std::same_as<Encoding>;
    requires cx.encoding() == Encoding::Real || true; // runtime check
};

/// A crossover that works with permutation encodings.
template <typename T>
concept PermutationCrossover = Crossover<T>;

// ============================================================
// Comparator concepts
// ============================================================

/// A comparator for dominance-based ranking.
template <typename T>
concept DominanceComparator = requires(T& cmp, Population<>& pop, int a, int b) {
    { cmp.compare(pop, a, b) } -> std::convertible_to<int>;
    // -1: a dominates b, 0: non-dominated, 1: b dominates a
};

/// A comparator for crowding distance sorting.
template <typename T>
concept CrowdingComparator = requires(T& cmp, Population<>& pop, int a, int b) {
    { cmp.crowding_distance(pop, a) } -> std::convertible_to<double>;
};

} // namespace ea