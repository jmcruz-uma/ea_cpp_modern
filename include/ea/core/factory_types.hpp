#pragma once
/// @file factory_types.hpp
/// @brief Polymorphic base types and factory aliases for runtime dispatch.
///
/// Provides type-erased wrappers (adapters) that bridge concept-based
/// concrete operators to virtual-dispatch base interfaces.
///
/// Design decision: Only operators with "standard" signatures are
/// registered in the factory. Special-sigature operators (DECrossover,
/// ParentCentricCrossover, selection with ranks/crowding) remain
/// compile-time only — they are algorithm-specific and don't benefit
/// from runtime dispatch.
///
/// Zero cost when using compile-time composition (existing template path).
/// Small virtual dispatch cost only when using factory-based runtime creation.

#include <ea/core/encoding.hpp>
#include <ea/core/factory.hpp>
#include <ea/core/population.hpp>
#include <memory>
#include <string_view>
#include <vector>

namespace ea {

// ============================================================
// Polymorphic bases
// ============================================================

/// Base for crossover operators with standard 2-parent signature.
struct CrossoverBase {
    virtual ~CrossoverBase() = default;
    virtual void apply(Population& pop, int parent_a, int parent_b, int child_start) = 0;
    virtual int arity() const = 0;
    virtual Encoding encoding() const = 0;
};

/// Base for mutation operators.
struct MutationBase {
    virtual ~MutationBase() = default;
    virtual void apply(Population& pop, int idx) = 0;
    virtual Encoding encoding() const = 0;
};

/// Base for problems.
struct ProblemBase {
    virtual ~ProblemBase() = default;
    virtual void evaluate(Population& pop, int idx) = 0;
    virtual int num_objectives() const = 0;
    virtual int dimension() const = 0;
    virtual std::string_view name() const = 0;
};

// ============================================================
// Factory aliases
// ============================================================

using CrossoverFactory = Factory<CrossoverBase, 32>;
using MutationFactory = Factory<MutationBase, 32>;
using ProblemFactory = Factory<ProblemBase, 64>;

// ============================================================
// Adapters: concrete → base (type erasure)
// ============================================================

/// Wraps a Crossover-concept type (2-parent signature) into CrossoverBase.
template <typename T> struct CrossoverAdapter : CrossoverBase {
    T op;
    explicit CrossoverAdapter(T t = T{})
        : op(std::move(t)) {}
    void apply(Population& pop, int a, int b, int child) override { op.apply(pop, a, b, child); }
    int arity() const override { return op.arity(); }
    Encoding encoding() const override { return op.encoding(); }
};

/// Wraps a Mutation-concept type into MutationBase.
template <typename T> struct MutationAdapter : MutationBase {
    T op;
    explicit MutationAdapter(T t = T{})
        : op(std::move(t)) {}
    void apply(Population& pop, int idx) override { op.apply(pop, idx); }
    Encoding encoding() const override { return op.encoding(); }
};

/// Wraps a Problem-concept type into ProblemBase.
/// Falls back to static name if the problem type doesn't have name().
template <typename T, typename = void> struct ProblemAdapter : ProblemBase {
    T prob;
    explicit ProblemAdapter(T t = T{})
        : prob(std::move(t)) {}
    void evaluate(Population& pop, int idx) override {
        const_cast<const T&>(prob).evaluate(pop, idx);
    }
    int num_objectives() const override { return prob.num_objectives(); }
    int dimension() const override { return prob.dimension(); }
    std::string_view name() const override { return "unknown"; }
};

// ============================================================
// Registration macros (adapter-aware)
// ============================================================

/// Register a 2-parent crossover type in CrossoverFactory.
#define EA_REGISTER_CROSSOVER(Name, Type)                                                      \
    static inline const bool ea_reg_cx_##Type =                                                \
        ea::CrossoverFactory::register_type(Name, []() -> std::unique_ptr<ea::CrossoverBase> { \
            return std::make_unique<ea::CrossoverAdapter<Type>>();                             \
        })

/// Register a mutation type in MutationFactory.
#define EA_REGISTER_MUTATION(Name, Type)                                                     \
    static inline const bool ea_reg_mt_##Type =                                              \
        ea::MutationFactory::register_type(Name, []() -> std::unique_ptr<ea::MutationBase> { \
            return std::make_unique<ea::MutationAdapter<Type>>();                            \
        })

/// Register a problem type in ProblemFactory.
#define EA_REGISTER_PROBLEM(Name, Type)                                                    \
    static inline const bool ea_reg_prob_##Type =                                          \
        ea::ProblemFactory::register_type(Name, []() -> std::unique_ptr<ea::ProblemBase> { \
            return std::make_unique<ea::ProblemAdapter<Type>>();                           \
        })

} // namespace ea