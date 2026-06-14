// tests/unit/test_concepts.cpp
//
// Compile-time concept verification — zero runtime cost.
// Every static_assert here fires at compile time; if any fails the whole build
// fails with a clear message pointing at the offending line.
//
// Strategy:
//   Positive checks: real operators must satisfy their concept.
//   Negative checks: unrelated types must NOT satisfy the concept.
//
// Reuse pattern: add a new operator/algorithm to this file to verify it
// satisfies the intended concept before integrating it into an algorithm.
//
// Build (standalone):
//   g++-14 -std=c++23 -O0 -I include tests/unit/test_concepts.cpp -o /dev/null

#include <ea/core/concepts.hpp>
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/crossover/blx_alpha.hpp>
#include <ea/operator/crossover/blx_alpha_beta.hpp>
#include <ea/operator/crossover/de.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/operator/mutation/gaussian.hpp>
#include <ea/operator/mutation/bit_flip.hpp>
#include <ea/operator/selection/binary_tournament.hpp>

// ============================================================
// Minimal structs for negative tests
// ============================================================

namespace {

// Satisfies EvalFunctor — the minimal callable pattern used with run()
struct MinimalEvalFn {
    void operator()(ea::Population<>& pop, int idx) { (void)pop; (void)idx; }
};

// Wrong return type for EvalFunctor (returns int instead of void)
struct BadEvalFnReturnType {
    int operator()(ea::Population<>& pop, int idx) { (void)pop; return idx; }
};

// Empty struct — satisfies nothing
struct Empty {};

// Struct with almost-right Crossover API — missing arity()
struct AlmostCrossover {
    static constexpr ea::Encoding encoding() { return ea::Encoding::Real; }
    void apply(ea::Population<>& pop, int a, int b, int c) { (void)pop; (void)a; (void)b; (void)c; }
};

// Struct with almost-right Mutation API — missing encoding()
struct AlmostMutation {
    void apply(ea::Population<>& pop, int idx) { (void)pop; (void)idx; }
};

} // namespace

// ============================================================
// EvalFunctor — used as run() parameter constraint
// ============================================================

// Positive: lambda-compatible struct
static_assert(ea::EvalFunctor<MinimalEvalFn>,
    "MinimalEvalFn must satisfy EvalFunctor");

// Positive: lambdas — verified via a type-erased functor stand-in
// (lambdas can't be named at namespace scope, but MinimalEvalFn covers the pattern)

// Negative: wrong return type
static_assert(!ea::EvalFunctor<BadEvalFnReturnType>,
    "EvalFunctor must require void return");

// Negative: empty struct
static_assert(!ea::EvalFunctor<Empty>,
    "Empty struct must not satisfy EvalFunctor");

// Negative: Crossover operators are not EvalFunctors
static_assert(!ea::EvalFunctor<ea::SBXCrossover>,
    "SBXCrossover must not satisfy EvalFunctor");

// ============================================================
// Crossover
// ============================================================

// Positive: standard crossover operators
static_assert(ea::Crossover<ea::SBXCrossover>,
    "SBXCrossover must satisfy Crossover");

static_assert(ea::Crossover<ea::BLXAlphaCrossover>,
    "BLXAlphaCrossover must satisfy Crossover");

static_assert(ea::Crossover<ea::BLXAlphaBetaCrossover>,
    "BLXAlphaBetaCrossover must satisfy Crossover");

// Note: DECrossover intentionally does NOT satisfy Crossover — its apply() takes
// (pop, const int* parent_indices, int num_parents, int child), a non-binary API.
// GDE3 and MOEAD_DE use it directly as a concrete member, bypassing the concept.
static_assert(!ea::Crossover<ea::DECrossover>,
    "DECrossover must not satisfy standard binary Crossover (uses n-ary API)");

// Negative: mutation operators are not crossovers
static_assert(!ea::Crossover<ea::PolynomialMutation>,
    "PolynomialMutation must not satisfy Crossover");

static_assert(!ea::Crossover<ea::GaussianMutation>,
    "GaussianMutation must not satisfy Crossover");

// Negative: missing arity()
static_assert(!ea::Crossover<AlmostCrossover>,
    "AlmostCrossover (missing arity()) must not satisfy Crossover");

// Negative: empty struct
static_assert(!ea::Crossover<Empty>,
    "Empty struct must not satisfy Crossover");

// ============================================================
// Mutation
// ============================================================

// Positive: standard mutation operators
static_assert(ea::Mutation<ea::PolynomialMutation>,
    "PolynomialMutation must satisfy Mutation");

static_assert(ea::Mutation<ea::GaussianMutation>,
    "GaussianMutation must satisfy Mutation");

static_assert(ea::Mutation<ea::BitFlipMutation>,
    "BitFlipMutation must satisfy Mutation");

// Negative: crossover operators are not mutations
static_assert(!ea::Mutation<ea::SBXCrossover>,
    "SBXCrossover must not satisfy Mutation");

static_assert(!ea::Mutation<ea::DECrossover>,
    "DECrossover must not satisfy Mutation");

// Negative: missing encoding()
static_assert(!ea::Mutation<AlmostMutation>,
    "AlmostMutation (missing encoding()) must not satisfy Mutation");

// Negative: empty struct
static_assert(!ea::Mutation<Empty>,
    "Empty struct must not satisfy Mutation");

// ============================================================
// Algorithm template instantiation smoke-test
// Verifying that constrained templates accept valid operator types.
// These instantiations would fail to compile if constraints are violated.
// ============================================================

#include <ea/algorithm/nsga2.hpp>
#include <ea/algorithm/nsga3.hpp>
#include <ea/algorithm/spea2.hpp>
#include <ea/algorithm/paes.hpp>
#include <ea/algorithm/smsemoa.hpp>
#include <ea/algorithm/agemoea.hpp>
#include <ea/algorithm/mocell.hpp>
#include <ea/algorithm/moead.hpp>
#include <ea/algorithm/ibea.hpp>

// Each of these template instantiations confirms that the algorithm accepts
// valid operator types under the new Crossover/Mutation constraints.
static_assert(sizeof(ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation>) > 0,
    "NSGAII<SBX, PM> must instantiate");

static_assert(sizeof(ea::NSGAIII<ea::SBXCrossover, ea::PolynomialMutation>) > 0,
    "NSGAIII<SBX, PM> must instantiate");

static_assert(sizeof(ea::SPEA2<ea::SBXCrossover, ea::PolynomialMutation>) > 0,
    "SPEA2<SBX, PM> must instantiate");

static_assert(sizeof(ea::PAES<ea::PolynomialMutation>) > 0,
    "PAES<PM> must instantiate");

static_assert(sizeof(ea::SMSEMOA<ea::SBXCrossover, ea::PolynomialMutation>) > 0,
    "SMSEMOA<SBX, PM> must instantiate");

static_assert(sizeof(ea::AGEMOEA<ea::SBXCrossover, ea::PolynomialMutation>) > 0,
    "AGEMOEA<SBX, PM> must instantiate");

static_assert(sizeof(ea::MOCell<ea::SBXCrossover, ea::PolynomialMutation>) > 0,
    "MOCell<SBX, PM> must instantiate");

static_assert(sizeof(ea::MOEAD<ea::SBXCrossover, ea::PolynomialMutation>) > 0,
    "MOEAD<SBX, PM> must instantiate");

static_assert(sizeof(ea::IBEA<ea::SBXCrossover, ea::PolynomialMutation>) > 0,
    "IBEA<SBX, PM> must instantiate");

int main() { return 0; }
