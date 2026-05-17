#pragma once
/// @file registry.hpp
/// @brief Self-registration of all built-in operators and problems.
///
/// Include this ONCE in your main.cpp (or a dedicated .cpp) to populate
/// the factories. After that, you can create operators by name:
///
///   #include <ea/core/registry.hpp>
///   int main() {
///       auto cx   = ea::CrossoverFactory::create("sbx");
///       auto mt   = ea::MutationFactory::create("polynomial");
///       auto prob = ea::ProblemFactory::create("zdt1");
///   }
///
/// Only operators with standard signatures are registered here.
/// Special-signature operators (DE, ParentCentric, rank-based selection)
/// remain compile-time only — they are algorithm-internal.

#include <ea/core/factory_types.hpp>

// Crossover operators (2-parent standard signature)
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/crossover/blx_alpha.hpp>
#include <ea/operator/crossover/blx_alpha_beta.hpp>
#include <ea/operator/crossover/arithmetic.hpp>
#include <ea/operator/crossover/whole_arithmetic.hpp>
#include <ea/operator/crossover/laplace.hpp>
#include <ea/operator/crossover/pmx.hpp>
#include <ea/operator/crossover/cycle.hpp>
#include <ea/operator/crossover/edge_recombination.hpp>
#include <ea/operator/crossover/oxd.hpp>
#include <ea/operator/crossover/position_based.hpp>
#include <ea/operator/crossover/uniform.hpp>

// Mutation operators
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/operator/mutation/bit_flip.hpp>
#include <ea/operator/mutation/swap.hpp>
#include <ea/operator/mutation/inversion.hpp>
#include <ea/operator/mutation/insert.hpp>
#include <ea/operator/mutation/scramble.hpp>
#include <ea/operator/mutation/uniform.hpp>
#include <ea/operator/mutation/non_uniform.hpp>
#include <ea/operator/mutation/power_law.hpp>
#include <ea/operator/mutation/levy_flight.hpp>
#include <ea/operator/mutation/simple_random.hpp>

// Problems
#include <ea/problem/zdt.hpp>
#include <ea/problem/dtlz.hpp>

namespace ea::registry {

// ============================================================
// Crossover (2-parent standard signature)
// ============================================================
EA_REGISTER_CROSSOVER("sbx",                SBXCrossover);
EA_REGISTER_CROSSOVER("blx_alpha",           BLXAlphaCrossover);
EA_REGISTER_CROSSOVER("blx_alpha_beta",      BLXAlphaBetaCrossover);
EA_REGISTER_CROSSOVER("arithmetic",          ArithmeticCrossover);
EA_REGISTER_CROSSOVER("whole_arithmetic",     WholeArithmeticCrossover);
EA_REGISTER_CROSSOVER("laplace",             LaplaceCrossover);
EA_REGISTER_CROSSOVER("pmx",                 PMXCrossover);
EA_REGISTER_CROSSOVER("cycle",               CycleCrossover);
EA_REGISTER_CROSSOVER("edge_recombination",  EdgeRecombinationCrossover);
EA_REGISTER_CROSSOVER("oxd",                 OXDCrossover);
EA_REGISTER_CROSSOVER("position_based",       PositionBasedCrossover);
EA_REGISTER_CROSSOVER("uniform_cx",           UniformCrossover);

// Note: DECrossover and ParentCentricCrossover have non-standard
// signatures (pointer-to-parents, 3-parent). They remain compile-time
// only. Users who need them via factory should provide their own adapter.

// ============================================================
// Mutation
// ============================================================
EA_REGISTER_MUTATION("polynomial",           PolynomialMutation);
EA_REGISTER_MUTATION("bit_flip",             BitFlipMutation);
EA_REGISTER_MUTATION("swap",                 SwapMutation);
EA_REGISTER_MUTATION("inversion",            InversionMutation);
EA_REGISTER_MUTATION("insert",               InsertMutation);
EA_REGISTER_MUTATION("scramble",            ScrambleMutation);
EA_REGISTER_MUTATION("uniform_mt",           UniformMutation);
EA_REGISTER_MUTATION("non_uniform",          NonUniformMutation);
EA_REGISTER_MUTATION("power_law",           PowerLawMutation);
EA_REGISTER_MUTATION("levy_flight",          LevyFlightMutation);
EA_REGISTER_MUTATION("simple_random",        SimpleRandomMutation);

// ============================================================
// Problems
// ============================================================
EA_REGISTER_PROBLEM("zdt1", ZDT1);
EA_REGISTER_PROBLEM("zdt2", ZDT2);
EA_REGISTER_PROBLEM("zdt3", ZDT3);
EA_REGISTER_PROBLEM("zdt4", ZDT4);
EA_REGISTER_PROBLEM("zdt6", ZDT6);
EA_REGISTER_PROBLEM("dtlz1", DTLZ1);
EA_REGISTER_PROBLEM("dtlz2", DTLZ2);
EA_REGISTER_PROBLEM("dtlz3", DTLZ3);
EA_REGISTER_PROBLEM("dtlz4", DTLZ4);
EA_REGISTER_PROBLEM("dtlz5", DTLZ5);
EA_REGISTER_PROBLEM("dtlz6", DTLZ6);
EA_REGISTER_PROBLEM("dtlz7", DTLZ7);

} // namespace ea::registry