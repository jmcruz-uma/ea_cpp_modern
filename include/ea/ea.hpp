#pragma once
/// @file ea.hpp
/// @brief Main include — pulls in all ea-cpp components.

// Core
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/core/concepts.hpp>
#include <ea/core/comparator.hpp>

// Util
#include <ea/util/random.hpp>

// Operators — Crossover
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/crossover/blx_alpha.hpp>
#include <ea/operator/crossover/blx_alpha_beta.hpp>
#include <ea/operator/crossover/arithmetic.hpp>
#include <ea/operator/crossover/whole_arithmetic.hpp>
#include <ea/operator/crossover/laplace.hpp>
#include <ea/operator/crossover/de.hpp>
#include <ea/operator/crossover/parent_centric.hpp>
#include <ea/operator/crossover/pmx.hpp>
#include <ea/operator/crossover/cycle.hpp>
#include <ea/operator/crossover/edge_recombination.hpp>
#include <ea/operator/crossover/oxd.hpp>
#include <ea/operator/crossover/position_based.hpp>
#include <ea/operator/crossover/uniform.hpp>

// Operators — Mutation
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

// Operators — Selection
#include <ea/operator/selection/binary_tournament.hpp>
#include <ea/operator/selection/nary_tournament.hpp>
#include <ea/operator/selection/random.hpp>

// Operators — Replacement
#include <ea/operator/replacement/mu_plus_lambda.hpp>
#include <ea/operator/replacement/mu_comma_lambda.hpp>

// Problems
#include <ea/problem/zdt.hpp>
#include <ea/problem/dtlz.hpp>

// Indicators
#include <ea/indicator/hypervolume.hpp>
#include <ea/indicator/igd.hpp>

// Algorithms
#include <ea/algorithm/nsga2.hpp>
#include <ea/algorithm/moead.hpp>
#include <ea/algorithm/spea2.hpp>
#include <ea/algorithm/smsemoa.hpp>
#include <ea/algorithm/gde3.hpp>

// Aggregation utilities
#include <ea/util/aggregation.hpp>