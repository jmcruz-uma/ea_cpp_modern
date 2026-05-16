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
#include <ea/operator/crossover/de.hpp>
#include <ea/operator/crossover/pmx.hpp>
#include <ea/operator/crossover/cycle.hpp>
#include <ea/operator/crossover/uniform.hpp>

// Operators — Mutation
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/operator/mutation/bit_flip.hpp>
#include <ea/operator/mutation/swap.hpp>
#include <ea/operator/mutation/inversion.hpp>
#include <ea/operator/mutation/uniform.hpp>
#include <ea/operator/mutation/non_uniform.hpp>

// Operators — Selection
#include <ea/operator/selection/binary_tournament.hpp>
#include <ea/operator/selection/nary_tournament.hpp>
#include <ea/operator/selection/random.hpp>

// Operators — Replacement
#include <ea/operator/replacement/mu_plus_lambda.hpp>
#include <ea/operator/replacement/nsga2_replacement.hpp>

// Problems
#include <ea/problem/zdt.hpp>
#include <ea/problem/dtlz.hpp>

// Algorithms
#include <ea/algorithm/nsga2.hpp>

// Indicators
// TODO: #include <ea/indicator/hypervolume.hpp>
// TODO: #include <ea/indicator/igd.hpp>