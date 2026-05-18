#pragma once
/// @file de.hpp
/// @brief Differential Evolution Crossover for real-valued encodings.
/// Reference: jMetal DifferentialEvolutionCrossover.java
///
/// Implements DE/rand/1/bin, DE/best/1/bin, DE/rand-to-best/1/bin,
/// DE/current-to-rand/1, and exponential (EXP) variants.
///
/// DE combines mutation and crossover: a mutant vector is created from
/// parent differences, then binomially or exponentially crossed with the
/// target vector.
///
/// Reference: Storn, R., & Price, K. (1997). Differential evolution—a simple
/// and efficient heuristic for global optimization over continuous spaces.
/// Journal of Global Optimization, 11(4), 341-359.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <stdexcept>
#include <string_view>

namespace ea {

/// DE variant types
enum class DEVariant : uint8_t {
    Rand1Bin,
    Rand1Exp,
    Rand2Bin,
    Rand2Exp,
    Best1Bin,
    Best1Exp,
    Best2Bin,
    Best2Exp,
    RandToBest1Bin,
    RandToBest1Exp,
    CurrentToRand1Bin,
    CurrentToRand1Exp,
};

/// Crossover type: binomial or exponential
enum class DECrossoverType : uint8_t { Bin, Exp };

/// Mutation type for DE
enum class DEMutationType : uint8_t { Rand, Best, RandToBest, CurrentToRand };

/// Parse DE variant from string
constexpr DEVariant de_variant_from_string(std::string_view s) {
    if (s == "RAND_1_BIN")
        return DEVariant::Rand1Bin;
    if (s == "RAND_1_EXP")
        return DEVariant::Rand1Exp;
    if (s == "RAND_2_BIN")
        return DEVariant::Rand2Bin;
    if (s == "RAND_2_EXP")
        return DEVariant::Rand2Exp;
    if (s == "BEST_1_BIN")
        return DEVariant::Best1Bin;
    if (s == "BEST_1_EXP")
        return DEVariant::Best1Exp;
    if (s == "BEST_2_BIN")
        return DEVariant::Best2Bin;
    if (s == "BEST_2_EXP")
        return DEVariant::Best2Exp;
    if (s == "RAND_TO_BEST_1_BIN")
        return DEVariant::RandToBest1Bin;
    if (s == "RAND_TO_BEST_1_EXP")
        return DEVariant::RandToBest1Exp;
    if (s == "CURRENT_TO_RAND_1_BIN")
        return DEVariant::CurrentToRand1Bin;
    if (s == "CURRENT_TO_RAND_1_EXP")
        return DEVariant::CurrentToRand1Exp;
    throw std::invalid_argument("Invalid DE variant string");
}

/// Differential Evolution Crossover for real-valued encodings.
/// Supports multiple mutation and crossover strategies.
struct DECrossover {
    double cr = 0.5; ///< Crossover rate [0, 1]
    double f = 0.5;  ///< Differential weight (scaling factor)
    DEVariant variant = DEVariant::Rand1Bin;

    static constexpr int arity() { return 4; } // Max: target + 3 difference vectors
    static constexpr Encoding encoding() { return Encoding::Real; }

    /// Number of difference vectors for this variant
    static constexpr int num_difference_vectors(DEVariant v) {
        switch (v) {
        case DEVariant::Rand1Bin:
        case DEVariant::Rand1Exp:
        case DEVariant::Best1Bin:
        case DEVariant::Best1Exp:
        case DEVariant::RandToBest1Bin:
        case DEVariant::RandToBest1Exp:
        case DEVariant::CurrentToRand1Bin:
        case DEVariant::CurrentToRand1Exp:
            return 1;
        case DEVariant::Rand2Bin:
        case DEVariant::Rand2Exp:
        case DEVariant::Best2Bin:
        case DEVariant::Best2Exp:
            return 2;
        }
        return 1;
    }

    /// Crossover type for this variant
    static constexpr DECrossoverType crossover_type(DEVariant v) {
        switch (v) {
        case DEVariant::Rand1Bin:
        case DEVariant::Best1Bin:
        case DEVariant::RandToBest1Bin:
        case DEVariant::CurrentToRand1Bin:
        case DEVariant::Rand2Bin:
        case DEVariant::Best2Bin:
            return DECrossoverType::Bin;
        case DEVariant::Rand1Exp:
        case DEVariant::Best1Exp:
        case DEVariant::RandToBest1Exp:
        case DEVariant::CurrentToRand1Exp:
        case DEVariant::Rand2Exp:
        case DEVariant::Best2Exp:
            return DECrossoverType::Exp;
        }
        return DECrossoverType::Bin;
    }

    /// Mutation type for this variant
    static constexpr DEMutationType mutation_type(DEVariant v) {
        switch (v) {
        case DEVariant::Rand1Bin:
        case DEVariant::Rand1Exp:
        case DEVariant::Rand2Bin:
        case DEVariant::Rand2Exp:
            return DEMutationType::Rand;
        case DEVariant::Best1Bin:
        case DEVariant::Best1Exp:
        case DEVariant::Best2Bin:
        case DEVariant::Best2Exp:
            return DEMutationType::Best;
        case DEVariant::CurrentToRand1Bin:
        case DEVariant::CurrentToRand1Exp:
            return DEMutationType::CurrentToRand;
        case DEVariant::RandToBest1Bin:
        case DEVariant::RandToBest1Exp:
            return DEMutationType::RandToBest;
        }
        return DEMutationType::Rand;
    }

    /// Number of required parents for this variant
    static constexpr int required_parents(DEVariant v) { return 1 + num_difference_vectors(v) * 2; }

    /// Apply DE crossover. Produces 1 child at child_start.
    /// Parents: [0] = target/current, [1..n] = difference vector parents.
    /// best_idx is the index of the best solution in the population (for BEST variants).
    void apply(this auto& self, Population& pop, const int* parent_indices, int num_parents,
               int child_start, int best_idx = -1) {
        auto& rng = Random::instance();
        (void)num_parents; // may be unused in some paths

        // Copy target vector as trial vector base
        int target = parent_indices[0];
        for (int j = 0; j < pop.dim; ++j) {
            pop.gene(child_start, j) = pop.gene(target, j);
        }

        auto cx_type = crossover_type(self.variant);
        auto mut_type = mutation_type(self.variant);
        int n_diff = num_difference_vectors(self.variant);

        if (cx_type == DECrossoverType::Bin) {
            // Binomial crossover: at least one gene from mutant
            int jrand = rng.uniform_int(0, pop.dim - 1);

            for (int j = 0; j < pop.dim; ++j) {
                if (rng.uniform() < self.cr || j == jrand) {
                    double mutant = self.mutate(pop, parent_indices, j, mut_type, best_idx, n_diff);
                    double lb = pop.lower_bounds[j];
                    double ub = pop.upper_bounds[j];
                    pop.gene(child_start, j) = std::clamp(mutant, lb, ub);
                }
            }
        } else {
            // Exponential crossover: starting at random position, copy contiguous genes
            int j = rng.uniform_int(0, pop.dim - 1);
            int l = 0;

            do {
                double mutant = self.mutate(pop, parent_indices, j, mut_type, best_idx, n_diff);
                double lb = pop.lower_bounds[j];
                double ub = pop.upper_bounds[j];
                pop.gene(child_start, j) = std::clamp(mutant, lb, ub);

                j = (j + 1) % pop.dim;
                l++;
            } while (rng.uniform() < self.cr && l < pop.dim);
        }

        pop.set_evaluated(child_start, false);
    }

private:
    /// Compute mutant value for a single gene
    double mutate(this auto& self, Population& pop, const int* parent_indices, int gene_idx,
                  DEMutationType mut_type, int best_idx, int n_diff) {
        switch (mut_type) {
        case DEMutationType::Rand: {
            // rand/1: p[2] + F * (p[0] - p[1])
            // rand/2: p[4] + F * (p[0] - p[1]) + F * (p[2] - p[3])
            int base = 2 * n_diff; // base vector index
            double value = pop.gene(parent_indices[base], gene_idx);
            for (int d = 0; d < n_diff; ++d) {
                value += self.f * (pop.gene(parent_indices[2 * d], gene_idx) -
                                   pop.gene(parent_indices[2 * d + 1], gene_idx));
            }
            return value;
        }
        case DEMutationType::Best: {
            // best/1: best + F * (p[0] - p[1])
            // best/2: best + F * (p[0] - p[1]) + F * (p[2] - p[3])
            double value = pop.gene(best_idx, gene_idx);
            for (int d = 0; d < n_diff; ++d) {
                value += self.f * (pop.gene(parent_indices[2 * d], gene_idx) -
                                   pop.gene(parent_indices[2 * d + 1], gene_idx));
            }
            return value;
        }
        case DEMutationType::RandToBest: {
            // rand-to-best/1: target + F * (best - target) + F * (p[0] - p[1])
            int target = parent_indices[0];
            return pop.gene(target, gene_idx) +
                   self.f * (pop.gene(best_idx, gene_idx) - pop.gene(target, gene_idx)) +
                   self.f * (pop.gene(parent_indices[0], gene_idx) -
                             pop.gene(parent_indices[1], gene_idx));
        }
        case DEMutationType::CurrentToRand: {
            // current-to-rand/1: target + K * (p[0] - target) + F * (p[1] - p[2])
            // where K = F (simplified as in jMetal)
            int target = parent_indices[0];
            double k = self.f;
            return pop.gene(target, gene_idx) +
                   k * (pop.gene(parent_indices[0], gene_idx) - pop.gene(target, gene_idx)) +
                   self.f * (pop.gene(parent_indices[1], gene_idx) -
                             pop.gene(parent_indices[2], gene_idx));
        }
        }
        return pop.gene(parent_indices[0], gene_idx);
    }
};

} // namespace ea
