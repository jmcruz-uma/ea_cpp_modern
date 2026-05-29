#pragma once
/// @file position_based.hpp
/// @brief Position-Based Crossover for permutation encodings.
/// Reference: jMetal PositionBasedCrossover — https://github.com/jMetal/jMetal
/// C++23 with deducing this, SoA Population<> access.

#include <algorithm>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <vector>

namespace ea {

/// Position-Based Crossover for permutation-based encodings.
/// Selects random positions to inherit genes from one parent.
/// Remaining genes are filled from the other parent in order of appearance.
/// Genes are integers stored as doubles in the Population<> array.
struct PositionBasedCrossover {
    double crossover_probability = 1.0; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Permutation; }

    /// Apply Position-Based crossover. Produces 2 children starting at child_start.
    void apply(this PositionBasedCrossover& self, Population<>& pop, int parent_a, int parent_b,
               int child_start) {
        auto& rng = Random::instance();
        int dim = pop.dim;

        // Copy parents to children (fallback)
        for (int j = 0; j < dim; ++j) {
            pop.gene(child_start, j) = pop.gene(parent_a, j);
            pop.gene(child_start + 1, j) = pop.gene(parent_b, j);
        }

        if (dim < 2 || !rng.coin_flip(self.crossover_probability)) {
            pop.set_evaluated(child_start, false);
            pop.set_evaluated(child_start + 1, false);
            return;
        }

        // Build child 1: selected positions from parent_a, rest from parent_b
        self.build_child(pop, parent_a, parent_b, child_start, rng);
        // Build child 2: selected positions from parent_b, rest from parent_a
        self.build_child(pop, parent_b, parent_a, child_start + 1, rng);

        pop.set_evaluated(child_start, false);
        pop.set_evaluated(child_start + 1, false);
    }

private:
    /// Build a child using position-based crossover from two parents.
    void build_child(this PositionBasedCrossover&, Population<>& pop, int parent_from, int parent_to,
                     int child_idx, Random& rng) {
        int dim = pop.dim;

        // Select random positions to keep from parent_from
        int num_positions = rng.uniform_int(1, dim / 2);
        std::vector<bool> selected(dim, false);
        int count = 0;
        while (count < num_positions) {
            int pos = rng.uniform_int(0, dim - 1);
            if (!selected[pos]) {
                selected[pos] = true;
                ++count;
            }
        }

        // Place selected positions from parent_from
        std::vector<int> child(dim, -1);
        for (int i = 0; i < dim; ++i) {
            if (selected[i]) {
                child[i] = static_cast<int>(pop.gene(parent_from, i));
            }
        }

        // Fill remaining positions from parent_to in order of appearance
        int child_idx_pos = 0;
        for (int j = 0; j < dim; ++j) {
            int val = static_cast<int>(pop.gene(parent_to, j));
            // Check if val is already in child
            bool already_present = false;
            for (int k = 0; k < dim; ++k) {
                if (child[k] == val) {
                    already_present = true;
                    break;
                }
            }
            if (!already_present) {
                // Find first empty position
                while (child_idx_pos < dim && child[child_idx_pos] != -1) {
                    ++child_idx_pos;
                }
                if (child_idx_pos < dim) {
                    child[child_idx_pos] = val;
                }
            }
        }

        // Write child to population
        for (int j = 0; j < dim; ++j) {
            pop.gene(child_idx, j) = static_cast<double>(child[j]);
        }
    }
};

} // namespace ea
