#pragma once
/// @file pmx.hpp
/// @brief Partially Mapped Crossover (PMX) for permutation encodings.
/// Reference: jMetal PMXCrossover — https://github.com/jMetal/jMetal
/// C++23 with deducing this, SoA Population access.

#include <algorithm>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <vector>

namespace ea {

/// Partially Mapped Crossover (PMX) for permutation-based encodings.
/// Produces two children from two parents while preserving valid permutations.
/// Genes are integers stored as doubles in the Population array.
struct PMXCrossover {
    double crossover_probability = 1.0; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Permutation; }

    /// Apply PMX crossover. Produces 2 children starting at child_start.
    void apply(this PMXCrossover& self, Population& pop, int parent_a, int parent_b,
               int child_start) {
        auto& rng = Random::instance();
        int dim = pop.dim;

        // Copy parents to children
        for (int j = 0; j < dim; ++j) {
            pop.gene(child_start, j) = pop.gene(parent_a, j);
            pop.gene(child_start + 1, j) = pop.gene(parent_b, j);
        }

        if (dim < 2 || !rng.coin_flip(self.crossover_probability)) {
            pop.set_evaluated(child_start, false);
            pop.set_evaluated(child_start + 1, false);
            return;
        }

        // STEP 1: Select two distinct cutting points
        int cp1 = rng.uniform_int(0, dim - 1);
        int cp2 = rng.uniform_int(0, dim - 1);
        while (cp2 == cp1)
            cp2 = rng.uniform_int(0, dim - 1);
        if (cp1 > cp2)
            std::swap(cp1, cp2);

        // STEP 2: Build replacement maps for repair
        // replacement maps from parent values to child values
        std::vector<int> repl1(dim, -1); // map parent_a value -> parent_b value
        std::vector<int> repl2(dim, -1); // map parent_b value -> parent_a value

        // STEP 3: Swap the segment between cutting points
        for (int i = cp1; i <= cp2; ++i) {
            int val_a = static_cast<int>(pop.gene(parent_a, i));
            int val_b = static_cast<int>(pop.gene(parent_b, i));

            pop.gene(child_start, i) = static_cast<double>(val_b);
            pop.gene(child_start + 1, i) = static_cast<double>(val_a);

            repl1[val_b] = val_a;
            repl2[val_a] = val_b;
        }

        // STEP 4: Repair positions outside the segment
        for (int i = 0; i < dim; ++i) {
            if (i >= cp1 && i <= cp2)
                continue;

            // Repair child 0 (from parent_a, resolving with repl1)
            int n1 = static_cast<int>(pop.gene(parent_a, i));
            int m1 = repl1[n1];
            while (m1 != -1) {
                n1 = m1;
                m1 = repl1[n1];
            }
            pop.gene(child_start, i) = static_cast<double>(n1);

            // Repair child 1 (from parent_b, resolving with repl2)
            int n2 = static_cast<int>(pop.gene(parent_b, i));
            int m2 = repl2[n2];
            while (m2 != -1) {
                n2 = m2;
                m2 = repl2[n2];
            }
            pop.gene(child_start + 1, i) = static_cast<double>(n2);
        }

        pop.set_evaluated(child_start, false);
        pop.set_evaluated(child_start + 1, false);
    }
};

} // namespace ea
