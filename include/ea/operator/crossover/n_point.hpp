#pragma once
/// @file n_point.hpp
/// @brief N-Point Crossover — generalization of single and two-point crossover.
/// Reference: jMetal NPointCrossover.java
///
/// Selects N distinct crossover points, sorts them, and alternates segments
/// between parents to create offspring.

#include <algorithm>
#include <vector>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// N-Point Crossover for binary and general encodings.
/// Selects `num_points` distinct crossover points and alternates segments.
struct NPointCrossover {
    double crossover_probability = 0.9; ///< Probability of applying crossover
    int num_points = 1;                ///< Number of crossover points (N >= 1)

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Binary; }
    // Also works for Real, Integer at runtime

    void apply(this NPointCrossover& self, Population& pop,
               int parent_a, int parent_b, int child_start) {
        auto& rng = Random::instance();

        if (!rng.coin_flip(self.crossover_probability)) {
            for (int j = 0; j < pop.dim; ++j) {
                pop.gene(child_start, j) = pop.gene(parent_a, j);
                pop.gene(child_start + 1, j) = pop.gene(parent_b, j);
            }
            pop.set_evaluated(child_start, false);
            pop.set_evaluated(child_start + 1, false);
            return;
        }

        int n = std::max(1, std::min(self.num_points, pop.dim - 1));

        // Generate n distinct crossover points
        std::vector<int> points;
        points.reserve(n);
        std::vector<bool> used(pop.dim, false);

        for (int i = 0; i < n; ++i) {
            int p;
            do {
                p = rng.uniform_int(0, pop.dim - 2); // points in [0, dim-2]
            } while (used[p]);
            used[p] = true;
            points.push_back(p);
        }

        std::sort(points.begin(), points.end());

        // Alternate segments: swap=false starts with parent_a for child1
        bool swap = false;
        int pt_idx = 0;

        for (int j = 0; j < pop.dim; ++j) {
            if (swap) {
                pop.gene(child_start, j) = pop.gene(parent_b, j);
                pop.gene(child_start + 1, j) = pop.gene(parent_a, j);
            } else {
                pop.gene(child_start, j) = pop.gene(parent_a, j);
                pop.gene(child_start + 1, j) = pop.gene(parent_b, j);
            }

            // Check if we've hit a crossover point
            if (pt_idx < static_cast<int>(points.size()) && j == points[pt_idx]) {
                swap = !swap;
                ++pt_idx;
            }
        }

        pop.set_evaluated(child_start, false);
        pop.set_evaluated(child_start + 1, false);
    }
};

} // namespace ea
