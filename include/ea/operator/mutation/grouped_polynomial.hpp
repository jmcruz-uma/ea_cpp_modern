#pragma once
/// @file grouped_polynomial.hpp
/// @brief Grouped Polynomial Mutation — applies polynomial mutation to a group of variables.
/// Reference: jMetal GroupedPolynomialMutation.java
///        https://doi.org/10.1109/SSCI.2016.7850214
///
/// Variables are divided into groups. A random group is selected and all variables
/// in that group are mutated using polynomial mutation. The grouping strategy
/// is configurable via a callable that maps group index -> vector of variable indices.

#include <cmath>
#include <algorithm>
#include <vector>
#include <functional>
#include <ea/core/population.hpp>
#include <ea/core/encoding.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// Grouped Polynomial Mutation for real-valued encodings.
/// Mutates all variables within a randomly selected group using polynomial mutation.
struct GroupedPolynomialMutation {
    double distribution_index = 20.0;                      ///< Distribution index (η)
    int num_groups = 1;                                       ///< Number of variable groups
    std::function<std::vector<int>(int)> group_fn;            ///< Maps group index -> variable indices

    static constexpr Encoding encoding() { return Encoding::Real; }

    void apply(this GroupedPolynomialMutation& self, Population& pop, int idx) {
        auto& rng = Random::instance();

        // Default grouping: contiguous blocks
        int num_groups = std::max(1, self.num_groups);
        int group_idx = rng.uniform_int(0, num_groups - 1);

        std::vector<int> variable_indices;
        if (self.group_fn) {
            variable_indices = self.group_fn(group_idx);
        } else {
            // Default: contiguous blocks
            int base = pop.dim / num_groups;
            int rem = pop.dim % num_groups;
            int start = 0;
            for (int g = 0; g < group_idx; ++g) {
                start += base + (g < rem ? 1 : 0);
            }
            int count = base + (group_idx < rem ? 1 : 0);
            variable_indices.reserve(count);
            for (int i = 0; i < count; ++i) {
                variable_indices.push_back(start + i);
            }
        }

        double mut_pow = 1.0 / (self.distribution_index + 1.0);

        for (int var_idx : variable_indices) {
            if (var_idx < 0 || var_idx >= pop.dim) continue;

            double y = pop.gene(idx, var_idx);
            double lb = pop.lower_bounds[var_idx];
            double ub = pop.upper_bounds[var_idx];

            if (lb >= ub) continue;

            double delta1 = (y - lb) / (ub - lb);
            double delta2 = (ub - y) / (ub - lb);
            double rnd = rng.uniform();
            double deltaq;

            if (rnd <= 0.5) {
                double xy = 1.0 - delta1;
                double val = 2.0 * rnd + (1.0 - 2.0 * rnd) * std::pow(xy, self.distribution_index + 1.0);
                deltaq = std::pow(val, mut_pow) - 1.0;
            } else {
                double xy = 1.0 - delta2;
                double val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * std::pow(xy, self.distribution_index + 1.0);
                deltaq = 1.0 - std::pow(val, mut_pow);
            }

            y += deltaq * (ub - lb);
            pop.gene(idx, var_idx) = std::clamp(y, lb, ub);
        }

        pop.set_evaluated(idx, false);
    }
};

} // namespace ea
