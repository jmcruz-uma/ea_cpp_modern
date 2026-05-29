#pragma once
/// @file crowding_replacement.hpp
/// @brief Crowding-based Replacement — general environmental selection with crowding.
///
/// AGE-MOEA style crowding replacement: uses non-dominated sort + crowding distance
/// to select individuals from combined parent+offspring population.
///
/// This is similar to NSGA-II replacement but as a standalone operator for use
/// in algorithms that need crowding-based environmental selection.

#include <algorithm>
#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <limits>
#include <vector>

namespace ea {

/// Crowding-based environmental selection (replacement operator).
/// Combines non-dominated sort with crowding distance for survival selection.
struct CrowdingReplacement {
    /// Perform crowding-based environmental selection.
    /// @param combined          Combined population (parents + offspring)
    /// @param offspring_indices  Indices of offspring (unused, for API compat)
    /// @param target_size       Target population size
    /// @return Indices of selected individuals
    std::vector<int> replace(this CrowdingReplacement&, Population<>& combined,
                           const std::vector<int>& offspring_indices, int target_size) {
        (void)offspring_indices;
        const int total = combined.pop_size;

        // Step 1: Fast non-dominated sort
        auto fronts = fast_non_dominated_sort(combined);

        // Step 2: Compute crowding distance for each front
        std::vector<double> crowding(total, 0.0);
        for (auto& front : fronts) {
            std::vector<double> front_cd;
            compute_crowding_distance(combined, front, front_cd);
            for (size_t i = 0; i < front.size(); ++i) {
                crowding[front[i]] = front_cd[i];
            }
        }

        // Step 3: Select individuals
        std::vector<int> selected;
        selected.reserve(target_size);

        for (auto& front : fronts) {
            if (static_cast<int>(selected.size() + front.size()) <= target_size) {
                // Full front fits
                selected.insert(selected.end(), front.begin(), front.end());
            } else {
                // Partial front — select by crowding distance (descending)
                int remaining = target_size - static_cast<int>(selected.size());

                std::vector<std::pair<double, int>> cd_index(front.size());
                for (size_t i = 0; i < front.size(); ++i) {
                    cd_index[i] = {crowding[front[i]], front[i]};
                }
                std::sort(cd_index.begin(), cd_index.end(),
                          [](const auto& a, const auto& b) { return a.first > b.first; });

                for (int i = 0; i < remaining; ++i) {
                    selected.push_back(cd_index[i].second);
                }
                break;
            }

            if (static_cast<int>(selected.size()) >= target_size)
                break;
        }

        return selected;
    }

    /// Compact selected individuals into the first target_size positions.
    void compact(this CrowdingReplacement& self, Population<>& combined, int target_size) {
        auto selected = self.replace(combined, {}, target_size);

        for (int i = 0; i < target_size; ++i) {
            if (selected[i] != i) {
                combined.copy_individual(selected[i], i);
            }
        }
    }
};

} // namespace ea
