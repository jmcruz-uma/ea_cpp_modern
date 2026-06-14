#pragma once
/// @file nsga2_replacement.hpp
/// @brief NSGA-II Environmental Selection — the core of NSGA-II's survival selection.
/// Combines parents and offspring, sorts by non-dominated rank then crowding distance,
/// and selects the best pop_size individuals.

#include <algorithm>
#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <numeric>
#include <vector>

namespace ea {

/// NSGA-II environmental selection (replacement operator).
/// Takes a combined population of 2*pop_size individuals,
/// performs fast non-dominated sort, and selects the best pop_size
/// using rank then crowding distance as tiebreaker.
struct NSGAIIReplacement {

    /// Perform NSGA-II environmental selection.
    /// @param combined          Population<> of size 2*pop_size (parents + offspring)
    /// @param offspring_indices Ignored — NSGA-II ranks all combined solutions jointly.
    /// @param pop_size          Target population size
    /// @return Vector of selected individual indices (sorted by rank, then crowding)
    std::vector<int> replace(this NSGAIIReplacement&, Population<>& combined,
                             const std::vector<int>& offspring_indices, int pop_size) {
        (void)offspring_indices;
        // Step 1: Fast non-dominated sort
        auto fronts = fast_non_dominated_sort(combined);

        // Step 2: Compute crowding distance for each front
        // Store per-individual crowding distance
        int total = combined.pop_size;
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
        selected.reserve(pop_size);

        for (auto& front : fronts) {
            if (static_cast<int>(selected.size() + front.size()) <= pop_size) {
                // Full front fits — add all
                selected.insert(selected.end(), front.begin(), front.end());
            } else {
                // Partial front — select by crowding distance (descending)
                int remaining = pop_size - static_cast<int>(selected.size());

                // Create index-crowding pairs for this front
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

            if (static_cast<int>(selected.size()) >= pop_size)
                break;
        }

        return selected;
    }

    /// Compact selected individuals into the first pop_size positions of the population.
    /// After calling this, indices [0, pop_size) contain the selected individuals.
    void compact(this NSGAIIReplacement& self, Population<>& combined, int pop_size) {
        auto selected = self.replace(combined, {}, pop_size);

        // Copy selected individuals to positions [0, pop_size)
        for (int i = 0; i < pop_size; ++i) {
            if (selected[i] != i) {
                combined.copy_individual(selected[i], i);
            }
        }
    }
};

} // namespace ea