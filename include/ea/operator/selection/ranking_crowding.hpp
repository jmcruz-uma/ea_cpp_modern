#pragma once
/// @file ranking_crowding.hpp
/// @brief RankingAndCrowdingSelection — NSGA-II environmental selection.
///        Selects individuals by Pareto rank (lower is better), then by crowding
///        distance (higher is better) within the last front that partially fits.
///
/// Reference: jMetal
///   jmetal-core/src/main/java/org/uma/jmetal/operator/selection/impl/RankingAndCrowdingSelection.java

#include <algorithm>
#include <ea/core/population.hpp>
#include <numeric>
#include <vector>

namespace ea {

/// RankingAndCrowdingSelection — NSGA-II environmental selection.
///
/// Given a population, ranks, and crowding distances, selects `pool_size`
/// individuals such that:
/// 1. All individuals from lower-ranked fronts are taken first.
/// 2. If a front partially fits, its members are sorted by crowding distance
///    (descending) and the top-k are taken.
struct RankingAndCrowdingSelection {
    int pool_size = 0; ///< Number of individuals to select

    /// Select `pool_size` individuals into `mating_pool`.
    /// @param pop          Population<> (used for bounds checking)
    /// @param mating_pool  Output: selected individual indices (resized)
    /// @param ranks        Rank of each individual (0 = best)
    /// @param crowding_dist Crowding distance of each individual
    void select(this RankingAndCrowdingSelection& self, Population<>& pop,
                std::vector<int>& mating_pool, const std::vector<int>& ranks,
                const std::vector<double>& crowding_dist) {
        const int n = pop.pop_size;
        mating_pool.clear();

        if (self.pool_size <= 0)
            return;
        if (self.pool_size > n)
            self.pool_size = n;

        // Group indices by rank
        int max_rank = 0;
        for (int r : ranks) {
            if (r > max_rank)
                max_rank = r;
        }

        std::vector<std::vector<int>> fronts(max_rank + 1);
        for (int i = 0; i < n; ++i) {
            fronts[ranks[i]].push_back(i);
        }

        int selected = 0;
        int rank_idx = 0;

        // Take whole fronts while they fit
        while (rank_idx <= max_rank &&
               selected + static_cast<int>(fronts[rank_idx].size()) <= self.pool_size) {
            mating_pool.insert(mating_pool.end(), fronts[rank_idx].begin(), fronts[rank_idx].end());
            selected += static_cast<int>(fronts[rank_idx].size());
            ++rank_idx;
        }

        // Partial front: sort by crowding distance descending, take top-k
        if (selected < self.pool_size && rank_idx <= max_rank) {
            auto& last_front = fronts[rank_idx];
            std::sort(last_front.begin(), last_front.end(), [&crowding_dist](int a, int b) {
                return crowding_dist[a] > crowding_dist[b];
            });
            int need = self.pool_size - selected;
            mating_pool.insert(mating_pool.end(), last_front.begin(), last_front.begin() + need);
        }
    }
};

} // namespace ea
