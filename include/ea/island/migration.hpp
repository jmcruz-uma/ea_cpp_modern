#pragma once
/// @file migration.hpp
/// @brief Migration operator — selects migrants and integrates them into populations.

#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <ea/island/migration_policy.hpp>
#include <ea/util/random.hpp>
#include <algorithm>
#include <vector>

namespace ea {

/// Policy for how incoming migrants are integrated.
enum class IntegrationPolicy : uint8_t {
    ReplaceWorst,  ///< Replace worst individuals (by rank/crowding)
    ReplaceRandom  ///< Replace random individuals
};

/// MigrationOperator — selects migrants from a population and integrates incoming migrants.
///
/// Usage:
///   MigrationOperator migrator;
///   auto migrants = migrator.select_migrants(pop, 5, MigrantSelection::Best);
///   migrator.integrate_migrants(pop, migrants, IntegrationPolicy::ReplaceWorst);
struct MigrationOperator {
    IntegrationPolicy integration_policy = IntegrationPolicy::ReplaceWorst;

    /// Select migrants from a population according to the given policy.
    /// @param pop Source population
    /// @param num_migrants Number of individuals to select
    /// @param selection How to select migrants
    /// @return Population containing the selected migrants
    Population select_migrants(const Population& pop, int num_migrants,
                                MigrantSelection selection) const {
        if (num_migrants <= 0 || pop.pop_size == 0) {
            return Population(0, pop.dim, pop.n_obj, pop.encoding, pop.n_const);
        }

        int nm = num_migrants > pop.pop_size ? pop.pop_size : num_migrants;
        Population migrants(nm, pop.dim, pop.n_obj, pop.encoding, pop.n_const);

        std::vector<int> selected_indices;
        selected_indices.reserve(nm);

        switch (selection) {
        case MigrantSelection::Best:
            selected_indices = select_best(pop, nm);
            break;
        case MigrantSelection::Random:
            selected_indices = select_random(pop, nm);
            break;
        case MigrantSelection::Tournament:
            selected_indices = select_tournament(pop, nm);
            break;
        }

        // Copy selected individuals to migrants population
        for (int i = 0; i < nm; ++i) {
            int src = selected_indices[i];
            std::copy_n(pop.genes_ptr(src), pop.dim, migrants.genes_ptr(i));
            std::copy_n(pop.objectives_ptr(src), pop.n_obj, migrants.objectives_ptr(i));
            if (pop.n_const > 0) {
                std::copy_n(pop.constraints.data() + src * pop.n_const, pop.n_const,
                            migrants.constraints.data() + i * pop.n_const);
            }
            migrants.flags[i] = pop.flags[src];
        }

        return migrants;
    }

    /// Integrate incoming migrants into a population.
    /// @param pop Target population (modified in place)
    /// @param migrants Incoming migrants to integrate
    void integrate_migrants(Population& pop, const Population& migrants) const {
        if (migrants.pop_size == 0) return;

        int nm = migrants.pop_size;
        if (nm > pop.pop_size) nm = pop.pop_size;

        std::vector<int> target_indices;
        target_indices.reserve(nm);

        switch (integration_policy) {
        case IntegrationPolicy::ReplaceWorst:
            target_indices = find_worst(pop, nm);
            break;
        case IntegrationPolicy::ReplaceRandom:
            target_indices = find_random(pop, nm);
            break;
        }

        // Replace target individuals with migrants
        for (int i = 0; i < nm; ++i) {
            int dst = target_indices[i];
            std::copy_n(migrants.genes_ptr(i), migrants.dim, pop.genes_ptr(dst));
            std::copy_n(migrants.objectives_ptr(i), migrants.n_obj, pop.objectives_ptr(dst));
            if (pop.n_const > 0 && migrants.n_const > 0) {
                std::copy_n(migrants.constraints.data() + i * migrants.n_const,
                            std::min(pop.n_const, migrants.n_const),
                            pop.constraints.data() + dst * pop.n_const);
            }
            pop.flags[dst] = migrants.flags[i];
        }
    }

private:
    /// Select best individuals using non-dominated sort + crowding.
    std::vector<int> select_best(const Population& pop, int nm) const {
        // Make a mutable copy for fast_non_dominated_sort
        Population pop_copy = pop;
        auto fronts = fast_non_dominated_sort(pop_copy);

        std::vector<int> selected;
        selected.reserve(nm);

        for (const auto& front : fronts) {
            if (static_cast<int>(selected.size() + front.size()) <= nm) {
                selected.insert(selected.end(), front.begin(), front.end());
            } else {
                int remaining = nm - static_cast<int>(selected.size());
                std::vector<double> front_cd;
                compute_crowding_distance(pop_copy, front, front_cd);

                std::vector<std::pair<double, int>> cd_index(front.size());
                for (size_t i = 0; i < front.size(); ++i) {
                    cd_index[i] = {front_cd[i], front[i]};
                }
                std::sort(cd_index.begin(), cd_index.end(),
                          [](const auto& a, const auto& b) { return a.first > b.first; });

                for (int i = 0; i < remaining; ++i) {
                    selected.push_back(cd_index[i].second);
                }
                break;
            }
        }

        return selected;
    }

    /// Random selection.
    std::vector<int> select_random(const Population& pop, int nm) const {
        auto& rng = Random::instance();
        std::vector<int> indices(pop.pop_size);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), rng.engine());
        indices.resize(nm);
        return indices;
    }

    /// Tournament selection on rank + crowding.
    std::vector<int> select_tournament(const Population& pop, int nm) const {
        Population pop_copy = pop;
        auto fronts = fast_non_dominated_sort(pop_copy);

        std::vector<int> ranks(pop.pop_size, 0);
        std::vector<double> crowding(pop.pop_size, 0.0);

        for (int r = 0; r < static_cast<int>(fronts.size()); ++r) {
            std::vector<double> front_cd;
            compute_crowding_distance(pop_copy, fronts[r], front_cd);
            for (size_t i = 0; i < fronts[r].size(); ++i) {
                ranks[fronts[r][i]] = r;
                crowding[fronts[r][i]] = front_cd[i];
            }
        }

        auto& rng = Random::instance();
        std::vector<int> selected;
        selected.reserve(nm);

        for (int i = 0; i < nm; ++i) {
            int a = rng.uniform_int(0, pop.pop_size - 1);
            int b = rng.uniform_int(0, pop.pop_size - 1);
            int winner = a;
            if (ranks[b] < ranks[a] || (ranks[b] == ranks[a] && crowding[b] > crowding[a])) {
                winner = b;
            }
            selected.push_back(winner);
        }

        return selected;
    }

    /// Find worst individuals to replace.
    std::vector<int> find_worst(const Population& pop, int nm) const {
        Population pop_copy = pop;
        auto fronts = fast_non_dominated_sort(pop_copy);

        // Reverse order: worst front first
        std::vector<int> worst;
        worst.reserve(nm);

        for (auto it = fronts.rbegin(); it != fronts.rend() && static_cast<int>(worst.size()) < nm; ++it) {
            const auto& front = *it;
            int remaining = nm - static_cast<int>(worst.size());

            if (static_cast<int>(front.size()) <= remaining) {
                worst.insert(worst.end(), front.begin(), front.end());
            } else {
                // Within the worst front, pick least crowded (replace boundary solutions last)
                std::vector<double> front_cd;
                compute_crowding_distance(pop_copy, front, front_cd);

                std::vector<std::pair<double, int>> cd_index(front.size());
                for (size_t i = 0; i < front.size(); ++i) {
                    cd_index[i] = {front_cd[i], front[i]};
                }
                // Sort by crowding ascending (least crowded = worst)
                std::sort(cd_index.begin(), cd_index.end(),
                          [](const auto& a, const auto& b) { return a.first < b.first; });

                for (int i = 0; i < remaining; ++i) {
                    worst.push_back(cd_index[i].second);
                }
            }
        }

        return worst;
    }

    /// Random individuals to replace.
    std::vector<int> find_random(const Population& pop, int nm) const {
        auto& rng = Random::instance();
        std::vector<int> indices(pop.pop_size);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), rng.engine());
        indices.resize(nm);
        return indices;
    }
};

} // namespace ea
