#pragma once
/// @file smpso_replacement.hpp
/// @brief SMPSO Replacement — leader-based archive update using crowding distance.
///
/// Reference: Nebro, A.J., et al. "SMPSO: A new PSO-based metaheuristic for
/// multi-objective optimization", IEEE MCDM 2009.
///
/// SMPSO maintains an external archive of non-dominated solutions (leaders).
/// Replacement updates the archive by:
/// 1. Adding non-dominated offspring
/// 2. Removing archive members dominated by offspring
/// 3. If archive overflows, removing most crowded individuals

#include <algorithm>
#include <cmath>
#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <limits>
#include <vector>

namespace ea {

/// SMPSO archive replacement (leader-based environmental selection).
/// Manages an external archive of non-dominated solutions.
struct SMPSOReplacement {
    int archive_size = 100; ///< Maximum archive capacity

    /// Add individual to archive (replacing dominated or most crowded).
    /// @param archive_pop       Archive population (will be modified)
    /// @param candidate_idx     Index of candidate individual in source population
    /// @param source            Source population containing the candidate
    /// @return true if candidate was added to archive
    bool add_to_archive(this SMPSOReplacement& self, Population& archive_pop,
                        int candidate_idx, const Population& source) {
        const int n_obj = source.n_obj;
        const int dim = source.dim;
        const int current_size = archive_pop.pop_size;

        // Check dominance against current archive members
        std::vector<int> to_remove;
        bool dominated = false;

        for (int a = 0; a < current_size; ++a) {
            auto dom = compare_dominance_source(source, candidate_idx, archive_pop, a);

            if (dom == Dominance::Dominated) {
                // Candidate is dominated by archive member
                dominated = true;
                break;
            } else if (dom == Dominance::Dominates) {
                // Candidate dominates archive member
                to_remove.push_back(a);
            }
        }

        if (dominated)
            return false;

        // Remove dominated archive members (in reverse order)
        std::sort(to_remove.begin(), to_remove.end(), std::greater<int>());
        for (int r : to_remove) {
            // Shift archive members down
            for (int i = r; i < archive_pop.pop_size - 1; ++i) {
                archive_pop.copy_individual(i + 1, i);
            }
            archive_pop.pop_size--;
        }

        if (archive_pop.pop_size < self.archive_size) {
            // Add directly
            int slot = archive_pop.pop_size;
            for (int j = 0; j < dim; ++j) {
                archive_pop.gene(slot, j) = source.gene(candidate_idx, j);
            }
            for (int o = 0; o < n_obj; ++o) {
                archive_pop.objective(slot, o) = source.objective(candidate_idx, o);
            }
            archive_pop.set_evaluated(slot, true);
            archive_pop.pop_size++;
            return true;
        } else if (archive_pop.pop_size > 0) {
            // Archive full — replace most crowded individual if candidate is better
            double min_crowding = std::numeric_limits<double>::max();
            int worst_idx = -1;

            std::vector<int> archive_indices(archive_pop.pop_size);
            std::iota(archive_indices.begin(), archive_indices.end(), 0);
            std::vector<double> cd;
            compute_crowding_distance(archive_pop, archive_indices, cd);

            for (int i = 0; i < archive_pop.pop_size; ++i) {
                if (cd[i] < min_crowding) {
                    min_crowding = cd[i];
                    worst_idx = i;
                }
            }

            if (worst_idx >= 0) {
                // Replace worst with candidate
                for (int j = 0; j < dim; ++j) {
                    archive_pop.gene(worst_idx, j) = source.gene(candidate_idx, j);
                }
                for (int o = 0; o < n_obj; ++o) {
                    archive_pop.objective(worst_idx, o) = source.objective(candidate_idx, o);
                }
                archive_pop.set_evaluated(worst_idx, true);
                return true;
            }
        }

        return false;
    }

    /// Compatibility: replace() signature for algorithm templates.
    /// Returns current archive indices (SMPSO updates in-place).
    std::vector<int> replace(this SMPSOReplacement& self, Population& archive_pop,
                             const std::vector<int>& offspring_indices, int target_size) {
        (void)target_size;
        // Add each offspring to archive
        for (int idx : offspring_indices) {
            // Note: this assumes offspring are in archive_pop itself
            // In practice, SMPSO calls add_to_archive with source=main pop
            self.add_to_archive(archive_pop, idx, archive_pop);
        }

        std::vector<int> selected(archive_pop.pop_size);
        std::iota(selected.begin(), selected.end(), 0);
        return selected;
    }

private:
    /// Compare dominance between individuals from different populations.
    static Dominance compare_dominance_source(const Population& pop_a, int idx_a,
                                               const Population& pop_b, int idx_b) {
        const int n_obj = pop_a.n_obj;
        bool a_dominates_b = false;
        bool b_dominates_a = false;

        for (int o = 0; o < n_obj; ++o) {
            double fa = pop_a.objective(idx_a, o);
            double fb = pop_b.objective(idx_b, o);
            if (fa < fb)
                a_dominates_b = true;
            else if (fb < fa)
                b_dominates_a = true;
            if (a_dominates_b && b_dominates_a)
                return Dominance::Equal;
        }

        if (a_dominates_b && !b_dominates_a)
            return Dominance::Dominates;
        if (b_dominates_a && !a_dominates_b)
            return Dominance::Dominated;
        return Dominance::Equal;
    }
};

} // namespace ea
