#pragma once
/// @file mu_plus_lambda.hpp
/// @brief (μ+λ) Replacement — elitist replacement merging parents and offspring.
///
/// Parents and offspring compete together. The best μ individuals survive.
/// This is the elitist strategy: good parents are never lost.
///
/// Reference: jMetal jmetal-core/src/main/java/org/uma/jmetal/operator/selection/impl/RankingAndCrowdingSelection.java

#include <ea/core/population.hpp>
#include <ea/core/comparator.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>

namespace ea {

/// (μ+λ) Replacement — keep best μ from merged parents + offspring.
/// Elitist: parents survive if they dominate offspring.
struct MuPlusLambdaReplacement {
    /// Replace population with best μ individuals from parents ∪ offspring.
    /// @param pop Population (after offspring have been appended)
    /// @param offspring_indices Indices of offspring in pop (typically [μ, μ+λ))
    /// @param mu Target population size (μ)
    /// @return Indices of the selected μ individuals
    std::vector<int> replace(Population& pop,
                             const std::vector<int>& offspring_indices,
                             int mu) {
        // Merge parent indices [0, mu) with offspring_indices
        std::vector<int> merged;
        merged.reserve(mu + offspring_indices.size());

        for (int i = 0; i < mu; ++i) {
            merged.push_back(i); // parents
        }
        for (int idx : offspring_indices) {
            merged.push_back(idx); // offspring
        }

        // Compute ranks and crowding distances for merged set
        auto fronts = fast_non_dominated_sort_subset(pop, merged);

        // Build selected set from fronts until we reach μ
        std::vector<int> selected;
        selected.reserve(mu);

        for (const auto& front : fronts) {
            if (static_cast<int>(selected.size()) + static_cast<int>(front.size()) <= mu) {
                // Take entire front
                for (int idx : front) {
                    selected.push_back(idx);
                }
            } else {
                // Need to select from this front by crowding distance
                std::vector<double> cdist;
                compute_crowding_distance_subset(pop, front, cdist);

                // Sort front by descending crowding distance
                std::vector<int> sorted_front = front;
                std::sort(sorted_front.begin(), sorted_front.end(),
                    [&cdist, &front](int a, int b) {
                        auto pos_a = std::find(front.begin(), front.end(), a) - front.begin();
                        auto pos_b = std::find(front.begin(), front.end(), b) - front.begin();
                        return cdist[pos_a] > cdist[pos_b];
                    });

                int needed = mu - static_cast<int>(selected.size());
                for (int i = 0; i < needed; ++i) {
                    selected.push_back(sorted_front[i]);
                }
                break;
            }
        }

        // Compact selected individuals into positions [0, mu)
        std::vector<int> new_order = selected;
        for (int idx : merged) {
            if (std::find(selected.begin(), selected.end(), idx) == selected.end()) {
                new_order.push_back(idx);
            }
        }

        // Reorder population according to new_order
        std::vector<double> temp_genes(pop.pop_size * pop.dim);
        std::vector<double> temp_objectives(pop.pop_size * pop.n_obj);
        std::vector<IndFlags> temp_flags(pop.pop_size);

        for (int new_pos = 0; new_pos < static_cast<int>(new_order.size()); ++new_pos) {
            int old_pos = new_order[new_pos];
            std::copy_n(pop.genes_ptr(old_pos), pop.dim,
                       temp_genes.data() + new_pos * pop.dim);
            std::copy_n(pop.objectives_ptr(old_pos), pop.n_obj,
                       temp_objectives.data() + new_pos * pop.n_obj);
            temp_flags[new_pos] = pop.flags[old_pos];
        }

        // Copy back
        std::copy(temp_genes.begin(), temp_genes.end(), pop.genes.begin());
        std::copy(temp_objectives.begin(), temp_objectives.end(), pop.objectives.begin());
        std::copy(temp_flags.begin(), temp_flags.end(), pop.flags.begin());

        // Resize pop back to mu
        pop.resize(mu);

        // Return indices of survivors
        std::vector<int> result(mu);
        std::iota(result.begin(), result.end(), 0);
        return result;
    }

private:
    /// Fast non-dominated sort for a subset of indices.
    std::vector<std::vector<int>> fast_non_dominated_sort_subset(
        Population& pop, const std::vector<int>& indices) {
        const int n = static_cast<int>(indices.size());
        std::vector<int> domination_count(n, 0);
        std::vector<std::vector<int>> dominated_by(n);
        std::vector<std::vector<int>> fronts;
        fronts.push_back({});

        for (int p = 0; p < n; ++p) {
            for (int q = 0; q < n; ++q) {
                if (p == q) continue;
                auto dom = compare_dominance(pop, indices[p], indices[q]);
                if (dom == Dominance::Dominates) {
                    dominated_by[p].push_back(q);
                } else if (dom == Dominance::Dominated) {
                    domination_count[p]++;
                }
            }
            if (domination_count[p] == 0) {
                fronts[0].push_back(indices[p]);
            }
        }

        int k = 0;
        while (!fronts[k].empty()) {
            std::vector<int> next_front;
            for (int p_idx : fronts[k]) {
                int p_local = static_cast<int>(std::find(indices.begin(), indices.end(), p_idx) - indices.begin());
                for (int q_local : dominated_by[p_local]) {
                    domination_count[q_local]--;
                    if (domination_count[q_local] == 0) {
                        next_front.push_back(indices[q_local]);
                    }
                }
            }
            k++;
            fronts.push_back(std::move(next_front));
        }

        if (!fronts.empty() && fronts.back().empty()) {
            fronts.pop_back();
        }

        return fronts;
    }

    /// Compute crowding distances for a subset of indices.
    void compute_crowding_distance_subset(const Population& pop,
                                          const std::vector<int>& indices,
                                          std::vector<double>& crowding_dist) {
        const int n = static_cast<int>(indices.size());
        crowding_dist.assign(n, 0.0);

        if (n <= 2) {
            for (int i = 0; i < n; ++i) {
                crowding_dist[i] = std::numeric_limits<double>::infinity();
            }
            return;
        }

        for (int o = 0; o < pop.n_obj; ++o) {
            std::vector<int> sorted_indices(indices);
            std::sort(sorted_indices.begin(), sorted_indices.end(),
                [&pop, o](int a, int b) {
                    return pop.objective(a, o) < pop.objective(b, o);
                });

            double f_min = pop.objective(sorted_indices[0], o);
            double f_max = pop.objective(sorted_indices[n - 1], o);
            double range = f_max - f_min;

            auto find_pos = [&indices](int idx) -> int {
                return static_cast<int>(std::find(indices.begin(), indices.end(), idx) - indices.begin());
            };

            int pos_min = find_pos(sorted_indices[0]);
            int pos_max = find_pos(sorted_indices[n - 1]);
            crowding_dist[pos_min] = std::numeric_limits<double>::infinity();
            crowding_dist[pos_max] = std::numeric_limits<double>::infinity();

            if (range < 1e-12) continue;

            for (int i = 1; i < n - 1; ++i) {
                double f_prev = pop.objective(sorted_indices[i - 1], o);
                double f_next = pop.objective(sorted_indices[i + 1], o);
                int pos = find_pos(sorted_indices[i]);
                crowding_dist[pos] += (f_next - f_prev) / range;
            }
        }
    }
};

} // namespace ea
