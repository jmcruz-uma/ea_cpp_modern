#pragma once
/// @file island_model.hpp
/// @brief Island Model — orchestrates multiple islands for parallel distributed evolution.
///
/// The Island Model distributes a population across multiple islands, each
/// evolving independently. Periodically, individuals migrate between islands
/// according to a topology and migration policy. This promotes diversity and
/// can improve convergence on complex multi-modal landscapes.
///
/// Usage:
///   NSGAII<SBXCrossover, PolynomialMutation> nsga;
///   RingTopology topo(4);
///   MigrationPolicy policy{10, 5, MigrantSelection::Best};
///   IslandModel model(nsga, topo, policy, 4, 50, 30, 2);
///   model.run(problem, 25000);
///   auto front = model.best_front();
///
/// Architecture:
///   - IslandModel owns N Island instances
///   - Each Island has its own sub-population and algorithm instance
///   - Synchronous migration: all islands evolve, then exchange migrants
///   - Global Pareto front is merged from all island non-dominated fronts

#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <ea/island/island.hpp>
#include <ea/island/migration.hpp>
#include <ea/island/migration_policy.hpp>
#include <ea/island/topology.hpp>
#include <ea/util/random.hpp>
#include <algorithm>
#include <span>
#include <vector>

namespace ea {

/// Island Model — high-level orchestrator for island-based evolution.
/// @tparam Algorithm Type of EA algorithm used on each island
template <typename Algorithm>
struct IslandModel {
    Algorithm base_algorithm;              /// Template algorithm (copied per island)
    MigrationPolicy migration_policy;      /// Global migration policy
    MigrationOperator migrator;           /// Migration operator

    int num_islands = 0;                   /// Number of islands
    int pop_per_island = 0;                /// Population<> per island
    int dim = 0;                           /// Number of decision variables
    int n_obj = 0;                         /// Number of objectives

    std::vector<Island<Algorithm>> islands; /// Island instances

    IslandModel() = default;

    /// Construct an island model with RingTopology.
    /// @param algo Base algorithm template (copied to each island)
    /// @param topo Ring topology defining island connectivity
    /// @param policy Migration policy
    /// @param num_islands Number of islands
    /// @param pop_per_island Population<> size per island
    /// @param dim Number of decision variables
    /// @param n_obj Number of objectives
    IslandModel(Algorithm algo, RingTopology topo, MigrationPolicy policy,
                int num_islands, int pop_per_island, int dim, int n_obj)
        : base_algorithm(std::move(algo))
        , migration_policy(policy)
        , num_islands(num_islands)
        , pop_per_island(pop_per_island)
        , dim(dim)
        , n_obj(n_obj) {
        islands.reserve(num_islands);
        for (int i = 0; i < num_islands; ++i) {
            Island<Algorithm> island(base_algorithm, migration_policy, pop_per_island, dim, n_obj);
            islands.push_back(std::move(island));
        }
    }

    /// Construct an island model with MeshTopology.
    /// @param algo Base algorithm template (copied to each island)
    /// @param topo Mesh topology defining island connectivity
    /// @param policy Migration policy
    /// @param num_islands Number of islands
    /// @param pop_per_island Population<> size per island
    /// @param dim Number of decision variables
    /// @param n_obj Number of objectives
    IslandModel(Algorithm algo, MeshTopology topo, MigrationPolicy policy,
                int num_islands, int pop_per_island, int dim, int n_obj)
        : base_algorithm(std::move(algo))
        , migration_policy(policy)
        , num_islands(num_islands)
        , pop_per_island(pop_per_island)
        , dim(dim)
        , n_obj(n_obj) {
        islands.reserve(num_islands);
        for (int i = 0; i < num_islands; ++i) {
            Island<Algorithm> island(base_algorithm, migration_policy, pop_per_island, dim, n_obj);
            islands.push_back(std::move(island));
        }
    }

    /// Construct an island model with FullyConnectedTopology.
    /// @param algo Base algorithm template (copied to each island)
    /// @param topo Fully connected topology
    /// @param policy Migration policy
    /// @param num_islands Number of islands
    /// @param pop_per_island Population<> size per island
    /// @param dim Number of decision variables
    /// @param n_obj Number of objectives
    IslandModel(Algorithm algo, FullyConnectedTopology topo, MigrationPolicy policy,
                int num_islands, int pop_per_island, int dim, int n_obj)
        : base_algorithm(std::move(algo))
        , migration_policy(policy)
        , num_islands(num_islands)
        , pop_per_island(pop_per_island)
        , dim(dim)
        , n_obj(n_obj) {
        islands.reserve(num_islands);
        for (int i = 0; i < num_islands; ++i) {
            Island<Algorithm> island(base_algorithm, migration_policy, pop_per_island, dim, n_obj);
            islands.push_back(std::move(island));
        }
    }

    /// Run the island model evolution.
    /// @param problem Evaluation function: void(Population<>&, int)
    /// @param max_evals Maximum total evaluations across all islands
    template <typename Problem>
    void run(this auto& self, Problem&& problem, int max_evals) {
        if (self.num_islands == 0 || self.pop_per_island == 0) return;

        // Determine bounds from problem if available
        std::vector<double> lb(self.dim, 0.0);
        std::vector<double> ub(self.dim, 1.0);

        // Try to get bounds from problem if it has the interface
        if constexpr (requires(Problem p) { p.lower_bounds(); }) {
            auto lb_span = problem.lower_bounds();
            auto ub_span = problem.upper_bounds();
            std::copy_n(lb_span.begin(), self.dim, lb.begin());
            std::copy_n(ub_span.begin(), self.dim, ub.begin());
        }

        // Initialize all islands
        for (int i = 0; i < self.num_islands; ++i) {
            self.islands[i].initialize(std::span<const double>(lb.data(), lb.size()),
                                        std::span<const double>(ub.data(), ub.size()),
                                        problem);
        }

        int total_evals = self.num_islands * self.pop_per_island;
        int generation = 0;

        // Main loop: evolve islands, then migrate
        while (total_evals < max_evals) {
            // 1. Evolve each island for migration_rate generations
            int epoch_generations = self.migration_policy.migration_rate;
            if (epoch_generations <= 0) epoch_generations = 10; // Default

            // Clamp to not exceed max_evals
            int epoch_evals = self.num_islands * self.pop_per_island * epoch_generations;
            if (total_evals + epoch_evals > max_evals) {
                int remaining_evals = max_evals - total_evals;
                epoch_generations = remaining_evals / (self.num_islands * self.pop_per_island);
                if (epoch_generations <= 0) epoch_generations = 1;
            }

            for (int i = 0; i < self.num_islands; ++i) {
                self.islands[i].run_epoch(problem, epoch_generations);
            }
            total_evals += self.num_islands * self.pop_per_island * epoch_generations;
            generation += epoch_generations;

            // 2. Migration phase
            if (self.migration_policy.should_migrate(generation)) {
                self.perform_migration(problem);
            }
        }
    }

    /// Perform one migration step: each island sends migrants to neighbors.
    /// @param problem Problem reference (for re-evaluation if needed)
    template <typename Problem>
    void perform_migration(this auto& self, Problem&& problem) {
        (void)problem; // May be used for re-evaluation in future extensions

        // Collect migrants from each island
        std::vector<Population<>> all_migrants;
        all_migrants.reserve(self.num_islands);
        for (int i = 0; i < self.num_islands; ++i) {
            all_migrants.push_back(self.islands[i].migrate_out());
        }

        // Distribute migrants to neighbors using ring topology
        for (int i = 0; i < self.num_islands; ++i) {
            auto neighbors = self.get_neighbors(i);
            for (int neighbor : neighbors) {
                if (neighbor >= 0 && neighbor < self.num_islands && all_migrants[neighbor].pop_size > 0) {
                    self.islands[i].migrate_in(all_migrants[neighbor]);
                }
            }
        }
    }

    /// Get the global best Pareto front across all islands.
    /// Merges all non-dominated individuals and returns the global non-dominated set.
    /// @return Population<> containing the global Pareto front
    [[nodiscard]] Population<> best_front(this auto& self) {
        // Count total non-dominated individuals
        int total_nd = 0;
        std::vector<std::vector<int>> all_fronts;
        all_fronts.reserve(self.num_islands);

        for (int i = 0; i < self.num_islands; ++i) {
            auto front = self.islands[i].non_dominated_front();
            all_fronts.push_back(front);
            total_nd += static_cast<int>(front.size());
        }

        if (total_nd == 0) {
            return Population<>(0, self.dim, self.n_obj);
        }

        // Merge all non-dominated individuals into a single population
        Population<> merged(total_nd, self.dim, self.n_obj);
        merged.lower_bounds = self.islands[0].population.lower_bounds;
        merged.upper_bounds = self.islands[0].population.upper_bounds;

        int idx = 0;
        for (int i = 0; i < self.num_islands; ++i) {
            const auto& front = all_fronts[i];
            const auto& pop = self.islands[i].population;
            for (int individual_idx : front) {
                std::copy_n(pop.genes_ptr(individual_idx), self.dim, merged.genes_ptr(idx));
                std::copy_n(pop.objectives_ptr(individual_idx), self.n_obj, merged.objectives_ptr(idx));
                merged.flags[idx] = pop.flags[individual_idx];
                ++idx;
            }
        }

        // Non-dominated sort the merged population and return only the first front
        Population<> merged_copy = merged;
        auto global_fronts = fast_non_dominated_sort(merged_copy);
        if (global_fronts.empty()) {
            return merged;
        }

        int nd_size = static_cast<int>(global_fronts[0].size());
        Population<> result(nd_size, self.dim, self.n_obj);
        result.lower_bounds = merged.lower_bounds;
        result.upper_bounds = merged.upper_bounds;

        for (int i = 0; i < nd_size; ++i) {
            int src = global_fronts[0][i];
            std::copy_n(merged.genes_ptr(src), self.dim, result.genes_ptr(i));
            std::copy_n(merged.objectives_ptr(src), self.n_obj, result.objectives_ptr(i));
            result.flags[i] = merged.flags[src];
        }

        return result;
    }

private:
    /// Get neighbors for a given island using ring topology.
    /// @param island_id Island index
    /// @return Vector of neighbor indices
    [[nodiscard]] std::vector<int> get_neighbors(int island_id) const {
        if (num_islands <= 1) return {};
        int prev = (island_id - 1 + num_islands) % num_islands;
        int next = (island_id + 1) % num_islands;
        return {prev, next};
    }
};

} // namespace ea
