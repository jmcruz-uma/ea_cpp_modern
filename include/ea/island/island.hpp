#pragma once
/// @file island.hpp
/// @brief Island — wraps a sub-population, algorithm, and migration state.
///
/// An Island is a self-contained evolutionary unit that evolves its own
/// sub-population for a number of generations (an "epoch"), then participates
/// in migration with neighboring islands.
///
/// Usage:
///   Island island(nsga, policy, 50, 30, 2);
///   island.initialize(problem);  // Random init + evaluate
///   island.run_epoch(problem, 10); // Evolve 10 generations
///   auto migrants = island.migrate_out(); // Export migrants
///   island.migrate_in(other_migrants);    // Import migrants

#include <ea/core/population.hpp>
#include <ea/island/migration.hpp>
#include <ea/island/migration_policy.hpp>
#include <ea/util/random.hpp>
#include <algorithm>
#include <span>

namespace ea {

/// Island — evolutionary unit with local population and algorithm.
/// @tparam Algorithm Type of the EA algorithm (e.g., NSGAII<CX, MT>)
template <typename Algorithm>
struct Island {
    Algorithm algorithm;                    ///< Local EA algorithm instance
    MigrationPolicy migration_policy;       ///< Migration policy for this island
    MigrationOperator migrator;            ///< Migration operator (selection + integration)

    Population population;                  ///< Local sub-population
    int pop_size = 0;                       ///< Population size for this island
    int dim = 0;                            ///< Number of decision variables
    int n_obj = 0;                          ///< Number of objectives

    // --- Stats ---
    int evaluations = 0;                    ///< Total evaluations performed on this island
    int generation = 0;                     ///< Current generation counter

    Island() = default;

    /// Construct an island with given algorithm and parameters.
    /// @param algo Algorithm instance (copied)
    /// @param policy Migration policy
    /// @param island_pop_size Population size for this island
    /// @param dim Number of decision variables
    /// @param n_obj Number of objectives
    Island(Algorithm algo, MigrationPolicy policy, int island_pop_size, int dim, int n_obj)
        : algorithm(std::move(algo))
        , migration_policy(policy)
        , pop_size(island_pop_size)
        , dim(dim)
        , n_obj(n_obj) {}

    /// Initialize the island population with random individuals within bounds.
    /// @param lb Lower bounds (span of size dim)
    /// @param ub Upper bounds (span of size dim)
    /// @param problem Evaluation function: void(Population&, int)
    template <typename Problem>
    void initialize(std::span<const double> lb, std::span<const double> ub, Problem&& problem) {
        population = Population(pop_size, dim, n_obj);
        population.lower_bounds.assign(lb.begin(), lb.end());
        population.upper_bounds.assign(ub.begin(), ub.end());

        auto& rng = Random::instance();
        for (int i = 0; i < pop_size; ++i) {
            for (int j = 0; j < dim; ++j) {
                population.gene(i, j) = rng.uniform(lb[j], ub[j]);
            }
        }

        // Evaluate initial population
        for (int i = 0; i < pop_size; ++i) {
            problem(population, i);
            population.set_evaluated(i, true);
        }
        evaluations = pop_size;
        generation = 0;
    }

    /// Run local evolution for a given number of generations.
    /// Uses the algorithm's run() method with adjusted max_evals.
    /// @param problem Evaluation function
    /// @param generations Number of generations to evolve
    template <typename Problem>
    void run_epoch(Problem&& problem, int generations) {
        if (generations <= 0) return;
        if (population.pop_size == 0) return;

        // Set max_evals to allow exactly `generations` more evaluations
        // Each generation evaluates pop_size individuals
        int epoch_evals = generations * pop_size;
        int prev_max_evals = algorithm.max_evals;
        algorithm.max_evals = evaluations + epoch_evals;
        int prev_pop_size = algorithm.pop_size;
        algorithm.pop_size = pop_size;

        // Run the algorithm
        algorithm.run(population, problem);
        evaluations = algorithm.max_evals; // Approximate

        // Restore
        algorithm.max_evals = prev_max_evals;
        algorithm.pop_size = prev_pop_size;

        generation += generations;
    }

    /// Select and export migrants from this island.
    /// @return Population containing the selected migrants
    [[nodiscard]] Population migrate_out() const {
        int nm = migration_policy.migrants_to_send(population.pop_size);
        return migrator.select_migrants(population, nm, migration_policy.migrant_selection);
    }

    /// Import and integrate incoming migrants.
    /// @param migrants Incoming migrants to integrate
    void migrate_in(const Population& migrants) {
        migrator.integrate_migrants(population, migrants);
    }

    /// Get the current non-dominated front of this island's population.
    /// @return Vector of indices of non-dominated individuals
    [[nodiscard]] std::vector<int> non_dominated_front() const {
        Population pop_copy = population;
        auto fronts = fast_non_dominated_sort(pop_copy);
        if (fronts.empty()) return {};
        return fronts[0];
    }
};

} // namespace ea
