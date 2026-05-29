#pragma once
/// @file migration_policy.hpp
/// @brief Migration policy for island model — controls when and how many individuals migrate.

#include <ea/core/population.hpp>
#include <cstdint>
#include <string_view>

namespace ea {

/// Migrant selection strategy.
enum class MigrantSelection : uint8_t {
    Best,      ///< Select best individuals by non-dominated rank and crowding
    Random,    ///< Random selection
    Tournament ///< Binary tournament selection
};

/// Migration policy — controls when and how individuals migrate between islands.
///
/// Usage:
///   MigrationPolicy policy;
///   policy.migration_rate = 10;      // Every 10 generations
///   policy.num_migrants = 5;           // Send 5 individuals per migration
///   policy.migrant_selection = MigrantSelection::Best;
struct MigrationPolicy {
    int migration_rate = 10;           ///< Migrate every N generations
    int num_migrants = 5;              ///< Number of individuals sent per migration
    MigrantSelection migrant_selection = MigrantSelection::Best; ///< How to select migrants

    static constexpr std::string_view name() { return "MigrationPolicy"; }

    /// Check if migration should occur at the given generation.
    /// @param generation Current generation number (0-indexed)
    /// @return true if migration should occur
    [[nodiscard]] bool should_migrate(int generation) const {
        return migration_rate > 0 && generation > 0 && (generation % migration_rate == 0);
    }

    /// Get the number of migrants to send, clamped to population size.
    /// @param pop_size Population<> size of the source island
    /// @return Clamped number of migrants
    [[nodiscard]] int migrants_to_send(int pop_size) const {
        if (num_migrants <= 0) return 0;
        return num_migrants > pop_size ? pop_size : num_migrants;
    }
};

} // namespace ea
