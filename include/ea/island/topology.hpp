#pragma once
/// @file topology.hpp
/// @brief Island connectivity topologies — define how islands are connected for migration.

#include <cmath>
#include <cstdint>
#include <string_view>
#include <vector>

namespace ea {

/// Ring topology — each island connected to its immediate neighbors.
/// Island i is connected to (i+1) % n and (i-1+n) % n.
struct RingTopology {
    int num_islands = 0;

    explicit RingTopology(int n = 0) : num_islands(n) {}

    static constexpr std::string_view name() { return "Ring"; }

    /// Get neighbors of island i in the ring.
    /// @param island_id Island index (0 to num_islands-1)
    /// @return Vector of neighbor island indices
    [[nodiscard]] std::vector<int> neighbors(int island_id) const {
        if (num_islands <= 1) return {};
        int prev = (island_id - 1 + num_islands) % num_islands;
        int next = (island_id + 1) % num_islands;
        return {prev, next};
    }
};

/// Mesh topology — 2D grid of islands, connected to North/South/East/West neighbors.
/// Islands are arranged in a rows x cols grid (row-major order).
///
/// Island index: i = row * cols + col
/// Neighbors: N=(row-1,col), S=(row+1,col), E=(row,col+1), W=(row,col-1)
/// With toroidal wrapping (periodic boundary conditions).
struct MeshTopology {
    int rows = 0;        ///< Number of rows in the grid
    int cols = 0;        ///< Number of columns in the grid
    int num_islands = 0; ///< Total islands (rows * cols)

    MeshTopology() = default;

    /// Create a mesh topology with given rows and columns.
    /// @param r Number of rows
    /// @param c Number of columns
    MeshTopology(int r, int c) : rows(r), cols(c), num_islands(r * c) {}

    static constexpr std::string_view name() { return "Mesh"; }

    /// Get neighbors of island i in the 2D mesh (with wrapping).
    /// @param island_id Island index (0 to num_islands-1)
    /// @return Vector of neighbor island indices (up to 4)
    [[nodiscard]] std::vector<int> neighbors(int island_id) const {
        if (num_islands <= 1 || cols == 0) return {};

        int row = island_id / cols;
        int col = island_id % cols;

        int north_row = (row - 1 + rows) % rows;
        int south_row = (row + 1) % rows;
        int east_col  = (col + 1) % cols;
        int west_col  = (col - 1 + cols) % cols;

        return {
            north_row * cols + col,  // North
            south_row * cols + col,  // South
            row * cols + east_col,   // East
            row * cols + west_col    // West
        };
    }
};

/// Fully connected topology — every island is connected to every other island.
/// Migration broadcasts to all islands (or sends to all).
struct FullyConnectedTopology {
    int num_islands = 0;

    explicit FullyConnectedTopology(int n = 0) : num_islands(n) {}

    static constexpr std::string_view name() { return "FullyConnected"; }

    /// Get neighbors of island i — all other islands.
    /// @param island_id Island index (0 to num_islands-1)
    /// @return Vector of all other island indices
    [[nodiscard]] std::vector<int> neighbors(int island_id) const {
        std::vector<int> result;
        result.reserve(num_islands - 1);
        for (int i = 0; i < num_islands; ++i) {
            if (i != island_id) result.push_back(i);
        }
        return result;
    }
};

} // namespace ea
