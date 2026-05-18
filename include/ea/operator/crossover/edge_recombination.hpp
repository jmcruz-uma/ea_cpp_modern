#pragma once
/// @file edge_recombination.hpp
/// @brief Edge Recombination Crossover (ERX) for permutation encodings.
/// Reference: jMetal EdgeRecombinationCrossover — https://github.com/jMetal/jMetal
/// C++23 with deducing this, SoA Population access.

#include <algorithm>
#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <limits>
#include <vector>

namespace ea {

/// Edge Recombination Crossover (ERX) for permutation-based encodings.
/// Constructs offspring by preserving adjacency (edge) information from parents.
/// Builds an edge map and constructs offspring by traversing it.
/// Genes are integers stored as doubles in the Population array.
struct EdgeRecombinationCrossover {
    double crossover_probability = 1.0; ///< Probability of applying crossover

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Permutation; }

    /// Apply Edge Recombination crossover. Produces 2 children starting at child_start.
    void apply(this EdgeRecombinationCrossover& self, Population& pop, int parent_a, int parent_b,
               int child_start) {
        auto& rng = Random::instance();
        int dim = pop.dim;

        // Copy parents to children (fallback if no crossover)
        for (int j = 0; j < dim; ++j) {
            pop.gene(child_start, j) = pop.gene(parent_a, j);
            pop.gene(child_start + 1, j) = pop.gene(parent_b, j);
        }

        if (dim < 2 || !rng.coin_flip(self.crossover_probability)) {
            pop.set_evaluated(child_start, false);
            pop.set_evaluated(child_start + 1, false);
            return;
        }

        // Build child 1: edge recombination starting from parent_a
        self.build_child(pop, parent_a, parent_b, child_start, rng);
        // Build child 2: edge recombination starting from parent_b
        self.build_child(pop, parent_b, parent_a, child_start + 1, rng);

        pop.set_evaluated(child_start, false);
        pop.set_evaluated(child_start + 1, false);
    }

private:
    /// Build a child using edge recombination from two parents.
    void build_child(this EdgeRecombinationCrossover&, Population& pop, int parent_from,
                     int parent_to, int child_idx, Random& rng) {
        int dim = pop.dim;

        // Build edge map: for each value, store adjacent values from both parents
        // Using fixed-size arrays for small-medium permutations; vector for generality
        std::vector<std::vector<int>> edge_map(dim);

        // Add edges from parent_from (circular adjacency)
        for (int i = 0; i < dim; ++i) {
            int val = static_cast<int>(pop.gene(parent_from, i));
            int next = static_cast<int>(pop.gene(parent_from, (i + 1) % dim));
            int prev = static_cast<int>(pop.gene(parent_from, (i - 1 + dim) % dim));
            edge_map[val].push_back(next);
            edge_map[val].push_back(prev);
        }

        // Add edges from parent_to (circular adjacency)
        for (int i = 0; i < dim; ++i) {
            int val = static_cast<int>(pop.gene(parent_to, i));
            int next = static_cast<int>(pop.gene(parent_to, (i + 1) % dim));
            int prev = static_cast<int>(pop.gene(parent_to, (i - 1 + dim) % dim));
            edge_map[val].push_back(next);
            edge_map[val].push_back(prev);
        }

        // Remove duplicate edges (keep unique neighbors)
        for (auto& neighbors : edge_map) {
            std::sort(neighbors.begin(), neighbors.end());
            neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
        }

        // Build child
        std::vector<int> child(dim);
        std::vector<bool> visited(dim, false);

        int current = static_cast<int>(pop.gene(parent_from, 0));
        int pos = 0;

        while (pos < dim) {
            child[pos] = current;
            visited[current] = true;
            ++pos;

            // Remove current from all adjacency lists
            for (auto& neighbors : edge_map) {
                auto it = std::find(neighbors.begin(), neighbors.end(), current);
                if (it != neighbors.end()) {
                    neighbors.erase(it);
                }
            }

            if (pos >= dim)
                break;

            // Find next candidate
            const auto& neighbors = edge_map[current];
            int next_val = -1;

            if (!neighbors.empty()) {
                // Select neighbor with fewest remaining neighbors
                int min_size = std::numeric_limits<int>::max();
                std::vector<int> candidates;
                for (int n : neighbors) {
                    if (visited[n])
                        continue;
                    int n_size = static_cast<int>(edge_map[n].size());
                    if (n_size < min_size) {
                        min_size = n_size;
                        candidates.clear();
                        candidates.push_back(n);
                    } else if (n_size == min_size) {
                        candidates.push_back(n);
                    }
                }

                if (!candidates.empty()) {
                    next_val =
                        candidates[rng.uniform_int(0, static_cast<int>(candidates.size()) - 1)];
                }
            }

            if (next_val == -1) {
                // No valid neighbors, pick any unvisited value
                for (int v = 0; v < dim; ++v) {
                    if (!visited[v]) {
                        next_val = v;
                        break;
                    }
                }
            }

            current = next_val;
        }

        // Write child to population
        for (int j = 0; j < dim; ++j) {
            pop.gene(child_idx, j) = static_cast<double>(child[j]);
        }
    }
};

} // namespace ea
