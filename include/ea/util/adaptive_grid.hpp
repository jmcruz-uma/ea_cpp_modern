#pragma once
/// @file adaptive_grid.hpp
/// @brief Adaptive grid archive for PAES (Pareto Archived Evolution Strategy).
/// Reference: Knowles, J.D. & Corne, D.W. "Approximating the nondominated front
/// using the Pareto Archived Evolution Strategy", Evolutionary Computation, 2000.
///
/// The adaptive grid divides objective space into a fixed number of hypercubes
/// (determined by bisections per dimension). When the archive is full, the most
/// crowded hypercube's solution is removed.

#include <algorithm>
#include <cmath>
#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <limits>
#include <vector>

namespace ea {

/// Grid-based archive for PAES.
/// Maintains a bounded set of non-dominated solutions using an adaptive grid.
/// Template param: max_size — maximum number of archived solutions.
/// For SoA, we store indices into an external population.
struct AdaptiveGridArchive {
    int max_size = 100; ///< Maximum archive size
    int bisections = 5; ///< Number of bisections per objective dimension
    int n_obj = 2;      ///< Number of objectives

    // --- Grid state ---
    /// Stored individuals: each is a vector of objectives + individual info
    struct ArchivedIndividual {
        std::vector<double> objectives; ///< [n_obj] objective values
        std::vector<double> genes;      ///< genes for reproduction (copied on add)
        int grid_location = -1;         ///< Which hypercube this individual occupies
        int original_index = -1;        ///< Original population index (if applicable)
    };

    std::vector<ArchivedIndividual> archive_; ///< Archived non-dominated solutions
    std::vector<double> grid_min_;            ///< [n_obj] minimum objective values
    std::vector<double> grid_max_;            ///< [n_obj] maximum objective values
    std::vector<int> grid_occupancy_;         ///< Number of individuals per grid cell
    std::vector<int> grid_bisection_counts_;  ///< Number of bisections per objective
    int total_grid_cells_ = 0;

    // --- Constructor ---
    AdaptiveGridArchive() = default;

    AdaptiveGridArchive(int max_size_, int bisections_, int n_obj_)
        : max_size(max_size_)
        , bisections(bisections_)
        , n_obj(n_obj_)
        , grid_min_(n_obj_, std::numeric_limits<double>::infinity())
        , grid_max_(n_obj_, -std::numeric_limits<double>::infinity())
        , grid_bisection_counts_(n_obj_, bisections_) {
        total_grid_cells_ = static_cast<int>(std::pow(2, bisections * n_obj));
        grid_occupancy_.assign(total_grid_cells_, 0);
    }

    /// Add an individual to the archive.
    /// Returns true if the individual was added (was non-dominated and improved diversity).
    bool add(this auto& self, const std::vector<double>& genes,
             const std::vector<double>& objectives) {
        // Check dominance against all archived individuals
        bool dominated = false;
        std::vector<size_t> dominates_indices;

        for (size_t i = 0; i < self.archive_.size(); ++i) {
            auto dom = self.compare_dominance(objectives, self.archive_[i].objectives);
            if (dom == Dominance::Dominated) {
                dominated = true;
                break;
            } else if (dom == Dominance::Dominates) {
                dominates_indices.push_back(i);
            }
        }

        if (dominated)
            return false;

        // Remove all dominated individuals
        // Remove in reverse order to preserve indices
        std::sort(dominates_indices.begin(), dominates_indices.end(), std::greater<size_t>());
        for (size_t idx : dominates_indices) {
            self.remove_individual(static_cast<int>(idx));
        }

        // If archive not full, add directly
        if (static_cast<int>(self.archive_.size()) < self.max_size) {
            self.insert_individual(genes, objectives);
            return true;
        }

        // Archive is full — check if we can replace someone in a crowded cell
        int loc = self.find_grid_location(objectives);
        int max_occupancy = 0;
        for (int occ : self.grid_occupancy_) {
            max_occupancy = std::max(max_occupancy, occ);
        }

        // If this location is NOT the most crowded, we can potentially add
        if (self.grid_occupancy_[loc] < max_occupancy) {
            // Remove one from the most crowded cell
            self.remove_from_most_crowded_cell();
            self.insert_individual(genes, objectives);
            return true;
        }

        // Archive full and this is in a crowded cell — refuse
        return false;
    }

    /// Add from a population using SoA access.
    bool add_from_population(this auto& self, const Population& pop, int idx) {
        std::vector<double> objs(pop.n_obj);
        for (int o = 0; o < pop.n_obj; ++o) {
            objs[o] = pop.objective(idx, o);
        }
        std::vector<double> genes(pop.dim);
        for (int j = 0; j < pop.dim; ++j) {
            genes[j] = pop.gene(idx, j);
        }
        return self.add(genes, objs);
    }

    /// Get the size of the archive.
    int size(this auto& self) { return static_cast<int>(self.archive_.size()); }

    /// Check if archive is empty.
    bool empty(this auto& self) { return self.archive_.empty(); }

    /// Get objectives of archived individual i.
    const std::vector<double>& objectives(this auto& self, int i) {
        return self.archive_[i].objectives;
    }

    /// Get genes of archived individual i.
    const std::vector<double>& genes(this auto& self, int i) { return self.archive_[i].genes; }

    /// Get a random individual from the archive (for selection).
    int random_index(this auto& self) {
        auto& rng = Random::instance();
        return rng.uniform_int(0, static_cast<int>(self.archive_.size()) - 1);
    }

    /// Copy genes from archive individual into population position.
    void copy_to_population(this auto& self, int archive_idx, Population& pop, int pop_idx) {
        const auto& ind = self.archive_[archive_idx];
        for (int j = 0; j < pop.dim; ++j) {
            pop.gene(pop_idx, j) = ind.genes[j];
        }
        for (int o = 0; o < pop.n_obj; ++o) {
            pop.objective(pop_idx, o) = ind.objectives[o];
        }
        pop.set_evaluated(pop_idx, true);
    }

    /// Find grid location for external objective vector (after compute_density_estimator).
    int find_grid_location_for(this auto& self, const std::vector<double>& objectives) {
        return self.find_grid_location(objectives);
    }

    /// Compute density estimator (update grid).
    void compute_density_estimator(this auto& self) {
        if (self.archive_.empty())
            return;

        // Update grid bounds
        for (int o = 0; o < self.n_obj; ++o) {
            self.grid_min_[o] = std::numeric_limits<double>::infinity();
            self.grid_max_[o] = -std::numeric_limits<double>::infinity();
        }

        for (const auto& ind : self.archive_) {
            for (int o = 0; o < self.n_obj; ++o) {
                self.grid_min_[o] = std::min(self.grid_min_[o], ind.objectives[o]);
                self.grid_max_[o] = std::max(self.grid_max_[o], ind.objectives[o]);
            }
        }

        // Add small epsilon to avoid zero range
        for (int o = 0; o < self.n_obj; ++o) {
            double range = self.grid_max_[o] - self.grid_min_[o];
            if (range < 1e-12) {
                self.grid_max_[o] += 1e-6;
                self.grid_min_[o] -= 1e-6;
            }
        }

        // Reset occupancy
        std::fill(self.grid_occupancy_.begin(), self.grid_occupancy_.end(), 0);

        // Recompute locations
        for (auto& ind : self.archive_) {
            ind.grid_location = self.find_grid_location(ind.objectives);
            self.grid_occupancy_[ind.grid_location]++;
        }
    }

private:
    /// Compare two objective vectors for dominance.
    Dominance compare_dominance(this auto&, const std::vector<double>& a,
                                const std::vector<double>& b) {
        bool a_dominates_b = false;
        bool b_dominates_a = false;
        const int n = static_cast<int>(a.size());
        for (int o = 0; o < n; ++o) {
            if (a[o] < b[o])
                a_dominates_b = true;
            else if (b[o] < a[o])
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

    /// Find grid location (cell index) for a set of objectives.
    int find_grid_location(this auto& self, const std::vector<double>& objectives) {
        if (self.archive_.empty())
            return 0;
        int location = 0;
        int multiplier = 1;
        const int cells_per_dim = static_cast<int>(std::pow(2, self.bisections));

        for (int o = 0; o < self.n_obj; ++o) {
            double range = self.grid_max_[o] - self.grid_min_[o];
            if (range < 1e-14)
                range = 1e-14;
            double normalized = (objectives[o] - self.grid_min_[o]) / range;
            int grid_coord = static_cast<int>(normalized * cells_per_dim);
            grid_coord = std::clamp(grid_coord, 0, cells_per_dim - 1);
            location += grid_coord * multiplier;
            multiplier *= cells_per_dim;
        }
        return location % self.total_grid_cells_;
    }

    /// Insert a new individual into the archive.
    void insert_individual(this auto& self, const std::vector<double>& genes,
                           const std::vector<double>& objectives) {
        ArchivedIndividual ind;
        ind.objectives = objectives;
        ind.genes = genes;
        ind.grid_location = self.find_grid_location(objectives);
        self.archive_.push_back(std::move(ind));
        if (self.grid_location_valid(ind.grid_location)) {
            self.grid_occupancy_[ind.grid_location]++;
        }
    }

    /// Remove an individual at the given archive index.
    void remove_individual(this auto& self, int idx) {
        if (idx < 0 || idx >= static_cast<int>(self.archive_.size()))
            return;
        int loc = self.archive_[idx].grid_location;
        if (self.grid_location_valid(loc) && self.grid_occupancy_[loc] > 0) {
            self.grid_occupancy_[loc]--;
        }
        self.archive_.erase(self.archive_.begin() + idx);
    }

    /// Remove one individual from the most crowded cell.
    void remove_from_most_crowded_cell(this auto& self) {
        if (self.archive_.empty())
            return;

        int max_occupancy = -1;
        int target_loc = -1;
        for (int i = 0; i < self.total_grid_cells_; ++i) {
            if (self.grid_occupancy_[i] > max_occupancy) {
                max_occupancy = self.grid_occupancy_[i];
                target_loc = i;
            }
        }

        if (target_loc < 0)
            return;

        // Find and remove one individual in that cell
        for (size_t i = 0; i < self.archive_.size(); ++i) {
            if (self.archive_[i].grid_location == target_loc) {
                self.remove_individual(static_cast<int>(i));
                return;
            }
        }
    }

    bool grid_location_valid(this auto& self, int loc) {
        return loc >= 0 && loc < self.total_grid_cells_;
    }
};

} // namespace ea
