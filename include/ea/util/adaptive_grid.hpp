#pragma once
/// @file adaptive_grid.hpp
/// @brief Adaptive grid archive for PAES (Pareto Archived Evolution Strategy).
/// Reference: Knowles, J.D. & Corne, D.W. "Approximating the nondominated front
/// using the Pareto Archived Evolution Strategy", Evolutionary Computation, 2000.
///
/// The adaptive grid divides objective space into a fixed number of hypercubes
/// (determined by bisections per dimension). When the archive is full, the most
/// crowded hypercube's solution is removed.
///
/// Performance note — incremental grid update:
///   When add() receives an individual whose objectives lie within the current
///   grid bounds (steady-state), grid_occupancy_ is maintained in O(1) via
///   swap-and-pop removal + direct cell-counter increment/decrement.
///   A full O(N) rebuild (compute_density_estimator) only occurs when bounds
///   must expand to accommodate a new extreme individual, which happens only
///   in the initial fill-up phase (~100 iterations for archive_size=100).
///   After calling add(), density_is_fresh() is always true; callers do NOT
///   need to call compute_density_estimator() separately.

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
struct AdaptiveGridArchive {
    int max_size  = 100; ///< Maximum archive size
    int bisections = 5;  ///< Number of bisections per objective dimension
    int n_obj      = 2;  ///< Number of objectives

    // Exposed for direct density reads from paes.hpp after add()
    std::vector<int> grid_occupancy_; ///< Count of archive members per grid cell

    // Debug counters — how often each path was taken in add()
    mutable int n_incremental_ = 0; ///< times O(1) incremental path was used
    mutable int n_rebuild_     = 0; ///< times full O(N) rebuild was needed

    struct ArchivedIndividual {
        std::vector<double> objectives;
        std::vector<double> genes;
        int grid_location  = -1;
        int original_index = -1;
    };

    std::vector<ArchivedIndividual> archive_;
    std::vector<double> grid_min_; ///< [n_obj] bounds used for current grid
    std::vector<double> grid_max_;
    std::vector<int>    grid_bisection_counts_;
    int  total_grid_cells_ = 0;
    bool density_fresh_    = false; ///< true when grid_occupancy_ is accurate

    // ── Constructor ────────────────────────────────────────────────────────────
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
        archive_.reserve(max_size_ + 1);
    }

    // ── Public API ─────────────────────────────────────────────────────────────

    /// Returns true iff grid_occupancy_ is accurate after the last add().
    /// After any add() call, this is always true.
    bool density_is_fresh() const { return density_fresh_; }

    /// Add an individual to the archive.
    /// Returns true if the individual was added (non-dominated w.r.t. archive).
    ///
    /// After returning, grid_occupancy_ is guaranteed accurate (density_is_fresh()
    /// == true). Callers do NOT need to call compute_density_estimator() afterward.
    ///
    /// Fast path (steady-state): O(1) incremental occupancy update.
    /// Slow path (bounds expansion): O(N) full grid rebuild.
    bool add(this auto& self, const std::vector<double>& genes,
             const std::vector<double>& objectives) {
        // ── 1. Dominance check ──────────────────────────────────────────────
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

        // ── 2. Remove dominated — reverse order preserves indices ───────────
        std::sort(dominates_indices.begin(), dominates_indices.end(), std::greater<size_t>());
        for (size_t idx : dominates_indices)
            self.remove_individual(static_cast<int>(idx));

        // ── 3. Choose update path ────────────────────────────────────────────
        bool within = self.is_within_bounds(objectives);

        // ── 4. Insert new individual ────────────────────────────────────────
        self.insert_individual(genes, objectives);

        if (within) {
            // FAST PATH: bounds stable → grid_occupancy_ is accurate.
            // insert_individual() already incremented the new cell; remove_individual()
            // above already decremented removed cells. Just prune if over capacity.
            if (static_cast<int>(self.archive_.size()) > self.max_size) {
                int worst = self.find_worst_by_occupancy();
                if (worst >= 0)
                    self.remove_individual(worst);
            }
            self.density_fresh_ = true;
            ++self.n_incremental_;
        } else {
            // SLOW PATH: bounds expanded (or archive was empty) → full rebuild.
            // compute_density_estimator() resets all occupancies and grid locations,
            // so the stale location assigned by insert_individual() is corrected.
            self.compute_density_estimator(); // sets density_fresh_ = true
            if (static_cast<int>(self.archive_.size()) > self.max_size) {
                int worst = self.find_worst_by_occupancy();
                if (worst >= 0)
                    self.remove_individual(worst);
                // remove_individual() decrements occupancy → still fresh
            }
            ++self.n_rebuild_;
        }

        return true;
    }

    /// Add from a population using SoA access.
    bool add_from_population(this auto& self, const Population<>& pop, int idx) {
        std::vector<double> objs(pop.n_obj);
        for (int o = 0; o < pop.n_obj; ++o)
            objs[o] = pop.objective(idx, o);
        std::vector<double> g(pop.dim);
        for (int j = 0; j < pop.dim; ++j)
            g[j] = pop.gene(idx, j);
        return self.add(g, objs);
    }

    int  size(this auto& self)  { return static_cast<int>(self.archive_.size()); }
    bool empty(this auto& self) { return self.archive_.empty(); }

    const std::vector<double>& objectives(this auto& self, int i) {
        return self.archive_[i].objectives;
    }
    const std::vector<double>& genes(this auto& self, int i) {
        return self.archive_[i].genes;
    }

    int random_index(this auto& self) {
        auto& rng = Random::instance();
        return rng.uniform_int(0, static_cast<int>(self.archive_.size()) - 1);
    }

    void copy_to_population(this auto& self, int archive_idx, Population<>& pop, int pop_idx) {
        const auto& ind = self.archive_[archive_idx];
        for (int j = 0; j < pop.dim; ++j)
            pop.gene(pop_idx, j) = ind.genes[j];
        for (int o = 0; o < pop.n_obj; ++o)
            pop.objective(pop_idx, o) = ind.objectives[o];
        pop.set_evaluated(pop_idx, true);
    }

    /// Compute grid location for an arbitrary objective vector using current bounds.
    int find_grid_location_for(this auto& self, const std::vector<double>& objectives) {
        return self.find_grid_location(objectives);
    }

    /// Full O(N) grid rebuild: recomputes bounds, all grid locations, and occupancy.
    /// After this call, density_is_fresh() == true.
    void compute_density_estimator(this auto& self) {
        if (self.archive_.empty())
            return;

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
        for (int o = 0; o < self.n_obj; ++o) {
            double range = self.grid_max_[o] - self.grid_min_[o];
            if (range < 1e-12) {
                self.grid_max_[o] += 1e-6;
                self.grid_min_[o] -= 1e-6;
            }
        }

        std::fill(self.grid_occupancy_.begin(), self.grid_occupancy_.end(), 0);
        for (auto& ind : self.archive_) {
            ind.grid_location = self.find_grid_location(ind.objectives);
            self.grid_occupancy_[ind.grid_location]++;
        }
        self.density_fresh_ = true;
    }

private:
    /// True if objectives lie within the current grid bounds (and archive is non-empty).
    bool is_within_bounds(this auto& self, const std::vector<double>& objectives) {
        if (self.archive_.empty())
            return false; // no bounds established yet
        for (int o = 0; o < self.n_obj; ++o) {
            if (objectives[o] < self.grid_min_[o] || objectives[o] > self.grid_max_[o])
                return false;
        }
        return true;
    }

    /// Return index of archive member in the most crowded grid cell.
    int find_worst_by_occupancy(this auto& self) {
        int worst_idx     = -1;
        int worst_density = -1;
        for (int i = 0; i < static_cast<int>(self.archive_.size()); ++i) {
            int loc = self.archive_[i].grid_location;
            int d   = self.grid_location_valid(loc) ? self.grid_occupancy_[loc] : 0;
            if (d > worst_density) {
                worst_density = d;
                worst_idx     = i;
            }
        }
        return worst_idx;
    }

    Dominance compare_dominance(this auto&, const std::vector<double>& a,
                                const std::vector<double>& b) {
        bool a_dom = false;
        bool b_dom = false;
        const int n = static_cast<int>(a.size());
        for (int o = 0; o < n; ++o) {
            if (a[o] < b[o])
                a_dom = true;
            else if (b[o] < a[o])
                b_dom = true;
            if (a_dom && b_dom)
                return Dominance::Equal;
        }
        if (a_dom && !b_dom)
            return Dominance::Dominates;
        if (b_dom && !a_dom)
            return Dominance::Dominated;
        return Dominance::Equal;
    }

    int find_grid_location(this auto& self, const std::vector<double>& objectives) {
        if (self.archive_.empty())
            return 0;
        int location      = 0;
        int multiplier    = 1;
        const int cpd     = static_cast<int>(std::pow(2, self.bisections)); // cells per dim
        for (int o = 0; o < self.n_obj; ++o) {
            double range = self.grid_max_[o] - self.grid_min_[o];
            if (range < 1e-14)
                range = 1e-14;
            double normalized = (objectives[o] - self.grid_min_[o]) / range;
            int coord = static_cast<int>(normalized * cpd);
            coord     = std::clamp(coord, 0, cpd - 1);
            location += coord * multiplier;
            multiplier *= cpd;
        }
        return location % self.total_grid_cells_;
    }

    void insert_individual(this auto& self, const std::vector<double>& genes,
                           const std::vector<double>& objectives) {
        ArchivedIndividual ind;
        ind.objectives    = objectives;
        ind.genes         = genes;
        ind.grid_location = self.find_grid_location(objectives);
        int loc           = ind.grid_location;
        self.archive_.push_back(std::move(ind));
        if (self.grid_location_valid(loc))
            self.grid_occupancy_[loc]++;
    }

    /// Swap-and-pop: O(1) removal (avoids O(N) erase shift).
    /// Decrements grid_occupancy_ for the removed individual's cell.
    void remove_individual(this auto& self, int idx) {
        if (idx < 0 || idx >= static_cast<int>(self.archive_.size()))
            return;
        int loc = self.archive_[idx].grid_location;
        if (self.grid_location_valid(loc) && self.grid_occupancy_[loc] > 0)
            self.grid_occupancy_[loc]--;
        int last = static_cast<int>(self.archive_.size()) - 1;
        if (idx != last)
            self.archive_[idx] = std::move(self.archive_[last]);
        self.archive_.pop_back();
    }

    bool grid_location_valid(this auto& self, int loc) {
        return loc >= 0 && loc < self.total_grid_cells_;
    }
};

} // namespace ea
