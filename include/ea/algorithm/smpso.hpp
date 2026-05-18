#pragma once
/// @file smpso.hpp
/// @brief SMPSO (Speed-constrained Multi-objective PSO)
/// Reference: Nebro, A.J., et al. "SMPSO: A new PSO-based metaheuristic for
/// multi-objective optimization", IEEE MCDM 2009.
///
/// Particle swarm optimization with velocity updates, personal best tracking,
/// leader selection via crowding distance from an external archive, and
/// polynomial mutation applied periodically.

#include <algorithm>
#include <cmath>
#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/util/random.hpp>
#include <limits>
#include <string_view>
#include <vector>

namespace ea {

/// SMPSO — Speed-constrained Multi-objective Particle Swarm Optimization.
/// Uses an internal PolynomialMutation for perturbation.
/// Template composition: no external operator templates (mutation is internal).
struct SMPSO {
    // --- PSO parameters ---
    double c1_min = 1.5;            ///< Cognitive coefficient min
    double c1_max = 2.5;            ///< Cognitive coefficient max
    double c2_min = 1.5;            ///< Social coefficient min
    double c2_max = 2.5;            ///< Social coefficient max
    double weight_min = 0.1;        ///< Inertia weight min
    double weight_max = 0.5;        ///< Inertia weight max
    double r1_min = 0.0;            ///< Random factor 1 min
    double r1_max = 1.0;            ///< Random factor 1 max
    double r2_min = 0.0;            ///< Random factor 2 min
    double r2_max = 1.0;            ///< Random factor 2 max
    double change_velocity1 = -1.0; ///< Velocity damping on lower bound hit (<0 = no change)
    double change_velocity2 = -1.0; ///< Velocity damping on upper bound hit (<0 = no change)

    // --- Algorithm parameters ---
    int pop_size = 100;
    int max_evals = 25000;
    int archive_size = 100; ///< External archive size (leaders)

    // --- Internal mutation ---
    PolynomialMutation mutation;

    // --- Internal state ---
    std::vector<std::vector<double>> speed_;         ///< [pop_size][dim] particle velocities
    std::vector<std::vector<double>> personal_best_; ///< [pop_size][dim] personal best genes
    std::vector<std::vector<double>>
        personal_best_obj_;                          ///< [pop_size][n_obj] personal best objectives
    std::vector<std::vector<double>> archive_genes_; ///< [archive_size][dim] archived genes
    std::vector<std::vector<double>> archive_obj_;   ///< [archive_size][n_obj] archived objectives
    int archive_count_ = 0;

    static constexpr std::string_view name() { return "SMPSO"; }

    /// Run SMPSO.
    template <typename Problem> void run(this auto& self, Population& pop, Problem&& problem) {
        const int n = self.pop_size;
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;

        if (pop.pop_size != n)
            pop.resize(n);

        // Initialize velocities and personal bests
        self.speed_.assign(n, std::vector<double>(dim, 0.0));
        self.personal_best_.resize(n, std::vector<double>(dim));
        self.personal_best_obj_.resize(n, std::vector<double>(n_obj));
        self.archive_genes_.resize(self.archive_size, std::vector<double>(dim));
        self.archive_obj_.resize(self.archive_size, std::vector<double>(n_obj));
        self.archive_count_ = 0;

        // Compute delta max/min for velocity constriction
        std::vector<double> delta_max(dim);
        std::vector<double> delta_min(dim);
        for (int j = 0; j < dim; ++j) {
            delta_max[j] = (pop.upper_bounds[j] - pop.lower_bounds[j]) / 2.0;
            delta_min[j] = -delta_max[j];
        }

        auto& rng = Random::instance();

        // === Evaluate initial population ===
        for (int i = 0; i < n; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
            // Set personal best = current
            for (int j = 0; j < dim; ++j) {
                self.personal_best_[i][j] = pop.gene(i, j);
            }
            for (int o = 0; o < n_obj; ++o) {
                self.personal_best_obj_[i][o] = pop.objective(i, o);
            }
        }
        int evals = n;

        // Initialize archive with all individuals
        for (int i = 0; i < n; ++i) {
            self.add_to_archive(pop, i);
        }

        // === Main loop ===
        int iteration = 1;
        while (evals < self.max_evals) {
            // --- Update velocities ---
            for (int i = 0; i < n; ++i) {
                if (self.archive_count_ == 0)
                    break;

                // Select global best from archive (binary tournament by crowding)
                int leader_idx = self.select_leader();

                double r1 = rng.uniform(self.r1_min, self.r1_max);
                double r2 = rng.uniform(self.r2_min, self.r2_max);
                double c1 = rng.uniform(self.c1_min, self.c1_max);
                double c2 = rng.uniform(self.c2_min, self.c2_max);
                double w = self.inertia_weight(iteration);

                double rho = c1 + c2;
                double chi = self.constriction_coefficient(c1, c2);

                for (int j = 0; j < dim; ++j) {
                    double v = self.speed_[i][j];
                    double pb = self.personal_best_[i][j];
                    double x = pop.gene(i, j);
                    double gb = self.archive_genes_[leader_idx][j];

                    v = chi * (w * v + c1 * r1 * (pb - x) + c2 * r2 * (gb - x));

                    // Velocity constriction
                    v = self.velocity_constriction(v, delta_max[j], delta_min[j]);
                    self.speed_[i][j] = v;
                }
            }

            // --- Update positions ---
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < dim; ++j) {
                    double new_x = pop.gene(i, j) + self.speed_[i][j];
                    double lb = pop.lower_bounds[j];
                    double ub = pop.upper_bounds[j];

                    if (new_x < lb) {
                        new_x = lb;
                        if (self.change_velocity1 > 0) {
                            self.speed_[i][j] *= self.change_velocity1;
                        }
                    } else if (new_x > ub) {
                        new_x = ub;
                        if (self.change_velocity2 > 0) {
                            self.speed_[i][j] *= self.change_velocity2;
                        }
                    }

                    pop.gene(i, j) = new_x;
                }
                pop.set_evaluated(i, false);
            }

            // --- Perturbation (mutation) ---
            // Every 6th particle gets mutated
            for (int i = 0; i < n; i += 6) {
                self.mutation.apply(pop, i);
            }

            // --- Evaluate ---
            for (int i = 0; i < n; ++i) {
                if (!pop.evaluated(i)) {
                    problem(pop, i);
                    pop.set_evaluated(i, true);
                    evals++;
                    if (evals >= self.max_evals)
                        break;
                }
            }
            if (evals >= self.max_evals)
                break;

            // --- Update personal bests ---
            for (int i = 0; i < n; ++i) {
                auto dom = self.compare_dominance_obj(pop, i, self.personal_best_obj_[i].data());

                if (dom == Dominance::Dominates || dom == Dominance::Equal) {
                    // Current dominates or equals personal best — update
                    for (int j = 0; j < dim; ++j) {
                        self.personal_best_[i][j] = pop.gene(i, j);
                    }
                    for (int o = 0; o < n_obj; ++o) {
                        self.personal_best_obj_[i][o] = pop.objective(i, o);
                    }
                }
            }

            // --- Update archive (leaders) ---
            for (int i = 0; i < n; ++i) {
                self.add_to_archive(pop, i);
            }

            iteration++;
        }
    }

private:
    /// Compare population individual against objective vector.
    Dominance compare_dominance_obj(this auto& self, const Population& pop, int idx,
                                    const double* other_obj) {
        const int n_obj = pop.n_obj;
        bool a_dominates_b = false;
        bool b_dominates_a = false;

        for (int o = 0; o < n_obj; ++o) {
            double fa = pop.objective(idx, o);
            double fb = other_obj[o];
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

    /// Add individual from population to archive (if non-dominated).
    void add_to_archive(this auto& self, const Population& pop, int idx) {
        const int n_obj = pop.n_obj;
        const int dim = pop.dim;

        // Check dominance against archive
        bool dominated = false;
        std::vector<int> to_remove;

        for (int a = 0; a < self.archive_count_; ++a) {
            bool a_dominates_new = false;
            bool new_dominates_a = false;

            for (int o = 0; o < n_obj; ++o) {
                double fa = self.archive_obj_[a][o];
                double fb = pop.objective(idx, o);
                if (fa < fb)
                    a_dominates_new = true;
                else if (fb < fa)
                    new_dominates_a = true;
            }

            if (a_dominates_new && !new_dominates_a) {
                dominated = true;
                break;
            }
            if (new_dominates_a && !a_dominates_new) {
                to_remove.push_back(a);
            }
        }

        if (dominated)
            return;

        // Remove dominated individuals (in reverse order)
        std::sort(to_remove.begin(), to_remove.end(), std::greater<int>());
        for (int r : to_remove) {
            self.remove_from_archive(r);
        }

        if (self.archive_count_ < self.archive_size) {
            // Add directly
            for (int j = 0; j < dim; ++j) {
                self.archive_genes_[self.archive_count_][j] = pop.gene(idx, j);
            }
            for (int o = 0; o < n_obj; ++o) {
                self.archive_obj_[self.archive_count_][o] = pop.objective(idx, o);
            }
            self.archive_count_++;
        } else {
            // Archive full — add if improves diversity (simpler: random replacement)
            // More sophisticated: crowding-based replacement
            auto& rng = Random::instance();
            int replace_idx = rng.uniform_int(0, self.archive_size - 1);
            for (int j = 0; j < dim; ++j) {
                self.archive_genes_[replace_idx][j] = pop.gene(idx, j);
            }
            for (int o = 0; o < n_obj; ++o) {
                self.archive_obj_[replace_idx][o] = pop.objective(idx, o);
            }
        }
    }

    /// Remove archived individual at index (shift remaining).
    void remove_from_archive(this auto& self, int idx) {
        if (idx < 0 || idx >= self.archive_count_)
            return;
        for (int i = idx; i < self.archive_count_ - 1; ++i) {
            self.archive_genes_[i] = self.archive_genes_[i + 1];
            self.archive_obj_[i] = self.archive_obj_[i + 1];
        }
        self.archive_count_--;
    }

    /// Select leader from archive using binary tournament (crowding distance).
    int select_leader(this auto& self) {
        auto& rng = Random::instance();
        if (self.archive_count_ <= 1)
            return 0;

        int a = rng.uniform_int(0, self.archive_count_ - 1);
        int b = rng.uniform_int(0, self.archive_count_ - 1);

        // Simple comparison: prefer one that is less crowded
        // For now, random selection (crowding distance would need recompute)
        // In jMetal, they use archive.comparator() which is crowding-based
        // Here we do random tournament
        if (a == b)
            return a;

        // Compute simple crowding: compare distance to nearest neighbor in archive
        double dist_a = self.archive_distance(a);
        double dist_b = self.archive_distance(b);

        return (dist_a > dist_b) ? a : b;
    }

    /// Compute approximate crowding distance for archive individual.
    double archive_distance(this auto& self, int idx) {
        if (self.archive_count_ <= 1)
            return std::numeric_limits<double>::infinity();

        const int n_obj = static_cast<int>(self.archive_obj_[0].size());
        double min_dist = std::numeric_limits<double>::infinity();

        for (int i = 0; i < self.archive_count_; ++i) {
            if (i == idx)
                continue;
            double dist = 0.0;
            for (int o = 0; o < n_obj; ++o) {
                double d = self.archive_obj_[idx][o] - self.archive_obj_[i][o];
                dist += d * d;
            }
            min_dist = std::min(min_dist, std::sqrt(dist));
        }
        return min_dist;
    }

    /// Velocity constriction — clamp velocity to bounds.
    double velocity_constriction(this auto&, double v, double dmax, double dmin) {
        if (v > dmax)
            return dmax;
        if (v < dmin)
            return dmin;
        return v;
    }

    /// Constriction coefficient (Clerc & Kennedy).
    double constriction_coefficient(this auto&, double c1, double c2) {
        double rho = c1 + c2;
        if (rho <= 4.0)
            return 1.0;
        return 2.0 / (2.0 - rho - std::sqrt(rho * rho - 4.0 * rho));
    }

    /// Inertia weight (constant in original SMPSO, could be linearly decreasing).
    double inertia_weight(this auto& self, int /*iteration*/) { return self.weight_max; }
};

} // namespace ea
