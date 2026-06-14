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
#include <ea/core/concepts.hpp>
#include <ea/core/population.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/util/random.hpp>
#include <limits>
#include <numeric>
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
    double weight_max = 0.1;        ///< Inertia weight max (jMetal default: constant 0.1)
    double r1_min = 0.0;            ///< Random factor 1 min
    double r1_max = 1.0;            ///< Random factor 1 max
    double r2_min = 0.0;            ///< Random factor 2 min
    double r2_max = 1.0;            ///< Random factor 2 max
    double change_velocity1 = -1.0; ///< Velocity multiplier on lower bound hit (-1 = bounce)
    double change_velocity2 = -1.0; ///< Velocity multiplier on upper bound hit (-1 = bounce)

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
    Population<> archive_;   ///< Non-dominated archive SoA (capacity = archive_size + 1)
    int archive_count_ = 0;

    static constexpr std::string_view name() { return "SMPSO"; }

    /// Run SMPSO.
    template <EvalFunctor F> void run(this auto& self, Population<>& pop, F&& problem) {
        const int n = self.pop_size;
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;

        if (pop.pop_size != n)
            pop.resize(n);

        // Initialize velocities and personal bests
        self.speed_.assign(n, std::vector<double>(dim, 0.0));
        self.personal_best_.resize(n, std::vector<double>(dim));
        self.personal_best_obj_.resize(n, std::vector<double>(n_obj));
        // +1 para el slot temporal que usa add_to_archive cuando el archivo está lleno
        self.archive_ = Population<>(self.archive_size + 1, dim, n_obj, pop.n_const);
        self.archive_.lower_bounds = pop.lower_bounds;
        self.archive_.upper_bounds = pop.upper_bounds;
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
                    double gb = self.archive_.gene(leader_idx, j);

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
                        self.speed_[i][j] *= self.change_velocity1;
                    } else if (new_x > ub) {
                        new_x = ub;
                        self.speed_[i][j] *= self.change_velocity2;
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
    Dominance compare_dominance_obj(this auto& self, const Population<>& pop, int idx,
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
    void add_to_archive(this auto& self, const Population<>& pop, int idx) {
        const int n_obj = pop.n_obj;
        const int dim = pop.dim;

        // Check dominance against archive
        bool dominated = false;
        std::vector<int> to_remove;

        for (int a = 0; a < self.archive_count_; ++a) {
            bool a_dominates_new = false;
            bool new_dominates_a = false;

            for (int o = 0; o < n_obj; ++o) {
                double fa = self.archive_.objective(a, o);
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
            std::copy_n(pop.genes_ptr(idx), dim, self.archive_.genes_ptr(self.archive_count_));
            std::copy_n(pop.objectives_ptr(idx), n_obj, self.archive_.objectives_ptr(self.archive_count_));
            self.archive_count_++;
        } else {
            // Archive full — protocolo jMetal CrowdingDistanceArchive.prune():
            // 1. Escribe el nuevo en el slot +1 (temporal).
            // 2. Recalcula crowding para archive_count_+1 soluciones (incluyendo la nueva).
            // 3. Elimina la de menor crowding distance — puede ser la nueva.
            int slot = self.archive_count_; // slot temporal pre-alojado
            std::copy_n(pop.genes_ptr(idx), dim, self.archive_.genes_ptr(slot));
            std::copy_n(pop.objectives_ptr(idx), n_obj, self.archive_.objectives_ptr(slot));

            auto crowding = self.compute_archive_crowding(self.archive_count_ + 1);
            int worst_idx = static_cast<int>(
                std::min_element(crowding.begin(), crowding.end()) - crowding.begin());

            if (worst_idx < self.archive_count_) {
                // Un miembro existente es el más crowded: sustituirlo por el nuevo
                self.archive_.copy_individual(slot, worst_idx);
            }
            // else worst_idx == slot: el nuevo es el más crowded, se descarta (no se toca archive_count_)
        }
    }

    /// Remove archived individual at index (shift remaining).
    void remove_from_archive(this auto& self, int idx) {
        if (idx < 0 || idx >= self.archive_count_)
            return;
        for (int i = idx; i < self.archive_count_ - 1; ++i) {
            self.archive_.copy_individual(i + 1, i);
        }
        self.archive_count_--;
    }

    /// Select leader from archive using binary tournament on NSGA-II crowding distance.
    /// Matches jMetal BinaryTournamentSelection with CrowdingDistanceComparator.
    int select_leader(this auto& self) {
        auto& rng = Random::instance();
        if (self.archive_count_ <= 1)
            return 0;

        int a = rng.uniform_int(0, self.archive_count_ - 1);
        int b = rng.uniform_int(0, self.archive_count_ - 1);

        if (a == b)
            return a;

        auto crowding = self.compute_archive_crowding();
        return (crowding[a] >= crowding[b]) ? a : b;
    }

    /// NSGA-II crowding distance para las primeras n entradas del archivo.
    /// n=-1 (default) → usa archive_count_. Acepta n=archive_count_+1 para prune.
    /// Matches jMetal CrowdingDistance.computeDensityEstimator().
    std::vector<double> compute_archive_crowding(this auto& self, int n = -1) {
        if (n < 0) n = self.archive_count_;
        const int m = self.archive_.n_obj;
        std::vector<double> crowding(n, 0.0);

        if (n <= 2) {
            for (int i = 0; i < n; ++i)
                crowding[i] = std::numeric_limits<double>::infinity();
            return crowding;
        }

        for (int o = 0; o < m; ++o) {
            std::vector<int> order(n);
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](int a, int b) {
                return self.archive_.objective(a, o) < self.archive_.objective(b, o);
            });

            crowding[order[0]]     = std::numeric_limits<double>::infinity();
            crowding[order[n - 1]] = std::numeric_limits<double>::infinity();

            double obj_min = self.archive_.objective(order[0], o);
            double obj_max = self.archive_.objective(order[n - 1], o);
            double range   = obj_max - obj_min;

            if (range < 1e-14)
                continue;

            for (int i = 1; i < n - 1; ++i) {
                if (crowding[order[i]] != std::numeric_limits<double>::infinity()) {
                    double prev_obj = self.archive_.objective(order[i - 1], o);
                    double next_obj = self.archive_.objective(order[i + 1], o);
                    crowding[order[i]] += (next_obj - prev_obj) / range;
                }
            }
        }

        return crowding;
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

    /// Inertia weight — constant. jMetal SMPSO ignores iter/wmin, always returns weightMax.
    double inertia_weight(this auto& self, int /*iteration*/) { return self.weight_max; }
};

} // namespace ea
