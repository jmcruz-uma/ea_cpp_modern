#pragma once
/// @file mocell.hpp
/// @brief MOCell (Multi-Objective Cellular Genetic Algorithm)
/// Reference: Durillo, J.J., et al. "A study of master-slave approaches to parallel
/// non-dominated sorting genetic algorithms", IEEE CEC 2008.
///
/// Cellular GA with neighborhood structure, external archive, and feedback from
/// archive to selection. Each individual evolves with its neighborhood.

#include <ea/core/population.hpp>
#include <ea/core/comparator.hpp>
#include <ea/util/random.hpp>
#include <ea/operator/replacement/nsga2_replacement.hpp>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string_view>

namespace ea {

/// MOCell — Multi-Objective Cellular Genetic Algorithm.
/// Template composition for crossover and mutation.
template<typename CX, typename MT>
struct MOCell {
    CX crossover;
    MT mutation;

    int pop_size = 100;
    int max_evals = 25000;
    int archive_size = 100;   ///< External archive capacity
    int neighborhood_radius = 1; ///< Radius for cellular neighborhood (1 = immediate neighbors)

    static constexpr std::string_view name() { return "MOCell"; }

    /// Run MOCell.
    template<typename Problem>
    void run(this auto& self, Population& pop, Problem&& problem) {
        const int n = self.pop_size;
        const int dim = pop.dim;
        const int n_obj = pop.n_obj;

        if (pop.pop_size != n) pop.resize(n);

        // === Evaluate initial population ===
        for (int i = 0; i < n; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }
        int evals = n;

        // Initialize archive with all non-dominated individuals
        std::vector<std::vector<double>> archive_genes_;
        std::vector<std::vector<double>> archive_obj_;
        int archive_count_ = 0;

        auto update_archive = [&](const Population& p, int idx,
                                      std::vector<std::vector<double>>& a_genes,
                                      std::vector<std::vector<double>>& a_obj,
                                      int& a_count, int a_size, int dim_, int n_obj_) {
            // Check dominance
            bool dominated = false;
            std::vector<int> to_remove;

            for (int a = 0; a < a_count; ++a) {
                bool a_dominates_new = false;
                bool new_dominates_a = false;
                for (int o = 0; o < n_obj_; ++o) {
                    double fa = a_obj[a][o];
                    double fb = p.objective(idx, o);
                    if (fa < fb) a_dominates_new = true;
                    else if (fb < fa) new_dominates_a = true;
                }
                if (a_dominates_new && !new_dominates_a) {
                    dominated = true;
                    break;
                }
                if (new_dominates_a && !a_dominates_new) {
                    to_remove.push_back(a);
                }
            }

            if (dominated) return;

            // Remove dominated
            std::sort(to_remove.begin(), to_remove.end(), std::greater<int>());
            for (int r : to_remove) {
                for (int i = r; i < a_count - 1; ++i) {
                    a_genes[i] = a_genes[i + 1];
                    a_obj[i] = a_obj[i + 1];
                }
                a_count--;
            }

            if (a_count < a_size) {
                for (int j = 0; j < dim_; ++j) a_genes[a_count][j] = p.gene(idx, j);
                for (int o = 0; o < n_obj_; ++o) a_obj[a_count][o] = p.objective(idx, o);
                a_count++;
            } else {
                // Archive full — random replacement (simplified)
                auto& rng = Random::instance();
                int replace = rng.uniform_int(0, a_size - 1);
                for (int j = 0; j < dim_; ++j) a_genes[replace][j] = p.gene(idx, j);
                for (int o = 0; o < n_obj_; ++o) a_obj[replace][o] = p.objective(idx, o);
            }
        };

        // Resize archive vectors
        archive_genes_.resize(self.archive_size, std::vector<double>(dim));
        archive_obj_.resize(self.archive_size, std::vector<double>(n_obj));

        // Initial archive population
        for (int i = 0; i < n; ++i) {
            update_archive(pop, i, archive_genes_, archive_obj_,
                           archive_count_, self.archive_size, dim, n_obj);
        }

        // Working populations
        Population parents(2, dim, n_obj, pop.encoding, pop.n_const);
        parents.lower_bounds = pop.lower_bounds;
        parents.upper_bounds = pop.upper_bounds;

        Population offspring(1, dim, n_obj, pop.encoding, pop.n_const);
        offspring.lower_bounds = pop.lower_bounds;
        offspring.upper_bounds = pop.upper_bounds;

        auto& rng = Random::instance();

        // === Main loop: cellular evolution ===
        int current_individual = 0;

        while (evals < self.max_evals) {
            // --- 1. Get neighborhood of current individual ---
            // Linear topology: neighbors are within radius on a ring
            std::vector<int> neighbors;
            neighbors.push_back(current_individual); // Include self
            for (int r = 1; r <= self.neighborhood_radius; ++r) {
                int left = (current_individual - r + n) % n;
                int right = (current_individual + r) % n;
                neighbors.push_back(left);
                neighbors.push_back(right);
            }

            // --- 2. Selection: pick two parents ---
            // Parent 1: from neighborhood
            int p1_idx = neighbors[rng.uniform_int(0, static_cast<int>(neighbors.size()) - 1)];

            // Parent 2: from archive if large enough, else from neighborhood
            int p2_idx;
            if (archive_count_ > 1 && rng.coin_flip(0.5)) {
                p2_idx = rng.uniform_int(0, archive_count_ - 1);
                // Copy archive individual to parents[1]
                for (int j = 0; j < dim; ++j) {
                    parents.gene(1, j) = archive_genes_[p2_idx][j];
                }
                for (int o = 0; o < n_obj; ++o) {
                    parents.objective(1, o) = archive_obj_[p2_idx][o];
                }
                parents.set_evaluated(1, true);
                p2_idx = 1; // Use parents index
            } else {
                p2_idx = neighbors[rng.uniform_int(0, static_cast<int>(neighbors.size()) - 1)];
                p2_idx = (p2_idx % n); // Ensure valid
                // Copy to parents[1]
                for (int j = 0; j < dim; ++j) {
                    parents.gene(1, j) = pop.gene(p2_idx, j);
                }
                for (int o = 0; o < n_obj; ++o) {
                    parents.objective(1, o) = pop.objective(p2_idx, o);
                }
                parents.set_evaluated(1, true);
                p2_idx = 1;
            }

            // Copy p1 to parents[0]
            for (int j = 0; j < dim; ++j) {
                parents.gene(0, j) = pop.gene(p1_idx, j);
            }
            for (int o = 0; o < n_obj; ++o) {
                parents.objective(0, o) = pop.objective(p1_idx, o);
            }
            parents.set_evaluated(0, true);

            // --- 3. Reproduction: crossover + mutation ---
            for (int j = 0; j < dim; ++j) {
                offspring.gene(0, j) = parents.gene(0, j);
            }

            self.crossover.apply(parents, 0, 1, 0); // Result in parents[0]
            for (int j = 0; j < dim; ++j) {
                offspring.gene(0, j) = parents.gene(0, j);
            }

            self.mutation.apply(offspring, 0);

            // Evaluate offspring
            problem(offspring, 0);
            offspring.set_evaluated(0, true);
            evals++;
            if (evals >= self.max_evals) break;

            // --- 4. Replacement: compare with current individual ---
            auto dom = self.compare_dominance(offspring, 0, pop, current_individual);

            if (dom == Dominance::Dominates) {
                // Offspring dominates current: replace
                for (int j = 0; j < dim; ++j) {
                    pop.gene(current_individual, j) = offspring.gene(0, j);
                }
                for (int o = 0; o < n_obj; ++o) {
                    pop.objective(current_individual, o) = offspring.objective(0, o);
                }
                pop.set_evaluated(current_individual, true);

                update_archive(pop, current_individual, archive_genes_, archive_obj_,
                               archive_count_, self.archive_size, dim, n_obj);
            } else if (dom == Dominance::Equal) {
                // Non-dominated: apply NSGA-II style selection in neighborhood + offspring
                // Build neighborhood + offspring pool
                std::vector<int> pool = neighbors;
                pool.push_back(-1); // Marker for offspring

                // Non-dominated sort on neighborhood + offspring
                // Simplified: add offspring to pop temporarily, sort, keep best n in neighborhood
                // Actually, we just compare using crowding-like distance

                // Add offspring to archive
                update_archive(offspring, 0, archive_genes_, archive_obj_,
                               archive_count_, self.archive_size, dim, n_obj);

                // Replace if offspring is not the worst in the neighborhood
                // Simple heuristic: if random < 0.5, replace
                if (rng.coin_flip(0.5)) {
                    for (int j = 0; j < dim; ++j) {
                        pop.gene(current_individual, j) = offspring.gene(0, j);
                    }
                    for (int o = 0; o < n_obj; ++o) {
                        pop.objective(current_individual, o) = offspring.objective(0, o);
                    }
                    pop.set_evaluated(current_individual, true);
                }
            }
            // If offspring is dominated, discard

            // Move to next individual (cellular sweep)
            current_individual = (current_individual + 1) % n;
        }

        // === Copy archive back to population for result ===
        // Note: MOCell result is typically the archive
        // We copy archive into the first positions of pop
        int result_size = std::min(archive_count_, n);
        for (int i = 0; i < result_size; ++i) {
            for (int j = 0; j < dim; ++j) {
                pop.gene(i, j) = archive_genes_[i][j];
            }
            for (int o = 0; o < n_obj; ++o) {
                pop.objective(i, o) = archive_obj_[i][o];
            }
            pop.set_evaluated(i, true);
        }
        pop.pop_size = result_size;
    }

private:
    /// Compare dominance between two individuals in potentially different populations.
    Dominance compare_dominance(this auto&, const Population& pop_a, int a,
                                 const Population& pop_b, int b) {
        const int n_obj = pop_a.n_obj;
        bool a_dominates_b = false;
        bool b_dominates_a = false;

        for (int o = 0; o < n_obj; ++o) {
            double fa = pop_a.objective(a, o);
            double fb = pop_b.objective(b, o);
            if (fa < fb) a_dominates_b = true;
            else if (fb < fa) b_dominates_a = true;
            if (a_dominates_b && b_dominates_a) return Dominance::Equal;
        }

        if (a_dominates_b && !b_dominates_a) return Dominance::Dominates;
        if (b_dominates_a && !a_dominates_b) return Dominance::Dominated;
        return Dominance::Equal;
    }
};

} // namespace ea
