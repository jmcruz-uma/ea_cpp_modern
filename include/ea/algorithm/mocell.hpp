#pragma once
/// @file mocell.hpp
/// @brief MOCell — Multi-Objective Cellular Genetic Algorithm.
/// Reference: Durillo, J.J. et al., CEC 2008.
///
/// Faithful port of jMetal 7.4 MOCell:
///   - C9 toroidal 2D mesh neighborhood (8 Moore neighbors + self = 9)
///   - Binary tournament selection (rank + crowding) from neighborhood
///   - Second parent: binary tournament from archive (if size > 1), else neighborhood
///   - Non-dominated replacement: NDS + crowding on (neighbors + offspring), find worst
///   - CrowdingDistanceArchive: maintains non-dominated set; prune by min crowding when full
///   - Per-individual rank/crowding updated from local NDS (stale values used in tournament)

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string_view>
#include <vector>

#include <ea/core/comparator.hpp>
#include <ea/core/concepts.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>

namespace ea {

/// MOCell — Multi-Objective Cellular Genetic Algorithm.
template <Crossover CX, Mutation MT> struct MOCell {
    CX crossover;
    MT mutation;

    int pop_size     = 100;
    int max_evals    = 25000;
    int archive_size = 100;
    int grid_rows    = 10;  // sqrt(pop_size); updated in run() if needed
    int grid_cols    = 10;

    static constexpr std::string_view name() { return "MOCell"; }

    template <EvalFunctor F> void run(this auto& self, Population<>& pop, F&& problem) {
        const int n     = self.pop_size;
        const int dim   = pop.dim;
        const int n_obj = pop.n_obj;
        auto& rng = Random::instance();

        if (pop.pop_size != n) pop.resize(n);

        // Compute grid dimensions from pop_size (default 10×10 for pop=100)
        const int rows = self.grid_rows;
        const int cols = self.grid_cols;

        // ── Evaluate initial population ───────────────────────────────────────
        for (int i = 0; i < n; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }
        int evals = 0;  // counts only offspring evaluations (matching jMetal)

        // ── Per-individual rank/crowding (updated during non-dominated replacement) ──
        // Initial values -1/0.0 match jMetal's unranked/uncomputed solution defaults.
        std::vector<int>    pop_rank(n, -1);
        std::vector<double> pop_crowding(n, 0.0);

        // ── Archive: non-dominated solution set with crowding-distance pruning ──
        // Pre-allocate archive_size+1 slots to temporarily accommodate one extra solution.
        const int ARC_CAP = self.archive_size;
        int arc_count = 0;
        bool arc_dirty = true;
        Population<> arc(ARC_CAP + 1, dim, n_obj, pop.n_const);
        arc.lower_bounds = pop.lower_bounds;
        arc.upper_bounds = pop.upper_bounds;
        std::vector<double> arc_crowding(ARC_CAP + 1, 0.0);

        // Compute crowding distances for arc_count archive members.
        auto compute_archive_crowding = [&]() {
            if (!arc_dirty || arc_count == 0) return;
            if (arc_count <= 2) {
                for (int a = 0; a < arc_count; ++a)
                    arc_crowding[a] = std::numeric_limits<double>::infinity();
                arc_dirty = false;
                return;
            }
            std::fill(arc_crowding.begin(), arc_crowding.begin() + arc_count, 0.0);
            std::vector<int> sorted(arc_count);
            std::iota(sorted.begin(), sorted.end(), 0);
            for (int o = 0; o < n_obj; ++o) {
                std::sort(sorted.begin(), sorted.end(),
                          [&](int a, int b) { return arc.objective(a, o) < arc.objective(b, o); });
                double fmin = arc.objective(sorted[0], o);
                double fmax = arc.objective(sorted[arc_count - 1], o);
                if (fmin == fmax) continue;
                double range = fmax - fmin;
                arc_crowding[sorted[0]] = std::numeric_limits<double>::infinity();
                arc_crowding[sorted[arc_count - 1]] = std::numeric_limits<double>::infinity();
                for (int i = 1; i < arc_count - 1; ++i) {
                    arc_crowding[sorted[i]] +=
                        (arc.objective(sorted[i + 1], o) - arc.objective(sorted[i - 1], o)) / range;
                }
            }
            arc_dirty = false;
        };

        // Add a solution to the non-dominated archive.
        // Rejects if dominated; removes solutions it dominates; prunes by min crowding if full.
        auto archive_add = [&](const std::vector<double>& g_in, const std::vector<double>& o_in) {
            bool dominated = false;
            std::vector<int> to_remove;
            for (int a = 0; a < arc_count; ++a) {
                bool arc_dom = false, new_dom = false;
                for (int o = 0; o < n_obj; ++o) {
                    if (arc.objective(a, o) < o_in[o]) arc_dom = true;
                    else if (o_in[o] < arc.objective(a, o)) new_dom = true;
                }
                if (arc_dom && !new_dom) { dominated = true; break; }
                if (new_dom && !arc_dom) to_remove.push_back(a);
            }
            if (dominated) return;

            // Remove dominated archive members (swap-with-last, descending order)
            std::sort(to_remove.rbegin(), to_remove.rend());
            for (int r : to_remove) {
                if (r != arc_count - 1)
                    arc.copy_individual(arc_count - 1, r);
                arc_count--;
            }

            // Add new solution (slot arc_count, within ARC_CAP+1 pre-allocation)
            std::copy(g_in.begin(), g_in.end(), arc.genes_ptr(arc_count));
            std::copy(o_in.begin(), o_in.end(), arc.objectives_ptr(arc_count));
            arc_count++;

            // Prune if over capacity: remove solution with minimum crowding distance
            if (arc_count > ARC_CAP) {
                arc_dirty = true;
                compute_archive_crowding();
                int worst = 0;
                for (int a = 1; a < arc_count; ++a)
                    if (arc_crowding[a] < arc_crowding[worst]) worst = a;
                if (worst != arc_count - 1)
                    arc.copy_individual(arc_count - 1, worst);
                arc_count--;
            }
            arc_dirty = true;
        };

        // Add initial population to archive (only non-dominated survive)
        {
            std::vector<double> g(dim), o(n_obj);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < dim; ++j) g[j] = pop.gene(i, j);
                for (int oo = 0; oo < n_obj; ++oo) o[oo] = pop.objective(i, oo);
                archive_add(g, o);
            }
        }

        // ── C9 neighborhood: 8-directional Moore on rows×cols toroidal grid ────
        // Returns 8 distinct neighbor indices; does not include self.
        auto get_neighbors = [&](int idx) -> std::vector<int> {
            int row = idx / cols;
            int col = idx % cols;
            static constexpr int DR[8] = {-1, -1, -1,  0,  0,  1,  1,  1};
            static constexpr int DC[8] = {-1,  0,  1, -1,  1, -1,  0,  1};
            std::vector<int> nbrs(8);
            for (int k = 0; k < 8; ++k)
                nbrs[k] = ((row + DR[k] + rows) % rows) * cols +
                           (col + DC[k] + cols) % cols;
            return nbrs;
        };

        // ── Binary tournament: lower rank wins; if tied, higher crowding wins ──
        // Pool uses population indices; rank/crowding from pop_rank/pop_crowding arrays.
        auto tournament_pop = [&](const std::vector<int>& pool) -> int {
            int a = pool[rng.uniform_int(0, (int)pool.size() - 1)];
            int b = pool[rng.uniform_int(0, (int)pool.size() - 1)];
            if (pop_rank[a] != pop_rank[b]) return pop_rank[a] < pop_rank[b] ? a : b;
            return pop_crowding[a] >= pop_crowding[b] ? a : b;
        };

        // Binary tournament from archive using crowding distance (all members non-dominated).
        auto tournament_archive = [&]() -> int {
            compute_archive_crowding();
            int a = rng.uniform_int(0, arc_count - 1);
            int b = rng.uniform_int(0, arc_count - 1);
            return arc_crowding[a] >= arc_crowding[b] ? a : b;
        };

        // ── Micro-NDS + crowding distance for small pools ────────────────────
        // pool_objs[i][o]: objectives of pool member i. Outputs rank and crowding per slot.
        auto micro_nds_crowding = [&](const std::vector<std::vector<double>>& pool_objs,
                                      int pool_size,
                                      std::vector<int>& rank_out,
                                      std::vector<double>& crowding_out) {
            rank_out.assign(pool_size, 0);
            crowding_out.assign(pool_size, 0.0);

            std::vector<int> dom_count(pool_size, 0);
            std::vector<std::vector<int>> dom_set(pool_size);
            for (int p = 0; p < pool_size; ++p) {
                for (int q = p + 1; q < pool_size; ++q) {
                    bool pd = false, qd = false;
                    for (int o = 0; o < n_obj; ++o) {
                        if (pool_objs[p][o] < pool_objs[q][o]) pd = true;
                        else if (pool_objs[q][o] < pool_objs[p][o]) qd = true;
                    }
                    if (pd && !qd) { dom_set[p].push_back(q); dom_count[q]++; }
                    else if (qd && !pd) { dom_set[q].push_back(p); dom_count[p]++; }
                }
            }

            std::vector<std::vector<int>> fronts;
            std::vector<int> cur;
            for (int p = 0; p < pool_size; ++p) if (dom_count[p] == 0) cur.push_back(p);
            for (int r = 0; !cur.empty(); ++r) {
                fronts.push_back(cur);
                for (int p : cur) rank_out[p] = r;
                std::vector<int> next;
                for (int p : cur)
                    for (int q : dom_set[p])
                        if (--dom_count[q] == 0) next.push_back(q);
                cur = std::move(next);
            }

            // Crowding distance per front
            for (const auto& front : fronts) {
                int fs = (int)front.size();
                if (fs <= 2) {
                    for (int i : front) crowding_out[i] = std::numeric_limits<double>::infinity();
                    continue;
                }
                // For each objective, sort front members and assign crowding
                std::vector<int> sf(fs);
                std::iota(sf.begin(), sf.end(), 0);
                for (int o = 0; o < n_obj; ++o) {
                    std::sort(sf.begin(), sf.end(), [&](int a, int b) {
                        return pool_objs[front[a]][o] < pool_objs[front[b]][o];
                    });
                    double fmin = pool_objs[front[sf[0]]][o];
                    double fmax = pool_objs[front[sf[fs - 1]]][o];
                    if (fmin == fmax) continue;
                    double range = fmax - fmin;
                    crowding_out[front[sf[0]]]      = std::numeric_limits<double>::infinity();
                    crowding_out[front[sf[fs - 1]]] = std::numeric_limits<double>::infinity();
                    for (int i = 1; i < fs - 1; ++i)
                        crowding_out[front[sf[i]]] +=
                            (pool_objs[front[sf[i + 1]]][o] - pool_objs[front[sf[i - 1]]][o]) / range;
                }
            }
        };

        // ── Cross-population dominance ────────────────────────────────────────
        auto cross_dom = [&](const Population<>& pa, int ia,
                              const Population<>& pb, int ib) -> Dominance {
            bool a_dom = false, b_dom = false;
            for (int o = 0; o < n_obj; ++o) {
                double fa = pa.objective(ia, o), fb = pb.objective(ib, o);
                if (fa < fb) a_dom = true;
                else if (fb < fa) b_dom = true;
                if (a_dom && b_dom) return Dominance::Equal;
            }
            if (a_dom && !b_dom) return Dominance::Dominates;
            if (b_dom && !a_dom) return Dominance::Dominated;
            return Dominance::Equal;
        };

        // ── Working populations ───────────────────────────────────────────────
        // slot 0 = child, slot 1 = second parent staging for crossover
        Population<> offspring(2, dim, n_obj, pop.n_const);
        offspring.lower_bounds = pop.lower_bounds;
        offspring.upper_bounds = pop.upper_bounds;

        // Reusable scratch vectors for archive_add calls
        std::vector<double> tmp_g(dim), tmp_o(n_obj);

        // ── Main loop ─────────────────────────────────────────────────────────
        int current = 0;
        while (evals < self.max_evals) {

            // 1. Neighborhood: 8 C9 neighbors + self = 9 total (jMetal selection pool)
            auto nbrs = get_neighbors(current);
            nbrs.push_back(current);  // index 8 = self

            // 2. Parent 1: binary tournament from neighborhood → offspring[0]
            int p1 = tournament_pop(nbrs);
            for (int j = 0; j < dim; ++j) offspring.gene(0, j) = pop.gene(p1, j);
            for (int o = 0; o < n_obj; ++o) offspring.objective(0, o) = pop.objective(p1, o);
            offspring.set_evaluated(0, true);

            // 3. Parent 2: binary tournament from archive (if size > 1), else neighborhood → offspring[1]
            if (arc_count > 1) {
                int ai = tournament_archive();
                std::copy_n(arc.genes_ptr(ai), dim, offspring.genes_ptr(1));
                std::copy_n(arc.objectives_ptr(ai), n_obj, offspring.objectives_ptr(1));
            } else {
                int p2 = tournament_pop(nbrs);
                for (int j = 0; j < dim; ++j) offspring.gene(1, j) = pop.gene(p2, j);
                for (int o = 0; o < n_obj; ++o) offspring.objective(1, o) = pop.objective(p2, o);
            }
            offspring.set_evaluated(1, true);

            // 4. Crossover + mutation → offspring[0]
            self.crossover.apply(offspring, 0, 1, 0);
            self.mutation.apply(offspring, 0);

            // 5. Evaluate and count
            problem(offspring, 0);
            offspring.set_evaluated(0, true);
            ++evals;

            // 6. Replacement
            auto dom = cross_dom(offspring, 0, pop, current);

            if (dom == Dominance::Dominates) {
                // Offspring dominates current → replace, add to archive
                for (int j = 0; j < dim; ++j) pop.gene(current, j) = offspring.gene(0, j);
                for (int o = 0; o < n_obj; ++o) pop.objective(current, o) = offspring.objective(0, o);
                pop.set_evaluated(current, true);
                // Rank/crowding from offspring: unranked (-1, 0.0), matching jMetal
                pop_rank[current]     = -1;
                pop_crowding[current] = 0.0;

                for (int j = 0; j < dim; ++j) tmp_g[j] = offspring.gene(0, j);
                for (int o = 0; o < n_obj; ++o) tmp_o[o] = offspring.objective(0, o);
                archive_add(tmp_g, tmp_o);

            } else if (dom == Dominance::Equal) {
                // Non-dominated: NDS + crowding on (8 neighbors + self + offspring) = 10
                const int pool_size = (int)nbrs.size() + 1;  // 9 + 1 = 10
                std::vector<std::vector<double>> pool_objs(pool_size, std::vector<double>(n_obj));
                for (int k = 0; k < (int)nbrs.size(); ++k)
                    for (int o = 0; o < n_obj; ++o)
                        pool_objs[k][o] = pop.objective(nbrs[k], o);
                for (int o = 0; o < n_obj; ++o)
                    pool_objs[pool_size - 1][o] = offspring.objective(0, o);

                std::vector<int>    pool_rank;
                std::vector<double> pool_crowding;
                micro_nds_crowding(pool_objs, pool_size, pool_rank, pool_crowding);

                // Update rank/crowding for neighborhood members (slots 0..8 → pop indices)
                for (int k = 0; k < (int)nbrs.size(); ++k) {
                    pop_rank[nbrs[k]]     = pool_rank[k];
                    pop_crowding[nbrs[k]] = pool_crowding[k];
                }

                // Find worst: highest rank; if tied, lowest crowding
                int worst_slot = 0;
                for (int k = 1; k < pool_size; ++k) {
                    if (pool_rank[k] > pool_rank[worst_slot] ||
                        (pool_rank[k] == pool_rank[worst_slot] &&
                         pool_crowding[k] < pool_crowding[worst_slot]))
                        worst_slot = k;
                }

                // Add offspring to archive
                for (int j = 0; j < dim; ++j) tmp_g[j] = offspring.gene(0, j);
                for (int o = 0; o < n_obj; ++o) tmp_o[o] = offspring.objective(0, o);
                archive_add(tmp_g, tmp_o);

                // Replace current if offspring is not the worst
                if (worst_slot != pool_size - 1) {
                    for (int j = 0; j < dim; ++j) pop.gene(current, j) = offspring.gene(0, j);
                    for (int o = 0; o < n_obj; ++o) pop.objective(current, o) = offspring.objective(0, o);
                    pop.set_evaluated(current, true);
                    // Offspring's pool rank/crowding (slot pool_size-1) assigned to current
                    pop_rank[current]     = pool_rank[pool_size - 1];
                    pop_crowding[current] = pool_crowding[pool_size - 1];
                }
            }
            // dom == Dominated: discard offspring, no archive update

            current = (current + 1) % n;
        }

        // ── Result: copy archive to population ───────────────────────────────
        int result_size = std::min(arc_count, n);
        for (int i = 0; i < result_size; ++i) {
            std::copy_n(arc.genes_ptr(i), dim, pop.genes_ptr(i));
            std::copy_n(arc.objectives_ptr(i), n_obj, pop.objectives_ptr(i));
            pop.set_evaluated(i, true);
        }
        pop.pop_size = result_size;
    }
};

} // namespace ea
