#pragma once
/// @file ibea.hpp
/// @brief IBEA — Indicator-Based Evolutionary Algorithm
///
/// Reference: E. Zitzler and S. Kunzli, "Indicator-Based Selection in
/// Multiobjective Search", PPSN 2004.
///
/// Algorithm logic mirrors jMetal 7.4 IBEA.java exactly:
///   - Binary hypervolume indicator (HVI) with rho = 2.0
///   - Separate archive starting empty; union = solutionSet ∪ archive each gen
///   - Fitness: sum_{i!=pos} exp(-I(i,pos) / (maxIV * kappa))
///     where I(j,i) < 0 when j dominates i, positive otherwise
///   - Environmental selection: remove argmax-fitness iteratively with
///     incremental fitness update (mirrors removeWorst())
///   - Parent selection: binary tournament — lower fitness wins
///     (mirrors BinaryTournamentSelection<FitnessComparator>)
///   - Offspring: SBX + PM on child 0 only (mirrors offspring.get(0))
///
/// Data-structure differences from jMetal for C++ performance:
///   Archive objs/genes : contiguous flat vectors (SoA layout, no pointer chains)
///   Indicator matrix   : flat row-major vector with fixed stride = 2·N
///   Removal            : swap-and-decrement O(N) instead of ArrayList.remove O(N²)
///
/// ── std::mdspan note (future improvement) ────────────────────────────────────
/// The flat indicator matrix ind[j * stride + i] would be written as:
///
///   // std::mdspan<double,
///   //     std::extents<int, std::dynamic_extent, std::dynamic_extent>>
///   //     ind_view(ind.data(), arch_sz, ind_stride);
///   // Access: ind_view[j, i]  instead of  ind[j * ind_stride + i]
///
/// std::mdspan is part of C++23 but not yet shipped in g++14.
/// Re-enable when migrating to g++16 (ships std::mdspan).
/// ─────────────────────────────────────────────────────────────────────────────

#include <algorithm>
#include <cmath>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <limits>
#include <vector>

namespace ea {

template <typename CX, typename MT>
struct IBEA {
    CX crossover;
    MT mutation;

    int    pop_size  = 100;
    int    max_evals = 25000;
    double kappa     = 0.05; ///< Fitness scaling factor (jMetal default)
    double rho       = 2.0;  ///< Reference-box expansion factor (jMetal default)

    static constexpr std::string_view name() { return "IBEA"; }

    /// Binary hypervolume indicator I_HV(a, b) — mirrors
    /// IBEA.java::calculateHypervolumeIndicator exactly.
    ///
    /// @param objs   Flat objective array with stride `stride`:
    ///               objs[individual * stride + objective_index]
    /// @param stride Number of objectives (= n_obj)
    /// @param a      Index of solution A
    /// @param b      Index of solution B (ignored when b_null = true)
    /// @param d      Current recursion depth, 1-indexed (starts at n_obj)
    /// @param b_null When true, use the reference box as B's boundary
    /// @param max_v  Per-objective maxima over the current archive
    /// @param min_v  Per-objective minima over the current archive
    double hvi(const double* objs, int stride,
               int a, int b, int d, bool b_null,
               const double* max_v, const double* min_v) const {
        const int k  = d - 1;
        double r     = rho * (max_v[k] - min_v[k]);
        double ref   = min_v[k] + r;
        double av    = objs[a * stride + k];
        double bv    = b_null ? ref : objs[b * stride + k];

        if (d == 1)
            return (av < bv) ? (bv - av) / r : 0.0;

        double vol;
        if (av < bv) {
            vol  = hvi(objs, stride, a, -1, d-1, true,  max_v, min_v) * (bv  - av) / r;
            vol += hvi(objs, stride, a,  b, d-1, false, max_v, min_v) * (ref - bv) / r;
        } else {
            vol  = hvi(objs, stride, a,  b, d-1, false, max_v, min_v) * (ref - av) / r;
        }
        return vol;
    }

    void run(this auto& self, Population<>& pop, auto&& problem) {
        const int N          = self.pop_size;
        const int dim        = pop.dim;
        const int n_obj      = pop.n_obj;
        const int cap        = 2 * N;   // max archive capacity (union size)
        const int ind_stride = cap;     // fixed column stride for flat ind matrix

        if (pop.pop_size != N)
            pop.resize(N);

        // Evaluate initial population (= jMetal solutionSet, initially)
        for (int i = 0; i < N; ++i) {
            if (!pop.evaluated(i)) {
                problem(pop, i);
                pop.set_evaluated(i, true);
            }
        }
        int evals = N;
        auto& rng = Random::instance();

        // ── Flat archive storage ──────────────────────────────────────────
        // arch_objs[i * n_obj + o]  = objective o of individual i
        // arch_genes[i * dim  + j]  = gene j of individual i
        //
        // See header comment for the pending std::mdspan upgrade.
        std::vector<double> arch_objs (cap * n_obj, 0.0);
        std::vector<double> arch_genes(cap * dim,   0.0);
        int arch_sz = 0;  // active individuals (0 ≤ arch_sz ≤ cap)

        // ── Flat indicator matrix (row-major, fixed stride = cap) ─────────
        // ind[j * ind_stride + i]  = I(j, i)
        // Active region: rows/cols 0 .. arch_sz-1.
        //
        // See header comment for the pending std::mdspan upgrade.
        std::vector<double> ind(cap * cap, 0.0);
        double max_ind = 1.0;
        std::vector<double> fitness(cap, 0.0);

        // Reusable bounds buffers (reallocating per call would be wasteful)
        std::vector<double> max_v(n_obj), min_v(n_obj);

        // ── calculate_fitness ─────────────────────────────────────────────
        // Rebuilds ind[][] and fitness[] for all arch_sz individuals.
        // Mirrors computeIndicatorValuesHD() + fitness() in jMetal.
        auto calculate_fitness = [&]() {
            const int sz          = arch_sz;
            const double* ao      = arch_objs.data();

            // Per-objective bounds over the current archive
            for (int o = 0; o < n_obj; ++o) {
                max_v[o] = -std::numeric_limits<double>::max();
                min_v[o] =  std::numeric_limits<double>::max();
            }
            for (int i = 0; i < sz; ++i)
                for (int o = 0; o < n_obj; ++o) {
                    double v = ao[i * n_obj + o];
                    if (v > max_v[o]) max_v[o] = v;
                    if (v < min_v[o]) min_v[o] = v;
                }

            // Build indicator matrix.
            // I(j,i) = -HVI(j,i)  if j strictly dominates i  [jMetal: flag==-1]
            // I(j,i) = +HVI(i,j)  otherwise                  [jMetal: else]
            double mi = -std::numeric_limits<double>::max();
            const double* mv = max_v.data();
            const double* nv = min_v.data();
            for (int j = 0; j < sz; ++j) {
                double* row_j = ind.data() + j * ind_stride;
                for (int i = 0; i < sz; ++i) {
                    // Strict dominance: j ≤ i in all objectives, j < i in at least one
                    bool leq = true, strict = false;
                    for (int o = 0; o < n_obj; ++o) {
                        double dj = ao[j * n_obj + o], di = ao[i * n_obj + o];
                        if (dj > di) { leq = false; break; }
                        if (dj < di) strict = true;
                    }
                    double value = (leq && strict)
                        ? -self.hvi(ao, n_obj, j, i, n_obj, false, mv, nv)
                        :  self.hvi(ao, n_obj, i, j, n_obj, false, mv, nv);
                    row_j[i] = value;
                    double av = std::abs(value);
                    if (av > mi) mi = av;
                }
            }
            max_ind = mi;

            // fitness[pos] = Σ_{i≠pos} exp(-ind[i*stride+pos] / (max_ind·κ))
            // Outer loop over row i → sequential access of row i → cache-friendly.
            const double denom = max_ind * self.kappa;
            std::fill_n(fitness.data(), sz, 0.0);
            for (int i = 0; i < sz; ++i) {
                const double* row_i = ind.data() + i * ind_stride;
                for (int pos = 0; pos < sz; ++pos)
                    if (pos != i)
                        fitness[pos] += std::exp(-row_i[pos] / denom);
            }
        };

        // ── remove_worst ──────────────────────────────────────────────────
        // Removes the individual with the highest fitness (most dominated).
        // Uses swap-and-decrement: O(N) per removal, all accesses sequential
        // or stride-N (L2-resident for N ≤ 200).
        //
        // Operation order is load-bearing:
        //   1. Fitness update uses ind[worst][*] → must precede row swap.
        //   2. Row/col swaps in ind after fitness update.
        //   3. Archive and fitness[] entry swaps.
        //   4. arch_sz--.
        auto remove_worst = [&]() {
            const int sz   = arch_sz;
            const int last = sz - 1;

            // 1. Find worst (argmax fitness)
            int worst = 0;
            for (int i = 1; i < sz; ++i)
                if (fitness[i] > fitness[worst]) worst = i;

            // 2. Incremental fitness update — subtract worst's contribution.
            //    Sequential read of ind row `worst` → prefetch-friendly.
            const double denom  = max_ind * self.kappa;
            const double* row_w = ind.data() + worst * ind_stride;
            for (int i = 0; i < sz; ++i)
                if (i != worst)
                    fitness[i] -= std::exp(-row_w[i] / denom);

            // Fast path: worst is already the last slot — just shrink.
            if (worst == last) { --arch_sz; return; }

            // 3. Swap ind rows worst ↔ last (sequential ranges, no overlap).
            std::swap_ranges(ind.data() + worst * ind_stride,
                             ind.data() + worst * ind_stride + sz,
                             ind.data() + last  * ind_stride);

            // 4. Swap ind cols worst ↔ last.
            //    Stride-access (ind_stride = 2N doubles apart).
            //    Both columns are L2-resident for typical N ≤ 200.
            for (int j = 0; j < sz; ++j)
                std::swap(ind[j * ind_stride + worst], ind[j * ind_stride + last]);

            // 5. Swap archive individuals worst ↔ last (flat SoA, no alloc).
            std::swap_ranges(arch_objs.data()  + worst * n_obj,
                             arch_objs.data()  + worst * n_obj + n_obj,
                             arch_objs.data()  + last  * n_obj);
            std::swap_ranges(arch_genes.data() + worst * dim,
                             arch_genes.data() + worst * dim + dim,
                             arch_genes.data() + last  * dim);

            // 6. Swap fitness entries.
            std::swap(fitness[worst], fitness[last]);

            --arch_sz;
        };

        // ── binary_tournament ─────────────────────────────────────────────
        // Returns the archive index with lower fitness from two random picks.
        // Mirrors BinaryTournamentSelection<FitnessComparator> in jMetal:
        //   FitnessComparator returns -1 when fitness_a < fitness_b,
        //   so the tournament returns the individual with lower fitness.
        auto binary_tournament = [&]() -> int {
            int a = rng.uniform_int(0, arch_sz - 1);
            int b = rng.uniform_int(0, arch_sz - 1);
            return (fitness[a] <= fitness[b]) ? a : b;
        };

        // Scratch Population<> for SBX — two slots (child 0 and child 1).
        // Only child 0 is evaluated and used (jMetal: offspring.get(0)).
        Population<> scratch(2, dim, n_obj);
        scratch.lower_bounds = pop.lower_bounds;
        scratch.upper_bounds = pop.upper_bounds;

        while (evals < self.max_evals) {
            // === 1. union = solutionSet ∪ archive (archive starts empty) ===
            // Append current solutionSet (pop[0..N-1]) to the archive.
            for (int i = 0; i < N; ++i) {
                const int dst = arch_sz + i;
                for (int o = 0; o < n_obj; ++o)
                    arch_objs[dst * n_obj + o] = pop.objective(i, o);
                for (int j = 0; j < dim; ++j)
                    arch_genes[dst * dim + j] = pop.gene(i, j);
            }
            arch_sz += N;

            // === 2. calculateFitness(union); archive = union ===
            calculate_fitness();

            // === 3. Prune archive to N ===
            while (arch_sz > N)
                remove_worst();

            // === 4. Generate N offspring from archive ===
            // Two independent binary tournaments → SBX → PM on child 0 → evaluate.
            // Offspring written directly into pop as the new solutionSet.
            for (int i = 0; i < N && evals < self.max_evals; ++i) {
                const int p1 = binary_tournament();
                const int p2 = binary_tournament();

                // Load parents into scratch slots 0 (p1) and 1 (p2)
                for (int j = 0; j < dim; ++j) {
                    scratch.gene(0, j) = arch_genes[p1 * dim + j];
                    scratch.gene(1, j) = arch_genes[p2 * dim + j];
                }
                scratch.set_evaluated(0, false);
                scratch.set_evaluated(1, false);

                // SBX — writes children into slots 0 and 1
                self.crossover.apply(scratch, 0, 1, 0);

                // PM on child 0 only (jMetal: offspring.get(0))
                self.mutation.apply(scratch, 0);

                // Evaluate child 0
                problem(scratch, 0);
                scratch.set_evaluated(0, true);
                ++evals;

                // Store in pop[i] (= next iteration's solutionSet entry)
                for (int j = 0; j < dim;   ++j) pop.gene(i, j)      = scratch.gene(0, j);
                for (int o = 0; o < n_obj; ++o) pop.objective(i, o) = scratch.objective(0, o);
                pop.set_evaluated(i, true);
            }
            // solutionSet for next iteration = pop[0..N-1] written above
        }

        // === Result: copy archive (arch_sz = N) to pop ===
        // jMetal: result() = getNonDominatedSolutions(archive).
        // We expose the full archive; caller may filter non-dominated if needed.
        for (int i = 0; i < arch_sz; ++i) {
            for (int j = 0; j < dim;   ++j) pop.gene(i, j)      = arch_genes[i * dim + j];
            for (int o = 0; o < n_obj; ++o) pop.objective(i, o) = arch_objs[i * n_obj + o];
            pop.set_evaluated(i, true);
        }
    }
};

} // namespace ea
