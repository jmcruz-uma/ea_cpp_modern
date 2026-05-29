#pragma once
/// @file stats_collector.hpp
/// @brief Per-generation statistics collector for single-objective, minimisation EAs.
///
/// Produces JSON output byte-compatible with Bio4Res/ea (ea_cpp_original):
///
///   [
///     {
///       "run": 0,  "seed": 1,  "time": 7.27,
///       "rundata": [{
///         "idata": { "evals": [...], "best": [...], "mean": [...], "diversity": [...] },
///         "isols": { "evals": [...], "fitness": [...], "genome": [[...], ...] }
///       }]
///     }, ...
///   ]
///
/// Usage (one StatsCollector across all runs of an experiment):
///
///   StatsCollector stats;
///   stats.new_run(seed);
///   // -- initialise population, evaluate --
///   stats.take(pop, evals);          // snapshot after init
///   while (evals < max_evals) {
///       // -- stepUp: select, vary, evaluate, replace --
///       stats.take(pop, evals);      // snapshot each generation
///   }
///   stats.close_run();
///   // ... repeat for next run ...
///   std::ofstream f("problem-n-seed-stats.json");
///   f << stats.to_json();
///
/// Diversity formula mirrors Bio4Res/ea VarianceDiversity::apply():
///   For each gene dimension j, compute the population std-dev (biased, N divisor),
///   then average across all dimensions.

#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include <ea/core/population.hpp>
#include <ea/util/json.hpp>

namespace ea {

struct StatsCollector {

    // ------------------------------------------------------------------
    // Internal records
    // ------------------------------------------------------------------

    struct Snapshot {
        uint64_t evals;
        double   best;
        double   mean;
        double   diversity;
    };

    struct Improvement {
        uint64_t             evals;
        double               fitness;
        std::vector<double>  genome;
    };

    struct RunRecord {
        uint64_t                  seed   = 0;
        double                    time_s = 0.0;
        std::vector<Snapshot>     snapshots;
        std::vector<Improvement>  improvements;
    };

    // ------------------------------------------------------------------
    // State
    // ------------------------------------------------------------------

    std::vector<RunRecord> runs;

private:
    RunRecord   current_;
    double      best_seen_ = std::numeric_limits<double>::max();
    bool        active_    = false;
    std::chrono::steady_clock::time_point tic_;

    // ------------------------------------------------------------------
    // Diversity — exact replica of Bio4Res/ea VarianceDiversity::apply()
    //
    //   For each dimension i:
    //     K     = pop.gene(0, i)                  (Welford shift)
    //     sum   = Σ_j (gene(j,i) − K)
    //     sum2  = Σ_j (gene(j,i) − K)²
    //     σ_i   = sqrt( (sum2 − sum²/mu) / mu )   (biased std-dev)
    //   return  Σ_i σ_i / n_dims
    // ------------------------------------------------------------------
    static double variance_diversity(const Population<>& pop) {
        const int mu = pop.pop_size;
        const int n  = pop.dim;
        if (mu < 2 || n == 0) return 0.0;

        double sigma = 0.0;
        for (int i = 0; i < n; ++i) {
            const double K = pop.gene(0, i);
            double sum = 0.0, sum2 = 0.0;
            for (int j = 0; j < mu; ++j) {
                const double v = pop.gene(j, i) - K;
                sum  += v;
                sum2 += v * v;
            }
            const double var = (sum2 - (sum * sum) / mu) / mu;
            sigma += std::sqrt(var > 0.0 ? var : 0.0);
        }
        return sigma / n;
    }

public:

    // ------------------------------------------------------------------
    // Run lifecycle
    // ------------------------------------------------------------------

    /// Start a new run.  Must be called before the first take().
    void new_run(uint64_t seed) {
        if (active_) close_run();
        current_   = RunRecord{};
        current_.seed = seed;
        best_seen_ = std::numeric_limits<double>::max();
        active_    = true;
        tic_       = std::chrono::steady_clock::now();
    }

    /// Snapshot after initialisation or after each generation.
    /// @param pop   Current population (post-replacement).
    /// @param evals Total objective-function evaluations so far.
    void take(const Population<>& pop, uint64_t evals) {
        // --- find best and compute mean (single-objective, minimisation) ---
        int    best_idx     = 0;
        double best_fitness = pop.objective(0, 0);
        double mean_fitness = best_fitness;

        for (int i = 1; i < pop.pop_size; ++i) {
            const double f = pop.objective(i, 0);
            mean_fitness += f;
            if (f < best_fitness) {
                best_fitness = f;
                best_idx     = i;
            }
        }
        mean_fitness /= pop.pop_size;

        // --- diversity ---
        const double div = variance_diversity(pop);

        // --- record snapshot ---
        current_.snapshots.push_back({evals, best_fitness, mean_fitness, div});

        // --- record improvement ---
        if (best_fitness < best_seen_) {
            best_seen_ = best_fitness;
            std::vector<double> genome(pop.dim);
            for (int j = 0; j < pop.dim; ++j)
                genome[j] = pop.gene(best_idx, j);
            current_.improvements.push_back({evals, best_fitness, std::move(genome)});
        }
    }

    /// Close the current run and append it to runs[].
    void close_run() {
        if (!active_) return;
        const auto toc = std::chrono::steady_clock::now();
        current_.time_s =
            std::chrono::duration<double>(toc - tic_).count();
        runs.push_back(std::move(current_));
        active_ = false;
    }

    // ------------------------------------------------------------------
    // JSON serialisation — matches ea_cpp_original EAStatistics::toJSON()
    // ------------------------------------------------------------------

    nlohmann::json to_json() const {
        auto out = nlohmann::json::array();

        for (int r = 0; r < static_cast<int>(runs.size()); ++r) {
            const auto& run = runs[r];

            // idata — parallel arrays
            auto evals_arr = nlohmann::json::array();
            auto best_arr  = nlohmann::json::array();
            auto mean_arr  = nlohmann::json::array();
            auto div_arr   = nlohmann::json::array();
            for (const auto& s : run.snapshots) {
                evals_arr.push_back(s.evals);
                best_arr.push_back(s.best);
                mean_arr.push_back(s.mean);
                div_arr.push_back(s.diversity);
            }
            nlohmann::json idata;
            idata["evals"]     = std::move(evals_arr);
            idata["best"]      = std::move(best_arr);
            idata["mean"]      = std::move(mean_arr);
            idata["diversity"] = std::move(div_arr);

            // isols — improvement history
            auto isols_evals   = nlohmann::json::array();
            auto isols_fitness = nlohmann::json::array();
            auto isols_genome  = nlohmann::json::array();
            for (const auto& imp : run.improvements) {
                isols_evals.push_back(imp.evals);
                isols_fitness.push_back(imp.fitness);
                auto g = nlohmann::json::array();
                for (double v : imp.genome) g.push_back(v);
                isols_genome.push_back(std::move(g));
            }
            nlohmann::json isols;
            isols["evals"]   = std::move(isols_evals);
            isols["fitness"] = std::move(isols_fitness);
            isols["genome"]  = std::move(isols_genome);

            // rundata (one island)
            nlohmann::json island_data;
            island_data["idata"] = std::move(idata);
            island_data["isols"] = std::move(isols);

            nlohmann::json run_obj;
            run_obj["run"]     = r;
            run_obj["seed"]    = run.seed;
            run_obj["time"]    = run.time_s;
            run_obj["rundata"] = nlohmann::json::array({std::move(island_data)});

            out.push_back(std::move(run_obj));
        }
        return out;
    }

    // ------------------------------------------------------------------
    // Convenience
    // ------------------------------------------------------------------

    /// Number of completed runs.
    int num_runs() const { return static_cast<int>(runs.size()); }

    /// Wall-clock time of run i.
    double time(int i) const { return runs.at(i).time_s; }

    /// Best fitness achieved across all completed runs.
    double overall_best() const {
        double b = std::numeric_limits<double>::max();
        for (const auto& r : runs)
            if (!r.snapshots.empty())
                b = std::min(b, r.snapshots.back().best);
        return b;
    }
};

} // namespace ea
