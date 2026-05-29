#pragma once
/// @file seed_manager.hpp
/// @brief File-based seed management compatible with Bio4Res/ea (ea_cpp_original).
///
/// Replicates the exact protocol from ea_cpp_original/ea_main.cpp:
///
///   Protocol
///   --------
///   LOAD — determine seed for this batch:
///     • If lastseed.txt does not exist → use config_seed as-is (first-ever run).
///     • If lastseed.txt exists         → seed = file_value + 1.
///
///   COMMIT — update the file after numruns runs complete:
///     • file_seed = seed + numruns - 1   (last seed consumed in this batch)
///     • If file_seed < maxruns → write file_seed to lastseed.txt.
///     • If file_seed >= maxruns → delete lastseed.txt (experiment finished).
///
///   Output filename convention (mirrors original):
///     {problem}-{numvars}-{seed}-stats.json
///   where seed is the value returned by load(), i.e. the FIRST seed of the batch.
///
///   Example — 30 runs, numruns=1 per invocation, config_seed=1:
///     inv 1:  seed=1   → runs with 1, writes 1  → file: "problem-n-1-stats.json"
///     inv 2:  seed=2   → runs with 2, writes 2  → file: "problem-n-2-stats.json"
///     ...
///     inv 30: seed=30  → runs with 30, 30 >= 30 → deletes lastseed.txt

#include <cstdio>   // std::remove
#include <fstream>
#include <stdexcept>
#include <string>

namespace ea {

struct SeedManager {
    /// The seed to use for the current batch of runs (and for the output filename).
    long seed = 1;

    /// Path to the state file (default matches the original).
    std::string filename = "lastseed.txt";

    // ------------------------------------------------------------------
    // Factory
    // ------------------------------------------------------------------

    /// Load state from disk.
    ///
    /// @param config_seed  Fallback seed from the algorithm config JSON
    ///                     (used only when no lastseed.txt exists).
    /// @param filename     Path to the state file.
    static SeedManager load(long config_seed, const std::string& filename = "lastseed.txt") {
        SeedManager mgr;
        mgr.filename = filename;

        std::ifstream f(filename);
        if (f.fail()) {
            // First run: no file → use config seed directly (original: seed = conf.seed)
            mgr.seed = config_seed;
        } else {
            long last = 0;
            f >> last;
            // Subsequent run: original does `conf.seed = ++seed` after reading
            mgr.seed = last + 1;
        }
        return mgr;
    }

    // ------------------------------------------------------------------
    // Commit
    // ------------------------------------------------------------------

    /// Update the seed file after completing @p numruns runs.
    ///
    /// @param numruns  Number of runs just completed in this invocation.
    /// @param maxruns  Total runs planned for the experiment (from experiment.json).
    void commit(int numruns, long maxruns) {
        // Mirrors original: seed += numruns - 1
        const long file_seed = seed + static_cast<long>(numruns) - 1L;

        if (file_seed < maxruns) {
            std::ofstream f(filename);
            if (!f) throw std::runtime_error("SeedManager: cannot write " + filename);
            f << file_seed;
        } else {
            std::remove(filename.c_str());
        }
    }

    // ------------------------------------------------------------------
    // Convenience
    // ------------------------------------------------------------------

    /// Returns true if lastseed.txt exists (experiment in progress).
    static bool in_progress(const std::string& filename = "lastseed.txt") {
        return std::ifstream(filename).good();
    }

    /// Reset: delete the file so the next load() starts from config_seed.
    static void reset(const std::string& filename = "lastseed.txt") {
        std::remove(filename.c_str());
    }
};

} // namespace ea
