#pragma once
/// @file config.hpp
/// @brief Runtime configuration for ea-cpp algorithms via JSON.
///
/// Allows creating and running algorithms without recompiling.
/// Uses nlohmann/json for parsing.

#include <ea/core/population.hpp>
#include <ea/algorithm/nsga2.hpp>
#include <ea/algorithm/rvea.hpp>
#include <ea/algorithm/ibea.hpp>
#include <ea/algorithm/random_search.hpp>
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/crossover/blx_alpha.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/operator/mutation/bit_flip.hpp>
#include <nlohmann/json.hpp>
#include <functional>
#include <stdexcept>
#include <fstream>
#include <string>

namespace ea {

using json = nlohmann::json;

/// Configuration for ea-cpp algorithms
struct Config {
    std::string algorithm = "NSGAII";
    std::string crossover_type = "sbx";
    std::string mutation_type = "polynomial";
    int pop_size = 100;
    int max_evals = 25000;
    int dim = 30;
    int n_obj = 2;

    // Crossover params
    double cx_probability = 0.9;
    double cx_distribution_index = 20.0;
    double cx_alpha = 0.5;

    // Mutation params
    double mt_rate = -1.0; // -1 = auto 1/dim
    double mt_distribution_index = 20.0;

    /// Load from JSON file
    static Config from_file(const std::string& filename) {
        std::ifstream f(filename);
        if (!f.is_open()) {
            throw std::runtime_error("Cannot open config file: " + filename);
        }
        json j;
        f >> j;
        return from_json(j);
    }

    /// Load from JSON object
    static Config from_json(const json& j) {
        Config cfg;
        cfg.algorithm = j.value("algorithm", "NSGAII");

        auto cx = j.value("crossover", json::object());
        cfg.crossover_type = cx.value("type", "sbx");
        cfg.cx_probability = cx.value("probability", 0.9);
        cfg.cx_distribution_index = cx.value("distribution_index", 20.0);
        cfg.cx_alpha = cx.value("alpha", 0.5);

        auto mt = j.value("mutation", json::object());
        cfg.mutation_type = mt.value("type", "polynomial");
        cfg.mt_rate = mt.value("rate", -1.0);
        cfg.mt_distribution_index = mt.value("distribution_index", 20.0);

        cfg.pop_size = j.value("pop_size", 100);
        cfg.max_evals = j.value("max_evals", 25000);
        cfg.dim = j.value("dim", 30);
        cfg.n_obj = j.value("n_obj", 2);

        return cfg;
    }
};

/// Example JSON config file:
/// {
///   "algorithm": "NSGAII",
///   "crossover": { "type": "sbx", "probability": 0.9, "distribution_index": 20.0 },
///   "mutation": { "type": "polynomial", "rate": -1.0, "distribution_index": 20.0 },
///   "pop_size": 100,
///   "max_evals": 25000
/// }

} // namespace ea
