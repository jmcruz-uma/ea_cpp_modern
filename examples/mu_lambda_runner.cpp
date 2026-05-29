/// @file mu_lambda_runner.cpp
/// @brief (µ,λ)-ES runner for single-objective continuous optimisation.
///
/// Reads two JSON files (same format as ea_cpp_original):
///   numeric.json   — algorithm config: pop_size, offspring, max_evals, operators
///   experiment.json — problem config: problem name, numvars, range, numruns
///
/// Outputs a stats JSON file byte-compatible with ea_cpp_original's output,
/// using StatsCollector and SeedManager.
///
/// Usage:
///   ./mu_lambda_runner [numeric.json] [experiment.json] [output.json]
/// Defaults: numeric.json  experiment.json  <problem>-<nvars>-<seed>-stats.json

#include <ea/algorithm/mu_lambda_es.hpp>
#include <ea/operator/crossover/blx_alpha_beta.hpp>
#include <ea/operator/mutation/gaussian.hpp>
#include <ea/problem/continuous_so.hpp>
#include <ea/util/json.hpp>
#include <ea/util/random.hpp>
#include <ea/util/seed_manager.hpp>
#include <ea/util/stats_collector.hpp>

#include <ea/core/population.hpp>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

// -----------------------------------------------------------------------
// Config structs
// -----------------------------------------------------------------------

struct OperatorSpec {
    std::string name;
    std::vector<std::string> parameters;
};

struct IslandConfig {
    int pop_size      = 100;
    int offspring     = 99;
    long max_evals    = 10'000'000;
    int tournament_k  = 2;
    std::vector<OperatorSpec> variation;
};

struct NumericConfig {
    int numruns     = 1;
    uint64_t seed   = 1;
    IslandConfig island;
};

struct ExperimentConfig {
    int numruns  = 30;
    std::string problem = "rastrigin";
    int numvars  = 20;
    double range = 5.12;
};

// -----------------------------------------------------------------------
// JSON parsing
// -----------------------------------------------------------------------

NumericConfig parse_numeric(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open " + path);
    auto j = nlohmann::json::parse(f);

    NumericConfig cfg;
    cfg.numruns = j.value("numruns", 1);
    cfg.seed    = j.value("seed", uint64_t{1});

    const auto& island_j = j["islands"][0];
    cfg.island.pop_size  = island_j.value("popsize",   100);
    cfg.island.offspring = island_j.value("offspring", 99);
    cfg.island.max_evals = island_j.value("maxevals",  10'000'000L);

    // Selection: tournament — read k from parameters[0]
    if (island_j.contains("selection")) {
        const auto& sel = island_j["selection"];
        const auto& pars = sel["parameters"];
        if (!pars.empty())
            cfg.island.tournament_k = std::stoi(pars[0].get<std::string>());
    }

    // Variation operators (skip "initialization" at index 0 in the original — not present here)
    if (island_j.contains("variation")) {
        for (const auto& op : island_j["variation"]) {
            OperatorSpec spec;
            spec.name = op["name"];
            for (const auto& p : op["parameters"])
                spec.parameters.push_back(p.get<std::string>());
            cfg.island.variation.push_back(spec);
        }
    }
    return cfg;
}

ExperimentConfig parse_experiment(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open " + path);
    auto j = nlohmann::json::parse(f);

    ExperimentConfig cfg;
    cfg.numruns = j.value("numruns",  30);
    cfg.problem = j.value("problem",  std::string("rastrigin"));
    cfg.numvars = j.value("numvars",  20);
    cfg.range   = j.value("range",    5.12);
    return cfg;
}

// -----------------------------------------------------------------------
// Algorithm factory — builds MuLambdaES from parsed specs
// -----------------------------------------------------------------------

using ESType = ea::MuLambdaES<ea::BLXAlphaBetaCrossover, ea::GaussianMutation>;

ESType build_es(const IslandConfig& cfg) {
    ESType es;
    es.pop_size      = cfg.pop_size;
    es.offspring_size = cfg.offspring;
    es.max_evals     = static_cast<int>(cfg.max_evals);
    es.tournament_k  = cfg.tournament_k;

    // Map variation operators by name (lowercase)
    for (const auto& op : cfg.variation) {
        std::string lname = op.name;
        for (auto& c : lname) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));

        if (lname == "blx") {
            // pars[0] = crossover_probability, pars[1] = alpha (symmetric BLX)
            if (op.parameters.size() >= 1)
                es.crossover.crossover_probability = std::stod(op.parameters[0]);
            if (op.parameters.size() >= 2) {
                es.crossover.alpha = std::stod(op.parameters[1]);
                es.crossover.beta  = es.crossover.alpha; // symmetric BLX-alpha
            }
        } else if (lname == "gaussian") {
            // pars[0] = mutation_rate, pars[1] = step_size
            if (op.parameters.size() >= 1)
                es.mutation.mutation_rate = std::stod(op.parameters[0]);
            if (op.parameters.size() >= 2)
                es.mutation.step_size = std::stod(op.parameters[1]);
        }
        // extend here for other operators as needed
    }
    return es;
}

// -----------------------------------------------------------------------
// Problem factory
// -----------------------------------------------------------------------

struct ProblemVariant {
    enum class Kind { Sphere, Rastrigin, Ackley, Rosenbrock } kind;
    int dim; double range;

    void evaluate(ea::Population<>& pop, int idx) const {
        switch (kind) {
            case Kind::Sphere:     { ea::Sphere{dim, range}.evaluate(pop, idx);     break; }
            case Kind::Rastrigin:  { ea::Rastrigin{dim, range}.evaluate(pop, idx);  break; }
            case Kind::Ackley:     { ea::Ackley{dim, range}.evaluate(pop, idx);     break; }
            case Kind::Rosenbrock: { ea::Rosenbrock{dim, range}.evaluate(pop, idx); break; }
        }
    }
    double lower_bound() const { return -range; }
    double upper_bound() const { return  range; }
};

ProblemVariant build_problem(const ExperimentConfig& cfg) {
    std::string name = cfg.problem;
    for (auto& c : name) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));

    ProblemVariant::Kind kind;
    if      (name == "sphere")     kind = ProblemVariant::Kind::Sphere;
    else if (name == "rastrigin")  kind = ProblemVariant::Kind::Rastrigin;
    else if (name == "ackley")     kind = ProblemVariant::Kind::Ackley;
    else if (name == "rosenbrock") kind = ProblemVariant::Kind::Rosenbrock;
    else throw std::runtime_error("Unknown problem: " + cfg.problem);

    return {kind, cfg.numvars, cfg.range};
}

// -----------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------

int main(int argc, char* argv[]) {
    const std::string numeric_path    = (argc > 1) ? argv[1] : "numeric.json";
    const std::string experiment_path = (argc > 2) ? argv[2] : "experiment.json";

    NumericConfig    ncfg = parse_numeric(numeric_path);
    ExperimentConfig ecfg = parse_experiment(experiment_path);

    // experiment.json numruns takes precedence
    const int total_runs = ecfg.numruns;

    ProblemVariant prob = build_problem(ecfg);
    ESType         es   = build_es(ncfg.island);

    std::string output_path = (argc > 3) ? argv[3]
        : (ecfg.problem + "-" + std::to_string(ecfg.numvars)
           + "-" + std::to_string(ncfg.seed) + "-stats.json");

    // SeedManager — protocol identical to ea_cpp_original (lastseed.txt)
    ea::SeedManager seed_mgr = ea::SeedManager::load(static_cast<long>(ncfg.seed), "lastseed.txt");
    uint64_t seed = static_cast<uint64_t>(seed_mgr.seed);

    ea::StatsCollector stats;
    auto& rng = ea::Random::instance();

    std::cout << "=== MuLambdaES (" << ecfg.problem << ", d=" << ecfg.numvars << ") ===\n"
              << "µ=" << es.pop_size << "  λ=" << es.offspring_size
              << "  max_evals=" << es.max_evals
              << "  runs=" << total_runs << "\n\n";

    for (int r = 0; r < total_runs; ++r) {
        uint64_t run_seed = seed + static_cast<uint64_t>(r);
        rng.seed(run_seed);

        // Initialise population
        ea::Population<> pop(es.pop_size, ecfg.numvars, 1, 0,
                             prob.lower_bound(), prob.upper_bound());

        for (int i = 0; i < es.pop_size; ++i) {
            for (int j = 0; j < ecfg.numvars; ++j)
                pop.gene(i, j) = rng.uniform(prob.lower_bound(), prob.upper_bound());
            pop.set_evaluated(i, false);
        }

        stats.new_run(run_seed);
        es.run(pop, [&prob](ea::Population<>& p, int idx) { prob.evaluate(p, idx); }, stats);
        stats.close_run();

        double best = std::numeric_limits<double>::max();
        for (int i = 0; i < pop.pop_size; ++i)
            best = std::min(best, pop.objective(i, 0));

        std::cout << "  Run " << (r + 1) << "/" << total_runs
                  << "  seed=" << run_seed
                  << "  best=" << std::scientific << best << "\n";
    }

    seed_mgr.commit(total_runs, static_cast<long>(total_runs));

    // Write stats JSON
    std::ofstream out(output_path);
    if (!out) {
        std::cerr << "Cannot write " << output_path << "\n";
        return 1;
    }
    out << stats.to_json().dump(2) << "\n";
    std::cout << "\nStats written to " << output_path << "\n";
    return 0;
}
