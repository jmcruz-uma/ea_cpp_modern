// Example: Using JSON configuration to run algorithms without recompiling
//
// Compile: g++-14 -std=c++23 -O2 -I include config_example.cpp -o config_example -lm
// Usage: ./config_example config.json

#include <ea/ea.hpp>
#include <ea/core/config.hpp>
#include <iostream>
#include <fstream>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config.json>\n";
        return 1;
    }

    // Load configuration
    ea::Config cfg = ea::Config::from_file(argv[1]);

    std::cout << "Algorithm: " << cfg.algorithm << "\n";
    std::cout << "Population: " << cfg.pop_size << "\n";
    std::cout << "Evaluations: " << cfg.max_evals << "\n";
    std::cout << "Dimensions: " << cfg.dim << "\n";
    std::cout << "Crossover: " << cfg.crossover_type << "\n";
    std::cout << "Mutation: " << cfg.mutation_type << "\n";

    // Initialize population
    ea::Population<> pop(cfg.pop_size, cfg.dim, cfg.n_obj);
    pop.lower_bounds = std::vector<double>(cfg.dim, 0.0);
    pop.upper_bounds = std::vector<double>(cfg.dim, 1.0);

    auto& rng = ea::Random::instance();
    for (int i = 0; i < cfg.pop_size; ++i) {
        for (int j = 0; j < cfg.dim; ++j) {
            pop.gene(i, j) = rng.uniform(0.0, 1.0);
        }
    }

    // Create and run algorithm based on config
    // Note: In production, use a factory pattern or switch statement
    // This example shows the Config struct usage
    ea::ZDT1 problem(cfg.dim);
    auto problem_fn = [&problem](ea::Population<>& p, int idx) {
        problem.evaluate(p, idx);
    };

    if (cfg.algorithm == "NSGAII") {
        ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> alg;
        alg.pop_size = cfg.pop_size;
        alg.max_evals = cfg.max_evals;
        alg.crossover.crossover_probability = cfg.cx_probability;
        alg.crossover.distribution_index = cfg.cx_distribution_index;
        alg.mutation.distribution_index = cfg.mt_distribution_index;
        alg.mutation.mutation_rate = cfg.mt_rate;
        alg.run(pop, problem_fn);
    } else if (cfg.algorithm == "IBEA") {
        ea::IBEA<ea::SBXCrossover, ea::PolynomialMutation> alg;
        alg.pop_size = cfg.pop_size;
        alg.max_evals = cfg.max_evals;
        alg.run(pop, problem_fn);
    } else if (cfg.algorithm == "RandomSearch") {
        ea::RandomSearch alg;
        alg.pop_size = cfg.pop_size;
        alg.max_evals = cfg.max_evals;
        alg.run(pop, problem_fn);
    } else {
        std::cerr << "Unknown algorithm: " << cfg.algorithm << "\n";
        return 1;
    }

    std::cout << "Done. Final population: " << pop.pop_size << "\n";

    return 0;
}
