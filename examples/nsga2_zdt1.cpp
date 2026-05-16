// End-to-end test: NSGA-II on ZDT1
// Compile: g++-14 -std=c++23 -O2 -I include test_nsga2_zdt1.cpp -o test_nsga2_zdt1

#include <ea/ea.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

int main() {
    std::cout << "=== NSGA-II on ZDT1 ===" << std::endl;

    // Problem setup
    const int pop_size = 100;
    const int dim = 30;
    const int max_evals = 25000;
    const int n_obj = 2;

    ea::ZDT1 zdt1(dim);

    // Create population
    ea::Population pop(pop_size, dim, n_obj);
    pop.lower_bounds = std::vector<double>(zdt1.lower_bounds().begin(), zdt1.lower_bounds().end());
    pop.upper_bounds = std::vector<double>(zdt1.upper_bounds().begin(), zdt1.upper_bounds().end());

    // Initialize random population
    auto& rng = ea::Random::instance();
    for (int i = 0; i < pop_size; ++i) {
        for (int j = 0; j < dim; ++j) {
            pop.gene(i, j) = rng.uniform(pop.lower_bounds[j], pop.upper_bounds[j]);
        }
        pop.set_evaluated(i, false);
    }

    // Create NSGA-II
    ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> nsga;
    nsga.crossover.distribution_index = 20.0;
    nsga.crossover.crossover_probability = 0.9;
    nsga.mutation.distribution_index = 20.0;
    nsga.pop_size = pop_size;
    nsga.max_evals = max_evals;

    // Define problem evaluation
    auto problem = [&zdt1](ea::Population& p, int idx) {
        zdt1.evaluate(p, idx);
    };

    // Run
    std::cout << "Running NSGA-II with " << pop_size << " individuals, "
              << dim << " variables, " << max_evals << " max evaluations..." << std::endl;

    nsga.run(pop, problem);

    std::cout << "Done. Final population size: " << pop.pop_size << std::endl;

    // Print first 10 individuals (f1, f2)
    std::cout << "\nFirst 10 individuals (f1, f2):" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    for (int i = 0; i < std::min(10, pop.pop_size); ++i) {
        std::cout << "  [" << i << "] f1=" << pop.objective(i, 0)
                  << " f2=" << pop.objective(i, 1) << std::endl;
    }

    // Compute IGD (approximation — using known Pareto front of ZDT1)
    // ZDT1 Pareto front: f2 = 1 - sqrt(f1), f1 in [0,1]
    std::vector<std::vector<double>> reference_front;
    for (int i = 0; i <= 100; ++i) {
        double f1 = static_cast<double>(i) / 100.0;
        double f2 = 1.0 - std::sqrt(f1);
        reference_front.push_back({f1, f2});
    }

    double igd = ea::compute_igd(pop, reference_front);
    std::cout << "\nIGD: " << igd << std::endl;

    // Write final front to file
    std::ofstream ofs("zdt1_front.csv");
    ofs << "f1,f2\n";
    for (int i = 0; i < pop.pop_size; ++i) {
        ofs << pop.objective(i, 0) << "," << pop.objective(i, 1) << "\n";
    }
    ofs.close();
    std::cout << "Front written to zdt1_front.csv" << std::endl;

    // Success criteria: IGD should be < 0.1 for a reasonable run
    bool success = igd < 0.1;
    std::cout << "\nTest " << (success ? "PASSED" : "FAILED") << " (IGD=" << igd << ", threshold=0.1)" << std::endl;

    return success ? 0 : 1;
}