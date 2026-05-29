// Test: RandomSearch on ZDT1
// Compile: g++-14 -std=c++23 -O2 -I include examples/test_random_search.cpp -o test_random_search -lm

#include <ea/ea.hpp>
#include <ea/algorithm/random_search.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    std::cout << "=== Random Search on ZDT1 ===" << std::endl;

    // Problem setup
    const int pop_size = 100;
    const int dim = 30;
    const int max_evals = 5000;  // Small number for a quick baseline run
    const int n_obj = 2;

    ea::ZDT1 zdt1(dim);

    // Create population
    ea::Population<> pop(pop_size, dim, n_obj);
    pop.lower_bounds = std::vector<double>(zdt1.lower_bounds().begin(), zdt1.lower_bounds().end());
    pop.upper_bounds = std::vector<double>(zdt1.upper_bounds().begin(), zdt1.upper_bounds().end());

    // Create Random Search
    ea::RandomSearch rs;
    rs.pop_size = pop_size;
    rs.max_evals = max_evals;

    // Define problem evaluation
    auto problem = [&zdt1](ea::Population<>& p, int idx) {
        zdt1.evaluate(p, idx);
    };

    // Run
    std::cout << "Running Random Search with " << pop_size << " individuals per batch, "
              << max_evals << " max evaluations..." << std::endl;

    rs.run(pop, problem);

    std::cout << "Done. Final non-dominated front size: " << pop.pop_size << std::endl;

    // Print first 5 individuals (f1, f2)
    std::cout << "\nFirst 5 individuals (f1, f2):" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    for (int i = 0; i < std::min(5, pop.pop_size); ++i) {
        std::cout << "  [" << i << "] f1=" << pop.objective(i, 0)
                  << " f2=" << pop.objective(i, 1) << std::endl;
    }

    std::cout << "\nAlgorithm: " << ea::RandomSearch::name() << std::endl;

    return 0;
}