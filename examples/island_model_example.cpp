// End-to-end test: Island Model with NSGA-II on ZDT1
// Compile: g++-14 -std=c++23 -O2 -I include island_model_example.cpp -o island_model_example

#include <ea/ea.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

int main() {
    std::cout << "=== Island Model (NSGA-II) on ZDT1 ===" << std::endl;

    // Problem setup
    const int num_islands = 4;
    const int pop_per_island = 26;   // 4 islands x 26 = 104 total (must be even)
    const int dim = 30;
    const int max_evals = 25000;
    const int n_obj = 2;

    ea::ZDT1 zdt1(dim);

    // Create algorithm template
    ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> nsga;
    nsga.crossover.distribution_index = 20.0;
    nsga.crossover.crossover_probability = 0.9;
    nsga.mutation.distribution_index = 20.0;
    nsga.pop_size = pop_per_island;

    // Create topology: ring of 4 islands
    ea::RingTopology topo(num_islands);

    // Migration policy: every 10 generations, send 3 best individuals
    ea::MigrationPolicy policy;
    policy.migration_rate = 10;
    policy.num_migrants = 3;
    policy.migrant_selection = ea::MigrantSelection::Best;

    // Create island model
    ea::IslandModel model(nsga, topo, policy, num_islands, pop_per_island, dim, n_obj);

    // Define problem evaluation
    auto problem = [&zdt1](ea::Population<>& p, int idx) {
        zdt1.evaluate(p, idx);
    };

    // Run
    std::cout << "Running Island Model with " << num_islands << " islands, "
              << pop_per_island << " individuals each, "
              << dim << " variables, " << max_evals << " max evaluations..." << std::endl;

    model.run(problem, max_evals);

    // Extract global best front
    auto front = model.best_front();
    std::cout << "Done. Global Pareto front size: " << front.pop_size << std::endl;

    // Print first 10 individuals (f1, f2)
    std::cout << "\nFirst 10 front individuals (f1, f2):" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    for (int i = 0; i < std::min(10, front.pop_size); ++i) {
        std::cout << "  [" << i << "] f1=" << front.objective(i, 0)
                  << " f2=" << front.objective(i, 1) << std::endl;
    }

    // Compute IGD (approximation — using known Pareto front of ZDT1)
    std::vector<std::vector<double>> reference_front;
    for (int i = 0; i <= 100; ++i) {
        double f1 = static_cast<double>(i) / 100.0;
        double f2 = 1.0 - std::sqrt(f1);
        reference_front.push_back({f1, f2});
    }

    std::vector<int> all_indices(front.pop_size);
    std::iota(all_indices.begin(), all_indices.end(), 0);
    double igd_val = ea::igd(front, all_indices, reference_front);
    std::cout << "\nIGD: " << igd_val << std::endl;

    // Write final front to file
    std::ofstream ofs("island_zdt1_front.csv");
    ofs << "f1,f2\n";
    for (int i = 0; i < front.pop_size; ++i) {
        ofs << front.objective(i, 0) << "," << front.objective(i, 1) << "\n";
    }
    ofs.close();
    std::cout << "Front written to island_zdt1_front.csv" << std::endl;

    // Success criteria: IGD should be < 0.1 for a reasonable run
    bool success = igd_val < 0.1;
    std::cout << "\nTest " << (success ? "PASSED" : "FAILED")
              << " (IGD=" << igd_val << ", threshold=0.1)" << std::endl;

    return success ? 0 : 1;
}
