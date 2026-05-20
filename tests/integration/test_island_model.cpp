// Test: Island Model with 2 islands on ZDT1 — verify convergence
// Compile: g++-14 -std=c++23 -O2 -I include test_island_model.cpp -o test_island_model

#include <ea/ea.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

int main() {
    std::cout << "=== Test: Island Model (2 islands, NSGA-II, ZDT1) ===" << std::endl;

    // Test parameters
    const int num_islands = 2;
    const int pop_per_island = 50;   // 2 x 50 = 100 total population
    const int dim = 30;
    const int max_evals = 25000;
    const int n_obj = 2;

    ea::ZDT1 zdt1(dim);

    // Algorithm template
    ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> nsga;
    nsga.crossover.distribution_index = 20.0;
    nsga.crossover.crossover_probability = 0.9;
    nsga.mutation.distribution_index = 20.0;
    nsga.pop_size = pop_per_island;

    // Topology: ring of 2 islands
    ea::RingTopology topo(num_islands);

    // Migration policy
    ea::MigrationPolicy policy;
    policy.migration_rate = 10;
    policy.num_migrants = 5;
    policy.migrant_selection = ea::MigrantSelection::Best;

    // Create island model
    ea::IslandModel model(nsga, topo, policy, num_islands, pop_per_island, dim, n_obj);

    // Evaluation function
    auto problem = [&zdt1](ea::Population& p, int idx) {
        zdt1.evaluate(p, idx);
    };

    // Run
    model.run(problem, max_evals);

    // Get global front
    auto front = model.best_front();
    std::cout << "Global Pareto front size: " << front.pop_size << std::endl;

    // Compute IGD
    std::vector<std::vector<double>> reference_front;
    for (int i = 0; i <= 100; ++i) {
        double f1 = static_cast<double>(i) / 100.0;
        double f2 = 1.0 - std::sqrt(f1);
        reference_front.push_back({f1, f2});
    }

    std::vector<int> all_indices(front.pop_size);
    std::iota(all_indices.begin(), all_indices.end(), 0);
    double igd_val = ea::igd(front, all_indices, reference_front);
    std::cout << "IGD: " << std::fixed << std::setprecision(6) << igd_val << std::endl;

    // Test 1: IGD should be < 0.1 (convergence)
    bool test1 = igd_val < 0.1;
    std::cout << "[Test 1] Convergence (IGD < 0.1): " << (test1 ? "PASS" : "FAIL")
              << " (IGD=" << igd_val << ")" << std::endl;

    // Test 2: Front should have reasonable size
    bool test2 = front.pop_size >= 10 && front.pop_size <= 100;
    std::cout << "[Test 2] Front size (10 <= " << front.pop_size << " <= 100): "
              << (test2 ? "PASS" : "FAIL") << std::endl;

    // Test 3: Front should span the objective space
    // Check that f1 ranges roughly [0, 1]
    double f1_min = 1e9, f1_max = -1e9;
    for (int i = 0; i < front.pop_size; ++i) {
        f1_min = std::min(f1_min, front.objective(i, 0));
        f1_max = std::max(f1_max, front.objective(i, 0));
    }
    bool test3 = f1_min < 0.05 && f1_max > 0.5;
    std::cout << "[Test 3] Front spread (f1_min=" << f1_min
              << ", f1_max=" << f1_max << "): "
              << (test3 ? "PASS" : "FAIL") << std::endl;

    // Test 4: Test MeshTopology constructor
    ea::MeshTopology mesh(2, 2);
    auto mesh_neighbors = mesh.neighbors(0);
    bool test4 = mesh_neighbors.size() == 4; // North, South, East, West
    std::cout << "[Test 4] MeshTopology neighbors count=" << mesh_neighbors.size()
              << ": " << (test4 ? "PASS" : "FAIL") << std::endl;

    // Test 5: Test FullyConnectedTopology
    ea::FullyConnectedTopology fct(4);
    auto fct_neighbors = fct.neighbors(0);
    bool test5 = fct_neighbors.size() == 3; // 3 other islands
    std::cout << "[Test 5] FullyConnected neighbors count=" << fct_neighbors.size()
              << ": " << (test5 ? "PASS" : "FAIL") << std::endl;

    bool all_passed = test1 && test2 && test3 && test4 && test5;
    std::cout << "\nOverall: " << (all_passed ? "ALL TESTS PASSED" : "SOME TESTS FAILED")
              << std::endl;

    return all_passed ? 0 : 1;
}
