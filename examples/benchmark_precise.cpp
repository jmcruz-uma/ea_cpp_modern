// Precise benchmark: NSGA-II on ZDT1 — measures evaluations, IGD, and wall-clock time
// Compile: g++-14 -std=c++23 -O2 -I include benchmark_precise.cpp -o benchmark_precise -lm

#include <ea/ea.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <vector>
#include <numeric>
#include <algorithm>

// ZDT1 reference front
std::vector<std::vector<double>> zdt1_reference_front() {
    std::vector<std::vector<double>> ref;
    for (int i = 0; i <= 1000; ++i) {
        double f1 = static_cast<double>(i) / 1000.0;
        ref.push_back({f1, 1.0 - std::sqrt(f1)});
    }
    return ref;
}

// Compute IGD manually for precise measurement
double compute_igd(const ea::Population& pop, const std::vector<std::vector<double>>& ref_front) {
    if (ref_front.empty() || pop.pop_size == 0) return 0.0;
    int n_obj = pop.n_obj;
    double sum = 0.0;
    for (const auto& ref_point : ref_front) {
        double min_dist = std::numeric_limits<double>::max();
        for (int i = 0; i < pop.pop_size; ++i) {
            double dist = 0.0;
            for (int o = 0; o < n_obj; ++o) {
                double diff = pop.objective(i, o) - ref_point[o];
                dist += diff * diff;
            }
            min_dist = std::min(min_dist, std::sqrt(dist));
        }
        sum += min_dist;
    }
    return sum / static_cast<double>(ref_front.size());
}

int main() {
    const int POP_SIZE = 100;
    const int MAX_EVALS = 25000;
    const int DIM = 30;
    const int N_OBJ = 2;
    const int RUNS = 10;

    std::cout << "=== ea-cpp vs jMetal Benchmark: NSGA-II on ZDT1 ===" << std::endl;
    std::cout << "Population: " << POP_SIZE << std::endl;
    std::cout << "Max evaluations: " << MAX_EVALS << std::endl;
    std::cout << "Dimensions: " << DIM << std::endl;
    std::cout << "Runs: " << RUNS << std::endl;
    std::cout << std::endl;

    auto ref_front = zdt1_reference_front();

    std::vector<double> igds;
    std::vector<double> times;

    for (int r = 0; r < RUNS; ++r) {
        ea::ZDT1 zdt1(DIM);
        ea::Population pop(POP_SIZE, DIM, N_OBJ);
        pop.lower_bounds = std::vector<double>(zdt1.lower_bounds().begin(), zdt1.lower_bounds().end());
        pop.upper_bounds = std::vector<double>(zdt1.upper_bounds().begin(), zdt1.upper_bounds().end());

        auto& rng = ea::Random::instance();
        for (int i = 0; i < POP_SIZE; ++i) {
            for (int j = 0; j < DIM; ++j) {
                pop.gene(i, j) = rng.uniform(pop.lower_bounds[j], pop.upper_bounds[j]);
            }
            pop.set_evaluated(i, false);
        }

        ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> nsga;
        nsga.crossover.distribution_index = 20.0;
        nsga.crossover.crossover_probability = 0.9;
        nsga.mutation.distribution_index = 20.0;
        nsga.pop_size = POP_SIZE;
        nsga.max_evals = MAX_EVALS;

        auto problem = [&zdt1](ea::Population& p, int idx) {
            zdt1.evaluate(p, idx);
        };

        auto start = std::chrono::high_resolution_clock::now();
        nsga.run(pop, problem);
        auto end = std::chrono::high_resolution_clock::now();

        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        double igd = compute_igd(pop, ref_front);

        igds.push_back(igd);
        times.push_back(time_ms);

        std::cout << "  Run " << (r+1) << ": IGD=" << std::fixed << std::setprecision(6)
                  << igd << "  Time=" << std::setprecision(2) << time_ms << "ms" << std::endl;
    }

    // Statistics
    double mean_igd = std::accumulate(igds.begin(), igds.end(), 0.0) / RUNS;
    double mean_time = std::accumulate(times.begin(), times.end(), 0.0) / RUNS;
    double std_igd = 0, std_time = 0;
    for (int i = 0; i < RUNS; ++i) {
        std_igd += (igds[i] - mean_igd) * (igds[i] - mean_igd);
        std_time += (times[i] - mean_time) * (times[i] - mean_time);
    }
    std_igd = std::sqrt(std_igd / RUNS);
    std_time = std::sqrt(std_time / RUNS);

    std::cout << std::endl;
    std::cout << "=== RESULTS (ea-cpp C++23) ===" << std::endl;
    std::cout << "Mean IGD:  " << std::fixed << std::setprecision(6) << mean_igd
              << " ± " << std::setprecision(6) << std_igd << std::endl;
    std::cout << "Mean Time: " << std::fixed << std::setprecision(2) << mean_time
              << " ± " << std::setprecision(2) << std_time << "ms" << std::endl;
    std::cout << "Best IGD:  " << std::setprecision(6) << *std::min_element(igds.begin(), igds.end()) << std::endl;
    std::cout << "Worst IGD: " << std::setprecision(6) << *std::max_element(igds.begin(), igds.end()) << std::endl;
    std::cout << "Best Time: " << std::setprecision(2) << *std::min_element(times.begin(), times.end()) << "ms" << std::endl;
    std::cout << "Worst Time:" << std::setprecision(2) << *std::max_element(times.begin(), times.end()) << "ms" << std::endl;

    return 0;
}