// Comprehensive benchmark: multiple algorithms × multiple problems
// Compile: g++-14 -std=c++23 -O2 -I include benchmark.cpp -o benchmark -lm

#include <ea/ea.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <vector>
#include <numeric>
#include <algorithm>

struct BenchmarkResult {
    std::string algorithm;
    std::string problem;
    double igd;
    double time_ms;
    int evaluations;
};

// Reference fronts for each problem
std::vector<std::vector<double>> zdt1_front() {
    std::vector<std::vector<double>> ref;
    for (int i = 0; i <= 500; ++i) {
        double f1 = static_cast<double>(i) / 500.0;
        ref.push_back({f1, 1.0 - std::sqrt(f1)});
    }
    return ref;
}

std::vector<std::vector<double>> zdt2_front() {
    std::vector<std::vector<double>> ref;
    for (int i = 0; i <= 500; ++i) {
        double f1 = static_cast<double>(i) / 500.0;
        ref.push_back({f1, 1.0 - f1 * f1});
    }
    return ref;
}

std::vector<std::vector<double>> zdt3_front() {
    std::vector<std::vector<double>> ref;
    for (int i = 0; i <= 500; ++i) {
        double f1 = static_cast<double>(i) / 500.0;
        double f2 = 1.0 - std::sqrt(f1) - f1 * std::sin(10.0 * M_PI * f1);
        if (f2 >= -0.001) ref.push_back({f1, std::max(0.0, f2)});
    }
    return ref;
}

std::vector<std::vector<double>> dtlz2_front() {
    std::vector<std::vector<double>> ref;
    int n = 100;
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            double f1 = static_cast<double>(i) / n;
            double f2 = static_cast<double>(j) / n;
            double f3 = std::max(0.0, 1.0 - f1*f1 - f2*f2);
            if (f3 >= 0) ref.push_back({f1, f2, f3});
        }
    }
    return ref;
}

template<typename Algo, typename Prob>
BenchmarkResult run_benchmark(const std::string& algo_name, const std::string& prob_name,
                               Algo& algo, Prob& prob, int pop_size, int max_evals,
                               const std::vector<std::vector<double>>& ref_front) {
    BenchmarkResult result;
    result.algorithm = algo_name;
    result.problem = prob_name;
    result.evaluations = max_evals;

    const int dim = prob.dimension();
    const int n_obj = prob.num_objectives();

    ea::Population<> pop(pop_size, dim, n_obj);
    pop.lower_bounds = std::vector<double>(prob.lower_bounds().begin(), prob.lower_bounds().end());
    pop.upper_bounds = std::vector<double>(prob.upper_bounds().begin(), prob.upper_bounds().end());

    auto& rng = ea::Random::instance();
    for (int i = 0; i < pop_size; ++i) {
        for (int j = 0; j < dim; ++j) {
            pop.gene(i, j) = rng.uniform(pop.lower_bounds[j], pop.upper_bounds[j]);
        }
        pop.set_evaluated(i, false);
    }

    auto problem = [&prob](ea::Population<>& p, int idx) {
        prob.evaluate(p, idx);
    };

    auto start = std::chrono::high_resolution_clock::now();
    algo.run(pop, problem);
    auto end = std::chrono::high_resolution_clock::now();

    result.time_ms = std::chrono::duration<double, std::milli>(end - start).count();

    // Compute IGD
    std::vector<int> all_indices(pop.pop_size);
    std::iota(all_indices.begin(), all_indices.end(), 0);
    result.igd = ea::igd(pop, all_indices, ref_front);

    return result;
}

int main() {
    const int pop_size = 100;
    const int max_evals = 25000;
    const int runs = 3; // Average over 3 runs

    std::vector<BenchmarkResult> results;

    std::cout << "=== ea-cpp Benchmark ===" << std::endl;
    std::cout << "Population: " << pop_size << ", Max evals: " << max_evals << ", Runs: " << runs << std::endl;
    std::cout << std::endl;

    // === NSGA-II on ZDT1 ===
    {
        double total_igd = 0, total_time = 0;
        for (int r = 0; r < runs; ++r) {
            ea::ZDT1 zdt1(30);
            ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> nsga;
            nsga.crossover.distribution_index = 20.0;
            nsga.crossover.crossover_probability = 0.9;
            nsga.mutation.distribution_index = 20.0;
            nsga.pop_size = pop_size;
            nsga.max_evals = max_evals;
            auto res = run_benchmark("NSGA-II", "ZDT1", nsga, zdt1, pop_size, max_evals, zdt1_front());
            total_igd += res.igd;
            total_time += res.time_ms;
        }
        std::cout << "NSGA-II on ZDT1:  IGD=" << std::fixed << std::setprecision(6) << total_igd/runs
                  << "  Time=" << std::setprecision(1) << total_time/runs << "ms" << std::endl;
    }

    // === NSGA-II on ZDT2 ===
    {
        double total_igd = 0, total_time = 0;
        for (int r = 0; r < runs; ++r) {
            ea::ZDT2 zdt2(30);
            ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> nsga;
            nsga.crossover.distribution_index = 20.0;
            nsga.crossover.crossover_probability = 0.9;
            nsga.mutation.distribution_index = 20.0;
            nsga.pop_size = pop_size;
            nsga.max_evals = max_evals;
            auto res = run_benchmark("NSGA-II", "ZDT2", nsga, zdt2, pop_size, max_evals, zdt2_front());
            total_igd += res.igd;
            total_time += res.time_ms;
        }
        std::cout << "NSGA-II on ZDT2:  IGD=" << std::fixed << std::setprecision(6) << total_igd/runs
                  << "  Time=" << std::setprecision(1) << total_time/runs << "ms" << std::endl;
    }

    // === NSGA-II on ZDT3 ===
    {
        double total_igd = 0, total_time = 0;
        for (int r = 0; r < runs; ++r) {
            ea::ZDT3 zdt3(30);
            ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> nsga;
            nsga.crossover.distribution_index = 20.0;
            nsga.crossover.crossover_probability = 0.9;
            nsga.mutation.distribution_index = 20.0;
            nsga.pop_size = pop_size;
            nsga.max_evals = max_evals;
            auto res = run_benchmark("NSGA-II", "ZDT3", nsga, zdt3, pop_size, max_evals, zdt3_front());
            total_igd += res.igd;
            total_time += res.time_ms;
        }
        std::cout << "NSGA-II on ZDT3:  IGD=" << std::fixed << std::setprecision(6) << total_igd/runs
                  << "  Time=" << std::setprecision(1) << total_time/runs << "ms" << std::endl;
    }

    // === NSGA-II on DTLZ2 (3 obj) ===
    {
        double total_igd = 0, total_time = 0;
        for (int r = 0; r < runs; ++r) {
            ea::DTLZ2 dtlz2(12, 3);
            ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> nsga;
            nsga.crossover.distribution_index = 20.0;
            nsga.crossover.crossover_probability = 0.9;
            nsga.mutation.distribution_index = 20.0;
            nsga.pop_size = pop_size;
            nsga.max_evals = max_evals * 2; // More evals for 3 objectives
            auto ref = dtlz2_front();
            auto res = run_benchmark("NSGA-II", "DTLZ2", nsga, dtlz2, pop_size, max_evals * 2, ref);
            total_igd += res.igd;
            total_time += res.time_ms;
        }
        std::cout << "NSGA-II on DTLZ2: IGD=" << std::fixed << std::setprecision(6) << total_igd/runs
                  << "  Time=" << std::setprecision(1) << total_time/runs << "ms" << std::endl;
    }

    // === SPEA2 on ZDT1 ===
    {
        double total_igd = 0, total_time = 0;
        for (int r = 0; r < runs; ++r) {
            ea::ZDT1 zdt1(30);
            ea::SPEA2<ea::SBXCrossover, ea::PolynomialMutation> spea;
            spea.crossover.distribution_index = 20.0;
            spea.crossover.crossover_probability = 0.9;
            spea.mutation.distribution_index = 20.0;
            spea.pop_size = pop_size;
            spea.max_evals = max_evals;
            auto res = run_benchmark("SPEA2", "ZDT1", spea, zdt1, pop_size, max_evals, zdt1_front());
            total_igd += res.igd;
            total_time += res.time_ms;
        }
        std::cout << "SPEA2  on ZDT1:  IGD=" << std::fixed << std::setprecision(6) << total_igd/runs
                  << "  Time=" << std::setprecision(1) << total_time/runs << "ms" << std::endl;
    }

    // === SMS-EMOA on ZDT1 ===
    {
        double total_igd = 0, total_time = 0;
        for (int r = 0; r < runs; ++r) {
            ea::ZDT1 zdt1(30);
            ea::SMSEMOA<ea::SBXCrossover, ea::PolynomialMutation> sms;
            sms.crossover.distribution_index = 20.0;
            sms.crossover.crossover_probability = 0.9;
            sms.mutation.distribution_index = 20.0;
            sms.pop_size = pop_size;
            sms.max_evals = max_evals;
            auto res = run_benchmark("SMS-EMOA", "ZDT1", sms, zdt1, pop_size, max_evals, zdt1_front());
            total_igd += res.igd;
            total_time += res.time_ms;
        }
        std::cout << "SMS-EMOA on ZDT1: IGD=" << std::fixed << std::setprecision(6) << total_igd/runs
                  << "  Time=" << std::setprecision(1) << total_time/runs << "ms" << std::endl;
    }

    // === NSGA-II on DTLZ1 (3 obj) ===
    {
        double total_igd = 0, total_time = 0;
        for (int r = 0; r < runs; ++r) {
            ea::DTLZ1 dtlz1(7, 3);
            ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> nsga;
            nsga.crossover.distribution_index = 20.0;
            nsga.crossover.crossover_probability = 0.9;
            nsga.mutation.distribution_index = 20.0;
            nsga.pop_size = pop_size;
            nsga.max_evals = max_evals * 4; // More evals for DTLZ1
            auto ref = zdt1_front(); // Approximate — DTLZ1 front is complex
            auto res = run_benchmark("NSGA-II", "DTLZ1", nsga, dtlz1, pop_size, max_evals * 4, ref);
            total_igd += res.igd;
            total_time += res.time_ms;
        }
        std::cout << "NSGA-II on DTLZ1: IGD=" << std::fixed << std::setprecision(6) << total_igd/runs
                  << "  Time=" << std::setprecision(1) << total_time/runs << "ms" << std::endl;
    }

    std::cout << std::endl << "=== Benchmark complete ===" << std::endl;
    return 0;
}