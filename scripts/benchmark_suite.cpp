// Comprehensive benchmark suite: all algorithms × all problems × 30 runs
// Compile: g++-14 -std=c++23 -O2 -I include scripts/benchmark_suite.cpp -o benchmark_suite -lm
// Usage: ./benchmark_suite --output results.csv

#include <ea/ea.hpp>
#include <ea/algorithm/nsga2.hpp>
#include <ea/algorithm/rvea.hpp>
#include <ea/algorithm/ibea.hpp>
#include <ea/algorithm/random_search.hpp>
#include <ea/indicator/igd.hpp>
#include <ea/indicator/gd.hpp>
#include <ea/indicator/spread.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <vector>
#include <string>
#include <memory>

using namespace ea;

struct BenchmarkConfig {
    int runs = 30;
    int pop_size = 100;
    int max_evals = 25000;
    int dim = 30;
};

// Problem factory
template<typename Problem>
std::vector<std::vector<double>> generate_reference_front(int dim) {
    std::vector<std::vector<double>> ref;
    for (int i = 0; i <= 1000; ++i) {
        double f1 = static_cast<double>(i) / 1000.0;
        if constexpr (std::is_same_v<Problem, ZDT1>) {
            ref.push_back({f1, 1.0 - std::sqrt(f1)});
        } else if constexpr (std::is_same_v<Problem, ZDT2>) {
            ref.push_back({f1, 1.0 - f1 * f1});
        } else if constexpr (std::is_same_v<Problem, ZDT3>) {
            ref.push_back({f1, 1.0 - std::sqrt(f1) - f1 * std::sin(10.0 * M_PI * f1)});
        } else if constexpr (std::is_same_v<Problem, ZDT6>) {
            ref.push_back({f1, 1.0 - std::exp(-4.0 * f1) * std::pow(std::sin(6.0 * M_PI * f1), 6)});
        }
    }
    return ref;
}

// Run single algorithm on single problem
template<typename Algorithm, typename Problem>
void run_benchmark(const BenchmarkConfig& cfg,
                    const std::string& alg_name,
                    const std::string& prob_name,
                    std::ofstream& out) {
    std::vector<double> igds, gds, spreads, times;

    for (int r = 0; r < cfg.runs; ++r) {
        Problem problem(cfg.dim);
        Population<> pop(cfg.pop_size, cfg.dim, 2);
        pop.lower_bounds = std::vector<double>(cfg.dim, 0.0);
        pop.upper_bounds = std::vector<double>(cfg.dim, 1.0);

        auto& rng = Random::instance();
        for (int i = 0; i < cfg.pop_size; ++i) {
            for (int j = 0; j < cfg.dim; ++j) {
                pop.gene(i, j) = rng.uniform(0.0, 1.0);
            }
        }

        auto problem_fn = [&problem](Population<>& p, int idx) {
            problem.evaluate(p, idx);
        };

        Algorithm alg;
        alg.pop_size = cfg.pop_size;
        alg.max_evals = cfg.max_evals;

        auto start = std::chrono::high_resolution_clock::now();
        alg.run(pop, problem_fn);
        auto end = std::chrono::high_resolution_clock::now();

        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();

        auto ref_front = generate_reference_front<Problem>(cfg.dim);
        std::vector<int> indices;
        for (int i = 0; i < cfg.pop_size; ++i) indices.push_back(i);

        igds.push_back(igd(pop, indices, ref_front));
        gds.push_back(gd(pop, indices, ref_front));
        spreads.push_back(spread(pop, indices, ref_front));
        times.push_back(time_ms);
    }

    auto mean = [](const std::vector<double>& v) {
        double sum = 0;
        for (double x : v) sum += x;
        return sum / v.size();
    };

    auto stddev = [](const std::vector<double>& v, double m) {
        double sum = 0;
        for (double x : v) sum += (x - m) * (x - m);
        return std::sqrt(sum / v.size());
    };

    double mean_igd = mean(igds);
    double mean_gd = mean(gds);
    double mean_spread = mean(spreads);
    double mean_time = mean(times);
    double std_igd = stddev(igds, mean_igd);

    out << alg_name << "," << prob_name << ","
        << std::fixed << std::setprecision(6)
        << mean_igd << "," << std_igd << ","
        << mean_gd << "," << mean_spread << ","
        << std::setprecision(2) << mean_time << "\n";

    std::cout << "  " << alg_name << " \u00d7 " << prob_name
              << ": IGD=" << mean_igd << " Time=" << mean_time << "ms\n";
}

int main(int argc, char** argv) {
    BenchmarkConfig cfg;
    std::string output_file = "benchmark_results.csv";

    for (int i = 1; i < argc; i += 2) {
        if (i + 1 >= argc) break;
        std::string key = argv[i];
        std::string val = argv[i + 1];
        if (key == "--runs") cfg.runs = std::stoi(val);
        else if (key == "--pop") cfg.pop_size = std::stoi(val);
        else if (key == "--evals") cfg.max_evals = std::stoi(val);
        else if (key == "--output") output_file = val;
    }

    std::cout << "=== ea-cpp Comprehensive Benchmark Suite ===\n";
    std::cout << "Runs: " << cfg.runs << ", Pop: " << cfg.pop_size
              << ", Evals: " << cfg.max_evals << "\n\n";

    std::ofstream out(output_file);
    out << "algorithm,problem,IGD_mean,IGD_std,GD_mean,Spread_mean,Time_ms\n";

    run_benchmark<NSGAII<SBXCrossover, PolynomialMutation>, ZDT1>(
        cfg, "NSGAII", "ZDT1", out);
    run_benchmark<NSGAII<SBXCrossover, PolynomialMutation>, ZDT2>(
        cfg, "NSGAII", "ZDT2", out);
    run_benchmark<NSGAII<SBXCrossover, PolynomialMutation>, ZDT3>(
        cfg, "NSGAII", "ZDT3", out);
    run_benchmark<NSGAII<SBXCrossover, PolynomialMutation>, ZDT6>(
        cfg, "NSGAII", "ZDT6", out);

    run_benchmark<RVEA<SBXCrossover, PolynomialMutation>, ZDT1>(
        cfg, "RVEA", "ZDT1", out);
    run_benchmark<IBEA<SBXCrossover, PolynomialMutation>, ZDT1>(
        cfg, "IBEA", "ZDT1", out);
    run_benchmark<RandomSearch, ZDT1>(
        cfg, "RandomSearch", "ZDT1", out);

    std::cout << "\nResults saved to: " << output_file << "\n";
    return 0;
}
