// Extended Benchmark Suite — algorithms × problems × N runs
// Compile: g++-14 -std=c++23 -O2 -I include benchmark_extended.cpp -o benchmark_extended -lm
// Usage: ./benchmark_extended [runs] > results.csv

#include <ea/ea.hpp>
#include <ea/algorithm/rvea.hpp>
#include <ea/algorithm/ibea.hpp>
#include <ea/algorithm/random_search.hpp>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

// === Reference fronts ===

std::vector<std::vector<double>> zdt1_ref() {
    std::vector<std::vector<double>> ref;
    for (int i = 0; i <= 500; ++i) {
        double f1 = static_cast<double>(i) / 500.0;
        ref.push_back({f1, 1.0 - std::sqrt(f1)});
    }
    return ref;
}

std::vector<std::vector<double>> zdt2_ref() {
    std::vector<std::vector<double>> ref;
    for (int i = 0; i <= 500; ++i) {
        double f1 = static_cast<double>(i) / 500.0;
        ref.push_back({f1, 1.0 - f1 * f1});
    }
    return ref;
}

std::vector<std::vector<double>> zdt3_ref() {
    std::vector<std::vector<double>> ref;
    for (int i = 0; i <= 10000; ++i) {
        double f1 = static_cast<double>(i) / 10000.0;
        double h = 1.0 - std::sqrt(f1) - f1 * std::sin(10.0 * M_PI * f1);
        if (h >= 0.0 && h <= 1.2)
            ref.push_back({f1, h});
    }
    return ref;
}

std::vector<std::vector<double>> zdt6_ref() {
    // ZDT6 Pareto front: g=1 on front, so f2 = 1 - f1^2
    std::vector<std::vector<double>> ref;
    for (int i = 0; i <= 500; ++i) {
        double t = static_cast<double>(i) / 500.0;
        double f1 = 1.0 - std::exp(-4.0 * t) * std::pow(std::sin(6.0 * M_PI * t), 6.0);
        double f2 = 1.0 - f1 * f1;
        if (f1 >= 0.0 && f1 <= 1.0 && f2 >= 0.0)
            ref.push_back({f1, f2});
    }
    return ref;
}

// === IGD ===
double compute_igd(const ea::Population<>& pop, const std::vector<std::vector<double>>& ref_front) {
    if (ref_front.empty() || pop.pop_size == 0) return 0.0;
    int n_obj = pop.n_obj;
    double sum = 0.0;
    for (const auto& rp : ref_front) {
        double min_d = std::numeric_limits<double>::max();
        for (int i = 0; i < pop.pop_size; ++i) {
            double d = 0.0;
            for (int o = 0; o < n_obj; ++o) d += std::pow(pop.objective(i, o) - rp[o], 2);
            min_d = std::min(min_d, std::sqrt(d));
        }
        sum += min_d;
    }
    return sum / static_cast<double>(ref_front.size());
}

// === Problem config ===
struct ProblemConfig {
    std::string name;
    int dim, n_obj, pop_size, max_evals;
    std::vector<std::vector<double>> ref_front;
    std::function<void(ea::Population<>&, int)> evaluate;
    std::vector<double> lower_bounds, upper_bounds;
};

std::vector<ProblemConfig> get_problems() {
    std::vector<ProblemConfig> problems;
    auto add = [&](auto& p, std::string name, int pop, int evals, std::vector<std::vector<double>> ref) {
        problems.push_back({name, p.dim_, 2, pop, evals, std::move(ref),
            [p](ea::Population<>& pop, int i) mutable { p.evaluate(pop, i); },
            std::vector<double>(p.lower_bounds().begin(), p.lower_bounds().end()),
            std::vector<double>(p.upper_bounds().begin(), p.upper_bounds().end())});
    };
    { ea::ZDT1 p(30); add(p, "ZDT1", 100, 25000, zdt1_ref()); }
    { ea::ZDT2 p(30); add(p, "ZDT2", 100, 25000, zdt2_ref()); }
    { ea::ZDT3 p(30); add(p, "ZDT3", 100, 25000, zdt3_ref()); }
    { ea::ZDT6 p(10); add(p, "ZDT6", 100, 25000, zdt6_ref()); }
    return problems;
}

ea::Population<> init_pop(const ProblemConfig& pc) {
    ea::Population<> pop(pc.pop_size, pc.dim, pc.n_obj);
    pop.lower_bounds = pc.lower_bounds;
    pop.upper_bounds = pc.upper_bounds;
    auto& rng = ea::Random::instance();
    for (int i = 0; i < pc.pop_size; ++i)
        for (int j = 0; j < pc.dim; ++j)
            pop.gene(i, j) = rng.uniform(pc.lower_bounds[j], pc.upper_bounds[j]);
    return pop;
}

struct RunResult { double igd; double time_ms; };

// === Runners ===
RunResult run_nsga2(const ProblemConfig& pc) {
    ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> algo;
    algo.pop_size = pc.pop_size; algo.max_evals = pc.max_evals;
    algo.crossover.crossover_probability = 0.9;
    algo.crossover.distribution_index = 20.0;
    algo.mutation.distribution_index = 20.0;
    auto pop = init_pop(pc);
    auto t0 = std::chrono::high_resolution_clock::now();
    algo.run(pop, pc.evaluate);
    return {compute_igd(pop, pc.ref_front), std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count()};
}

RunResult run_spea2(const ProblemConfig& pc) {
    ea::SPEA2<ea::SBXCrossover, ea::PolynomialMutation> algo;
    algo.pop_size = pc.pop_size; algo.max_evals = pc.max_evals;
    algo.crossover.crossover_probability = 0.9;
    algo.crossover.distribution_index = 20.0;
    algo.mutation.distribution_index = 20.0;
    auto pop = init_pop(pc);
    auto t0 = std::chrono::high_resolution_clock::now();
    algo.run(pop, pc.evaluate);
    return {compute_igd(pop, pc.ref_front), std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count()};
}

RunResult run_smpso(const ProblemConfig& pc) {
    ea::SMPSO algo;
    algo.pop_size = pc.pop_size; algo.max_evals = pc.max_evals;
    auto pop = init_pop(pc);
    auto t0 = std::chrono::high_resolution_clock::now();
    algo.run(pop, pc.evaluate);
    return {compute_igd(pop, pc.ref_front), std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count()};
}

RunResult run_agemoea(const ProblemConfig& pc) {
    ea::AGEMOEA<ea::SBXCrossover, ea::PolynomialMutation> algo;
    algo.pop_size = pc.pop_size; algo.max_evals = pc.max_evals;
    algo.crossover.crossover_probability = 0.9;
    algo.crossover.distribution_index = 20.0;
    algo.mutation.distribution_index = 20.0;
    auto pop = init_pop(pc);
    auto t0 = std::chrono::high_resolution_clock::now();
    algo.run(pop, pc.evaluate);
    return {compute_igd(pop, pc.ref_front), std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count()};
}

RunResult run_paes(const ProblemConfig& pc) {
    ea::PAES<ea::PolynomialMutation> algo;
    algo.max_evals = pc.max_evals; algo.archive_size = pc.pop_size;
    auto pop = init_pop(pc);
    auto t0 = std::chrono::high_resolution_clock::now();
    algo.run(pop, pc.evaluate);
    return {compute_igd(pop, pc.ref_front), std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count()};
}

RunResult run_mocell(const ProblemConfig& pc) {
    ea::MOCell<ea::SBXCrossover, ea::PolynomialMutation> algo;
    algo.pop_size = pc.pop_size; algo.max_evals = pc.max_evals;
    algo.crossover.crossover_probability = 0.9;
    algo.crossover.distribution_index = 20.0;
    algo.mutation.distribution_index = 20.0;
    auto pop = init_pop(pc);
    auto t0 = std::chrono::high_resolution_clock::now();
    algo.run(pop, pc.evaluate);
    return {compute_igd(pop, pc.ref_front), std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count()};
}

RunResult run_random_search(const ProblemConfig& pc) {
    ea::RandomSearch algo;
    algo.pop_size = pc.pop_size; algo.max_evals = pc.max_evals;
    auto pop = init_pop(pc);
    auto t0 = std::chrono::high_resolution_clock::now();
    algo.run(pop, pc.evaluate);
    return {compute_igd(pop, pc.ref_front), std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count()};
}

// === Registry ===
struct AlgorithmEntry {
    std::string name;
    std::function<RunResult(const ProblemConfig&)> runner;
};

std::vector<AlgorithmEntry> get_algorithms() {
    return {
        {"NSGA-II",      run_nsga2},
        {"SPEA2",        run_spea2},
        {"SMPSO",        run_smpso},
        {"AGE-MOEA",     run_agemoea},
        {"PAES",         run_paes},
        {"MOCell",       run_mocell},
        {"RandomSearch", run_random_search},
    };
}

// === Stats ===
struct Stats { double mean, stddev, median, min_val, max_val; };
Stats compute_stats(std::vector<double>& v) {
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / n;
    double var = 0.0;
    for (auto x : v) var += (x - mean) * (x - mean);
    return {mean, n > 1 ? std::sqrt(var / (n - 1)) : 0.0,
            n % 2 == 0 ? (v[n/2-1] + v[n/2]) / 2.0 : v[n/2],
            v.front(), v.back()};
}

int main(int argc, char** argv) {
    int runs = 10;
    if (argc > 1) runs = std::atoi(argv[1]);
    if (runs < 1) runs = 10;
    auto problems = get_problems();
    auto algorithms = get_algorithms();
    std::cout << "algorithm,problem,runs,igd_mean,igd_std,igd_median,igd_min,igd_max,time_mean_ms,time_std_ms" << std::endl;
    for (const auto& algo : algorithms) {
        for (const auto& prob : problems) {
            std::vector<double> igds, times;
            for (int r = 0; r < runs; ++r) {
                auto res = algo.runner(prob);
                igds.push_back(res.igd);
                times.push_back(res.time_ms);
            }
            auto is = compute_stats(igds);
            auto ts = compute_stats(times);
            std::cout << std::fixed << std::setprecision(6);
            std::cout << algo.name << "," << prob.name << "," << runs << ","
                      << is.mean << "," << is.stddev << "," << is.median << "," << is.min_val << "," << is.max_val << ","
                      << ts.mean << "," << ts.stddev << std::endl;
            std::cerr << "[" << algo.name << " x " << prob.name << "] IGD=" << is.mean << " +/- " << is.stddev
                      << "  Time=" << ts.mean << " +/- " << ts.stddev << "ms" << std::endl;
        }
    }
    return 0;
}