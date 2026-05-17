/// @file bench_v1_vs_current.cpp
/// @brief Benchmark: current ea-cpp NSGA-II on ZDT1-6, DTLZ1-2
///
/// Measures wall time, CPU time, and IGD quality.
/// Run multiple times for statistical significance.
///
/// Usage:
///   clang++-18 -std=c++23 -O3 -march=native -I../include bench_v1_vs_current.cpp -o bench
///   ./bench [--runs 30] [--csv] [--verbose]

#include "bench_timer.hpp"

#include <ea/ea.hpp>
#include <cmath>
#include <numeric>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstring>

// ============================================================
// Known Pareto-optimal fronts for IGD calculation
// ============================================================

namespace ea::bench {

std::vector<std::pair<double,double>> zdt1_front(int n = 500) {
    std::vector<std::pair<double,double>> pf(n);
    for (int i = 0; i < n; ++i) {
        double f1 = static_cast<double>(i) / (n - 1);
        pf[i] = {f1, 1.0 - std::sqrt(f1)};
    }
    return pf;
}

std::vector<std::pair<double,double>> zdt2_front(int n = 500) {
    std::vector<std::pair<double,double>> pf(n);
    for (int i = 0; i < n; ++i) {
        double f1 = static_cast<double>(i) / (n - 1);
        pf[i] = {f1, 1.0 - f1 * f1};
    }
    return pf;
}

std::vector<std::pair<double,double>> zdt3_front() {
    std::vector<std::pair<double,double>> pf;
    for (double f1 = 0; f1 <= 0.0830; f1 += 0.0830/50)
        pf.push_back({f1, 1 - std::sqrt(f1) - f1 * std::sin(10*M_PI*f1)});
    for (double f1 = 0.1822; f1 <= 0.2578; f1 += 0.0756/50)
        pf.push_back({f1, 1 - std::sqrt(f1) - f1 * std::sin(10*M_PI*f1)});
    for (double f1 = 0.4099; f1 <= 0.4539; f1 += 0.0440/50)
        pf.push_back({f1, 1 - std::sqrt(f1) - f1 * std::sin(10*M_PI*f1)});
    for (double f1 = 0.6184; f1 <= 0.6525; f1 += 0.0341/50)
        pf.push_back({f1, 1 - std::sqrt(f1) - f1 * std::sin(10*M_PI*f1)});
    for (double f1 = 0.8233; f1 <= 0.8518; f1 += 0.0285/50)
        pf.push_back({f1, 1 - std::sqrt(f1) - f1 * std::sin(10*M_PI*f1)});
    return pf;
}

double compute_igd(const std::vector<std::pair<double,double>>& obtained,
                   const std::vector<std::pair<double,double>>& reference) {
    double sum = 0.0;
    for (auto& [f1, f2] : reference) {
        double min_dist = 1e30;
        for (auto& [g1, g2] : obtained) {
            double d = std::pow(f1 - g1, 2) + std::pow(f2 - g2, 2);
            min_dist = std::min(min_dist, d);
        }
        sum += std::sqrt(min_dist);
    }
    return sum / reference.size();
}

struct BenchConfig {
    int pop_size = 100;
    int max_evals = 25000;
    int runs = 30;
    int dim = 30;
};

// ============================================================
// Generic runner — templated on Problem type for zero-overhead
// ============================================================

template<typename Problem>
TimingResult run_nsga2_zdt_typed(const std::string& label,
                                  const std::vector<std::pair<double,double>>& ref_front,
                                  Problem& prob,
                                  const BenchConfig& cfg) {
    TimingResult result;
    result.label = label;
    result.population_size = cfg.pop_size;
    result.evaluations = cfg.max_evals;
    result.dimensions = cfg.dim;
    result.n_objectives = 2;
    result.runs = cfg.runs;

    std::vector<double> igd_values(cfg.runs);
    double total_wall = 0, total_cpu = 0;

    for (int r = 0; r < cfg.runs; ++r) {
        Population pop(cfg.pop_size, cfg.dim, 2,
                       Encoding::Real, 0, 0.0, 1.0);

        NSGAII<SBXCrossover, PolynomialMutation> nsga;
        nsga.crossover.distribution_index = 20.0;
        nsga.crossover.crossover_probability = 0.9;
        nsga.mutation.distribution_index = 20.0;
        nsga.pop_size = cfg.pop_size;
        nsga.max_evals = cfg.max_evals;

        // Wrap evaluate() as callable for NSGA-II
        auto evaluator = [&prob](Population& p, int i) { prob.evaluate(p, i); };

        Timer timer;
        timer.start();
        nsga.run(pop, evaluator);
        timer.stop();

        total_wall += timer.wall_ms();
        total_cpu += timer.cpu_ms();

        std::vector<std::pair<double,double>> obtained;
        for (int i = 0; i < pop.pop_size; ++i) {
            obtained.push_back({pop.objective(i, 0), pop.objective(i, 1)});
        }
        igd_values[r] = compute_igd(obtained, ref_front);
    }

    result.wall_ms = total_wall / cfg.runs;
    result.cpu_ms = total_cpu / cfg.runs;

    double mean = std::accumulate(igd_values.begin(), igd_values.end(), 0.0) / cfg.runs;
    double variance = 0;
    for (double v : igd_values) variance += (v - mean) * (v - mean);
    variance /= cfg.runs;
    result.mean_igd = mean;
    result.std_igd = std::sqrt(variance);
    return result;
}

// DTLZ runner
template<typename Problem>
TimingResult run_nsga2_dtlz_typed(const std::string& label,
                                   int n_obj, int dim,
                                   Problem& prob,
                                   const BenchConfig& cfg) {
    TimingResult result;
    result.label = label;
    result.population_size = cfg.pop_size;
    result.evaluations = cfg.max_evals;
    result.dimensions = dim;
    result.n_objectives = n_obj;
    result.runs = cfg.runs;

    std::vector<double> igd_values(cfg.runs);
    double total_wall = 0, total_cpu = 0;

    for (int r = 0; r < cfg.runs; ++r) {
        Population pop(cfg.pop_size, dim, n_obj,
                       Encoding::Real, 0, 0.0, 1.0);

        NSGAII<SBXCrossover, PolynomialMutation> nsga;
        nsga.crossover.distribution_index = 20.0;
        nsga.crossover.crossover_probability = 0.9;
        nsga.mutation.distribution_index = 20.0;
        nsga.pop_size = cfg.pop_size;
        nsga.max_evals = cfg.max_evals;

        // wrap evaluate() as callable for NSGA-II
        auto evaluator = [&prob](Population& p, int i) { prob.evaluate(p, i); };

        Timer timer;
        timer.start();
        nsga.run(pop, evaluator);
        timer.stop();

        total_wall += timer.wall_ms();
        total_cpu += timer.cpu_ms();

        // Proxy metric: avg distance to origin
        double sum_dist = 0;
        for (int i = 0; i < pop.pop_size; ++i) {
            double d = 0;
            for (int j = 0; j < n_obj; ++j) d += pop.objective(i, j) * pop.objective(i, j);
            sum_dist += std::sqrt(d);
        }
        igd_values[r] = sum_dist / pop.pop_size;
    }

    result.wall_ms = total_wall / cfg.runs;
    result.cpu_ms = total_cpu / cfg.runs;

    double mean = std::accumulate(igd_values.begin(), igd_values.end(), 0.0) / cfg.runs;
    double variance = 0;
    for (double v : igd_values) variance += (v - mean) * (v - mean);
    variance /= cfg.runs;
    result.mean_igd = mean;
    result.std_igd = std::sqrt(variance);
    return result;
}

} // namespace ea::bench

int main(int argc, char* argv[]) {
    int runs = 30;
    bool csv = false;
    bool verbose = false;

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--runs") == 0 && i + 1 < argc) runs = std::atoi(argv[++i]);
        if (std::strcmp(argv[i], "--csv") == 0) csv = true;
        if (std::strcmp(argv[i], "--verbose") == 0) verbose = true;
    }

    ea::bench::BenchConfig cfg;
    cfg.runs = runs;

    std::vector<ea::bench::TimingResult> results;

    // ZDT benchmarks
    auto zdt1_front = ea::bench::zdt1_front();
    auto zdt2_front = ea::bench::zdt2_front();
    auto zdt3_front = ea::bench::zdt3_front();

    {
        ea::ZDT1 prob(cfg.dim);
        results.push_back(ea::bench::run_nsga2_zdt_typed("NSGA-II-ZDT1", zdt1_front, prob, cfg));
    }
    {
        ea::ZDT2 prob(cfg.dim);
        results.push_back(ea::bench::run_nsga2_zdt_typed("NSGA-II-ZDT2", zdt2_front, prob, cfg));
    }
    {
        ea::ZDT3 prob(cfg.dim);
        results.push_back(ea::bench::run_nsga2_zdt_typed("NSGA-II-ZDT3", zdt3_front, prob, cfg));
    }
    {
        // ZDT4: same optimal front shape as ZDT1
        ea::ZDT4 prob(cfg.dim);
        results.push_back(ea::bench::run_nsga2_zdt_typed("NSGA-II-ZDT4", zdt1_front, prob, cfg));
    }
    {
        // ZDT6: concave front
        auto zdt6_front = ea::bench::zdt2_front();  // approximate
        ea::ZDT6 prob(cfg.dim);
        results.push_back(ea::bench::run_nsga2_zdt_typed("NSGA-II-ZDT6", zdt6_front, prob, cfg));
    }

    // DTLZ benchmarks (3 objectives)
    {
        ea::DTLZ1 prob(7, 3);
        results.push_back(ea::bench::run_nsga2_dtlz_typed("NSGA-II-DTLZ1", 3, 7, prob, cfg));
    }
    {
        ea::DTLZ2 prob(12, 3);
        results.push_back(ea::bench::run_nsga2_dtlz_typed("NSGA-II-DTLZ2", 3, 12, prob, cfg));
    }

    // Output
    if (csv) {
        ea::bench::print_csv_header();
        for (auto& r : results) ea::bench::print_csv_row(r);
    }

    for (auto& r : results) {
        ea::bench::print_result(r);
    }

    if (verbose) {
        std::cerr << "\nNOTE: Factory layer not used in this benchmark (compile-time dispatch)\n";
    }

    return 0;
}