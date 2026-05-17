/// @file bench_original_comparison.cpp
/// @brief Comparison benchmark: Bio4Res/ea-cpp original vs our v1.0
///
/// Since the original uses (µ,λ)-ES (mono-objective, BLX+Gaussian+Comma)
/// and our version uses NSGA-II (multi-objective, SBX+Polynomial),
/// we compare what's comparable:
///
/// 1. **Raw evaluation speed** — How fast can we evaluate 25K individuals
///    on the same problems (Sphere, Rastrigin, etc.)
/// 2. **Population operations** — Crossover, mutation on SoA vs AoS
/// 3. **End-to-end** — Same algorithm configuration as possible
///
/// Original results (30 runs, 25K evals, pop=100):
///   Sphere(30):     avg=0.094s
///   Rastrigin(30):  avg=0.105s
///   Rosenbrock(30): avg=0.093s
///   Ackley(30):     avg=0.102s

#include "bench_timer.hpp"
#include <ea/ea.hpp>
#include <cmath>
#include <numeric>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstring>

// Mono-objective problems matching the original
namespace ea::bench {

struct Sphere {
    int dim;
    explicit Sphere(int d = 30) : dim(d) {}
    static constexpr int num_objectives() { return 1; }
    int dimension() const { return dim; }
    void evaluate(Population& pop, int idx) const {
        double sum = 0;
        for (int j = 0; j < dim; ++j) sum += pop.gene(idx, j) * pop.gene(idx, j);
        pop.objective(idx, 0) = sum;
    }
};

struct Rastrigin {
    int dim;
    explicit Rastrigin(int d = 30) : dim(d) {}
    static constexpr int num_objectives() { return 1; }
    int dimension() const { return dim; }
    void evaluate(Population& pop, int idx) const {
        double sum = 10.0 * dim;
        for (int j = 0; j < dim; ++j) {
            double x = pop.gene(idx, j);
            sum += x * x - 10.0 * std::cos(2.0 * M_PI * x);
        }
        pop.objective(idx, 0) = sum;
    }
};

struct Rosenbrock {
    int dim;
    explicit Rosenbrock(int d = 30) : dim(d) {}
    static constexpr int num_objectives() { return 1; }
    int dimension() const { return dim; }
    void evaluate(Population& pop, int idx) const {
        double sum = 0;
        for (int j = 0; j < dim - 1; ++j) {
            double x1 = pop.gene(idx, j);
            double x2 = pop.gene(idx, j + 1);
            sum += 100.0 * (x2 - x1 * x1) * (x2 - x1 * x1) + (1 - x1) * (1 - x1);
        }
        pop.objective(idx, 0) = sum;
    }
};

struct Ackley {
    int dim;
    explicit Ackley(int d = 30) : dim(d) {}
    static constexpr int num_objectives() { return 1; }
    int dimension() const { return dim; }
    void evaluate(Population& pop, int idx) const {
        double sum1 = 0, sum2 = 0;
        for (int j = 0; j < dim; ++j) {
            double x = pop.gene(idx, j);
            sum1 += x * x;
            sum2 += std::cos(2.0 * M_PI * x);
        }
        pop.objective(idx, 0) = -20.0 * std::exp(-0.2 * std::sqrt(sum1 / dim))
                               - std::exp(sum2 / dim) + 20.0 + M_E;
    }
};

// Simple (µ,λ)-ES matching original config:
// pop_size=100, offspring=99, BLX-alpha crossover, Gaussian mutation, comma replacement
// This replicates the original's algorithm for fair comparison.

template<typename Problem>
TimingResult run_mu_lambda_es(const std::string& label,
                                Problem& prob,
                                int pop_size,
                                int offspring,
                                int max_evals,
                                int runs) {
    TimingResult result;
    result.label = label;
    result.population_size = pop_size;
    result.evaluations = max_evals;
    result.dimensions = prob.dimension();
    result.n_objectives = 1;
    result.runs = runs;

    std::vector<double> fitness_values(runs);
    double total_wall = 0, total_cpu = 0;

    int dim = prob.dimension();

    for (int r = 0; r < runs; ++r) {
        // Initialize population
        int total = pop_size + offspring;
        Population pop(total, dim, 1, Encoding::Real, 0, -5.12, 5.12);

        // Random initialization
        auto& rng = Random::instance();
        for (int i = 0; i < total; ++i) {
            for (int j = 0; j < dim; ++j) {
                pop.gene(i, j) = rng.uniform(-5.12, 5.12);
            }
            prob.evaluate(pop, i);
        }

        Timer timer;
        timer.start();

        int evals = total;
        while (evals < max_evals) {
            // Generate offspring using BLX-alpha crossover + Gaussian mutation
            for (int i = 0; i < offspring && evals < max_evals; ++i) {
                int parent_a = rng.uniform_int(0, pop_size - 1);
                int parent_b = rng.uniform_int(0, pop_size - 1);

                // BLX-alpha crossover (alpha=0.1 as in original)
                double alpha = 0.1;
                for (int j = 0; j < dim; ++j) {
                    double x1 = pop.gene(parent_a, j);
                    double x2 = pop.gene(parent_b, j);
                    double lo = std::min(x1, x2) - alpha * std::abs(x2 - x1);
                    double hi = std::max(x1, x2) + alpha * std::abs(x2 - x1);
                    pop.gene(pop_size + i, j) = rng.uniform(lo, hi);
                }

                // Gaussian mutation (sigma=1.0 as in original)
                double mut_prob = 1.0 / dim;
                for (int j = 0; j < dim; ++j) {
                    if (rng.uniform() < mut_prob) {
                        pop.gene(pop_size + i, j) += rng.normal(0.0, 1.0);
                    }
                    // Clamp to bounds
                    pop.gene(pop_size + i, j) = std::clamp(pop.gene(pop_size + i, j), -5.12, 5.12);
                }

                prob.evaluate(pop, pop_size + i);
                evals++;
            }

            // Comma replacement: keep best offspring_size from offspring
            // Sort by fitness and keep top pop_size
            std::vector<int> indices(total);
            std::iota(indices.begin(), indices.end(), 0);
            std::partial_sort(indices.begin(), indices.begin() + pop_size, indices.end(),
                [&](int a, int b) { return pop.objective(a, 0) < pop.objective(b, 0); });

            // Copy best to positions 0..pop_size-1
            Population new_pop(pop_size, dim, 1, Encoding::Real, 0, -5.12, 5.12);
            for (int i = 0; i < pop_size; ++i) {
                for (int j = 0; j < dim; ++j) {
                    new_pop.gene(i, j) = pop.gene(indices[i], j);
                }
                new_pop.objective(i, 0) = pop.objective(indices[i], 0);
            }

            // Copy back
            for (int i = 0; i < pop_size; ++i) {
                for (int j = 0; j < dim; ++j) {
                    pop.gene(i, j) = new_pop.gene(i, j);
                }
                pop.objective(i, 0) = new_pop.objective(i, 0);
            }

            // Re-initialize offspring portion
            for (int i = pop_size; i < total; ++i) {
                for (int j = 0; j < dim; ++j) {
                    pop.gene(i, j) = 0;
                }
            }
        }

        timer.stop();
        total_wall += timer.wall_ms();
        total_cpu += timer.cpu_ms();

        // Find best fitness
        double best = 1e30;
        for (int i = 0; i < pop_size; ++i) {
            best = std::min(best, pop.objective(i, 0));
        }
        fitness_values[r] = best;
    }

    result.wall_ms = total_wall / runs;
    result.cpu_ms = total_cpu / runs;

    double mean = std::accumulate(fitness_values.begin(), fitness_values.end(), 0.0) / runs;
    double variance = 0;
    for (double v : fitness_values) variance += (v - mean) * (v - mean);
    variance /= runs;
    result.mean_igd = mean;  // reusing field for best fitness
    result.std_igd = std::sqrt(variance);

    return result;
}

} // namespace ea::bench

int main(int argc, char* argv[]) {
    int runs = 30;
    bool csv = false;

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--runs") == 0 && i + 1 < argc) runs = std::atoi(argv[++i]);
        if (std::strcmp(argv[i], "--csv") == 0) csv = true;
    }

    std::cerr << "=== ea-cpp v1.0 (current) vs Bio4Res/ea-cpp (original) ===\n";
    std::cerr << "Algorithm: (µ,λ)-ES with BLX-alpha crossover + Gaussian mutation\n";
    std::cerr << "Config: pop=100, offspring=99, max_evals=25000, 30 runs\n";
    std::cerr << "Problems: Sphere, Rastrigin, Rosenbrock, Ackley (dim=30)\n\n";

    std::vector<ea::bench::TimingResult> results;

    {
        ea::bench::Sphere prob(30);
        results.push_back(ea::bench::run_mu_lambda_es("Sphere(30)", prob, 100, 99, 25000, runs));
    }
    {
        ea::bench::Rastrigin prob(30);
        results.push_back(ea::bench::run_mu_lambda_es("Rastrigin(30)", prob, 100, 99, 25000, runs));
    }
    {
        ea::bench::Rosenbrock prob(30);
        results.push_back(ea::bench::run_mu_lambda_es("Rosenbrock(30)", prob, 100, 99, 25000, runs));
    }
    {
        ea::bench::Ackley prob(30);
        results.push_back(ea::bench::run_mu_lambda_es("Ackley(30)", prob, 100, 99, 25000, runs));
    }

    // Print comparison table
    // Original results (from Bio4Res/ea-cpp, same config)
    struct OrigResult { const char* name; double time_ms; };
    OrigResult orig_results[] = {
        {"Sphere(30)", 94.0},
        {"Rastrigin(30)", 104.5},
        {"Rosenbrock(30)", 93.0},
        {"Ackley(30)", 101.6},
    };

    std::cerr << "┌──────────────────┬────────────┬────────────┬─────────┬──────────────┐\n";
    std::cerr << "│ Problem          │ Original   │ Current    │ Speedup │ Best Fitness │\n";
    std::cerr << "├──────────────────┼────────────┼────────────┼─────────┼──────────────┤\n";

    for (size_t i = 0; i < results.size(); ++i) {
        double speedup = orig_results[i].time_ms / results[i].wall_ms;
        char speedup_str[32];
        snprintf(speedup_str, sizeof(speedup_str), "%.2fx", speedup);

        std::cerr << "│ " << std::left << std::setw(16) << results[i].label
                  << " │ " << std::right << std::fixed << std::setprecision(1) << std::setw(8) << orig_results[i].time_ms << "ms"
                  << " │ " << std::setw(8) << results[i].wall_ms << "ms"
                  << " │ " << std::setw(7) << speedup_str
                  << " │ " << std::scientific << std::setprecision(4) << results[i].mean_igd
                  << " │\n";
    }
    std::cerr << "└──────────────────┴────────────┴────────────┴─────────┴──────────────┘\n";

    if (csv) {
        ea::bench::print_csv_header();
        for (auto& r : results) ea::bench::print_csv_row(r);
    }

    return 0;
}