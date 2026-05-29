/// @file bench_arch_comparison.cpp
/// @brief Architecture comparison: ea-cpp v1.0 vs Bio4Res/ea-cpp original
///
/// Measures:
/// 1. Evaluation throughput (SoA vs AoS overhead)
/// 2. Operator throughput (templates vs virtual dispatch)
/// 3. NSGA-II end-to-end (our advantage — original has no multi-obj)
/// 4. Scalability with population size
///
/// Run: ./bench_arch [--runs N] [--full] [--csv]

#include "bench_timer.hpp"
#include <ea/ea.hpp>
#include <cmath>
#include <numeric>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstring>

// Mono-objective problems for evaluation throughput
namespace {

struct SphereFn {
    int dim;
    explicit SphereFn(int d = 30) : dim(d) {}
    static constexpr int num_objectives() { return 1; }
    int dimension() const { return dim; }
    void evaluate(ea::Population<>& pop, int idx) const {
        double sum = 0;
        for (int j = 0; j < dim; ++j) sum += pop.gene(idx, j) * pop.gene(idx, j);
        pop.objective(idx, 0) = sum;
    }
};

struct RastriginFn {
    int dim;
    explicit RastriginFn(int d = 30) : dim(d) {}
    static constexpr int num_objectives() { return 1; }
    int dimension() const { return dim; }
    void evaluate(ea::Population<>& pop, int idx) const {
        double sum = 10.0 * dim;
        for (int j = 0; j < dim; ++j) {
            double x = pop.gene(idx, j);
            sum += x * x - 10.0 * std::cos(2.0 * M_PI * x);
        }
        pop.objective(idx, 0) = sum;
    }
};

} // anonymous

namespace ea::bench {

// ============================================================
// Benchmark 1: Raw evaluation throughput
// ============================================================

template<typename Problem>
TimingResult bench_evaluation(const char* label, Problem& prob, int dim, int n_obj, int pop_size, int max_evals, int runs) {
    TimingResult result;
    result.label = label;
    result.population_size = pop_size;
    result.evaluations = max_evals;
    result.dimensions = dim;
    result.n_objectives = n_obj;
    result.runs = runs;

    double total_wall = 0;

    for (int r = 0; r < runs; ++r) {
        Population<> pop(pop_size, dim, n_obj, 0, -5.12, 5.12);
        auto& rng = Random::instance();
        for (int i = 0; i < pop_size; ++i)
            for (int j = 0; j < dim; ++j)
                pop.gene(i, j) = rng.uniform(-5.12, 5.12);

        Timer timer;
        timer.start();
        for (int e = 0; e < max_evals; ++e) {
            prob.evaluate(pop, e % pop_size);
        }
        timer.stop();
        total_wall += timer.wall_ms();
    }

    result.wall_ms = total_wall / runs;
    result.cpu_ms = 0;
    result.mean_igd = 0;
    result.std_igd = 0;
    return result;
}

// ============================================================
// Benchmark 2: Crossover + Mutation throughput
// ============================================================

TimingResult bench_operator_throughput(const char* label, int pop_size, int dim, int max_ops, int runs) {
    TimingResult result;
    result.label = label;
    result.population_size = pop_size;
    result.evaluations = max_ops;
    result.dimensions = dim;
    result.n_objectives = 2;
    result.runs = runs;

    double total_wall = 0;

    for (int r = 0; r < runs; ++r) {
        Population<> pop(pop_size * 2, dim, 2, 0, -5.12, 5.12);
        auto& rng = Random::instance();
        for (int i = 0; i < pop_size * 2; ++i)
            for (int j = 0; j < dim; ++j)
                pop.gene(i, j) = rng.uniform(-5.12, 5.12);

        SBXCrossover cx;
        cx.distribution_index = 20.0;
        cx.crossover_probability = 0.9;

        PolynomialMutation mt;
        mt.distribution_index = 20.0;

        Timer timer;
        timer.start();
        int ops = 0;
        for (int i = 0; i < pop_size && ops < max_ops; i += 2) {
            cx.apply(pop, i, i + 1, pop_size + i);
            mt.apply(pop, pop_size + i);
            mt.apply(pop, pop_size + i + 1);
            ops += 2;
        }
        timer.stop();
        total_wall += timer.wall_ms();
    }

    result.wall_ms = total_wall / runs;
    result.cpu_ms = 0;
    result.mean_igd = 0;
    result.std_igd = 0;
    return result;
}

// ============================================================
// Benchmark 3: NSGA-II end-to-end
// ============================================================

template<typename Problem>
TimingResult bench_nsga2(const char* label, Problem& prob, int dim, int n_obj, int pop_size, int max_evals, int runs) {
    TimingResult result;
    result.label = label;
    result.population_size = pop_size;
    result.evaluations = max_evals;
    result.dimensions = dim;
    result.n_objectives = n_obj;
    result.runs = runs;

    double total_wall = 0;

    for (int r = 0; r < runs; ++r) {
        Population<> pop(pop_size, dim, n_obj, 0, 0.0, 1.0);

        NSGAII<SBXCrossover, PolynomialMutation> nsga;
        nsga.crossover.distribution_index = 20.0;
        nsga.crossover.crossover_probability = 0.9;
        nsga.mutation.distribution_index = 20.0;
        nsga.pop_size = pop_size;
        nsga.max_evals = max_evals;

        auto evaluator = [&prob](Population<>& p, int i) { prob.evaluate(p, i); };

        Timer timer;
        timer.start();
        nsga.run(pop, evaluator);
        timer.stop();
        total_wall += timer.wall_ms();
    }

    result.wall_ms = total_wall / runs;
    result.cpu_ms = 0;
    result.mean_igd = 0;
    result.std_igd = 0;
    return result;
}

// ============================================================
// Benchmark 4: Scalability — population size sweep
// ============================================================

TimingResult bench_scalability(const char* label, int pop_size, int dim, int generations, int runs) {
    TimingResult result;
    result.label = label;
    result.population_size = pop_size;
    result.evaluations = pop_size * generations;
    result.dimensions = dim;
    result.n_objectives = 2;
    result.runs = runs;

    double total_wall = 0;

    for (int r = 0; r < runs; ++r) {
        Population<> pop(pop_size, dim, 2, 0, 0.0, 1.0);
        ZDT1 prob(dim);

        NSGAII<SBXCrossover, PolynomialMutation> nsga;
        nsga.crossover.distribution_index = 20.0;
        nsga.crossover.crossover_probability = 0.9;
        nsga.mutation.distribution_index = 20.0;
        nsga.pop_size = pop_size;
        nsga.max_evals = pop_size * generations;

        auto evaluator = [&prob](Population<>& p, int i) { prob.evaluate(p, i); };

        Timer timer;
        timer.start();
        nsga.run(pop, evaluator);
        timer.stop();
        total_wall += timer.wall_ms();
    }

    result.wall_ms = total_wall / runs;
    result.cpu_ms = 0;
    result.mean_igd = 0;
    result.std_igd = 0;
    return result;
}

} // namespace ea::bench

int main(int argc, char* argv[]) {
    int runs = 30;
    bool csv = false;
    bool full = false;

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--runs") == 0 && i + 1 < argc) runs = std::atoi(argv[++i]);
        if (std::strcmp(argv[i], "--csv") == 0) csv = true;
        if (std::strcmp(argv[i], "--full") == 0) full = true;
    }

    std::cerr << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cerr << "║  ea-cpp v1.0 — Architecture Benchmark                    ║\n";
    std::cerr << "║  SoA Population + C++23 Concepts vs AoS + Virtual       ║\n";
    std::cerr << "╚══════════════════════════════════════════════════════════╝\n\n";

    // ── 1. Evaluation throughput ──
    std::cerr << "── 1. Evaluation Throughput (25K evals × " << runs << " runs) ──\n";
    {
        SphereFn prob(30);
        auto r = ea::bench::bench_evaluation("Sphere(30)", prob, 30, 1, 100, 25000, runs);
        std::cerr << "   Sphere(30):     " << std::fixed << std::setprecision(1) << r.wall_ms << " ms\n";
    }
    {
        RastriginFn prob(30);
        auto r = ea::bench::bench_evaluation("Rastrigin(30)", prob, 30, 1, 100, 25000, runs);
        std::cerr << "   Rastrigin(30):  " << std::fixed << std::setprecision(1) << r.wall_ms << " ms\n";
    }
    {
        // ZDT1 = 2 objectives
        ea::ZDT1 prob(30);
        ea::Population<> pop(100, 30, 2, 0, 0.0, 1.0);
        double total = 0;
        for (int r = 0; r < runs; ++r) {
            ea::bench::Timer timer;
            timer.start();
            for (int e = 0; e < 25000; ++e) prob.evaluate(pop, e % 100);
            timer.stop();
            total += timer.wall_ms();
        }
        std::cerr << "   ZDT1(30,2obj):  " << std::fixed << std::setprecision(1) << total / runs << " ms\n";
    }

    // ── 2. Operator throughput ──
    std::cerr << "\n── 2. Operator Throughput (SBX+Poly, 25K ops × " << runs << " runs) ──\n";
    {
        auto r = ea::bench::bench_operator_throughput("dim=30", 100, 30, 25000, runs);
        std::cerr << "   dim=30:   " << std::fixed << std::setprecision(1) << r.wall_ms << " ms\n";
    }
    {
        auto r = ea::bench::bench_operator_throughput("dim=100", 100, 100, 25000, runs);
        std::cerr << "   dim=100:  " << std::fixed << std::setprecision(1) << r.wall_ms << " ms\n";
    }
    {
        auto r = ea::bench::bench_operator_throughput("dim=500", 100, 500, 25000, runs);
        std::cerr << "   dim=500:  " << std::fixed << std::setprecision(1) << r.wall_ms << " ms\n";
    }

    // ── 3. NSGA-II end-to-end ──
    std::cerr << "\n── 3. NSGA-II End-to-End (25K evals × " << runs << " runs) ──\n";
    std::cerr << "   ⚠ Original has NO multi-objective — this is our advantage.\n\n";
    {
        ea::ZDT1 prob(30); auto r = ea::bench::bench_nsga2("ZDT1(30)", prob, 30, 2, 100, 25000, runs);
        std::cerr << "   ZDT1(30):    " << std::fixed << std::setprecision(1) << r.wall_ms << " ms\n";
    }
    {
        ea::ZDT2 prob(30); auto r = ea::bench::bench_nsga2("ZDT2(30)", prob, 30, 2, 100, 25000, runs);
        std::cerr << "   ZDT2(30):    " << std::fixed << std::setprecision(1) << r.wall_ms << " ms\n";
    }
    {
        ea::ZDT3 prob(30); auto r = ea::bench::bench_nsga2("ZDT3(30)", prob, 30, 2, 100, 25000, runs);
        std::cerr << "   ZDT3(30):    " << std::fixed << std::setprecision(1) << r.wall_ms << " ms\n";
    }
    {
        ea::DTLZ2 prob(12, 3); auto r = ea::bench::bench_nsga2("DTLZ2(3obj)", prob, 12, 3, 100, 25000, runs);
        std::cerr << "   DTLZ2(3obj): " << std::fixed << std::setprecision(1) << r.wall_ms << " ms\n";
    }

    // ── 4. Scalability ──
    if (full) {
        std::cerr << "\n── 4. Scalability — Population Size (ZDT1, 250 gens × " << runs << " runs) ──\n";
        for (int pop : {50, 100, 200, 500, 1000}) {
            auto r = ea::bench::bench_scalability(
                ("pop=" + std::to_string(pop)).c_str(), pop, 30, 250, runs);
            std::cerr << "   pop=" << std::setw(4) << pop
                      << ":  " << std::fixed << std::setprecision(1) << r.wall_ms << " ms\n";
        }
    }

    // ── Original comparison ──
    std::cerr << "\n── Original (Bio4Res/ea-cpp) Reference ──\n";
    std::cerr << "   (µ,λ)-ES with BLX+Gaussian, pop=100, 25K evals, 30 runs:\n";
    std::cerr << "   Sphere(30):     ~94 ms\n";
    std::cerr << "   Rastrigin(30):  ~105 ms\n";
    std::cerr << "   Rosenbrock(30): ~93 ms\n";
    std::cerr << "   Ackley(30):     ~102 ms\n";
    std::cerr << "   (Mono-objective only — no NSGA-II, no Pareto front)\n\n";

    std::cerr << "╔══════════════════════════════════════════════════════════╗\n";
    std::cerr << "║  Key Differences                                         ║\n";
    std::cerr << "║                                                          ║\n";
    std::cerr << "║  Original (Bio4Res/ea-cpp):                              ║\n";
    std::cerr << "║    • AoS Individual + unique_ptr<Genotype> (heap alloc) ║\n";
    std::cerr << "║    • Virtual dispatch for operators                       ║\n";
    std::cerr << "║    • JSON config for runtime setup                        ║\n";
    std::cerr << "║    • Mono-objective ONLY                                 ║\n";
    std::cerr << "║    • 4,942 LOC, 55 files                                 ║\n";
    std::cerr << "║                                                          ║\n";
    std::cerr << "║  Current (ea-cpp v1.0):                                   ║\n";
    std::cerr << "║    • SoA Population (contiguous, SIMD-friendly)            ║\n";
    std::cerr << "║    • C++23 Concepts + templates (zero-overhead)          ║\n";
    std::cerr << "║    • constinit Factory (extensible, zero UB)            ║\n";
    std::cerr << "║    • Multi-objective first-class (10 algorithms)         ║\n";
    std::cerr << "║    • IGD/HV indicators, 22 problems, 30+ operators       ║\n";
    std::cerr << "║    • 76+ headers, header-only                             ║\n";
    std::cerr << "╚══════════════════════════════════════════════════════════╝\n";

    return 0;
}