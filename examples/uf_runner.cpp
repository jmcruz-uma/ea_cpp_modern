/// @file uf_runner.cpp
/// @brief NSGA-II benchmark runner — ea_cpp_modern vs jMetal — UF1-UF7 (CEC2009).
///
/// Corre NSGA-II sobre UF1-UF7, NUM_RUNS veces por problema.
/// Salida CSV por stdout: system,problem,run,time_ms,igd
///
/// Compilar (benchmark):
///   g++-14 -std=c++23 -O3 -march=native -DNDEBUG \
///     -I include examples/uf_runner.cpp -o build/uf_runner -lm
///
/// Uso:
///   ./build/uf_runner [NUM_RUNS [BASE_SEED]]
///
/// Parámetros (deben coincidir con JMetalUFBenchmark.java):
///   Population : 100
///   Max evals  : 25000
///   Crossover  : SBX, pc=0.9, eta_c=20
///   Mutation   : PM,  pm=1/dim, eta_m=20
///   Variables  : 30 (todos UF1-7)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include <ea/algorithm/nsga2.hpp>
#include <ea/core/comparator.hpp>
#include <ea/indicator/igd.hpp>
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/problem/uf.hpp>
#include <ea/util/random.hpp>

static constexpr int      POP_SIZE  = 100;
static constexpr int      MAX_EVALS = 25000;
static constexpr double   CX_PROB   = 0.9;
static constexpr double   CX_ETA    = 20.0;
static constexpr double   MUT_ETA   = 20.0;
static constexpr int      NUM_RUNS  = 30;
static constexpr uint64_t BASE_SEED = 42;
static constexpr int      REF_PTS   = 500;

// ── Frentes de referencia analíticos ─────────────────────────────────────────
//
// En el PS (yj = 0 para todo j) los objetivos se reducen a:
//   UF1-3  : f1 = x0,       f2 = 1 - sqrt(x0)   → PF: f2 = 1 - sqrt(f1)
//   UF4    : f1 = x0,       f2 = 1 - x0^2        → PF: f2 = 1 - f1^2
//   UF5 (N=10): hterm = 0 ↔ x0 = k/20            → 21 puntos en f1+f2=1
//   UF6 (N=2) : hj = 0    ↔ f1 ∈ [1/4,1/2]∪[3/4,1] → 2 segmentos de f1+f2=1
//   UF7    : f1 = x0^0.2,  f2 = 1 - x0^0.2       → PF: f1+f2=1

static std::vector<std::vector<double>> make_ref_front(const std::string& name) {
    std::vector<std::vector<double>> ref;

    if (name == "UF1" || name == "UF2" || name == "UF3") {
        // f2 = 1 - sqrt(f1), f1 in [0,1]
        ref.reserve(REF_PTS);
        for (int i = 0; i < REF_PTS; ++i) {
            double f1 = static_cast<double>(i) / (REF_PTS - 1);
            ref.push_back({f1, 1.0 - std::sqrt(f1)});
        }
    } else if (name == "UF4") {
        // f2 = 1 - f1^2, f1 in [0,1]
        ref.reserve(REF_PTS);
        for (int i = 0; i < REF_PTS; ++i) {
            double f1 = static_cast<double>(i) / (REF_PTS - 1);
            ref.push_back({f1, 1.0 - f1 * f1});
        }
    } else if (name == "UF5") {
        // 21 puntos discretos en f1+f2=1 con f1 = k/20, k=0..20
        ref.reserve(21);
        for (int k = 0; k <= 20; ++k) {
            double f1 = k / 20.0;
            ref.push_back({f1, 1.0 - f1});
        }
    } else if (name == "UF6") {
        // Dos segmentos continuos de f1+f2=1:
        //   f1 in [1/4, 1/2]  (250 puntos)
        //   f1 in [3/4, 1]    (250 puntos)
        ref.reserve(REF_PTS);
        int half = REF_PTS / 2;
        for (int i = 0; i < half; ++i) {
            double f1 = 0.25 + 0.25 * i / (half - 1);
            ref.push_back({f1, 1.0 - f1});
        }
        for (int i = 0; i < half; ++i) {
            double f1 = 0.75 + 0.25 * i / (half - 1);
            ref.push_back({f1, 1.0 - f1});
        }
    } else if (name == "UF7") {
        // f1+f2 = 1, f1 in [0,1]
        ref.reserve(REF_PTS);
        for (int i = 0; i < REF_PTS; ++i) {
            double f1 = static_cast<double>(i) / (REF_PTS - 1);
            ref.push_back({f1, 1.0 - f1});
        }
    }
    return ref;
}

// ── Un run — devuelve {time_ms, igd} ─────────────────────────────────────────
template <typename Problem>
std::pair<double, double> run_once(Problem& prob, const std::string& prob_name, uint64_t seed) {
    ea::Random::instance().seed(seed);

    const int dim   = prob.dimension();
    const int n_obj = Problem::num_objectives();

    ea::Population<> pop(POP_SIZE, dim, n_obj);
    pop.lower_bounds.assign(prob.lower_bounds().begin(), prob.lower_bounds().end());
    pop.upper_bounds.assign(prob.upper_bounds().begin(), prob.upper_bounds().end());

    auto& rng = ea::Random::instance();
    for (int i = 0; i < POP_SIZE; ++i)
        for (int j = 0; j < dim; ++j)
            pop.gene(i, j) = rng.uniform(pop.lower_bounds[j], pop.upper_bounds[j]);

    ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> nsga;
    nsga.crossover.crossover_probability = CX_PROB;
    nsga.crossover.distribution_index   = CX_ETA;
    nsga.mutation.distribution_index    = MUT_ETA;
    nsga.pop_size  = POP_SIZE;
    nsga.max_evals = MAX_EVALS;

    auto t0 = std::chrono::high_resolution_clock::now();
    nsga.run(pop, [&](ea::Population<>& p, int idx) { prob.evaluate(p, idx); });
    auto t1 = std::chrono::high_resolution_clock::now();

    double time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    auto fronts  = ea::fast_non_dominated_sort(pop);
    auto ref     = make_ref_front(prob_name);
    double igd_v = ea::igd(pop, fronts[0], ref, 1.0);

    return {time_ms, igd_v};
}

// ── Runner por problema ───────────────────────────────────────────────────────
template <typename Problem>
void benchmark_problem(const std::string& name, int num_runs, uint64_t base_seed) {
    std::cerr << "  " << name << " (dim=30)...\n";
    Problem prob;
    for (int r = 0; r < num_runs; ++r) {
        auto [t, igd] = run_once(prob, name, base_seed + static_cast<uint64_t>(r));
        std::cout << "ea-cpp," << name << "," << (r + 1) << ","
                  << std::fixed << std::setprecision(3) << t << ","
                  << std::setprecision(6) << igd << "\n";
        std::cout.flush();
    }
}

int main(int argc, char* argv[]) {
    int      num_runs  = (argc > 1) ? std::stoi(argv[1])   : NUM_RUNS;
    uint64_t base_seed = (argc > 2) ? std::stoull(argv[2]) : BASE_SEED;

    std::cerr << "=== ea-cpp NSGA-II UF Benchmark ===\n";
    std::cerr << "pop=" << POP_SIZE << "  max_evals=" << MAX_EVALS
              << "  runs=" << num_runs << "  seed=" << base_seed << "\n";

    std::cout << "system,problem,run,time_ms,igd\n";

    benchmark_problem<ea::UF1>("UF1", num_runs, base_seed);
    benchmark_problem<ea::UF2>("UF2", num_runs, base_seed);
    benchmark_problem<ea::UF3>("UF3", num_runs, base_seed);
    benchmark_problem<ea::UF4>("UF4", num_runs, base_seed);
    benchmark_problem<ea::UF5>("UF5", num_runs, base_seed);
    benchmark_problem<ea::UF6>("UF6", num_runs, base_seed);
    benchmark_problem<ea::UF7>("UF7", num_runs, base_seed);

    std::cerr << "Done.\n";
    return 0;
}
