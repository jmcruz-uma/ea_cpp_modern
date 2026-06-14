/// @file spea2_runner.cpp
/// @brief SPEA2 MO benchmark runner — ea_cpp_modern vs jMetal.
///
/// Corre SPEA2 sobre ZDT1-ZDT4, NUM_RUNS veces por problema.
/// Salida CSV por stdout: system,problem,run,time_ms,igd
///
/// Compilar (benchmark):
///   g++-14 -std=c++23 -O3 -march=native -DNDEBUG -I include examples/spea2_runner.cpp -o build/spea2_runner -lm
///
/// Uso:
///   ./build/spea2_runner [NUM_RUNS [BASE_SEED]]

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include <ea/algorithm/spea2.hpp>
#include <ea/indicator/igd.hpp>
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/problem/zdt.hpp>
#include <ea/util/random.hpp>

// ── Parámetros del benchmark (deben coincidir con JMetalSPEA2Benchmark.java) ──
static constexpr int      POP_SIZE      = 100;
static constexpr int      MAX_EVALS     = 25000;  // 100 init + 249 iters × 100
static constexpr int      ARCHIVE_SIZE  = 100;
static constexpr double   CX_PROB       = 0.9;
static constexpr double   CX_ETA        = 20.0;
static constexpr double   MUT_ETA       = 20.0;
static constexpr int      REF_PTS       = 500;
static constexpr uint64_t BASE_SEED     = 42;
static constexpr int      NUM_RUNS      = 30;

// ── Frente de referencia analítico (500 puntos) ───────────────────────────────
static std::vector<std::vector<double>> make_ref_front(const std::string& name) {
    std::vector<std::vector<double>> ref;
    ref.reserve(REF_PTS);

    if (name == "ZDT1" || name == "ZDT4") {
        for (int i = 0; i < REF_PTS; ++i) {
            double f1 = static_cast<double>(i) / (REF_PTS - 1);
            ref.push_back({f1, 1.0 - std::sqrt(f1)});
        }
    } else if (name == "ZDT2") {
        for (int i = 0; i < REF_PTS; ++i) {
            double f1 = static_cast<double>(i) / (REF_PTS - 1);
            ref.push_back({f1, 1.0 - f1 * f1});
        }
    } else if (name == "ZDT3") {
        const double segs[5][2] = {
            {0.0,         0.0830015349},
            {0.182228780, 0.257762881},
            {0.409525530, 0.453861524},
            {0.617876091, 0.652585810},
            {0.823295506, 0.851832239}
        };
        int pps = REF_PTS / 5;
        for (const auto& seg : segs) {
            double lo = seg[0], hi = seg[1];
            for (int j = 0; j < pps; ++j) {
                double f1 = lo + (hi - lo) * j / (pps - 1);
                double f2 = 1.0 - std::sqrt(f1) - f1 * std::sin(10.0 * M_PI * f1);
                ref.push_back({f1, f2});
            }
        }
    }
    return ref;
}

// ── Un run de SPEA2 — devuelve {time_ms, igd} ─────────────────────────────────
template <typename Problem>
std::pair<double, double> run_once(Problem& prob, const std::string& prob_name, uint64_t seed) {
    ea::Random::instance().seed(seed);

    const int dim   = prob.dimension();
    const int n_obj = prob.num_objectives();

    ea::Population<> pop(POP_SIZE, dim, n_obj);
    pop.lower_bounds.assign(prob.lower_bounds().begin(), prob.lower_bounds().end());
    pop.upper_bounds.assign(prob.upper_bounds().begin(), prob.upper_bounds().end());

    auto& rng = ea::Random::instance();
    for (int i = 0; i < POP_SIZE; ++i)
        for (int j = 0; j < dim; ++j)
            pop.gene(i, j) = rng.uniform(pop.lower_bounds[j], pop.upper_bounds[j]);

    ea::SPEA2<ea::SBXCrossover, ea::PolynomialMutation> spea2;
    spea2.pop_size                       = POP_SIZE;
    spea2.archive_size                   = ARCHIVE_SIZE;
    spea2.max_evals                      = MAX_EVALS;
    spea2.crossover.crossover_probability = CX_PROB;
    spea2.crossover.distribution_index   = CX_ETA;
    spea2.mutation.distribution_index    = MUT_ETA;
    // mutation_rate = -1 → auto 1/dim (matches jMetal PolynomialMutation(1/D, 20))

    auto t0 = std::chrono::high_resolution_clock::now();
    spea2.run(pop, [&](ea::Population<>& p, int idx) { prob.evaluate(p, idx); });
    auto t1 = std::chrono::high_resolution_clock::now();

    double time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Resultado: pop contiene el archivo externo (coincide con jMetal algorithm.result())
    std::vector<std::vector<double>> archive_front;
    archive_front.reserve(POP_SIZE);
    for (int i = 0; i < POP_SIZE; ++i) {
        archive_front.push_back(
            std::vector<double>(pop.objectives_ptr(i), pop.objectives_ptr(i) + n_obj));
    }

    auto ref = make_ref_front(prob_name);
    double igd_val = ea::igd(archive_front, ref, 1.0);

    return {time_ms, igd_val};
}

// ── Runner por problema ───────────────────────────────────────────────────────
template <typename Problem>
void benchmark_problem(const std::string& name, int dim, int num_runs, uint64_t base_seed) {
    std::cerr << "  " << name << " (dim=" << dim << ")...\n";
    Problem prob(dim);
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

    std::cerr << "=== ea-cpp SPEA2 MO Benchmark ===\n";
    std::cerr << "pop=" << POP_SIZE << "  archive=" << ARCHIVE_SIZE
              << "  max_evals=" << MAX_EVALS
              << "  cx_prob=" << CX_PROB << "  cx_eta=" << CX_ETA
              << "  mut_eta=" << MUT_ETA
              << "  runs=" << num_runs << "  seed=" << base_seed << "\n";

    std::cout << "system,problem,run,time_ms,igd\n";

    benchmark_problem<ea::ZDT1>("ZDT1", 30, num_runs, base_seed);
    benchmark_problem<ea::ZDT2>("ZDT2", 30, num_runs, base_seed);
    benchmark_problem<ea::ZDT3>("ZDT3", 30, num_runs, base_seed);
    benchmark_problem<ea::ZDT4>("ZDT4", 10, num_runs, base_seed);

    std::cerr << "Done.\n";
    return 0;
}
