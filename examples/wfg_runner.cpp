/// @file wfg_runner.cpp
/// @brief NSGA-II benchmark runner — ea_cpp_modern vs jMetal — WFG1-9.
///
/// Corre NSGA-II sobre WFG1-9, NUM_RUNS veces por problema.
/// Salida CSV por stdout: system,problem,run,time_ms,igd
///
/// Compilar (benchmark):
///   g++-14 -std=c++23 -O3 -march=native -DNDEBUG \
///     -I include examples/wfg_runner.cpp -o build/wfg_runner -lm
///
/// Uso:
///   ./build/wfg_runner [NUM_RUNS [BASE_SEED [REF_FRONT_DIR]]]
///
/// Parámetros (deben coincidir con JMetalWFGBenchmark.java):
///   Population : 100
///   Max evals  : 25000
///   Crossover  : SBX, pc=0.9, eta_c=20
///   Mutation   : PM,  pm=1/dim, eta_m=20
///   k=2, l=4, M=2, dim=6  (jMetal DefaultWFGSettings)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <ea/algorithm/nsga2.hpp>
#include <ea/core/comparator.hpp>
#include <ea/indicator/igd.hpp>
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/problem/wfg.hpp>
#include <ea/util/random.hpp>

static constexpr int      POP_SIZE  = 100;
static constexpr int      MAX_EVALS = 25000;
static constexpr double   CX_PROB   = 0.9;
static constexpr double   CX_ETA    = 20.0;
static constexpr double   MUT_ETA   = 20.0;
static constexpr int      NUM_RUNS  = 30;
static constexpr uint64_t BASE_SEED = 42;

// jMetal DefaultWFGSettings: k=2, l=4, M=2
static constexpr int WFG_M = 2;
static constexpr int WFG_K = 2;
static constexpr int WFG_L = 4;

static const char* DEFAULT_REF_DIR =
    "/home/alumno/jMetal/resources/referenceFrontsCSV";

// ── Cargar frente de referencia desde CSV ─────────────────────────────────────
//
// Formato: f1,f2 (sin cabecera), ~1113 puntos por archivo.
// Ejemplo: /home/alumno/jMetal/resources/referenceFrontsCSV/WFG1.2D.csv

static std::vector<std::vector<double>> load_ref_front(const std::string& path) {
    std::vector<std::vector<double>> ref;
    std::ifstream f(path);
    if (!f) {
        std::cerr << "ERROR: no se puede abrir " << path << "\n";
        std::exit(1);
    }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string s1, s2;
        if (!std::getline(ss, s1, ',') || !std::getline(ss, s2)) continue;
        ref.push_back({std::stod(s1), std::stod(s2)});
    }
    return ref;
}

// ── Un run — devuelve {time_ms, igd} ─────────────────────────────────────────
template <typename Problem>
std::pair<double, double> run_once(Problem& prob,
                                   const std::vector<std::vector<double>>& ref,
                                   uint64_t seed) {
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
    double igd_v = ea::igd(pop, fronts[0], ref, 1.0);

    return {time_ms, igd_v};
}

// ── Runner por problema ───────────────────────────────────────────────────────
template <typename Problem>
void benchmark_problem(const std::string& name,
                       const std::vector<std::vector<double>>& ref,
                       int num_runs, uint64_t base_seed) {
    Problem prob(WFG_M, WFG_K, WFG_L);
    std::cerr << "  " << name << " (dim=" << prob.dimension() << ")...\n";
    for (int r = 0; r < num_runs; ++r) {
        auto [t, igd] = run_once(prob, ref, base_seed + static_cast<uint64_t>(r));
        std::cout << "ea-cpp," << name << "," << (r + 1) << ","
                  << std::fixed << std::setprecision(3) << t << ","
                  << std::setprecision(6) << igd << "\n";
        std::cout.flush();
    }
}

int main(int argc, char* argv[]) {
    int         num_runs  = (argc > 1) ? std::stoi(argv[1])   : NUM_RUNS;
    uint64_t    base_seed = (argc > 2) ? std::stoull(argv[2]) : BASE_SEED;
    std::string ref_dir   = (argc > 3) ? argv[3]              : DEFAULT_REF_DIR;

    std::cerr << "=== ea-cpp NSGA-II WFG Benchmark ===\n";
    std::cerr << "pop=" << POP_SIZE << "  max_evals=" << MAX_EVALS
              << "  M=" << WFG_M << "  k=" << WFG_K << "  l=" << WFG_L
              << "  dim=" << (WFG_K + WFG_L)
              << "  runs=" << num_runs << "  seed=" << base_seed << "\n";
    std::cerr << "ref_dir=" << ref_dir << "\n";

    std::cout << "system,problem,run,time_ms,igd\n";

    for (int n = 1; n <= 9; ++n) {
        std::string name = "WFG" + std::to_string(n);
        std::string csv  = ref_dir + "/" + name + ".2D.csv";
        auto ref = load_ref_front(csv);
        std::cerr << "  [" << name << "] ref front: " << ref.size() << " pts\n";

        switch (n) {
            case 1: benchmark_problem<ea::WFG1>(name, ref, num_runs, base_seed); break;
            case 2: benchmark_problem<ea::WFG2>(name, ref, num_runs, base_seed); break;
            case 3: benchmark_problem<ea::WFG3>(name, ref, num_runs, base_seed); break;
            case 4: benchmark_problem<ea::WFG4>(name, ref, num_runs, base_seed); break;
            case 5: benchmark_problem<ea::WFG5>(name, ref, num_runs, base_seed); break;
            case 6: benchmark_problem<ea::WFG6>(name, ref, num_runs, base_seed); break;
            case 7: benchmark_problem<ea::WFG7>(name, ref, num_runs, base_seed); break;
            case 8: benchmark_problem<ea::WFG8>(name, ref, num_runs, base_seed); break;
            case 9: benchmark_problem<ea::WFG9>(name, ref, num_runs, base_seed); break;
        }
    }

    std::cerr << "Done.\n";
    return 0;
}
