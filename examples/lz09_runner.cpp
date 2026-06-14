/// @file lz09_runner.cpp
/// @brief LZ09 benchmark runner — ea_cpp_modern vs jMetal — LZ09F1-F9 (Li & Zhang 2009).
///
/// Corre NSGA-II sobre LZ09F1-F5,F7-F9 (2 obj) y NSGA-III sobre LZ09F6 (3 obj).
/// Salida CSV por stdout: system,problem,run,time_ms,igd
///
/// Compilar (benchmark):
///   g++-14 -std=c++23 -O3 -march=native -DNDEBUG \
///     -I include examples/lz09_runner.cpp -o build/lz09_runner -lm
///
/// Uso:
///   ./build/lz09_runner [NUM_RUNS [BASE_SEED]]
///
/// Parámetros (deben coincidir con JMetalLZ09Benchmark.java):
///   NSGA-II  : pop=100, max_evals=25000, SBX(pc=0.9, η=20), PM(η=20, p=1/dim)
///   NSGA-III : divisions=12 → 91 ref pts → pop=92, max_evals=25000
///
/// Dimensiones por problema:
///   LZ09F1/F6/F7/F8 → dim=10;  LZ09F2/F3/F4/F5/F9 → dim=30

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include <ea/algorithm/nsga2.hpp>
#include <ea/algorithm/nsga3.hpp>
#include <ea/core/comparator.hpp>
#include <ea/indicator/igd.hpp>
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/problem/lz09.hpp>
#include <ea/util/random.hpp>
#include <ea/util/reference_point.hpp>

// ── Parámetros NSGA-II ─────────────────────────────────────────────────────────
static constexpr int      POP_SIZE  = 100;
static constexpr int      MAX_EVALS = 25000;
static constexpr double   CX_PROB   = 0.9;
static constexpr double   CX_ETA    = 20.0;
static constexpr double   MUT_ETA   = 20.0;
static constexpr int      NUM_RUNS  = 30;
static constexpr uint64_t BASE_SEED = 42;
static constexpr int      REF_PTS   = 500;

// ── Parámetros NSGA-III (LZ09F6, 3-obj) ───────────────────────────────────────
static constexpr int N_OBJ_F6  = 3;
static constexpr int DIVISIONS = 12;   // 91 ref pts → pop_size=92
static constexpr int REF_DIV   = 30;   // Das-Dennis divisions para frente de referencia

// ── Frentes de referencia 2-obj ───────────────────────────────────────────────
//
// Con yj=0 en el PS (betaFunction=0), los objetivos se reducen a alphaFunction:
//   ptype=21  → f2 = 1 - sqrt(f1),  f1∈[0,1]   (LZ09F1-F5, F7-F8)
//   ptype=22  → f2 = 1 - f1²,       f1∈[0,1]   (LZ09F9)

static std::vector<std::vector<double>> make_ref_front_2obj(int ptype) {
    std::vector<std::vector<double>> ref;
    ref.reserve(REF_PTS);
    for (int i = 0; i < REF_PTS; ++i) {
        double f1 = static_cast<double>(i) / (REF_PTS - 1);
        if (ptype == 21)
            ref.push_back({f1, 1.0 - std::sqrt(f1)});
        else // ptype == 22
            ref.push_back({f1, 1.0 - f1 * f1});
    }
    return ref;
}

// ── Frente de referencia 3-obj (cuarto de esfera, LZ09F6) ─────────────────────
//
// ptype=31 → f1²+f2²+f3²=1, fi≥0  (igual que DTLZ2)
// Generamos puntos Das-Dennis y normalizamos a la esfera unitaria.

static std::vector<std::vector<double>> make_ref_front_sphere() {
    std::vector<ea::ReferencePoint> pts;
    ea::generate_reference_points_das_dennis(pts, N_OBJ_F6, REF_DIV);
    std::vector<std::vector<double>> ref;
    ref.reserve(pts.size());
    for (const auto& rp : pts) {
        double norm = std::sqrt(rp.position[0] * rp.position[0] +
                                rp.position[1] * rp.position[1] +
                                rp.position[2] * rp.position[2]);
        if (norm < 1e-14) norm = 1.0;
        ref.push_back({rp.position[0] / norm,
                       rp.position[1] / norm,
                       rp.position[2] / norm});
    }
    return ref;
}

// ── Un run NSGA-II (2-obj) ─────────────────────────────────────────────────────
template <typename Problem>
std::pair<double, double> run_once_2obj(Problem& prob,
                                        const std::vector<std::vector<double>>& ref,
                                        uint64_t seed) {
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
    double igd_v = ea::igd(pop, fronts[0], ref, 1.0);

    return {time_ms, igd_v};
}

// ── Un run NSGA-III (3-obj, LZ09F6) ───────────────────────────────────────────
std::pair<double, double> run_once_3obj(ea::LZ09F6& prob,
                                        const std::vector<std::vector<double>>& ref,
                                        uint64_t seed) {
    ea::Random::instance().seed(seed);

    const int dim   = prob.dimension();
    const int n_obj = ea::LZ09F6::num_objectives();

    std::vector<ea::ReferencePoint> tmp_refs;
    ea::generate_reference_points_das_dennis(tmp_refs, n_obj, DIVISIONS);
    int n = ea::compute_population_size(static_cast<int>(tmp_refs.size()));

    ea::Population<> pop(n, dim, n_obj);
    pop.lower_bounds.assign(prob.lower_bounds().begin(), prob.lower_bounds().end());
    pop.upper_bounds.assign(prob.upper_bounds().begin(), prob.upper_bounds().end());

    auto& rng = ea::Random::instance();
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < dim; ++j)
            pop.gene(i, j) = rng.uniform(pop.lower_bounds[j], pop.upper_bounds[j]);

    ea::NSGAIII<ea::SBXCrossover, ea::PolynomialMutation> alg;
    alg.pop_size  = n;
    alg.max_evals = MAX_EVALS;
    alg.divisions = DIVISIONS;
    alg.crossover.crossover_probability = CX_PROB;
    alg.crossover.distribution_index    = CX_ETA;
    alg.mutation.distribution_index     = MUT_ETA;

    auto t0 = std::chrono::high_resolution_clock::now();
    alg.run(pop, [&](ea::Population<>& p, int idx) { prob.evaluate(p, idx); });
    auto t1 = std::chrono::high_resolution_clock::now();

    double time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    auto fronts  = ea::fast_non_dominated_sort(pop);
    double igd_v = ea::igd(pop, fronts[0], ref, 1.0);

    return {time_ms, igd_v};
}

// ── Runner por problema 2-obj ─────────────────────────────────────────────────
template <typename Problem>
void benchmark_2obj(const std::string& name, int ptype,
                    int num_runs, uint64_t base_seed) {
    Problem prob;
    std::cerr << "  " << name << " (dim=" << prob.dimension() << ")...\n";
    auto ref = make_ref_front_2obj(ptype);
    for (int r = 0; r < num_runs; ++r) {
        auto [t, igd] = run_once_2obj(prob, ref, base_seed + static_cast<uint64_t>(r));
        std::cout << "ea-cpp," << name << "," << (r + 1) << ","
                  << std::fixed << std::setprecision(3) << t << ","
                  << std::setprecision(6) << igd << "\n";
        std::cout.flush();
    }
}

// ── Runner LZ09F6 (NSGA-III, 3-obj) ──────────────────────────────────────────
void benchmark_lz09f6(int num_runs, uint64_t base_seed) {
    ea::LZ09F6 prob;
    std::cerr << "  LZ09F6 (dim=" << prob.dimension() << ", 3-obj, NSGA-III)...\n";
    auto ref = make_ref_front_sphere();
    for (int r = 0; r < num_runs; ++r) {
        auto [t, igd] = run_once_3obj(prob, ref, base_seed + static_cast<uint64_t>(r));
        std::cout << "ea-cpp,LZ09F6," << (r + 1) << ","
                  << std::fixed << std::setprecision(3) << t << ","
                  << std::setprecision(6) << igd << "\n";
        std::cout.flush();
    }
}

int main(int argc, char* argv[]) {
    int      num_runs  = (argc > 1) ? std::stoi(argv[1])   : NUM_RUNS;
    uint64_t base_seed = (argc > 2) ? std::stoull(argv[2]) : BASE_SEED;

    std::cerr << "=== ea-cpp LZ09 Benchmark ===\n";
    std::cerr << "NSGA-II:  pop=" << POP_SIZE << "  max_evals=" << MAX_EVALS
              << "  runs=" << num_runs << "  seed=" << base_seed << "\n";
    {
        std::vector<ea::ReferencePoint> pts;
        ea::generate_reference_points_das_dennis(pts, N_OBJ_F6, DIVISIONS);
        int pop_size = ea::compute_population_size(static_cast<int>(pts.size()));
        std::cerr << "NSGA-III (F6): divisions=" << DIVISIONS
                  << "  ref_pts=" << pts.size()
                  << "  pop_size=" << pop_size << "\n";
    }

    std::cout << "system,problem,run,time_ms,igd\n";

    benchmark_2obj<ea::LZ09F1>("LZ09F1", 21, num_runs, base_seed);  // ptype=21: f2=1-sqrt(f1)
    benchmark_2obj<ea::LZ09F2>("LZ09F2", 21, num_runs, base_seed);
    benchmark_2obj<ea::LZ09F3>("LZ09F3", 21, num_runs, base_seed);
    benchmark_2obj<ea::LZ09F4>("LZ09F4", 21, num_runs, base_seed);
    benchmark_2obj<ea::LZ09F5>("LZ09F5", 21, num_runs, base_seed);
    benchmark_lz09f6(num_runs, base_seed);                           // ptype=31: quarter sphere
    benchmark_2obj<ea::LZ09F7>("LZ09F7", 21, num_runs, base_seed);
    benchmark_2obj<ea::LZ09F8>("LZ09F8", 21, num_runs, base_seed);
    benchmark_2obj<ea::LZ09F9>("LZ09F9", 22, num_runs, base_seed);  // ptype=22: f2=1-f1²

    std::cerr << "Done.\n";
    return 0;
}
