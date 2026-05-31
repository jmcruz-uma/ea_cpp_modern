/// @file nsga3_runner.cpp
/// @brief NSGA-III benchmark runner — ea_cpp_modern vs jMetal 7.4 (DTLZ1-4).
///
/// Corre NSGA-III sobre DTLZ1-4, NUM_RUNS veces por problema.
/// Salida CSV por stdout: system,problem,run,time_ms,igd
///
/// Compilar (benchmark):
///   g++-14 -std=c++23 -O3 -march=native -DNDEBUG \
///     -I include examples/nsga3_runner.cpp -o build/nsga3_runner -lm
///
/// Uso:
///   ./build/nsga3_runner [NUM_RUNS [BASE_SEED]]
///
/// Parámetros (deben coincidir con JMetalNSGA3Benchmark.java):
///   3 objectives, divisions=12 → 91 ref points → pop_size=92
///   max_evals=25000, SBX(p=0.9, η=20), PolMut(p=1/dim, η=20)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include <ea/algorithm/nsga3.hpp>
#include <ea/core/comparator.hpp>
#include <ea/indicator/igd.hpp>
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/problem/dtlz.hpp>
#include <ea/util/random.hpp>
#include <ea/util/reference_point.hpp>

// ── Parámetros del benchmark ───────────────────────────────────────────────────
static constexpr int      N_OBJ         = 3;
static constexpr int      DIVISIONS     = 12;   // 91 ref points → pop_size=92
static constexpr int      MAX_EVALS     = 25000;
static constexpr double   SBX_PROB      = 0.9;
static constexpr double   SBX_ETA       = 20.0;
static constexpr double   MUT_ETA       = 20.0; // mutation_rate auto = 1/dim
static constexpr uint64_t BASE_SEED     = 42;
static constexpr int      NUM_RUNS      = 30;
static constexpr int      REF_DIV       = 30;   // Das-Dennis divisions for reference front
// C(32,2) = 496 ≈ 500 reference points for 3 objectives

// ── Frente de referencia para DTLZ (3 objetivos) ─────────────────────────────
static std::vector<std::vector<double>> make_dtlz1_ref() {
    // Pareto front: f1+f2+f3 = 0.5, f_i >= 0
    // Generate Das-Dennis simplex points (sum=1), scale by 0.5
    std::vector<ea::ReferencePoint> pts;
    ea::generate_reference_points_das_dennis(pts, N_OBJ, REF_DIV);
    std::vector<std::vector<double>> ref;
    ref.reserve(pts.size());
    for (const auto& rp : pts)
        ref.push_back({rp.position[0] * 0.5, rp.position[1] * 0.5, rp.position[2] * 0.5});
    return ref;
}

static std::vector<std::vector<double>> make_dtlz_sphere_ref() {
    // Pareto front: f1^2+f2^2+f3^2 = 1, f_i >= 0 (quarter sphere)
    // Generate Das-Dennis simplex points, normalize each to unit sphere
    std::vector<ea::ReferencePoint> pts;
    ea::generate_reference_points_das_dennis(pts, N_OBJ, REF_DIV);
    std::vector<std::vector<double>> ref;
    ref.reserve(pts.size());
    for (const auto& rp : pts) {
        double norm = std::sqrt(rp.position[0] * rp.position[0] +
                                rp.position[1] * rp.position[1] +
                                rp.position[2] * rp.position[2]);
        if (norm < 1e-14)
            norm = 1.0;
        ref.push_back({rp.position[0] / norm, rp.position[1] / norm, rp.position[2] / norm});
    }
    return ref;
}

// ── Un run de NSGA-III — devuelve {time_ms, igd} ──────────────────────────────
template <typename Problem>
std::pair<double, double> run_once(Problem& prob,
                                   const std::vector<std::vector<double>>& ref_front,
                                   uint64_t seed) {
    ea::Random::instance().seed(seed);

    const int dim   = prob.dimension();
    const int n_obj = prob.num_objectives();

    // Derive population size from reference points (matching jMetal derivation)
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
    alg.crossover.crossover_probability = SBX_PROB;
    alg.crossover.distribution_index    = SBX_ETA;
    alg.mutation.distribution_index     = MUT_ETA;
    // mutation.mutation_rate = -1 (default) → auto 1/dim

    auto t0 = std::chrono::high_resolution_clock::now();
    alg.run(pop, [&](ea::Population<>& p, int idx) { prob.evaluate(p, idx); });
    auto t1 = std::chrono::high_resolution_clock::now();

    double time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Use non-dominated front (matching jMetal result())
    auto fronts = ea::fast_non_dominated_sort(pop);
    std::vector<std::vector<double>> approx;
    approx.reserve(fronts[0].size());
    for (int idx : fronts[0])
        approx.push_back({pop.objective(idx, 0), pop.objective(idx, 1), pop.objective(idx, 2)});

    double igd_val = ea::igd(pop, fronts[0], ref_front, 1.0);

    return {time_ms, igd_val};
}

// ── Runner por problema ───────────────────────────────────────────────────────
template <typename Problem>
void benchmark_problem(const std::string& name, int dim,
                       const std::vector<std::vector<double>>& ref_front,
                       int num_runs, uint64_t base_seed) {
    std::cerr << "  " << name << " (dim=" << dim << ")...\n";
    Problem prob(dim, N_OBJ);
    for (int r = 0; r < num_runs; ++r) {
        auto [t, igd] = run_once(prob, ref_front, base_seed + static_cast<uint64_t>(r));
        std::cout << "ea-cpp," << name << "," << (r + 1) << ","
                  << std::fixed << std::setprecision(3) << t << ","
                  << std::setprecision(6) << igd << "\n";
        std::cout.flush();
    }
}

int main(int argc, char* argv[]) {
    int      num_runs  = (argc > 1) ? std::stoi(argv[1])   : NUM_RUNS;
    uint64_t base_seed = (argc > 2) ? std::stoull(argv[2]) : BASE_SEED;

    std::cerr << "=== ea-cpp NSGA-III Benchmark ===\n";
    int n_ref_pts;
    {
        std::vector<ea::ReferencePoint> pts;
        ea::generate_reference_points_das_dennis(pts, N_OBJ, DIVISIONS);
        n_ref_pts = static_cast<int>(pts.size());
    }
    int pop_size = ea::compute_population_size(n_ref_pts);
    std::cerr << "n_obj=" << N_OBJ << "  divisions=" << DIVISIONS
              << "  ref_pts=" << n_ref_pts << "  pop_size=" << pop_size
              << "  max_evals=" << MAX_EVALS << "\n";

    std::cout << "system,problem,run,time_ms,igd\n";

    auto ref1  = make_dtlz1_ref();
    auto ref24 = make_dtlz_sphere_ref();

    benchmark_problem<ea::DTLZ1>("DTLZ1", 7,  ref1,  num_runs, base_seed);
    benchmark_problem<ea::DTLZ2>("DTLZ2", 12, ref24, num_runs, base_seed);
    benchmark_problem<ea::DTLZ3>("DTLZ3", 12, ref24, num_runs, base_seed);
    benchmark_problem<ea::DTLZ4>("DTLZ4", 12, ref24, num_runs, base_seed);

    return 0;
}
