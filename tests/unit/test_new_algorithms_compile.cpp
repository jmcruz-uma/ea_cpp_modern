#include <ea/algorithm/smpso.hpp>
#include <ea/algorithm/agemoea.hpp>
#include <ea/algorithm/paes.hpp>
#include <ea/algorithm/mocell.hpp>
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/problem/zdt.hpp>
#include <ea/util/adaptive_grid.hpp>
#include <vector>
#include <cmath>

using namespace ea;

// ============================================================
// Standalone compilation test — does not link with gtest
// Compile: g++-14 -std=c++23 -I include -c tests/unit/test_new_algorithms_compile.cpp
// ============================================================

static void test_smpso_compiles() {
    SMPSO smpso;
    assert(std::string(smpso.name()) == "SMPSO");
    assert(smpso.pop_size == 100);
    assert(smpso.max_evals == 25000);
    assert(smpso.archive_size == 100);
    assert(smpso.c1_min == 1.5);
    assert(smpso.c1_max == 2.5);
    (void)smpso;
}

static void test_agemoea_compiles() {
    AGEMOEA<SBXCrossover, PolynomialMutation> agemoea;
    assert(std::string(agemoea.name()) == "AGE-MOEA");
    assert(agemoea.pop_size == 100);
    assert(agemoea.max_evals == 25000);
    (void)agemoea;
}

static void test_paes_compiles() {
    PAES<PolynomialMutation> paes;
    assert(std::string(paes.name()) == "PAES");
    assert(paes.max_evals == 25000);
    assert(paes.archive_size == 100);
    assert(paes.bisections == 5);
    (void)paes;
}

static void test_mocell_compiles() {
    MOCell<SBXCrossover, PolynomialMutation> mocell;
    assert(std::string(mocell.name()) == "MOCell");
    assert(mocell.pop_size == 100);
    assert(mocell.max_evals == 25000);
    assert(mocell.archive_size == 100);
    (void)mocell;
}

static void test_archive_compiles() {
    AdaptiveGridArchive archive(50, 3, 2);
    assert(archive.max_size == 50);
    assert(archive.bisections == 3);
    assert(archive.n_obj == 2);
    assert(archive.size() == 0);
    assert(archive.empty());
    
    std::vector<double> genes = {0.5, 0.5};
    std::vector<double> objs1 = {0.0, 1.0};
    std::vector<double> objs2 = {1.0, 0.0};
    
    assert(archive.add(genes, objs1));
    assert(archive.size() == 1);
    
    assert(archive.add(genes, objs2));
    assert(archive.size() == 2);
    
    // Dominated should not be added
    std::vector<double> objs3 = {0.5, 0.5};
    archive.add(genes, objs3); // dominated by objs1 and objs2 combined? no, but between
    (void)archive;
}

static void test_smpso_runs() {
    Population pop(30, 30, 2);
    ZDT1 zdt1;
    
    auto& rng = Random::instance();
    for (int i = 0; i < pop.pop_size; ++i) {
        for (int j = 0; j < pop.dim; ++j) {
            pop.gene(i, j) = rng.uniform();
        }
    }
    
    SMPSO smpso;
    smpso.pop_size = 30;
    smpso.max_evals = 500;
    smpso.archive_size = 30;
    
    auto problem = [&zdt1](Population& p, int idx) { zdt1.evaluate(p, idx); };
    smpso.run(pop, problem);
    
    bool any_evaluated = false;
    for (int i = 0; i < pop.pop_size; ++i) {
        if (pop.evaluated(i)) {
            any_evaluated = true;
            assert(pop.objective(i, 0) >= 0.0);
            assert(pop.objective(i, 1) >= 0.0);
        }
    }
    assert(any_evaluated);
}

static void test_agemoea_runs() {
    Population pop(30, 30, 2);
    ZDT1 zdt1;
    
    auto& rng = Random::instance();
    for (int i = 0; i < pop.pop_size; ++i) {
        for (int j = 0; j < pop.dim; ++j) {
            pop.gene(i, j) = rng.uniform();
        }
    }
    
    AGEMOEA<SBXCrossover, PolynomialMutation> agemoea;
    agemoea.pop_size = 30;
    agemoea.max_evals = 500;
    agemoea.crossover.distribution_index = 20.0;
    agemoea.mutation.distribution_index = 20.0;
    
    auto problem = [&zdt1](Population& p, int idx) { zdt1.evaluate(p, idx); };
    agemoea.run(pop, problem);
    
    for (int i = 0; i < pop.pop_size; ++i) {
        assert(pop.evaluated(i));
        assert(pop.objective(i, 0) >= 0.0);
        assert(pop.objective(i, 1) >= 0.0);
    }
}

static void test_paes_runs() {
    Population pop(1, 30, 2);
    ZDT1 zdt1;
    
    auto& rng = Random::instance();
    for (int j = 0; j < pop.dim; ++j) {
        pop.gene(0, j) = rng.uniform();
    }
    
    PAES<PolynomialMutation> paes;
    paes.max_evals = 500;
    paes.archive_size = 30;
    paes.mutation.distribution_index = 20.0;
    
    auto problem = [&zdt1](Population& p, int idx) { zdt1.evaluate(p, idx); };
    paes.run(pop, problem);
    
    assert(pop.pop_size > 0);
    for (int i = 0; i < pop.pop_size; ++i) {
        assert(pop.evaluated(i));
        assert(pop.objective(i, 0) >= 0.0);
        assert(pop.objective(i, 1) >= 0.0);
    }
}

static void test_mocell_runs() {
    Population pop(30, 30, 2);
    ZDT1 zdt1;
    
    auto& rng = Random::instance();
    for (int i = 0; i < pop.pop_size; ++i) {
        for (int j = 0; j < pop.dim; ++j) {
            pop.gene(i, j) = rng.uniform();
        }
    }
    
    MOCell<SBXCrossover, PolynomialMutation> mocell;
    mocell.pop_size = 30;
    mocell.max_evals = 500;
    mocell.archive_size = 30;
    mocell.crossover.distribution_index = 20.0;
    mocell.mutation.distribution_index = 20.0;
    
    auto problem = [&zdt1](Population& p, int idx) { zdt1.evaluate(p, idx); };
    mocell.run(pop, problem);
    
    assert(pop.pop_size > 0);
    for (int i = 0; i < pop.pop_size; ++i) {
        assert(pop.evaluated(i));
        assert(pop.objective(i, 0) >= 0.0);
        assert(pop.objective(i, 1) >= 0.0);
    }
}

int main() {
    test_smpso_compiles();
    test_agemoea_compiles();
    test_paes_compiles();
    test_mocell_compiles();
    test_archive_compiles();
    
    test_smpso_runs();
    test_agemoea_runs();
    test_paes_runs();
    test_mocell_runs();
    
    return 0;
}
