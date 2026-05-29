#include <ea/algorithm/smpso.hpp>
#include <ea/algorithm/agemoea.hpp>
#include <ea/algorithm/paes.hpp>
#include <ea/algorithm/mocell.hpp>
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/operator/mutation/uniform.hpp>
#include <ea/problem/zdt.hpp>
#include <ea/util/random.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <cmath>

using namespace ea;

// ============================================================
// Compilation tests
// ============================================================

TEST(NewAlgorithms, SMPSOCompiles) {
    SMPSO smpso;
    EXPECT_EQ(smpso.name(), "SMPSO");
    EXPECT_EQ(smpso.pop_size, 100);
    EXPECT_EQ(smpso.max_evals, 25000);
    EXPECT_EQ(smpso.archive_size, 100);
    EXPECT_DOUBLE_EQ(smpso.c1_min, 1.5);
    EXPECT_DOUBLE_EQ(smpso.c1_max, 2.5);
}

TEST(NewAlgorithms, AGEMOEACompiles) {
    AGEMOEA<SBXCrossover, PolynomialMutation> agemoea;
    EXPECT_EQ(agemoea.name(), "AGE-MOEA");
    EXPECT_EQ(agemoea.pop_size, 100);
    EXPECT_EQ(agemoea.max_evals, 25000);
}

TEST(NewAlgorithms, PAESCompiles) {
    PAES<PolynomialMutation> paes;
    EXPECT_EQ(paes.name(), "PAES");
    EXPECT_EQ(paes.max_evals, 25000);
    EXPECT_EQ(paes.archive_size, 100);
    EXPECT_EQ(paes.bisections, 5);
}

TEST(NewAlgorithms, MOCellCompiles) {
    MOCell<SBXCrossover, PolynomialMutation> mocell;
    EXPECT_EQ(mocell.name(), "MOCell");
    EXPECT_EQ(mocell.pop_size, 100);
    EXPECT_EQ(mocell.max_evals, 25000);
    EXPECT_EQ(mocell.archive_size, 100);
}

// ============================================================
// Run tests — ZDT1
// ============================================================

TEST(NewAlgorithms, SMPSORunsOnZDT1) {
    Population<> pop(30, 30, 2);
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
    
    auto problem = [&zdt1](Population<>& p, int idx) { zdt1.evaluate(p, idx); };
    smpso.run(pop, problem);
    
    // Archive or population should have evaluated individuals
    bool any_evaluated = false;
    for (int i = 0; i < pop.pop_size; ++i) {
        if (pop.evaluated(i)) {
            any_evaluated = true;
            EXPECT_GE(pop.objective(i, 0), 0.0);
            EXPECT_GE(pop.objective(i, 1), 0.0);
        }
    }
    EXPECT_TRUE(any_evaluated);
}

TEST(NewAlgorithms, AGEMOEARunsOnZDT1) {
    Population<> pop(30, 30, 2);
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
    
    auto problem = [&zdt1](Population<>& p, int idx) { zdt1.evaluate(p, idx); };
    agemoea.run(pop, problem);
    
    for (int i = 0; i < pop.pop_size; ++i) {
        EXPECT_TRUE(pop.evaluated(i));
        EXPECT_GE(pop.objective(i, 0), 0.0);
        EXPECT_GE(pop.objective(i, 1), 0.0);
    }
}

TEST(NewAlgorithms, PAESRunsOnZDT1) {
    Population<> pop(1, 30, 2); // PAES uses single individual
    ZDT1 zdt1;
    
    auto& rng = Random::instance();
    for (int j = 0; j < pop.dim; ++j) {
        pop.gene(0, j) = rng.uniform();
    }
    
    PAES<PolynomialMutation> paes;
    paes.max_evals = 500;
    paes.archive_size = 30;
    paes.mutation.distribution_index = 20.0;
    
    auto problem = [&zdt1](Population<>& p, int idx) { zdt1.evaluate(p, idx); };
    paes.run(pop, problem);
    
    // Result should be in archive (copied to pop)
    EXPECT_GT(pop.pop_size, 0);
    for (int i = 0; i < pop.pop_size; ++i) {
        EXPECT_TRUE(pop.evaluated(i));
        EXPECT_GE(pop.objective(i, 0), 0.0);
        EXPECT_GE(pop.objective(i, 1), 0.0);
    }
}

TEST(NewAlgorithms, MOCellRunsOnZDT1) {
    Population<> pop(30, 30, 2);
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
    
    auto problem = [&zdt1](Population<>& p, int idx) { zdt1.evaluate(p, idx); };
    mocell.run(pop, problem);
    
    EXPECT_GT(pop.pop_size, 0);
    for (int i = 0; i < pop.pop_size; ++i) {
        EXPECT_TRUE(pop.evaluated(i));
        EXPECT_GE(pop.objective(i, 0), 0.0);
        EXPECT_GE(pop.objective(i, 1), 0.0);
    }
}

// ============================================================
// AdaptiveGridArchive tests
// ============================================================

TEST(AdaptiveGridArchive, Construction) {
    AdaptiveGridArchive archive(50, 3, 2);
    EXPECT_EQ(archive.max_size, 50);
    EXPECT_EQ(archive.bisections, 3);
    EXPECT_EQ(archive.n_obj, 2);
    EXPECT_EQ(archive.size(), 0);
    EXPECT_TRUE(archive.empty());
}

TEST(AdaptiveGridArchive, AddNonDominated) {
    AdaptiveGridArchive archive(10, 2, 2);
    
    std::vector<double> genes = {0.5, 0.5};
    std::vector<double> objs1 = {0.0, 1.0};
    std::vector<double> objs2 = {1.0, 0.0};
    
    EXPECT_TRUE(archive.add(genes, objs1));
    EXPECT_EQ(archive.size(), 1);
    
    EXPECT_TRUE(archive.add(genes, objs2));
    EXPECT_EQ(archive.size(), 2);
}

TEST(AdaptiveGridArchive, DominatedNotAdded) {
    AdaptiveGridArchive archive(10, 2, 2);
    
    std::vector<double> genes = {0.5, 0.5};
    std::vector<double> objs1 = {0.0, 0.0}; // Ideal point
    std::vector<double> objs2 = {1.0, 1.0}; // Dominated by objs1
    
    EXPECT_TRUE(archive.add(genes, objs1));
    EXPECT_FALSE(archive.add(genes, objs2));
    EXPECT_EQ(archive.size(), 1);
}

TEST(AdaptiveGridArchive, ArchiveFullReplacesCrowded) {
    AdaptiveGridArchive archive(2, 2, 2);
    
    std::vector<double> genes = {0.5, 0.5};
    std::vector<double> objs1 = {0.0, 1.0};
    std::vector<double> objs2 = {1.0, 0.0};
    std::vector<double> objs3 = {0.5, 0.5};
    
    EXPECT_TRUE(archive.add(genes, objs1));
    EXPECT_TRUE(archive.add(genes, objs2));
    EXPECT_EQ(archive.size(), 2);
    
    // objs3 is non-dominated but archive is full
    // May or may not be added depending on crowding
    archive.add(genes, objs3);
    EXPECT_LE(archive.size(), 2);
}
