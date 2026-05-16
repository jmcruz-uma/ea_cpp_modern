#include <ea/algorithm/moead.hpp>
#include <ea/util/aggregation.hpp>
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/problem/zdt.hpp>
#include <ea/util/random.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <cmath>

using namespace ea;

// ============================================================
// Compilation tests
// ============================================================

TEST(MOEAD, CompilesWithSBXAndPolynomial) {
    MOEAD<SBXCrossover, PolynomialMutation> moead;
    EXPECT_EQ(moead.name(), "MOEA/D");
    EXPECT_EQ(moead.pop_size, 100);
    EXPECT_EQ(moead.max_evals, 25000);
    EXPECT_EQ(moead.neighbor_size, -1);  // default: auto
    EXPECT_EQ(moead.max_replaced_solutions, 2);
    EXPECT_DOUBLE_EQ(moead.neighborhood_prob, 0.9);
    EXPECT_EQ(moead.aggregation_type, AggregationType::Tchebycheff);
    EXPECT_FALSE(moead.normalize);
}

// ============================================================
// Aggregation function tests
// ============================================================

TEST(Aggregation, WeightedSumBasic) {
    WeightedSum ws;
    double obj[] = {2.0, 3.0};
    double weight[] = {0.5, 0.5};
    double ideal[] = {0.0, 0.0};
    double nadir[] = {5.0, 5.0};
    
    double result = ws.compute(obj, weight, ideal, nadir, 2);
    // (2-0)*0.5 + (3-0)*0.5 = 1.0 + 1.5 = 2.5
    EXPECT_DOUBLE_EQ(result, 2.5);
}

TEST(Aggregation, WeightedSumNormalized) {
    WeightedSum ws;
    ws.normalize = true;
    double obj[] = {2.0, 3.0};
    double weight[] = {0.5, 0.5};
    double ideal[] = {0.0, 0.0};
    double nadir[] = {5.0, 5.0};
    
    double result = ws.compute(obj, weight, ideal, nadir, 2);
    // (2/5)*0.5 + (3/5)*0.5 = 0.2 + 0.3 = 0.5
    EXPECT_DOUBLE_EQ(result, 0.5);
}

TEST(Aggregation, TchebycheffBasic) {
    Tchebycheff tch;
    double obj[] = {2.0, 3.0};
    double weight[] = {0.5, 0.5};
    double ideal[] = {0.0, 0.0};
    double nadir[] = {5.0, 5.0};
    
    double result = tch.compute(obj, weight, ideal, nadir, 2);
    // max(|2|*0.5, |3|*0.5) = max(1.0, 1.5) = 1.5
    EXPECT_DOUBLE_EQ(result, 1.5);
}

TEST(Aggregation, TchebycheffZeroWeight) {
    Tchebycheff tch;
    double obj[] = {2.0, 3.0};
    double weight[] = {0.0, 1.0};
    double ideal[] = {0.0, 0.0};
    double nadir[] = {5.0, 5.0};
    
    double result = tch.compute(obj, weight, ideal, nadir, 2);
    // max(2*0.0001, 3*1.0) = max(0.0002, 3.0) = 3.0
    EXPECT_DOUBLE_EQ(result, 3.0);
}

TEST(Aggregation, PBI_Basic) {
    PenaltyBoundaryIntersection pbi;
    pbi.theta = 5.0;
    double obj[] = {2.0, 3.0};
    double weight[] = {0.5, 0.5};
    double ideal[] = {0.0, 0.0};
    double nadir[] = {5.0, 5.0};
    
    double result = pbi.compute(obj, weight, ideal, nadir, 2);
    // d1 + theta * d2 (should be positive)
    EXPECT_GT(result, 0.0);
}

TEST(Aggregation, AggregationFunctionDispatch) {
    AggregationFunction aggr;
    aggr.type = AggregationType::WeightedSum;
    
    double obj[] = {2.0, 3.0};
    double weight[] = {0.5, 0.5};
    double ideal[] = {0.0, 0.0};
    double nadir[] = {5.0, 5.0};
    
    double result = aggr.compute(obj, weight, ideal, nadir, 2);
    EXPECT_DOUBLE_EQ(result, 2.5);
}

// ============================================================
// MOEA/D run tests
// ============================================================

TEST(MOEAD, RunsOnZDT1WithTchebycheff) {
    Population pop(50, 30, 2);
    ZDT1 zdt1;
    
    auto& rng = Random::instance();
    for (int i = 0; i < pop.pop_size; ++i) {
        for (int j = 0; j < pop.dim; ++j) {
            pop.gene(i, j) = rng.uniform();
        }
    }
    
    MOEAD<SBXCrossover, PolynomialMutation> moead;
    moead.pop_size = 50;
    moead.max_evals = 1000;
    moead.neighbor_size = 5;
    moead.aggregation_type = AggregationType::Tchebycheff;
    moead.crossover.distribution_index = 20.0;
    moead.mutation.distribution_index = 20.0;
    
    moead.run(pop, zdt1);
    
    // Should have performed evaluations
    EXPECT_GE(moead.evals_, 50);  // at least initial population
    EXPECT_LE(moead.evals_, 1000); // at most max_evals
    
    // Population should have objectives set
    for (int i = 0; i < pop.pop_size; ++i) {
        EXPECT_TRUE(pop.evaluated(i));
        EXPECT_GE(pop.objective(i, 0), 0.0);
        EXPECT_GE(pop.objective(i, 1), 0.0);
    }
}

TEST(MOEAD, RunsOnZDT1WithWeightedSum) {
    Population pop(50, 30, 2);
    ZDT1 zdt1;
    
    auto& rng = Random::instance();
    for (int i = 0; i < pop.pop_size; ++i) {
        for (int j = 0; j < pop.dim; ++j) {
            pop.gene(i, j) = rng.uniform();
        }
    }
    
    MOEAD<SBXCrossover, PolynomialMutation> moead;
    moead.pop_size = 50;
    moead.max_evals = 1000;
    moead.neighbor_size = 5;
    moead.aggregation_type = AggregationType::WeightedSum;
    moead.crossover.distribution_index = 20.0;
    moead.mutation.distribution_index = 20.0;
    
    moead.run(pop, zdt1);
    
    EXPECT_GE(moead.evals_, 50);
    EXPECT_LE(moead.evals_, 1000);
    
    for (int i = 0; i < pop.pop_size; ++i) {
        EXPECT_TRUE(pop.evaluated(i));
    }
}

TEST(MOEAD, WeightVectorsInitializedCorrectly) {
    // 2-objective case: weights should be uniformly distributed
    Population pop(10, 5, 2);
    ZDT1 zdt1;
    
    auto& rng = Random::instance();
    for (int i = 0; i < pop.pop_size; ++i) {
        for (int j = 0; j < pop.dim; ++j) {
            pop.gene(i, j) = rng.uniform();
        }
    }
    
    MOEAD<SBXCrossover, PolynomialMutation> moead;
    moead.pop_size = 10;
    moead.max_evals = 10;  // Just init
    moead.neighbor_size = 2;
    
    moead.run(pop, zdt1);
    
    // Check weights sum to ~1
    for (int i = 0; i < 10; ++i) {
        double sum = moead.lambda_[i][0] + moead.lambda_[i][1];
        EXPECT_NEAR(sum, 1.0, 1e-10);
    }
    
    // First weight should be 0, last should be 1
    EXPECT_DOUBLE_EQ(moead.lambda_[0][0], 0.0);
    EXPECT_DOUBLE_EQ(moead.lambda_[0][1], 1.0);
    EXPECT_DOUBLE_EQ(moead.lambda_[9][0], 1.0);
    EXPECT_DOUBLE_EQ(moead.lambda_[9][1], 0.0);
}

TEST(MOEAD, NeighborhoodInitializedCorrectly) {
    Population pop(20, 5, 2);
    ZDT1 zdt1;
    
    auto& rng = Random::instance();
    for (int i = 0; i < pop.pop_size; ++i) {
        for (int j = 0; j < pop.dim; ++j) {
            pop.gene(i, j) = rng.uniform();
        }
    }
    
    MOEAD<SBXCrossover, PolynomialMutation> moead;
    moead.pop_size = 20;
    moead.max_evals = 20;  // Just init
    moead.neighbor_size = 3;
    
    moead.run(pop, zdt1);
    
    // Each subproblem should have T neighbors
    EXPECT_EQ(moead.neighborhood_.size(), 20u);
    for (int i = 0; i < 20; ++i) {
        EXPECT_EQ(moead.neighborhood_[i].size(), 3u);
        // First neighbor is always self (distance = 0)
        EXPECT_EQ(moead.neighborhood_[i][0], i);
    }
}

TEST(MOEAD, IdealPointUpdatesCorrectly) {
    Population pop(10, 5, 2);
    ZDT1 zdt1;
    
    auto& rng = Random::instance();
    for (int i = 0; i < pop.pop_size; ++i) {
        for (int j = 0; j < pop.dim; ++j) {
            pop.gene(i, j) = rng.uniform();
        }
    }
    
    MOEAD<SBXCrossover, PolynomialMutation> moead;
    moead.pop_size = 10;
    moead.max_evals = 10;  // Just init
    moead.neighbor_size = 2;
    
    moead.run(pop, zdt1);
    
    // Ideal point should be <= all objectives
    for (int o = 0; o < 2; ++o) {
        for (int i = 0; i < 10; ++i) {
            EXPECT_LE(moead.ideal_point_[o], pop.objective(i, o));
        }
    }
}
