#include <gtest/gtest.h>
#include <ea/operator/selection/nary_tournament.hpp>
#include <ea/core/population.hpp>

using namespace ea;

TEST(NaryTournamentSelection, BasicSelect) {
    NaryTournamentSelection sel;
    sel.tournament_size = 2;
    Population pop(10, 2, 2);
    std::vector<int> mating_pool;
    std::vector<int> ranks(10, 0);
    std::vector<double> cd(10, 0.0);
    sel.select(pop, mating_pool, ranks, cd);
    EXPECT_EQ(mating_pool.size(), 20u);
    for (int idx : mating_pool) {
        EXPECT_GE(idx, 0);
        EXPECT_LT(idx, 10);
    }
}

TEST(NaryTournamentSelection, LargerTournament) {
    NaryTournamentSelection sel;
    sel.tournament_size = 4;
    Population pop(20, 2, 2);
    std::vector<int> mating_pool;
    std::vector<int> ranks(20, 0);
    std::vector<double> cd(20, 0.0);
    sel.select(pop, mating_pool, ranks, cd);
    EXPECT_EQ(mating_pool.size(), 40u);
}
