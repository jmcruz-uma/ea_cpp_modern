#include <ea/operator/selection/nary_tournament.hpp>
#include <ea/operator/selection/random.hpp>
#include <ea/core/population.hpp>
#include <vector>
#include <iostream>

using namespace ea;

int main() {
    // Test NaryTournamentSelection
    {
        NaryTournamentSelection sel;
        sel.tournament_size = 3;
        Population<> pop(10, 2, 2);
        std::vector<int> mating_pool;
        std::vector<int> ranks(10, 0);
        std::vector<double> cd(10, 0.0);
        sel.select(pop, mating_pool, ranks, cd);
        std::cout << "NaryTournamentSelection: " << mating_pool.size() << " indices\n";
        if (mating_pool.size() != 20) return 1;
    }

    // Test RandomSelection
    {
        RandomSelection sel;
        Population<> pop(10, 2, 2);
        std::vector<int> mating_pool;
        sel.select(pop, mating_pool);
        std::cout << "RandomSelection: " << mating_pool.size() << " indices\n";
        if (mating_pool.size() != 20) return 1;
    }

    std::cout << "Selection tests passed!\n";
    return 0;
}
