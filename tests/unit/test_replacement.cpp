#include <ea/operator/replacement/mu_plus_lambda.hpp>
#include <ea/operator/replacement/mu_comma_lambda.hpp>
#include <ea/core/population.hpp>
#include <ea/problem/zdt.hpp>
#include <iostream>

using namespace ea;

int main() {
    // Test MuPlusLambdaReplacement
    {
        MuPlusLambdaReplacement repl;
        ZDT1 prob(5);  // small dimension for test
        Population<> pop(4, prob.dimension(), prob.num_objectives());

        // Fill with known values
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < prob.dimension(); ++j) {
                pop.gene(i, j) = 0.5;
            }
            prob.evaluate(pop, i);
        }

        // Simulate offspring appended (in real use, offspring would be at indices 4,5,6,7)
        // For this test, we'll just use indices [2,3] as "offspring"
        std::vector<int> offspring = {2, 3};
        auto survivors = repl.replace(pop, offspring, 2);

        std::cout << "MuPlusLambda: " << survivors.size() << " survivors\n";
        if (survivors.size() != 2) return 1;
    }

    // Test MuCommaLambdaReplacement
    {
        MuCommaLambdaReplacement repl;
        ZDT1 prob(5);
        Population<> pop(6, prob.dimension(), prob.num_objectives());

        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < prob.dimension(); ++j) {
                pop.gene(i, j) = 0.5;
            }
            prob.evaluate(pop, i);
        }

        std::vector<int> offspring = {2, 3, 4, 5};
        auto survivors = repl.replace(pop, offspring, 2);

        std::cout << "MuCommaLambda: " << survivors.size() << " survivors\n";
        if (survivors.size() != 2) return 1;
    }

    std::cout << "Replacement tests passed!\n";
    return 0;
}
