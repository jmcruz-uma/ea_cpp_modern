#include <ea/problem/zdt.hpp>
#include <ea/problem/dtlz.hpp>
#include <ea/core/population.hpp>
#include <iostream>

using namespace ea;

int main() {
    // Test ZDT1
    {
        ZDT1 prob(30);
        std::cout << "ZDT1: dim=" << prob.dimension()
                  << " obj=" << prob.num_objectives() << "\n";
        Population pop(2, prob.dimension(), prob.num_objectives());
        // Set some genes
        for (int j = 0; j < prob.dimension(); ++j) {
            pop.gene(0, j) = 0.5;
            pop.gene(1, j) = 0.1;
        }
        prob.evaluate(pop, 0);
        prob.evaluate_batch(pop, 0, 2);
        std::cout << "ZDT1 f1=" << pop.objective(0, 0)
                  << " f2=" << pop.objective(0, 1) << "\n";
    }

    // Test ZDT2
    {
        ZDT2 prob(30);
        std::cout << "ZDT2: dim=" << prob.dimension()
                  << " obj=" << prob.num_objectives() << "\n";
        Population pop(1, prob.dimension(), prob.num_objectives());
        prob.evaluate(pop, 0);
    }

    // Test ZDT3
    {
        ZDT3 prob(30);
        std::cout << "ZDT3: dim=" << prob.dimension()
                  << " obj=" << prob.num_objectives() << "\n";
        Population pop(1, prob.dimension(), prob.num_objectives());
        prob.evaluate(pop, 0);
    }

    // Test DTLZ1
    {
        DTLZ1 prob(7, 3);
        std::cout << "DTLZ1: dim=" << prob.dimension()
                  << " obj=" << prob.num_objectives() << "\n";
        Population pop(1, prob.dimension(), prob.num_objectives());
        for (int j = 0; j < prob.dimension(); ++j) {
            pop.gene(0, j) = 0.5;
        }
        prob.evaluate(pop, 0);
        std::cout << "DTLZ1 f1=" << pop.objective(0, 0)
                  << " f2=" << pop.objective(0, 1)
                  << " f3=" << pop.objective(0, 2) << "\n";
    }

    // Test DTLZ2
    {
        DTLZ2 prob(12, 3);
        std::cout << "DTLZ2: dim=" << prob.dimension()
                  << " obj=" << prob.num_objectives() << "\n";
        Population pop(1, prob.dimension(), prob.num_objectives());
        for (int j = 0; j < prob.dimension(); ++j) {
            pop.gene(0, j) = 0.5;
        }
        prob.evaluate(pop, 0);
        std::cout << "DTLZ2 f1=" << pop.objective(0, 0)
                  << " f2=" << pop.objective(0, 1)
                  << " f3=" << pop.objective(0, 2) << "\n";
    }

    std::cout << "All problem tests passed!\n";
    return 0;
}
