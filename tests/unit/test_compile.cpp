// Quick compilation test — include all headers
#include <ea/ea.hpp>

int main() {
    // Create a population
    ea::Population pop(100, 30, 2);

    // Create operators
    ea::SBXCrossover sbx;
    ea::PolynomialMutation pmut;
    ea::BinaryTournamentSelection bts;

    // Create ZDT1 problem
    ea::ZDT1 zdt1;

    // Create NSGA-II (template takes CX and MT)
    ea::NSGAII<ea::SBXCrossover, ea::PolynomialMutation> nsga;

    return 0;
}