// Functional test for NSGA-III — runs a few iterations on a small problem
#include <ea/algorithm/nsga3.hpp>
#include <ea/operator/crossover/sbx.hpp>
#include <ea/operator/mutation/polynomial.hpp>
#include <ea/core/population.hpp>
#include <ea/util/reference_point.hpp>
#include <cmath>
#include <iostream>
#include <cassert>

// Simple 3-objective test problem (scaled DTLZ1)
struct TestProblem3D {
    static constexpr int num_objectives() { return 3; }
    static constexpr int dimension() { return 7; }

    std::vector<double> lower_bounds() const { return std::vector<double>(7, 0.0); }
    std::vector<double> upper_bounds() const { return std::vector<double>(7, 1.0); }

    // Callable operator for NSGA-III algorithm
    void operator()(ea::Population<>& pop, int idx) const {
        const int k = 5;
        double g = 0.0;
        for (int i = num_objectives() - 1; i < dimension(); ++i) {
            double xi = pop.gene(idx, i);
            g += (xi - 0.5) * (xi - 0.5) - std::cos(20.0 * M_PI * (xi - 0.5));
        }
        g = 100.0 * (k + g);

        double f1 = 0.5 * (1.0 + g);
        for (int i = 0; i < num_objectives() - 1; ++i) {
            f1 *= pop.gene(idx, i);
        }
        pop.objective(idx, 0) = f1;

        for (int m = 1; m < num_objectives(); ++m) {
            double fm = 0.5 * (1.0 + g);
            for (int i = 0; i < num_objectives() - m - 1; ++i) {
                fm *= pop.gene(idx, i);
            }
            if (num_objectives() - m - 1 >= 0) {
                fm *= (1.0 - pop.gene(idx, num_objectives() - m - 1));
            }
            pop.objective(idx, m) = fm;
        }
    }
};

int main() {
    std::cout << "=== NSGA-III Functional Test ===\n";

    // Test 1: Reference point generation
    {
        std::vector<ea::ReferencePoint> ref_points;
        ea::generate_reference_points_das_dennis(ref_points, 3, 4);
        std::cout << "Reference points for 3 objectives, 4 divisions: " << ref_points.size() << "\n";
        assert(ref_points.size() == 15);  // C(4+3-1, 3-1) = C(6,2) = 15

        // Verify they sum to 1
        for (const auto& rp : ref_points) {
            double sum = 0.0;
            for (double v : rp.position) sum += v;
            assert(std::abs(sum - 1.0) < 1e-6);
        }
        std::cout << "  All reference points sum to 1.0 ✓\n";
    }

    // Test 2: Population<> size computation
    {
        assert(ea::compute_population_size(15) == 16);
        assert(ea::compute_population_size(91) == 92);
        assert(ea::compute_population_size(92) == 92);
        std::cout << "Population size computation ✓\n";
    }

    // Test 3: Run NSGA-III for a few iterations
    {
        ea::NSGAIII<ea::SBXCrossover, ea::PolynomialMutation> nsga3;
        nsga3.divisions = 4;  // Small for fast test
        nsga3.max_evals = 200;

        // Create population with correct size
        int pop_size = ea::compute_population_size(15);  // 16
        ea::Population<> pop(pop_size, 7, 3);

        // Initialize with random values
        auto& rng = ea::Random::instance();
        for (int i = 0; i < pop_size; ++i) {
            for (int j = 0; j < 7; ++j) {
                pop.gene(i, j) = rng.uniform(0.0, 1.0);
            }
        }

        TestProblem3D problem;
        nsga3.run(pop, problem);

        std::cout << "NSGA-III completed " << nsga3.max_evals << " evaluations\n";
        std::cout << "Final population size: " << pop.pop_size << "\n";

        // Verify population is evaluated
        for (int i = 0; i < pop.pop_size; ++i) {
            assert(pop.evaluated(i));
        }
        std::cout << "All individuals evaluated ✓\n";
    }

    // Test 4: Two-layer reference points
    {
        std::vector<ea::ReferencePoint> ref_points;
        ea::generate_reference_points_two_layer(ref_points, 5, 3, 2);
        std::cout << "Two-layer reference points for 5 objectives: " << ref_points.size() << "\n";

        // Verify all sum to 1
        for (const auto& rp : ref_points) {
            double sum = 0.0;
            for (double v : rp.position) sum += v;
            assert(std::abs(sum - 1.0) < 1e-6);
        }
        std::cout << "All two-layer reference points sum to 1.0 ✓\n";
    }

    // Test 5: ReferencePoint operations
    {
        ea::ReferencePoint rp(3);
        rp.position = {0.3, 0.4, 0.3};
        rp.member_count = 0;

        rp.add_potential_member(1, 0.5);
        rp.add_potential_member(2, 0.3);
        rp.add_potential_member(3, 0.7);

        rp.sort_potential_members();

        assert(rp.potential_members[0].distance == 0.3);
        assert(rp.potential_members[0].individual_index == 2);

        int closest = rp.remove_closest_member();
        assert(closest == 2);
        assert(rp.potential_members.size() == 2);

        std::cout << "ReferencePoint member management ✓\n";
    }

    std::cout << "\nAll NSGA-III tests passed!\n";
    return 0;
}
