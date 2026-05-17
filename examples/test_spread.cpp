/// @file test_spread.cpp
/// @brief Test for the Spread (Δ) indicator.
///
/// Tests:
///   1. Perfect distribution on y = 1 - x → Spread ≈ 0
///   2. Clustered front → Spread > 0

#include <ea/indicator/spread.hpp>
#include <vector>
#include <cmath>
#include <iostream>

using namespace ea;

int main() {
    std::cout << "=== Spread (Δ) Indicator Tests ===" << std::endl;
    std::cout << std::endl;

    // --- Test 1: Perfect distribution on y = 1 - x ---
    std::cout << "Test 1: Perfect distribution on y = 1 - x" << std::endl;
    
    std::vector<std::vector<double>> perfect_front;
    std::vector<std::vector<double>> reference_front;
    
    const int N = 11;  // 11 points from x=0 to x=1
    for (int i = 0; i < N; ++i) {
        double x = i / 10.0;
        double y = 1.0 - x;
        perfect_front.push_back({x, y});
        reference_front.push_back({x, y});
    }
    
    double spread_perfect = spread(perfect_front, reference_front);
    std::cout << "  Spread (perfect front, same as reference): " << spread_perfect << std::endl;
    
    // With perfect distribution matching reference extremes, Spread should be very close to 0
    if (std::abs(spread_perfect) < 1e-10) {
        std::cout << "  ✓ PASS: Spread ≈ 0 for perfect distribution" << std::endl;
    } else {
        std::cout << "  ✗ FAIL: Expected Spread ≈ 0, got " << spread_perfect << std::endl;
        return 1;
    }
    std::cout << std::endl;

    // --- Test 2: Clustered front (solutions bunched together) ---
    std::cout << "Test 2: Clustered front" << std::endl;
    
    std::vector<std::vector<double>> clustered_front;
    // 8 points clustered near (0.5, 0.5), 1 near each extreme
    clustered_front.push_back({0.0, 1.0});   // extreme
    clustered_front.push_back({0.48, 0.52});
    clustered_front.push_back({0.49, 0.51});
    clustered_front.push_back({0.50, 0.50});
    clustered_front.push_back({0.50, 0.50});
    clustered_front.push_back({0.51, 0.49});
    clustered_front.push_back({0.52, 0.48});
    clustered_front.push_back({0.53, 0.47});
    clustered_front.push_back({0.54, 0.46});
    clustered_front.push_back({1.0, 0.0});   // extreme
    
    // Reference: evenly distributed 11 points
    std::vector<std::vector<double>> ref_for_clustered;
    for (int i = 0; i < N; ++i) {
        double x = i / 10.0;
        ref_for_clustered.push_back({x, 1.0 - x});
    }
    
    double spread_clustered = spread(clustered_front, ref_for_clustered);
    std::cout << "  Spread (clustered front): " << spread_clustered << std::endl;
    
    // Clustered front should have higher Spread than perfect
    if (spread_clustered > spread_perfect) {
        std::cout << "  ✓ PASS: Clustered front has higher Spread (" << spread_clustered 
                  << " > " << spread_perfect << ")" << std::endl;
    } else {
        std::cout << "  ✗ FAIL: Expected clustered Spread > " << spread_perfect 
                  << ", got " << spread_clustered << std::endl;
        return 1;
    }
    std::cout << std::endl;

    // --- Test 3: SpreadIndicator struct ---
    std::cout << "Test 3: SpreadIndicator struct" << std::endl;
    
    SpreadIndicator indicator;
    double spread_via_struct = indicator.compute(perfect_front, reference_front);
    std::cout << "  Spread via SpreadIndicator: " << spread_via_struct << std::endl;
    
    if (std::abs(spread_via_struct - spread_perfect) < 1e-15) {
        std::cout << "  ✓ PASS: SpreadIndicator matches free function" << std::endl;
    } else {
        std::cout << "  ✗ FAIL: SpreadIndicator result mismatch" << std::endl;
        return 1;
    }
    std::cout << std::endl;

    // --- Test 4: Edge cases ---
    std::cout << "Test 4: Edge cases" << std::endl;
    
    // Empty front
    std::vector<std::vector<double>> empty_front;
    double spread_empty = spread(empty_front, reference_front);
    std::cout << "  Spread (empty front): " << spread_empty << " (expected 0)" << std::endl;
    if (spread_empty != 0.0) {
        std::cout << "  ✗ FAIL: Expected 0 for empty front" << std::endl;
        return 1;
    }
    
    // Single point front
    std::vector<std::vector<double>> single_front = {{0.5, 0.5}};
    double spread_single = spread(single_front, reference_front);
    std::cout << "  Spread (single point): " << spread_single << " (expected 0)" << std::endl;
    if (spread_single != 0.0) {
        std::cout << "  ✗ FAIL: Expected 0 for single-point front" << std::endl;
        return 1;
    }
    
    // Empty reference
    std::vector<std::vector<double>> empty_ref;
    double spread_empty_ref = spread(perfect_front, empty_ref);
    std::cout << "  Spread (empty reference): " << spread_empty_ref << " (expected 0)" << std::endl;
    if (spread_empty_ref != 0.0) {
        std::cout << "  ✗ FAIL: Expected 0 for empty reference" << std::endl;
        return 1;
    }
    
    std::cout << "  ✓ PASS: All edge cases handled" << std::endl;
    std::cout << std::endl;

    std::cout << "=== All Tests Passed ===" << std::endl;
    return 0;
}
