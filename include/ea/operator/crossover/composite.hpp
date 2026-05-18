#pragma once
/// @file composite.hpp
/// @brief Composite Crossover — applies multiple crossovers with configurable probability.
/// Reference: jMetal CompositeCrossover.java
///
/// A meta-crossover that randomly selects one of several child crossover operators
/// and applies it. Useful for ensemble strategies or when different parts of the
/// solution space benefit from different crossover operators.

#include <ea/core/encoding.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <functional>
#include <optional>
#include <variant>
#include <vector>

namespace ea {

/// Composite Crossover — applies one of several crossover operators randomly.
///
/// Stores crossover operators as type-erased function objects and randomly
/// selects one at application time. Each operator has an associated weight
/// that determines its selection probability.
struct CompositeCrossover {
    /// Type-erased crossover application function
    using CrossoverFn = std::function<void(Population&, int, int, int)>;

    std::vector<CrossoverFn> operators; ///< Child crossover operators
    std::vector<double> weights;        ///< Selection weights (proportional)

    static constexpr int arity() { return 2; }
    static constexpr Encoding encoding() { return Encoding::Composite; }

    /// Add a crossover operator with given weight.
    template <typename CX> void add_crossover(CX&& cx, double weight = 1.0) {
        operators.emplace_back(
            [cx = std::forward<CX>(cx)](Population& pop, int a, int b, int child) mutable {
                cx.apply(pop, a, b, child);
            });
        weights.push_back(weight > 0.0 ? weight : 1.0);
    }

    /// Apply one randomly selected crossover operator.
    void apply(this CompositeCrossover& self, Population& pop, int parent_a, int parent_b,
               int child_start) {
        if (self.operators.empty()) {
            // No operators registered — just copy parents
            for (int j = 0; j < pop.dim; ++j) {
                pop.gene(child_start, j) = pop.gene(parent_a, j);
                pop.gene(child_start + 1, j) = pop.gene(parent_b, j);
            }
            pop.set_evaluated(child_start, false);
            pop.set_evaluated(child_start + 1, false);
            return;
        }

        auto& rng = Random::instance();

        // Weighted random selection
        double total = 0.0;
        for (double w : self.weights)
            total += w;

        double r = rng.uniform(0.0, total);
        double accum = 0.0;
        size_t selected = 0;
        for (size_t i = 0; i < self.weights.size(); ++i) {
            accum += self.weights[i];
            if (r <= accum) {
                selected = i;
                break;
            }
        }

        self.operators[selected](pop, parent_a, parent_b, child_start);
    }
};

} // namespace ea
