#pragma once
/// @file encoding.hpp
/// @brief Solution encoding types for evolutionary algorithms.

#include <cstdint>
#include <string_view>
#include <type_traits>

namespace ea {

/// Encoding type for solution representations.
/// Used at compile time for Concept constraints and at runtime for variant dispatch.
enum class Encoding : uint8_t {
    Real,        ///< Continuous real-valued optimization
    Binary,      ///< Binary string optimization
    Permutation, ///< Permutation-based optimization (TSP, scheduling)
    Integer,     ///< Integer-valued optimization
    Composite    ///< Mixed encoding (multiple sub-encodings)
};

/// Compile-time gene value type for a given encoding.
/// Real → double (8 B), Binary → uint8_t (1 B), Integer/Permutation → int32_t (4 B).
/// For Encoding::Real this resolves to exactly `double`, so Population<Real>
/// generates identical machine code to a hardcoded std::vector<double> version.
template <Encoding E>
using gene_t =
    std::conditional_t<E == Encoding::Real, double,
    std::conditional_t<E == Encoding::Binary, uint8_t,
    std::conditional_t<E == Encoding::Integer, int32_t,
    std::conditional_t<E == Encoding::Permutation, int32_t,
    double>>>>;

/// Get encoding name for logging/debugging.
constexpr std::string_view encoding_name(Encoding enc) {
    switch (enc) {
    case Encoding::Real:
        return "Real";
    case Encoding::Binary:
        return "Binary";
    case Encoding::Permutation:
        return "Permutation";
    case Encoding::Integer:
        return "Integer";
    case Encoding::Composite:
        return "Composite";
    }
    return "Unknown";
}

} // namespace ea