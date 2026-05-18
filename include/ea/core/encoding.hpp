#pragma once
/// @file encoding.hpp
/// @brief Solution encoding types for evolutionary algorithms.

#include <cstdint>
#include <string_view>

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