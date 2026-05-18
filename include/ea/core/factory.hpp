#pragma once
/// @file factory.hpp
/// @brief Self-registering factory with C++20 constinit.
/// Eliminates static initialization order UB.
/// Pattern: https://www.cppstories.com/2023/ub-factory-constinit/
///
/// The factory stores creator functions in a compile-time-initialized map.
/// Types self-register via static init — order-independent thanks to constinit.

#include <array>
#include <cassert>
#include <functional>
#include <memory>
#include <string_view>

namespace ea {

// ============================================================
// Compile-time fixed-capacity map (no heap, constinit-safe)
// ============================================================

/// A simple flat map that satisfies constinit requirements:
/// - No dynamic allocation (fixed-size array)
/// - Trivial default construction (all zeroed)
/// - constexpr insert/lookup
template <typename Key, typename Value, std::size_t Capacity> class ConstinitMap {
public:
    struct Entry {
        Key key;
        Value value;
    };

    constexpr bool insert(Key key, Value value) {
        for (std::size_t i = 0; i < size_; ++i) {
            if (entries_[i].key == key)
                return false; // duplicate
        }
        if (size_ >= Capacity)
            return false; // full
        entries_[size_++] = {key, std::move(value)};
        return true;
    }

    [[nodiscard]] constexpr const Value* find(Key key) const {
        for (std::size_t i = 0; i < size_; ++i) {
            if (entries_[i].key == key)
                return &entries_[i].value;
        }
        return nullptr;
    }

    [[nodiscard]] constexpr std::size_t size() const { return size_; }

    [[nodiscard]] constexpr const Entry* begin() const { return entries_.data(); }
    [[nodiscard]] constexpr const Entry* end() const { return entries_.data() + size_; }

private:
    std::array<Entry, Capacity> entries_{};
    std::size_t size_{0};
};

// ============================================================
// Generic Factory
// ============================================================

/// Factory<T, Capacity> manages runtime creation of T-derived objects by name.
///
/// Usage:
///   using MyFactory = Factory<MyBase, 32>;
///   MyFactory::register_type("name", []{ return std::make_unique<MyImpl>(); });
///   auto obj = MyFactory::create("name");
template <typename T, std::size_t Capacity = 64> class Factory {
public:
    using CreatorFn = std::unique_ptr<T> (*)();

    /// Register a creator function under a name.
    /// Returns true on success, false on duplicate or full.
    static constexpr bool register_type(std::string_view name, CreatorFn fn) {
        return registry_.insert(name, fn);
    }

    /// Create an instance by name. Returns nullptr if not found.
    static std::unique_ptr<T> create(std::string_view name) {
        auto* fn = registry_.find(name);
        return fn ? (*fn)() : nullptr;
    }

    /// Check if a type is registered.
    static constexpr bool has(std::string_view name) { return registry_.find(name) != nullptr; }

    /// Number of registered types.
    static constexpr std::size_t size() { return registry_.size(); }

    /// Iterate over registered entries.
    static constexpr const auto& entries() { return registry_; }

private:
    /// constinit guarantees the map is zero-initialized at compile time,
    /// before any dynamic initialization runs. This eliminates the
    /// static initialization order fiasco for the registry itself.
    static inline constinit ConstinitMap<std::string_view, CreatorFn, Capacity> registry_{};
};

} // namespace ea