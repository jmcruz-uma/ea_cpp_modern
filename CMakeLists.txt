cmake_minimum_required(VERSION 3.25)
project(ea-cpp VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compile flags for performance
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
    add_compile_options(-march=native -O3)
endif()

# Header-only library
add_library(ea INTERFACE)
target_include_directories(ea INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# nlohmann/json (header-only, already present in original ea-cpp)
find_package(nlohmann_json 3.11.3 QUIET)
if(nlohmann_json_FOUND)
    message(STATUS "Found nlohmann_json via find_package")
else()
    # Fallback: look for local json.hpp
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/include/ea/util/json.hpp)
        message(STATUS "Using bundled nlohmann/json.hpp")
    else()
        message(STATUS "Downloading nlohmann/json.hpp")
        file(DOWNLOAD
            https://github.com/nlohmann/json/releases/download/v3.11.3/json.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/ea/util/json.hpp
            SHOW_PROGRESS
        )
    endif()
endif()

# Source files (problems with explicit instantiation, benchmarks, examples)
file(GLOB EA_SOURCES src/*.cpp)

# Main executable
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/src/ea_main.cpp OR EA_SOURCES)
    add_executable(ea_main ${EA_SOURCES})
    target_link_libraries(ea_main PRIVATE ea)
    if(nlohmann_json_FOUND)
        target_link_libraries(ea_main PRIVATE nlohmann_json::nlohmann_json)
    endif()
endif()

# --- Testing ---
option(EA_BUILD_TESTS "Build tests" ON)
if(EA_BUILD_TESTS)
    enable_testing()
    find_package(GTest QUIET)
    if(NOT GTest_FOUND)
        message(STATUS "Fetching GoogleTest...")
        include(FetchContent)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG v1.14.0
        )
        FetchContent_MakeAvailable(googletest)
    endif()

    file(GLOB TEST_SOURCES tests/unit/*.cpp)
    if(TEST_SOURCES)
        add_executable(ea_tests ${TEST_SOURCES})
        target_link_libraries(ea_tests PRIVATE ea GTest::gtest_main)
        if(nlohmann_json_FOUND)
            target_link_libraries(ea_tests PRIVATE nlohmann_json::nlohmann_json)
        endif()
        include(GoogleTest)
        gtest_discover_tests(ea_tests)
    endif()
endif()

# --- Benchmarks ---
option(EA_BUILD_BENCHMARKS "Build benchmarks" ON)
if(EA_BUILD_BENCHMARKS)
    file(GLOB BENCH_SOURCES tests/benchmark/*.cpp)
    if(BENCH_SOURCES)
        add_executable(ea_bench ${BENCH_SOURCES})
        target_link_libraries(ea_bench PRIVATE ea)
        if(nlohmann_json_FOUND)
            target_link_libraries(ea_bench PRIVATE nlohmann_json::nlohmann_json)
        endif()
    endif()
endif()