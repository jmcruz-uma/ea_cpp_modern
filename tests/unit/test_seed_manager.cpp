// tests/unit/test_seed_manager.cpp
// Unit tests for ea::SeedManager — verifies exact protocol parity with
// Bio4Res/ea (ea_cpp_original) ea_main.cpp seed management.

#include <ea/util/seed_manager.hpp>

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

using ea::SeedManager;

// Use a throwaway filename so tests don't interfere with real experiments.
static const std::string TEST_FILE = "test_lastseed_tmp.txt";

static void remove_test_file() {
    std::remove(TEST_FILE.c_str());
}

static void write_test_file(long value) {
    std::ofstream f(TEST_FILE);
    f << value;
}

static long read_test_file() {
    std::ifstream f(TEST_FILE);
    long v = -999;
    f >> v;
    return v;
}

static bool file_exists() {
    return std::ifstream(TEST_FILE).good();
}

// ------------------------------------------------------------------
// Test A: no file → seed equals config_seed (first ever run)
// ------------------------------------------------------------------
static bool test_load_no_file() {
    std::cout << "TEST: load — no file → config_seed ... " << std::flush;
    remove_test_file();

    auto mgr = SeedManager::load(1L, TEST_FILE);

    bool pass = (mgr.seed == 1L);
    std::cout << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

// ------------------------------------------------------------------
// Test B: file present → seed = file_value + 1
// ------------------------------------------------------------------
static bool test_load_with_file() {
    std::cout << "TEST: load — file present → last+1 ... " << std::flush;

    // Simulate file left by a previous invocation
    write_test_file(7L);
    auto mgr = SeedManager::load(1L, TEST_FILE);  // config_seed ignored

    bool pass = (mgr.seed == 8L);
    remove_test_file();
    std::cout << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

// ------------------------------------------------------------------
// Test C: commit — seed < maxruns → file updated with seed + numruns - 1
// ------------------------------------------------------------------
static bool test_commit_updates_file() {
    std::cout << "TEST: commit — seed < maxruns → file written ... " << std::flush;
    remove_test_file();

    // Mirrors: first invocation, numruns=1, maxruns=30
    auto mgr = SeedManager::load(1L, TEST_FILE);  // seed = 1
    mgr.commit(1, 30L);

    // file_seed = 1 + 1 - 1 = 1; 1 < 30 → write 1
    bool pass = file_exists() && (read_test_file() == 1L);
    remove_test_file();
    std::cout << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

// ------------------------------------------------------------------
// Test D: commit — seed >= maxruns → file deleted
// ------------------------------------------------------------------
static bool test_commit_deletes_file() {
    std::cout << "TEST: commit — seed >= maxruns → file deleted ... " << std::flush;

    // Simulate: file contains 29, load gives seed=30, numruns=1, maxruns=30
    write_test_file(29L);
    auto mgr = SeedManager::load(1L, TEST_FILE);  // seed = 30
    mgr.commit(1, 30L);

    // file_seed = 30 + 0 = 30; 30 < 30 is false → delete
    bool pass = !file_exists();
    remove_test_file();
    std::cout << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

// ------------------------------------------------------------------
// Test E: full sequence — 3 invocations, numruns=1, maxruns=3
//   inv1: no file  → seed=1, commit → writes 1
//   inv2: file=1   → seed=2, commit → writes 2
//   inv3: file=2   → seed=3, commit → deletes (3 >= 3)
// ------------------------------------------------------------------
static bool test_full_sequence() {
    std::cout << "TEST: full sequence (3 runs, maxruns=3) ... " << std::flush;
    remove_test_file();
    bool pass = true;

    // Invocation 1
    {
        auto mgr = SeedManager::load(1L, TEST_FILE);
        if (mgr.seed != 1L) pass = false;
        mgr.commit(1, 3L);
        if (!file_exists() || read_test_file() != 1L) pass = false;
    }
    // Invocation 2
    {
        auto mgr = SeedManager::load(1L, TEST_FILE);
        if (mgr.seed != 2L) pass = false;
        mgr.commit(1, 3L);
        if (!file_exists() || read_test_file() != 2L) pass = false;
    }
    // Invocation 3
    {
        auto mgr = SeedManager::load(1L, TEST_FILE);
        if (mgr.seed != 3L) pass = false;
        mgr.commit(1, 3L);
        if (file_exists()) pass = false;   // must be deleted
    }

    remove_test_file();
    std::cout << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

// ------------------------------------------------------------------
// Test F: batch — numruns=5, maxruns=30
//   inv1: seed=1, commit → file_seed = 1+5-1 = 5; 5<30 → write 5
//   inv2: file=5 → seed=6, commit → file_seed=10; 10<30 → write 10
// ------------------------------------------------------------------
static bool test_batch_numruns() {
    std::cout << "TEST: batch mode (numruns=5, maxruns=30) ... " << std::flush;
    remove_test_file();
    bool pass = true;

    {
        auto mgr = SeedManager::load(1L, TEST_FILE);
        if (mgr.seed != 1L) pass = false;
        mgr.commit(5, 30L);
        if (!file_exists() || read_test_file() != 5L) pass = false;
    }
    {
        auto mgr = SeedManager::load(1L, TEST_FILE);
        if (mgr.seed != 6L) pass = false;
        mgr.commit(5, 30L);
        if (!file_exists() || read_test_file() != 10L) pass = false;
    }

    remove_test_file();
    std::cout << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

// ------------------------------------------------------------------
// Test G: in_progress() and reset()
// ------------------------------------------------------------------
static bool test_helpers() {
    std::cout << "TEST: in_progress() and reset() ... " << std::flush;
    remove_test_file();
    bool pass = true;

    if (SeedManager::in_progress(TEST_FILE)) pass = false;  // no file yet

    write_test_file(5L);
    if (!SeedManager::in_progress(TEST_FILE)) pass = false; // file exists

    SeedManager::reset(TEST_FILE);
    if (SeedManager::in_progress(TEST_FILE)) pass = false;  // deleted

    std::cout << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

// ------------------------------------------------------------------
// Main
// ------------------------------------------------------------------
int main() {
    std::cout << "========================================\n";
    std::cout << "ea-cpp SeedManager Unit Tests\n";
    std::cout << "========================================\n";

    int passed = 0, failed = 0;

    auto run = [&](bool result) { if (result) ++passed; else ++failed; };

    run(test_load_no_file());
    run(test_load_with_file());
    run(test_commit_updates_file());
    run(test_commit_deletes_file());
    run(test_full_sequence());
    run(test_batch_numruns());
    run(test_helpers());

    std::cout << "----------------------------------------\n";
    std::cout << "Results: " << passed << " passed, " << failed << " failed\n";
    std::cout << "========================================\n";

    return failed > 0 ? 1 : 0;
}
