#pragma once
/// @file bench_timer.hpp
/// @brief High-resolution timer for benchmarking.
/// Uses std::chrono::steady_clock for wall-clock time
/// and platform-specific APIs for CPU time.

#include <chrono>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <unistd.h>

namespace ea::bench {

struct TimingResult {
    std::string label;
    int evaluations = 0;
    int population_size = 0;
    int dimensions = 0;
    int n_objectives = 0;
    int runs = 0;
    double wall_ms = 0;
    double cpu_ms = 0;
    double mean_igd = 0;
    double std_igd = 0;
};

class Timer {
public:
    void start() {
        wall_start_ = std::chrono::steady_clock::now();
        cpu_start_ = cpu_time_ms();
    }

    void stop() {
        auto wall_end = std::chrono::steady_clock::now();
        wall_ms_ = std::chrono::duration<double, std::milli>(wall_end - wall_start_).count();
        cpu_ms_ = cpu_time_ms() - cpu_start_;
    }

    double wall_ms() const { return wall_ms_; }
    double cpu_ms() const { return cpu_ms_; }

private:
    std::chrono::steady_clock::time_point wall_start_;
    double cpu_start_ = 0;
    double wall_ms_ = 0;
    double cpu_ms_ = 0;

    static double cpu_time_ms() {
        std::ifstream proc("/proc/self/stat");
        std::string ignore;
        for (int i = 0; i < 13; ++i) proc >> ignore;
        long utime, stime;
        proc >> utime >> stime;
        static const long clk_tck = ::sysconf(_SC_CLK_TCK);
        return (utime + stime) * 1000.0 / clk_tck;
    }
};

/// Print a CSV header to stdout.
inline void print_csv_header() {
    std::cout << "label,algorithm,problem,dim,n_obj,pop_size,max_evals,runs,"
              << "wall_ms,cpu_ms,mean_igd,std_igd\n";
}

/// Print one CSV row.
inline void print_csv_row(const TimingResult& r) {
    std::cout << r.label << ","
              << r.evaluations << ","
              << r.population_size << ","
              << r.dimensions << ","
              << r.n_objectives << ","
              << r.runs << ","
              << std::fixed << std::setprecision(3)
              << r.wall_ms << ","
              << r.cpu_ms << ","
              << std::setprecision(8)
              << r.mean_igd << ","
              << r.std_igd << "\n";
}

/// Print a formatted table to stderr.
inline void print_result(const TimingResult& r) {
    std::cerr << "=== " << r.label << " ===\n"
              << "  Dimensions: " << r.dimensions << "\n"
              << "  Objectives:  " << r.n_objectives << "\n"
              << "  Pop size:    " << r.population_size << "\n"
              << "  Max evals:   " << r.evaluations << "\n"
              << "  Runs:        " << r.runs << "\n"
              << "  Wall time:   " << std::fixed << std::setprecision(3) << r.wall_ms << " ms\n"
              << "  CPU time:    " << r.cpu_ms << " ms\n"
              << "  Mean IGD:    " << std::setprecision(8) << r.mean_igd << "\n"
              << "  Std IGD:     " << std::setprecision(8) << r.std_igd << "\n";
}

} // namespace ea::bench