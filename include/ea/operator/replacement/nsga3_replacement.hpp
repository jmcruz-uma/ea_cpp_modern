#pragma once
/// @file nsga3_replacement.hpp
/// @brief NSGA-III Environmental Selection — reference-point-based niching.
///
/// Full implementation of Algorithm 2 (normalize) + Algorithm 3 (associate) +
/// Algorithm 4 (niching) from Deb & Jain, IEEE TEVC 2014.
///
/// jMetal reference: EnvironmentalSelection.java (nsgaiii/util/)
///
/// Bug history (corrected vs previous implementation):
///   - find_closest_ref used Euclidean distance to ref position; must use
///     perpendicular distance to the LINE from origin through ref direction.
///   - No hyperplane fitting: missing ASF extreme-point search, Gaussian
///     elimination, and normalization by (intercept - ideal).
///   - Niching always picked closest member; must pick RANDOM when
///     member_count > 0 (Algorithm 4, Section IV-E in the paper).

#include <algorithm>
#include <cmath>
#include <ea/core/comparator.hpp>
#include <ea/core/population.hpp>
#include <ea/util/random.hpp>
#include <ea/util/reference_point.hpp>
#include <limits>
#include <vector>

namespace ea {

/// NSGA-III environmental selection (replacement operator).
struct NSGAIIIReplacement {
    std::vector<ReferencePoint> reference_points;

    NSGAIIIReplacement() = default;
    explicit NSGAIIIReplacement(std::vector<ReferencePoint> refs)
        : reference_points(std::move(refs)) {}

    /// Environmental selection on a combined (2N) population.
    /// Returns target_size indices from combined to keep as the next population.
    std::vector<int> replace(this NSGAIIIReplacement& self, Population<>& combined,
                             int target_size);

private:
    /// Perpendicular distance from `point` to the LINE from origin through `dir`.
    static double perp_dist(const double* point, const double* dir, int m) {
        double dot = 0.0, dir_sq = 0.0;
        for (int i = 0; i < m; ++i) {
            dot += point[i] * dir[i];
            dir_sq += dir[i] * dir[i];
        }
        double dist_sq = 0.0;
        for (int i = 0; i < m; ++i) {
            double proj = (dir_sq > 1e-14) ? dot / dir_sq * dir[i] : 0.0;
            double diff = point[i] - proj;
            dist_sq += diff * diff;
        }
        return std::sqrt(dist_sq);
    }

    /// ASF (Achievement Scalarization Function) on translated objectives.
    /// axis f: weight_f = 1.0, weight_i = 1e-6 for i != f.
    /// jMetal: EnvironmentalSelection.ASF()
    static double asf(const double* obj, int f, int m) {
        double max_val = -std::numeric_limits<double>::max();
        for (int i = 0; i < m; ++i) {
            double w = (i == f) ? 1.0 : 1e-6;
            max_val = std::max(max_val, obj[i] / w);
        }
        return max_val;
    }

    /// Gaussian elimination (no partial pivoting, matching jMetal exactly).
    /// A is an n x (n+1) augmented matrix [M | b] modified in-place.
    /// jMetal: EnvironmentalSelection.guassianElimination()
    static std::vector<double> gaussian_elimination(std::vector<std::vector<double>> A, int n) {
        for (int base = 0; base < n - 1; ++base) {
            for (int target = base + 1; target < n; ++target) {
                if (std::abs(A[base][base]) < 1e-14)
                    continue;
                double ratio = A[target][base] / A[base][base];
                for (int col = 0; col <= n; ++col)
                    A[target][col] -= A[base][col] * ratio;
            }
        }
        std::vector<double> x(n, 0.0);
        for (int i = n - 1; i >= 0; --i) {
            for (int k = i + 1; k < n; ++k)
                A[i][n] -= A[i][k] * x[k];
            x[i] = (std::abs(A[i][i]) > 1e-14) ? A[i][n] / A[i][i] : 0.0;
        }
        return x;
    }
};

inline std::vector<int> NSGAIIIReplacement::replace(this NSGAIIIReplacement& self,
                                                    Population<>& combined, int target_size) {
    const int n_obj = combined.n_obj;

    // ── Step 1: Non-dominated sort, fill complete fronts ──────────────────
    auto fronts = fast_non_dominated_sort(combined);

    std::vector<int> selected;
    int last_f = -1;
    for (int f = 0; f < (int)fronts.size(); ++f) {
        if ((int)(selected.size() + fronts[f].size()) <= target_size) {
            selected.insert(selected.end(), fronts[f].begin(), fronts[f].end());
        } else {
            last_f = f;
            break;
        }
    }

    if ((int)selected.size() >= target_size || last_f < 0) {
        selected.resize(target_size);
        return selected;
    }

    const auto& last_front = fronts[last_f];
    int remaining = target_size - (int)selected.size();

    // ── Step 2: Ideal point from first front (Algorithm 2, line 1) ────────
    // jMetal: translateObjectives() — "min values must appear in the first front"
    std::vector<double> ideal(n_obj, std::numeric_limits<double>::max());
    for (int idx : fronts[0]) {
        for (int o = 0; o < n_obj; ++o)
            ideal[o] = std::min(ideal[o], combined.objective(idx, o));
    }

    // Helper: translated objective (obj - ideal)
    auto trans = [&](int idx, int o) { return combined.objective(idx, o) - ideal[o]; };

    // ── Step 3: Extreme points via ASF, from first front ──────────────────
    // jMetal: findExtremePoints()
    std::vector<int> extreme(n_obj);
    {
        std::vector<double> buf(n_obj);
        for (int f = 0; f < n_obj; ++f) {
            double best = std::numeric_limits<double>::max();
            extreme[f] = fronts[0][0];
            for (int idx : fronts[0]) {
                for (int o = 0; o < n_obj; ++o)
                    buf[o] = trans(idx, o);
                double v = asf(buf.data(), f, n_obj);
                if (v < best) {
                    best = v;
                    extreme[f] = idx;
                }
            }
        }
    }

    // ── Step 4: Hyperplane — construct intercepts ──────────────────────────
    // jMetal: constructHyperplane() — uses ORIGINAL (not translated) objectives
    bool duplicate = false;
    for (int i = 0; i < n_obj && !duplicate; ++i)
        for (int j = i + 1; j < n_obj && !duplicate; ++j)
            duplicate = (extreme[i] == extreme[j]);

    std::vector<double> intercept(n_obj);
    if (!duplicate) {
        // Solve A x = b where A[i][j] = extreme[i].obj[j], b = [1,...,1]
        std::vector<std::vector<double>> A(n_obj, std::vector<double>(n_obj + 1, 0.0));
        for (int i = 0; i < n_obj; ++i) {
            for (int j = 0; j < n_obj; ++j)
                A[i][j] = combined.objective(extreme[i], j);
            A[i][n_obj] = 1.0;
        }
        auto x = gaussian_elimination(std::move(A), n_obj);
        for (int f = 0; f < n_obj; ++f)
            intercept[f] = (std::abs(x[f]) > 1e-14) ? 1.0 / x[f] : 1e14;
    } else {
        // Fallback: original objective value of the extreme point for each axis
        for (int f = 0; f < n_obj; ++f)
            intercept[f] = combined.objective(extreme[f], f);
    }

    // Helper: normalized objective (translated then divided by intercept span)
    // jMetal: normalizeObjectives() — denom = intercept[f] - ideal[f]
    auto norm_obj = [&](int idx, int o) {
        double denom = intercept[o] - ideal[o];
        return trans(idx, o) / (std::abs(denom) > 1e-9 ? denom : 1e-9);
    };

    // Helper: associate individual with nearest reference point (perpendicular distance)
    // Returns {ref_index, distance}. jMetal: associate()
    auto associate = [&](int idx) -> std::pair<int, double> {
        std::vector<double> nrm(n_obj);
        for (int o = 0; o < n_obj; ++o)
            nrm[o] = norm_obj(idx, o);
        int best_r = 0;
        double best_d = std::numeric_limits<double>::max();
        for (int r = 0; r < (int)self.reference_points.size(); ++r) {
            double d = perp_dist(nrm.data(), self.reference_points[r].position.data(), n_obj);
            if (d < best_d) {
                best_d = d;
                best_r = r;
            }
        }
        return {best_r, best_d};
    };

    // ── Step 5: Associate all individuals with reference points ────────────
    for (auto& rp : self.reference_points)
        rp.clear();

    // Non-critical fronts → member counts
    for (int idx : selected) {
        auto [r, d] = associate(idx);
        (void)d;
        self.reference_points[r].add_member();
    }

    // Critical front → potential members with distance
    for (int idx : last_front) {
        auto [r, d] = associate(idx);
        self.reference_points[r].add_potential_member(idx, d);
    }

    // Sort potential members ascending by distance (closest first)
    for (auto& rp : self.reference_points)
        rp.sort_potential_members();

    // ── Step 6: Niching selection (Algorithm 4) ────────────────────────────
    // - Pick reference point with fewest members (random among ties)
    // - member_count == 0 → FindClosestMember; member_count >= 1 → RandomMember
    // jMetal: SelectClusterMember()
    auto& rng = Random::instance();
    std::vector<int> from_last;
    from_last.reserve(remaining);

    while ((int)from_last.size() < remaining) {
        // Minimum member count among refs that still have potential members
        int min_cnt = std::numeric_limits<int>::max();
        for (const auto& rp : self.reference_points)
            if (rp.has_potential_members() && rp.member_count < min_cnt)
                min_cnt = rp.member_count;

        if (min_cnt == std::numeric_limits<int>::max())
            break; // no potential members left (shouldn't happen)

        std::vector<int> cands;
        for (int r = 0; r < (int)self.reference_points.size(); ++r)
            if (self.reference_points[r].has_potential_members() &&
                self.reference_points[r].member_count == min_cnt)
                cands.push_back(r);

        int cr = cands[rng.uniform_int(0, (int)cands.size() - 1)];
        auto& rp = self.reference_points[cr];

        int sel = (rp.member_count == 0) ? rp.remove_closest_member()
                                         : rp.remove_random_member(rng.engine());
        rp.add_member();
        from_last.push_back(sel);
    }

    selected.insert(selected.end(), from_last.begin(), from_last.end());
    return selected;
}

} // namespace ea
