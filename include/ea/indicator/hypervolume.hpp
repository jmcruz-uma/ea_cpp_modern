#pragma once
/// @file hypervolume.hpp
/// @brief WFG hypervolume indicator — the most important quality indicator for MOEAs.
///
/// Reference: jMetal
///   jmetal-core/src/main/java/org/uma/jmetal/qualityindicator/impl/hypervolume/impl/WFGHypervolume.java
///
/// Algorithm: WFG "Walking Fish Group" hypervolume by While, Bradstreet, Barone.
/// Uses inclusive/exclusive recursive computation with base cases for 1-4 points.
/// For 2D: uses O(n log n) dimension sweep.
/// For 3D: uses WFG recursive method.
/// For 4D+: uses WFG recursive method.

#include <algorithm>
#include <cmath>
#include <ea/core/population.hpp>
#include <limits>
#include <numeric>
#include <vector>

namespace ea {

namespace detail {

struct Point {
    std::vector<double> objectives;
};

struct Front {
    int n_points = 0;
    std::vector<Point> points;

    explicit Front(int max_points, int n_obj) {
        points.resize(max_points);
        for (auto& p : points)
            p.objectives.resize(n_obj);
    }
};

inline bool beats(double x, double y) {
    return x > y;
}
inline double worse(double x, double y) {
    return x > y ? y : x;
}

/// -1: p dominates q, 1: q dominates p, 2: equal, 0: incomparable
inline int dominates2way(const Point& p, const Point& q, int k) {
    for (int i = k; i >= 0; --i) {
        if (beats(p.objectives[i], q.objectives[i])) {
            for (int j = i - 1; j >= 0; --j)
                if (beats(q.objectives[j], p.objectives[j]))
                    return 0;
            return -1;
        } else if (beats(q.objectives[i], p.objectives[i])) {
            for (int j = i - 1; j >= 0; --j)
                if (beats(p.objectives[j], q.objectives[j]))
                    return 0;
            return 1;
        }
    }
    return 2;
}

inline bool dominates1way(const Point& p, const Point& q, int k) {
    for (int i = k; i >= 0; --i)
        if (beats(q.objectives[i], p.objectives[i]))
            return false;
    return true;
}

inline double inclusive_hv(const Point& p, int n_obj) {
    double volume = 1.0;
    for (int i = 0; i < n_obj; ++i)
        volume *= p.objectives[i];
    return volume;
}

inline double inclusive_hv2(const Point& p, const Point& q, int n_obj) {
    double vp = 1.0, vq = 1.0, vpq = 1.0;
    for (int i = 0; i < n_obj; ++i) {
        vp *= p.objectives[i];
        vq *= q.objectives[i];
        vpq *= worse(p.objectives[i], q.objectives[i]);
    }
    return vp + vq - vpq;
}

inline double inclusive_hv3(const Point& p, const Point& q, const Point& r, int n_obj) {
    double vp = 1.0, vq = 1.0, vr = 1.0;
    double vpq = 1.0, vpr = 1.0, vqr = 1.0, vpqr = 1.0;
    for (int i = 0; i < n_obj; ++i) {
        vp *= p.objectives[i];
        vq *= q.objectives[i];
        vr *= r.objectives[i];
        if (beats(p.objectives[i], q.objectives[i])) {
            if (beats(q.objectives[i], r.objectives[i])) {
                vpq *= q.objectives[i];
                vpr *= r.objectives[i];
                vqr *= r.objectives[i];
                vpqr *= r.objectives[i];
            } else {
                vpq *= q.objectives[i];
                vpr *= worse(p.objectives[i], r.objectives[i]);
                vqr *= q.objectives[i];
                vpqr *= q.objectives[i];
            }
        } else if (beats(p.objectives[i], r.objectives[i])) {
            vpq *= p.objectives[i];
            vpr *= r.objectives[i];
            vqr *= r.objectives[i];
            vpqr *= r.objectives[i];
        } else {
            vpq *= p.objectives[i];
            vpr *= p.objectives[i];
            vqr *= worse(q.objectives[i], r.objectives[i]);
            vpqr *= p.objectives[i];
        }
    }
    return vp + vq + vr - vpq - vpr - vqr + vpqr;
}

inline double inclusive_hv4(const Point& p, const Point& q, const Point& r, const Point& s,
                            int n_obj) {
    double vp = 1.0, vq = 1.0, vr = 1.0, vs = 1.0;
    double vpq = 1.0, vpr = 1.0, vps = 1.0, vqr = 1.0, vqs = 1.0, vrs = 1.0;
    double vpqr = 1.0, vpqs = 1.0, vprs = 1.0, vqrs = 1.0, vpqrs = 1.0;
    for (int i = 0; i < n_obj; ++i) {
        double po = p.objectives[i], qo = q.objectives[i], ro = r.objectives[i],
               so = s.objectives[i];
        vp *= po;
        vq *= qo;
        vr *= ro;
        vs *= so;
        if (beats(po, qo)) {
            if (beats(qo, ro)) {
                if (beats(ro, so)) {
                    vpq *= qo;
                    vpr *= ro;
                    vps *= so;
                    vqr *= ro;
                    vqs *= so;
                    vrs *= so;
                    vpqr *= ro;
                    vpqs *= so;
                    vprs *= so;
                    vqrs *= so;
                    vpqrs *= so;
                } else {
                    double z1 = worse(qo, so);
                    vpq *= qo;
                    vpr *= ro;
                    vps *= worse(po, so);
                    vqr *= ro;
                    vqs *= z1;
                    vrs *= ro;
                    vpqr *= ro;
                    vpqs *= z1;
                    vprs *= ro;
                    vqrs *= ro;
                    vpqrs *= ro;
                }
            } else if (beats(qo, so)) {
                vpq *= qo;
                vpr *= worse(po, ro);
                vps *= so;
                vqr *= qo;
                vqs *= so;
                vrs *= so;
                vpqr *= qo;
                vpqs *= so;
                vprs *= so;
                vqrs *= so;
                vpqrs *= so;
            } else {
                double z1 = worse(po, ro);
                vpq *= qo;
                vpr *= z1;
                vps *= worse(po, so);
                vqr *= qo;
                vqs *= qo;
                vrs *= worse(ro, so);
                vpqr *= qo;
                vpqs *= qo;
                vprs *= worse(z1, so);
                vqrs *= qo;
                vpqrs *= qo;
            }
        } else if (beats(qo, ro)) {
            if (beats(po, so)) {
                double z1 = worse(po, ro), z2 = worse(ro, so);
                vpq *= po;
                vpr *= z1;
                vps *= so;
                vqr *= ro;
                vqs *= so;
                vrs *= z2;
                vpqr *= z1;
                vpqs *= so;
                vprs *= z2;
                vqrs *= z2;
                vpqrs *= z2;
            } else {
                double z1 = worse(po, ro), z2 = worse(ro, so);
                vpq *= po;
                vpr *= z1;
                vps *= po;
                vqr *= ro;
                vqs *= worse(qo, so);
                vrs *= z2;
                vpqr *= z1;
                vpqs *= po;
                vprs *= z1;
                vqrs *= z2;
                vpqrs *= z1;
            }
        } else if (beats(po, so)) {
            vpq *= po;
            vpr *= po;
            vps *= so;
            vqr *= qo;
            vqs *= so;
            vrs *= so;
            vpqr *= po;
            vpqs *= so;
            vprs *= so;
            vqrs *= so;
            vpqrs *= so;
        } else {
            double z1 = worse(qo, so);
            vpq *= po;
            vpr *= po;
            vps *= po;
            vqr *= qo;
            vqs *= z1;
            vrs *= worse(ro, so);
            vpqr *= po;
            vpqs *= po;
            vprs *= po;
            vqrs *= z1;
            vpqrs *= po;
        }
    }
    return vp + vq + vr + vs - vpq - vpr - vps - vqr - vqs - vrs + vpqr + vpqs + vprs + vqrs -
           vpqrs;
}

// Forward declarations
void make_dominated_bit(Front& ps, int p, int n_obj, std::vector<Front>& fs, int& fr, int& safe);
double hv(Front& ps, int n_obj, std::vector<Front>& fs, int& fr, int& safe);

inline double excl_hv(Front& ps, int p, int n_obj, std::vector<Front>& fs, int& fr, int& safe) {
    make_dominated_bit(ps, p, n_obj, fs, fr, safe);
    double a = inclusive_hv(ps.points[p], n_obj);
    double b = hv(fs[fr - 1], n_obj - 1, fs, fr, safe);
    fr--;
    return a - b;
}

inline double hv2(const Front& ps, int k, int /*n_obj*/) {
    // 2D hypervolume (sorted improving)
    double volume = ps.points[0].objectives[0] * ps.points[0].objectives[1];
    for (int i = 1; i < k; ++i) {
        volume += ps.points[i].objectives[1] *
                  (ps.points[i].objectives[0] - ps.points[i - 1].objectives[0]);
    }
    return volume;
}

inline void make_dominated_bit(Front& ps, int p, int n_obj, std::vector<Front>& fs, int& fr,
                               int& safe) {
    int l = 0;
    int u = p - 1;
    for (int i = p - 1; i >= 0; --i) {
        if (beats(ps.points[p].objectives[n_obj - 1], ps.points[i].objectives[n_obj - 1])) {
            fs[fr].points[u].objectives[n_obj - 1] = ps.points[i].objectives[n_obj - 1];
            for (int j = 0; j < n_obj - 1; ++j)
                fs[fr].points[u].objectives[j] =
                    worse(ps.points[p].objectives[j], ps.points[i].objectives[j]);
            --u;
        } else {
            fs[fr].points[l].objectives[n_obj - 1] = ps.points[p].objectives[n_obj - 1];
            for (int j = 0; j < n_obj - 1; ++j)
                fs[fr].points[l].objectives[j] =
                    worse(ps.points[p].objectives[j], ps.points[i].objectives[j]);
            ++l;
        }
    }

    // Points below l: equal in last objective; above l: worse in last objective
    fs[fr].n_points = 1;
    for (int i = 1; i < l; ++i) {
        int j = 0;
        while (j < fs[fr].n_points) {
            int dom = dominates2way(fs[fr].points[i], fs[fr].points[j], n_obj - 2);
            if (dom == 0) {
                ++j;
            } else if (dom == -1) {
                std::swap(fs[fr].points[j], fs[fr].points[i]);
                while (
                    j < fs[fr].n_points - 1 &&
                    dominates1way(fs[fr].points[j], fs[fr].points[fs[fr].n_points - 1], n_obj - 1))
                    fs[fr].n_points--;
                int k = j + 1;
                while (k < fs[fr].n_points) {
                    if (dominates1way(fs[fr].points[j], fs[fr].points[k], n_obj - 2)) {
                        std::swap(fs[fr].points[k], fs[fr].points[fs[fr].n_points - 1]);
                        fs[fr].n_points--;
                    } else {
                        ++k;
                    }
                }
            } else {
                j = fs[fr].n_points + 1;
            }
        }
        if (j == fs[fr].n_points) {
            std::swap(fs[fr].points[fs[fr].n_points], fs[fr].points[i]);
            fs[fr].n_points++;
        }
    }

    safe = static_cast<int>(worse(static_cast<double>(l), static_cast<double>(fs[fr].n_points)));
    for (int i = l; i < p; ++i) {
        int j = 0;
        while (j < safe) {
            if (dominates1way(fs[fr].points[j], fs[fr].points[i], n_obj - 2))
                j = fs[fr].n_points + 1;
            else
                ++j;
        }
        while (j < fs[fr].n_points) {
            int dom = dominates2way(fs[fr].points[i], fs[fr].points[j], n_obj - 1);
            if (dom == 0) {
                ++j;
            } else if (dom == -1) {
                std::swap(fs[fr].points[j], fs[fr].points[i]);
                while (
                    j < fs[fr].n_points - 1 &&
                    dominates1way(fs[fr].points[j], fs[fr].points[fs[fr].n_points - 1], n_obj - 1))
                    fs[fr].n_points--;
                int k = j + 1;
                while (k < fs[fr].n_points) {
                    if (dominates1way(fs[fr].points[j], fs[fr].points[k], n_obj - 1)) {
                        std::swap(fs[fr].points[k], fs[fr].points[fs[fr].n_points - 1]);
                        fs[fr].n_points--;
                    } else {
                        ++k;
                    }
                }
            } else {
                j = fs[fr].n_points + 1;
            }
        }
        if (j == fs[fr].n_points) {
            std::swap(fs[fr].points[fs[fr].n_points], fs[fr].points[i]);
            fs[fr].n_points++;
        }
    }
    ++fr;
}

inline double hv(Front& ps, int n_obj, std::vector<Front>& fs, int& fr, int& safe) {
    switch (ps.n_points) {
    case 1:
        return inclusive_hv(ps.points[0], n_obj);
    case 2:
        return inclusive_hv2(ps.points[0], ps.points[1], n_obj);
    case 3:
        return inclusive_hv3(ps.points[0], ps.points[1], ps.points[2], n_obj);
    case 4:
        return inclusive_hv4(ps.points[0], ps.points[1], ps.points[2], ps.points[3], n_obj);
    default:
        break;
    }

    // Sort by last objective descending
    std::vector<int> order(ps.n_points);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return ps.points[a].objectives[n_obj - 1] > ps.points[b].objectives[n_obj - 1];
    });

    // Reorder points in-place using the sorted order
    std::vector<Point> sorted_points(ps.n_points);
    for (int i = 0; i < ps.n_points; ++i)
        sorted_points[i] = ps.points[order[i]];
    for (int i = 0; i < ps.n_points; ++i)
        ps.points[i] = sorted_points[i];

    if (n_obj == 2)
        return hv2(ps, ps.n_points, n_obj);

    if (n_obj == 3 && safe > 0) {
        double volume = ps.points[0].objectives[2] * hv2(ps, safe, 2);
        for (int i = safe; i < ps.n_points; ++i) {
            volume += ps.points[i].objectives[n_obj - 1] * excl_hv(ps, i, n_obj, fs, fr, safe);
        }
        return volume;
    } else {
        double volume =
            inclusive_hv4(ps.points[0], ps.points[1], ps.points[2], ps.points[3], n_obj);
        for (int i = 4; i < ps.n_points; ++i) {
            volume += ps.points[i].objectives[n_obj - 1] * excl_hv(ps, i, n_obj, fs, fr, safe);
        }
        return volume;
    }
}

/// 2D hypervolume using efficient O(n log n) dimension sweep
inline double hypervolume_2d(std::vector<std::vector<double>> points,
                             const std::vector<double>& ref_point) {
    // Translate: reference point becomes origin
    for (auto& p : points) {
        for (size_t i = 0; i < p.size(); ++i) {
            p[i] = ref_point[i] - p[i];
            if (p[i] < 0)
                p[i] = 0;
        }
    }

    // Remove points with any zero coordinate (dominated by ref)
    points.erase(std::remove_if(points.begin(), points.end(),
                                [](const auto& p) {
                                    return std::any_of(p.begin(), p.end(),
                                                       [](double v) { return v <= 0; });
                                }),
                 points.end());

    if (points.empty())
        return 0.0;

    // Sort by first objective descending
    std::sort(points.begin(), points.end(),
              [](const auto& a, const auto& b) { return a[0] > b[0]; });

    double volume = points[0][0] * points[0][1];
    for (size_t i = 1; i < points.size(); ++i) {
        volume += points[i][1] * (points[i][0] - points[i - 1][0]);
    }
    return volume;
}

/// WFG hypervolume for 3D+
inline double hypervolume_wfg(const std::vector<std::vector<double>>& points, int n_obj,
                              const std::vector<double>& ref_point) {
    // Translate points
    std::vector<std::vector<double>> translated;
    translated.reserve(points.size());
    for (const auto& p : points) {
        std::vector<double> tp(n_obj);
        bool dominated = false;
        for (int i = 0; i < n_obj; ++i) {
            tp[i] = ref_point[i] - p[i];
            if (tp[i] < 0)
                dominated = true;
        }
        if (!dominated)
            translated.push_back(tp);
    }

    if (translated.empty())
        return 0.0;
    if (translated.size() == 1) {
        double vol = 1.0;
        for (double v : translated[0])
            vol *= v;
        return vol;
    }

    int n_points = static_cast<int>(translated.size());

    // Build front
    Front front(n_points, n_obj);
    front.n_points = n_points;
    for (int i = 0; i < n_points; ++i) {
        front.points[i].objectives = translated[i];
    }

    int max_depth = std::max(0, n_obj - 2);
    std::vector<Front> fs;
    fs.reserve(max_depth);
    for (int i = 0; i < max_depth; ++i) {
        fs.emplace_back(n_points, n_obj);
    }

    int fr = 0;
    int safe = 0;
    return hv(front, n_obj, fs, fr, safe);
}

} // namespace detail

/// Hypervolume indicator using WFG algorithm.
///
/// @param pop           Population<> containing the Pareto front approximation
/// @param indices       Indices of individuals in the front
/// @param ref_point     Reference point (must dominate all front points)
/// @return Hypervolume value
inline double hypervolume(const Population<>& pop, const std::vector<int>& indices,
                          const std::vector<double>& ref_point) {
    if (indices.empty())
        return 0.0;

    int n_obj = pop.n_obj;
    std::vector<std::vector<double>> points;
    points.reserve(indices.size());
    for (int idx : indices) {
        std::vector<double> p(n_obj);
        for (int o = 0; o < n_obj; ++o)
            p[o] = pop.objective(idx, o);
        points.push_back(p);
    }

    if (n_obj == 2) {
        return detail::hypervolume_2d(points, ref_point);
    }
    return detail::hypervolume_wfg(points, n_obj, ref_point);
}

/// Hypervolume indicator from a vector of objective vectors.
inline double hypervolume(const std::vector<std::vector<double>>& points,
                          const std::vector<double>& ref_point) {
    if (points.empty())
        return 0.0;
    int n_obj = static_cast<int>(points[0].size());
    if (n_obj == 2) {
        return detail::hypervolume_2d(points, ref_point);
    }
    return detail::hypervolume_wfg(points, n_obj, ref_point);
}

struct HypervolumeIndicator {
    std::vector<double> ref_point;

    double compute(this HypervolumeIndicator& self, const std::vector<std::vector<double>>& front) {
        if (self.ref_point.empty() || front.empty())
            return 0.0;
        return hypervolume(front, self.ref_point);
    }

    double compute(this HypervolumeIndicator& self, const Population<>& pop,
                   const std::vector<int>& indices) {
        if (self.ref_point.empty() || indices.empty())
            return 0.0;
        return hypervolume(pop, indices, self.ref_point);
    }
};

} // namespace ea
