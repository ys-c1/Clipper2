// Standalone self-test for the path statistics / self-intersection /
// area-conservation helpers added to clipper.core.h.
//
// No third-party dependencies: it only links against the Clipper2 sources
// that ship with this repository. Build and run it with
//   Tests/run_path_stats_selftest.sh
//
// Exit code is 0 only when every check passes.

#include "clipper2/clipper.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

using namespace Clipper2Lib;

namespace
{
  int g_failures = 0;
  int g_checks = 0;

  void Check(bool cond, const std::string& what)
  {
    ++g_checks;
    if (!cond)
    {
      ++g_failures;
      std::printf("  [FAIL] %s\n", what.c_str());
    }
  }

  // ---------------------------------------------------------------------
  // Brute-force reference implementations (deliberately written as plain
  // per-point loops, independent of the library's optimised pairing).
  // ---------------------------------------------------------------------

  template <typename T>
  long double RefShoelace(const Path<T>& p)
  {
    if (p.size() < 3) return 0.0L;
    long double a = 0.0L;
    for (size_t i = 0; i < p.size(); ++i)
    {
      const auto& cur = p[i];
      const auto& nxt = p[(i + 1) % p.size()];
      a += static_cast<long double>(cur.x) * static_cast<long double>(nxt.y) -
        static_cast<long double>(nxt.x) * static_cast<long double>(cur.y);
    }
    return a / 2.0L;
  }

  template <typename T>
  long double RefPerimeter(const Path<T>& p)
  {
    if (p.size() < 2) return 0.0L;
    long double total = 0.0L;
    for (size_t i = 0; i < p.size(); ++i)
    {
      const auto& a = p[i];
      const auto& b = p[(i + 1) % p.size()];
      const long double dx = static_cast<long double>(b.x) - a.x;
      const long double dy = static_cast<long double>(b.y) - a.y;
      total += std::sqrt(dx * dx + dy * dy);
    }
    return total;
  }

  // Exact orientation, independent of Clipper2Lib.
  __extension__ typedef __int128 i128;
  static int RefOrient(const Point64& a, const Point64& b, const Point64& c)
  {
    const i128 ab = static_cast<i128>(b.x - a.x) * (c.y - a.y);
    const i128 ac = static_cast<i128>(b.y - a.y) * (c.x - a.x);
    return (ab > ac) ? 1 : ((ab < ac) ? -1 : 0);
  }

  static bool RefInBox(const Point64& p, const Point64& a, const Point64& b)
  {
    return p.x >= (std::min)(a.x, b.x) && p.x <= (std::max)(a.x, b.x) &&
      p.y >= (std::min)(a.y, b.y) && p.y <= (std::max)(a.y, b.y);
  }

  static bool RefSegTouch(const Point64& a, const Point64& b,
    const Point64& c, const Point64& d)
  {
    const int o1 = RefOrient(a, b, c), o2 = RefOrient(a, b, d);
    const int o3 = RefOrient(c, d, a), o4 = RefOrient(c, d, b);
    if (o1 != o2 && o3 != o4) return true;
    if (o1 == 0 && RefInBox(c, a, b)) return true;
    if (o2 == 0 && RefInBox(d, a, b)) return true;
    if (o3 == 0 && RefInBox(a, c, d)) return true;
    if (o4 == 0 && RefInBox(b, c, d)) return true;
    return false;
  }

  static Path64 DedupPath(const Path64& path)
  {
    Path64 pts;
    for (const Point64& pt : path)
      if (pts.empty() || pts.back() != pt) pts.push_back(pt);
    if (pts.size() > 1 && pts.front() == pts.back()) pts.pop_back();
    return pts;
  }

  // O(n^2) reference. Same documented semantics: touching/crossing
  // non-adjacent edges count, adjacent edges count only when they retrace.
  static size_t RefCountSelfIntersections(const Path64& path)
  {
    const Path64 pts = DedupPath(path);
    const size_t n = pts.size();
    if (n < 3) return 0;
    size_t count = 0;
    for (size_t i = 0; i < n; ++i)
    {
      const size_t i2 = (i + 1) % n;
      for (size_t j = i + 1; j < n; ++j)
      {
        const size_t j2 = (j + 1) % n;
        const bool adjacent = (i == j) || (i2 == j) || (j2 == i) ||
          (i2 == j2);
        if (!adjacent)
        {
          if (RefSegTouch(pts[i], pts[i2], pts[j], pts[j2])) ++count;
        }
        else if (i2 == j || j2 == i) // consecutive, share exactly one vertex
        {
          // retrace when collinear and pointing the same way from the
          // shared vertex
          size_t shared, other0, other1;
          if (i2 == j) { shared = j; other0 = i; other1 = j2; }
          else { shared = i; other0 = j; other1 = i2; }
          if (RefOrient(pts[other0], pts[shared], pts[other1]) == 0)
          {
            const long long dx1 = pts[other0].x - pts[shared].x;
            const long long dy1 = pts[other0].y - pts[shared].y;
            const long long dx2 = pts[other1].x - pts[shared].x;
            const long long dy2 = pts[other1].y - pts[shared].y;
            if (dx1 * dx2 + dy1 * dy2 > 0) ++count;
          }
        }
      }
    }
    return count;
  }

  // ---------------------------------------------------------------------
  // Printing helpers
  // ---------------------------------------------------------------------

  double MsSince(std::chrono::steady_clock::time_point start)
  {
    return std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - start).count();
  }

  // ---------------------------------------------------------------------
  // 1. Statistics against per-point reference
  // ---------------------------------------------------------------------

  struct Sample
  {
    std::string name;
    Paths64 paths;
  };

  Path64 RegularPolygon(int64_t radius, int count, double start_deg)
  {
    Path64 p;
    p.reserve(count);
    for (int i = 0; i < count; ++i)
    {
      const double a = (start_deg + 360.0 * i / count) * PI / 180.0;
      p.emplace_back(static_cast<int64_t>(std::llround(radius * std::cos(a))),
        static_cast<int64_t>(std::llround(radius * std::sin(a))));
    }
    return p;
  }

  void TestStatistics()
  {
    std::printf("\n=== 1. Area / perimeter / vertex stats vs per-point loops "
      "===\n");

    std::vector<Sample> samples;

    {
      Sample s;
      s.name = "square 100x100 (1 path)";
      s.paths.push_back(MakePath({ 0,0, 100,0, 100,100, 0,100 }));
      samples.push_back(s);
    }
    {
      Sample s;
      s.name = "L-shape (1 path)";
      s.paths.push_back(MakePath({ 0,0, 60,0, 60,20, 20,20, 20,60, 0,60 }));
      samples.push_back(s);
    }
    {
      // pentagram: {5/2} star polygon drawn as 5 vertices in step order
      Sample s;
      s.name = "pentagram (1 path, self-intersecting)";
      Path64 star = RegularPolygon(100, 5, 90.0);
      Path64 ordered = { star[0], star[2], star[4], star[1], star[3] };
      s.paths.push_back(ordered);
      samples.push_back(s);
    }
    {
      // outer CCW square + inner CW square (hole)
      Sample s;
      s.name = "square with square hole (2 paths)";
      s.paths.push_back(MakePath({ 0,0, 100,0, 100,100, 0,100 }));
      s.paths.push_back(MakePath({ 25,25, 25,75, 75,75, 75,25 }));
      samples.push_back(s);
    }

    const double tol = 1e-9;
    for (const Sample& s : samples)
    {
      long double ref_area = 0.0L, ref_perim = 0.0L;
      size_t ref_verts = 0;
      for (const Path64& p : s.paths)
      {
        ref_area += RefShoelace(p);
        ref_perim += RefPerimeter(p);
        ref_verts += p.size();
      }
      const PathStats st = GetPathStats(s.paths);
      const double d_area = std::fabs((double)st.area - (double)ref_area);
      const double d_perim =
        std::fabs((double)st.perimeter - (double)ref_perim);

      std::printf("\n%-8s %s\n", "sample:", s.name.c_str());
      std::printf("  paths=%zu vertices=%zu (ref %zu)\n",
        st.path_count, st.vertex_count, ref_verts);
      std::printf("  area      = %.12f   ref = %.12Lf   |diff| = %.3e\n",
        st.area, ref_area, d_area);
      std::printf("  perimeter = %.12f   ref = %.12Lf   |diff| = %.3e\n",
        st.perimeter, ref_perim, d_perim);
      std::printf("  abs_area  = %.12f   signed area >= 0: %s\n",
        st.abs_area, st.is_positive ? "yes" : "no");

      Check(st.path_count == s.paths.size(), "path_count");
      Check(st.vertex_count == ref_verts, "vertex_count");
      Check(d_area < tol, "area within 1e-9");
      Check(d_perim < tol, "perimeter within 1e-9");
      Check(std::fabs(st.area - Area(s.paths)) == 0.0,
        "GetPathStats.area == Area(paths)");
      Check(std::fabs(st.perimeter - Perimeter(s.paths)) == 0.0,
        "GetPathStats.perimeter == Perimeter(paths)");
      Check(st.is_positive == (st.area >= 0.0), "is_positive");
    }

    // Single path overload and explicit open-path perimeter.
    const Path64 line = MakePath({ 0,0, 3,4, 3,8 });
    std::printf("\nopen path (0,0)-(3,4)-(3,8): perimeter open=%f closed=%f\n",
      Perimeter(line, false), Perimeter(line, true));
    Check(std::fabs(Perimeter(line, false) - 9.0) < 1e-12, "open perimeter");
    // closed path adds the 3-4-5 straight-line distance back to the start
    Check(std::fabs(Perimeter(line, true) - (9.0 + std::sqrt(73.0))) < 1e-12,
      "closed perimeter");
    Check(VertexCount(line) == 3, "VertexCount(Path)");
    Check(GetPathStats(line, false).perimeter == Perimeter(line, false),
      "GetPathStats(path, is_closed_path=false)");

    // Floating point (PathD) stats.
    const PathD dsq = { {0.0,0.0}, {2.5,0.0}, {2.5,4.0}, {0.0,4.0} };
    const PathStats dst = GetPathStats(PathsD{ dsq });
    std::printf("PathD square 2.5x4: area=%.6f perimeter=%.6f vertices=%zu\n",
      dst.area, dst.perimeter, dst.vertex_count);
    Check(std::fabs(dst.area - 10.0) < 1e-12, "PathD area");
    Check(std::fabs(dst.perimeter - 13.0) < 1e-12, "PathD perimeter");
  }

  // ---------------------------------------------------------------------
  // 2. Self-intersection detection: hand-built cases + random agreement
  // ---------------------------------------------------------------------

  struct SICase
  {
    std::string name;
    Path64 path;
    bool expect_self;
    int expect_count; // -1 when only the boolean is asserted
  };

  void TestSelfIntersectionCases()
  {
    std::printf("\n=== 2a. Self-intersection on hand-built polygons ===\n");

    std::vector<SICase> cases;

    cases.push_back({ "square (closing point absent)", MakePath(
      { 0,0, 10,0, 10,10, 0,10 }), false, 0 });
    cases.push_back({ "square (first point repeated at end)", MakePath(
      { 0,0, 10,0, 10,10, 0,10, 0,0 }), false, 0 });
    {
      SICase c;
      c.name = "convex 8-gon";
      c.path = RegularPolygon(100, 8, 0.0);
      c.expect_self = false;
      c.expect_count = 0;
      cases.push_back(c);
    }
    cases.push_back({ "L-shape (concave)", MakePath(
      { 0,0, 60,0, 60,20, 20,20, 20,60, 0,60 }), false, 0 });
    cases.push_back({ "U-shape (concave)", MakePath(
      { 0,0, 90,0, 90,60, 60,60, 60,20, 30,20, 30,60, 0,60 }),
      false, 0 });
    cases.push_back({ "rectangle with collinear midpoints", MakePath(
      { 0,0, 50,0, 100,0, 100,50, 100,100, 50,100, 0,100, 0,50 }),
      false, 0 });
    cases.push_back({ "clockwise square", MakePath(
      { 0,0, 0,10, 10,10, 10,0 }), false, 0 });
    cases.push_back({ "bow tie", MakePath(
      { 0,0, 10,10, 10,0, 0,10 }), true, 1 });
    {
      Path64 star = RegularPolygon(100, 5, 90.0);
      cases.push_back({ "pentagram", { star[0], star[2], star[4],
        star[1], star[3] }, true, 5 });
    }
    // the spike retraces one edge, and because the spike tip is 10,10 the
    // path also revisits that vertex, which touches two earlier edges
    cases.push_back({ "square with spike", MakePath(
      { 0,0, 10,0, 10,10, 15,10, 10,10, 0,10 }), true, 4 });
    cases.push_back({ "edge crossing long chord", MakePath(
      { 0,0, 10,0, 10,10, 5,10, 5,-5, 4,-5, 4,10, 0,10 }), true, -1 });
    cases.push_back({ "collinear edges overlapping", MakePath(
      { 0,0, 10,0, 10,10, 10,5, 0,5, 0,0 }), true, -1 });

    size_t tp = 0, fp = 0, fn = 0, tn = 0;
    std::printf("%-42s %8s %8s %8s\n", "case", "expect", "actual", "ok");
    for (const SICase& c : cases)
    {
      const size_t actual = CountSelfIntersections(c.path);
      const size_t reference = RefCountSelfIntersections(c.path);
      const bool actual_bool = actual > 0;
      if (c.expect_self && actual_bool) ++tp;
      else if (!c.expect_self && actual_bool) ++fp;
      else if (c.expect_self && !actual_bool) ++fn;
      else ++tn;

      std::string expected = c.expect_count >= 0 ?
        std::to_string(c.expect_count) : (c.expect_self ? ">0" : "0");
      std::printf("%-42s %8s %8zu %8s\n", c.name.c_str(), expected.c_str(),
        actual, (actual_bool == c.expect_self &&
          (c.expect_count < 0 || actual == (size_t)c.expect_count)) ?
          "ok" : "MISMATCH");
      std::printf("%44s reference O(n^2) count = %zu\n", "", reference);

      Check(actual_bool == c.expect_self, c.name + ": boolean");
      if (c.expect_count >= 0)
        Check(actual == (size_t)c.expect_count, c.name + ": count");
      Check(actual == reference, c.name + ": matches O(n^2) reference");
      Check(IsSelfIntersecting(c.path) == actual_bool, c.name + ": Is*");
    }
    std::printf("\nhand-built cases: TP=%zu FP=%zu FN=%zu TN=%zu\n",
      tp, fp, fn, tn);

    // double precision branch
    {
      const PathD square = { {0.0,0.0}, {1.0,0.0}, {1.0,1.0}, {0.0,1.0} };
      const PathD bowtie = { {0.0,0.0}, {1.0,1.0}, {1.0,0.0}, {0.0,1.0} };
      std::printf("PathD square count=%zu, PathD bow tie count=%zu\n",
        CountSelfIntersections(square), CountSelfIntersections(bowtie));
      Check(CountSelfIntersections(square) == 0, "PathD square");
      Check(CountSelfIntersections(bowtie) == 1, "PathD bow tie");
    }

    // Paths overloads: per-path counts add up, cross-path intersections
    // are (by design) not reported.
    {
      Paths64 two;
      two.push_back(MakePath({ 0,0, 10,0, 10,10, 0,10 }));
      two.push_back(MakePath({ 20,20, 30,20, 30,30, 20,30 }));
      Check(CountSelfIntersections(two) == 0, "Paths overload clean");
      two.push_back(MakePath({ 0,0, 30,30, 30,0, 0,30 }));
      Check(CountSelfIntersections(two) == 1, "Paths overload sums per path");
      Check(IsSelfIntersecting(two), "IsSelfIntersecting(Paths)");
    }
  }

  void TestSelfIntersectionRandom()
  {
    std::printf("\n=== 2b. Random polygons: agreement with O(n^2) reference "
      "===\n");

    uint64_t state = 0x9E3779B97F4A7C15ULL;
    const auto next = [&state]()
    {
      state = state * 6364136223846793005ULL + 1442695040888963407ULL;
      return (state >> 33);
    };

    size_t tp = 0, fp = 0, fn = 0, tn = 0, count_mismatch = 0;
    const int kSimple = 40, kRandom = 40;
    for (int i = 0; i < kSimple; ++i)
    {
      // star-shaped about the origin => guaranteed simple
      Path64 p;
      const int n = 12;
      for (int k = 0; k < n; ++k)
      {
        const double a = (360.0 * k / n) * PI / 180.0;
        const double r = 50.0 + static_cast<double>(next() % 5000) / 100.0;
        p.emplace_back(static_cast<int64_t>(std::llround(r * std::cos(a))),
          static_cast<int64_t>(std::llround(r * std::sin(a))));
      }
      const size_t actual = CountSelfIntersections(p);
      const size_t reference = RefCountSelfIntersections(p);
      if (reference == 0 && actual == 0) ++tn;
      else if (reference == 0 && actual > 0) ++fp;
      else if (reference > 0 && actual == 0) ++fn;
      else ++tp;
      if (actual != reference) ++count_mismatch;
      Check(actual == 0, "star-shaped polygon must be simple");
    }
    for (int i = 0; i < kRandom; ++i)
    {
      Path64 p;
      const int n = 12;
      for (int k = 0; k < n; ++k)
        p.emplace_back(static_cast<int64_t>(next() % 1001) - 500,
          static_cast<int64_t>(next() % 1001) - 500);
      const size_t actual = CountSelfIntersections(p);
      const size_t reference = RefCountSelfIntersections(p);
      const bool a = actual > 0, r = reference > 0;
      if (a && r) ++tp; else if (a && !r) ++fp;
      else if (!a && r) ++fn; else ++tn;
      if (actual != reference) ++count_mismatch;
    }
    const size_t pos = tp + fn, neg = fp + tn;
    std::printf("polygons: %d simple + %d random, 12 vertices each\n",
      kSimple, kRandom);
    std::printf("confusion: TP=%zu FP=%zu FN=%zu TN=%zu\n", tp, fp, fn, tn);
    std::printf("recall    = TP/(TP+FN) = %.4f\n",
      pos ? (double)tp / pos : 1.0);
    std::printf("precision = TP/(TP+FP) = %.4f\n",
      (tp + fp) ? (double)tp / (tp + fp) : 1.0);
    std::printf("false-positive rate = FP/(FP+TN) = %.4f\n",
      neg ? (double)fp / neg : 0.0);
    std::printf("exact-count mismatches vs O(n^2) reference: %zu / %d\n",
      count_mismatch, kSimple + kRandom);
    Check(count_mismatch == 0, "random polygons agree with reference");
    Check(fp == 0, "no false positives on simple polygons");
    Check(fn == 0, "no false negatives on self-intersecting polygons");
  }

  // ---------------------------------------------------------------------
  // 3. Performance at 100k vertices
  // ---------------------------------------------------------------------

  double TimeMs(const std::function<void()>& fn, int reps)
  {
    std::vector<double> times;
    times.reserve(reps);
    for (int i = 0; i < reps; ++i)
    {
      const auto t0 = std::chrono::steady_clock::now();
      fn();
      times.push_back(MsSince(t0));
    }
    std::sort(times.begin(), times.end());
    return times[times.size() / 2]; // median
  }

  void TestPerformance()
  {
    std::printf("\n=== 3. Performance with 100,000 vertices ===\n");
    Path64 big = RegularPolygon(1000000, 100000, 0.0);
    std::printf("polygon vertices = %zu\n", big.size());

    volatile double sink = 0.0;
    const double stats_ms = TimeMs([&]()
    {
      const PathStats s = GetPathStats(big);
      sink = s.area + s.perimeter + (double)s.vertex_count;
    }, 21);
    std::printf("GetPathStats (area + perimeter + vertex count): "
      "median %.3f ms\n", stats_ms);

    const double area_ms = TimeMs([&]()
    {
      sink = Area(big);
    }, 21);
    std::printf("Area only:                                        "
      "median %.3f ms\n", area_ms);

    const double perim_ms = TimeMs([&]()
    {
      sink = Perimeter(big);
    }, 21);
    std::printf("Perimeter only:                                   "
      "median %.3f ms\n", perim_ms);

    const size_t si = CountSelfIntersections(big);
    const double si_ms = TimeMs([&]()
    {
      sink = static_cast<double>(CountSelfIntersections(big));
    }, 5);
    std::printf("CountSelfIntersections (simple 100k-gon, result=%zu): "
      "median %.3f ms\n", si, si_ms);

    std::printf("reference area (long double) = %.6Lf\n", RefShoelace(big));
    std::printf("stats area                   = %.6f\n", GetPathStats(big).area);

    Check(stats_ms < 10.0, "GetPathStats(100k vertices) < 10 ms");
    Check(si == 0, "100k-gon is not self-intersecting");
    (void)sink;
  }

  // ---------------------------------------------------------------------
  // 4. Area conservation against clipping output
  // ---------------------------------------------------------------------

  void PrintAreaCheck(const char* name, const AreaCheckResult& r)
  {
    std::printf("%-46s expected=%14.6f actual=%14.6f "
      "abs_err=%.3e rel_err=%.3e tol=%.1e -> %s\n",
      name, r.expected_area, r.actual_area, r.abs_error, r.rel_error,
      r.tolerance, r.passed ? "PASS" : "FAIL");
  }

  void TestAreaConservation()
  {
    std::printf("\n=== 4. Area conservation with configurable tolerance ===\n");

    // (a) A rectangle split into two tiles must have the same area after
    // a union.
    const Paths64 rect = { MakePath({ 0,0, 100,0, 100,60, 0,60 }) };
    const Paths64 left = { MakePath({ 0,0, 40,0, 40,60, 0,60 }) };
    const Paths64 right = { MakePath({ 40,0, 100,0, 100,60, 40,60 }) };
    const Paths64 joined = Union(left, right, FillRule::NonZero);
    const AreaCheckResult r1 =
      CheckAreaConservation(rect, joined, 1e-9);
    PrintAreaCheck("(a) union of two tiles vs rectangle", r1);
    Check(r1.passed, "area conserved by union");

    // (b) Boolean identity: intersect(P, P) == P.
    const Paths64 lshape = {
      MakePath({ 0,0, 60,0, 60,20, 20,20, 20,60, 0,60 }) };
    const Paths64 self_hit = Intersect(lshape, lshape, FillRule::NonZero);
    const AreaCheckResult r2 =
      CheckAreaConservation(Area(lshape), self_hit, 1e-9);
    PrintAreaCheck("(b) intersect(P,P) vs P", r2);
    Check(r2.passed, "area conserved by intersect");

    // (c) XOR of a shape with itself is empty: expected area 0.
    const Paths64 empty = Xor(lshape, lshape, FillRule::NonZero);
    const AreaCheckResult r3 =
      CheckAreaConservation(Area(lshape), empty, 1e-9);
    PrintAreaCheck("(c) xor(P,P) vs P (expected failure)", r3);
    Check(!r3.passed, "xor(P,P) is not area conserving");

    // (d) Same comparison as (a) but with one vertex moved by 10 units.
    const Paths64 moved = { MakePath({ 0,0, 110,0, 110,60, 0,60 }) };
    const AreaCheckResult tight =
      CheckAreaConservation(rect, moved, 1e-9);
    const AreaCheckResult loose =
      CheckAreaConservation(rect, moved, 1000.0);
    PrintAreaCheck("(d) moved vertex, tol = 1e-9", tight);
    PrintAreaCheck("(d) moved vertex, tol = 1e+3", loose);
    Check(!tight.passed, "tight tolerance rejects changed area");
    Check(loose.passed, "loose tolerance accepts changed area");
    Check(std::fabs(tight.abs_error - 600.0) < 1e-9,
      "moved-vertex area delta is 600");

    // (e) Double precision pipeline.
    const PathsD d1 = { { {0.0,0.0}, {2.5,0.0}, {2.5,4.0}, {0.0,4.0} } };
    const PathsD d2 = { { {2.5,0.0}, {5.0,0.0}, {5.0,4.0}, {2.5,4.0} } };
    const PathsD dj = Union(d1, d2, FillRule::NonZero);
    const AreaCheckResult r5 = CheckAreaConservation(
      PathsD{ { {0.0,0.0}, {5.0,0.0}, {5.0,4.0}, {0.0,4.0} } },
      dj, 1e-9);
    PrintAreaCheck("(e) PathD union of halves vs 5x4 rectangle", r5);
    Check(r5.passed, "PathD area conserved");

    // (f) Perform a clipping pipeline and check every output is simple.
    const Paths64 a = { MakePath({ 0,0, 100,0, 100,100, 0,100 }) };
    const Paths64 b = { MakePath({ 50,50, 150,50, 150,150, 50,150 }) };
    const Paths64 u = Union(a, b, FillRule::NonZero);
    const Paths64 x = Xor(a, b, FillRule::NonZero);
    const Paths64 d = Difference(a, b, FillRule::NonZero);
    const Paths64 n = Intersect(a, b, FillRule::NonZero);
    const size_t bad = CountSelfIntersections(u) +
      CountSelfIntersections(x) + CountSelfIntersections(d) +
      CountSelfIntersections(n);
    std::printf("clipping outputs: union paths=%zu xor=%zu diff=%zu "
      "intersection=%zu, self-intersections=%zu\n",
      u.size(), x.size(), d.size(), n.size(), bad);
    Check(bad == 0, "clipping outputs are self-intersection free");

    // Conservation across a clip: area(A) = area(A∩B) + area(A\\B)
    const AreaCheckResult r6 = CheckAreaConservation(
      Area(a), Intersect(a, b, FillRule::NonZero), 1e-9);
    const AreaCheckResult r7 = CheckAreaConservation(
      Area(a), Difference(a, b, FillRule::NonZero), 1e-9);
    std::printf("(g) area(A)=%f, area(A&B)=%f + area(A-B)=%f = %f\n",
      Area(a), r6.actual_area, r7.actual_area,
      r6.actual_area + r7.actual_area);
    Check(std::fabs(r6.actual_area + r7.actual_area - Area(a)) < 1e-9,
      "area partition identity");
  }
} // namespace

int main()
{
  std::printf("Clipper2 path statistics / validation self-test\n");
  std::printf("(standalone: links only the Clipper2 sources in this repo)\n");

  TestStatistics();
  TestSelfIntersectionCases();
  TestSelfIntersectionRandom();
  TestPerformance();
  TestAreaConservation();

  std::printf("\n----------------------------------------------\n");
  std::printf("checks: %d, failures: %d -> %s\n", g_checks, g_failures,
    g_failures == 0 ? "ALL PASSED" : "FAILED");
  return g_failures == 0 ? 0 : 1;
}
