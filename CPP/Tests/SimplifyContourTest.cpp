// SimplifyContourTest.cpp
//
// Self-contained acceptance test for contour simplification:
//   1. Ramer-Douglas-Peucker thinning by max deviation (Clipper2's RDP),
//   2. area-constraint correction: |dArea| / area must stay within a given ratio,
//   3. the simplified contour is guaranteed non-self-intersecting
//      (epsilon is reduced / refined until both constraints hold; epsilon -> 0
//      is always feasible for a simple input contour).
//
// Builds standalone (no GoogleTest, no third-party deps):
//   g++ -std=c++17 -O2 -I ../Clipper2Lib/include SimplifyContourTest.cpp -o SimplifyContourTest
//
// Reports, for 4 contours (circle / zigzag / glyph / near-straight line):
// vertices before/after, max deviation, area error, self-intersection count,
// and the time to simplify 10k points (must be < 20 ms).

#include "clipper2/clipper.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace Clipper2Lib;

// ---------------------------------------------------------------------------
// geometry helpers (integer-exact orientation, coords are small enough for i64)
// ---------------------------------------------------------------------------

static long long Cross(const Point64& a, const Point64& b, const Point64& c)
{
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

static bool OnSegment(const Point64& a, const Point64& b, const Point64& p)
{
  return Cross(a, b, p) == 0 &&
    std::min(a.x, b.x) <= p.x && p.x <= std::max(a.x, b.x) &&
    std::min(a.y, b.y) <= p.y && p.y <= std::max(a.y, b.y);
}

static bool EdgesIntersect(const Point64& a, const Point64& b,
  const Point64& c, const Point64& d)
{
  // quick bbox reject
  if (std::max(a.x, b.x) < std::min(c.x, d.x) ||
      std::max(c.x, d.x) < std::min(a.x, b.x) ||
      std::max(a.y, b.y) < std::min(c.y, d.y) ||
      std::max(c.y, d.y) < std::min(a.y, b.y)) return false;
  long long c1 = Cross(a, b, c), c2 = Cross(a, b, d);
  long long c3 = Cross(c, d, a), c4 = Cross(c, d, b);
  if (((c1 > 0 && c2 < 0) || (c1 < 0 && c2 > 0)) &&
      ((c3 > 0 && c4 < 0) || (c3 < 0 && c4 > 0))) return true;
  if (c1 == 0 && OnSegment(a, b, c)) return true;
  if (c2 == 0 && OnSegment(a, b, d)) return true;
  if (c3 == 0 && OnSegment(c, d, a)) return true;
  if (c4 == 0 && OnSegment(c, d, b)) return true;
  return false;
}

// counts intersections between non-adjacent edges of a closed contour
static size_t SelfIntersectionCount(const Path64& path)
{
  const size_t n = path.size();
  size_t count = 0;
  for (size_t i = 0; i < n; ++i)
  {
    const Point64& a = path[i];
    const Point64& b = path[(i + 1) % n];
    for (size_t j = i + 1; j < n; ++j)
    {
      if (j == i + 1) continue;                 // adjacent
      if (i == 0 && j == n - 1) continue;       // adjacent (wrap)
      if (EdgesIntersect(a, b, path[j], path[(j + 1) % n])) ++count;
    }
  }
  return count;
}

static double PointSegDist(const Point64& p, const Point64& a, const Point64& b)
{
  double dx = double(b.x - a.x), dy = double(b.y - a.y);
  double len2 = dx * dx + dy * dy;
  if (len2 == 0) return std::hypot(double(p.x - a.x), double(p.y - a.y));
  double t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / len2;
  t = std::max(0.0, std::min(1.0, t));
  return std::hypot(p.x - a.x - t * dx, p.y - a.y - t * dy);
}

// max distance from any original vertex to the simplified contour
static double MaxDeviation(const Path64& orig, const Path64& simp)
{
  const size_t m = simp.size();
  double maxd = 0;
  for (const Point64& p : orig)
  {
    double d = HUGE_VAL;
    for (size_t i = 0; i < m; ++i)
      d = std::min(d, PointSegDist(p, simp[i], simp[(i + 1) % m]));
    maxd = std::max(maxd, d);
  }
  return maxd;
}

// ---------------------------------------------------------------------------
// simplify: RDP + area-constraint correction + self-intersection guarantee
// ---------------------------------------------------------------------------

// Douglas-Peucker for a CLOSED contour. The library's RamerDouglasPeucker
// treats its input as an open polyline (it anchors only the first/last
// vertices), so the loop is split at two anchor vertices -- vertex 0 and the
// vertex farthest from it -- and each open half is simplified separately.
// Both anchors survive RDP by construction, so the stitched result is closed.
static Path64 ClosedRDP(const Path64& path, double epsilon)
{
  const size_t n = path.size();
  if (n < 5) return Path64(path);

  size_t mid = 0;
  double max_d2 = -1;
  for (size_t i = 1; i < n; ++i)
  {
    double dx = double(path[i].x - path[0].x), dy = double(path[i].y - path[0].y);
    double d2 = dx * dx + dy * dy;
    if (d2 > max_d2) { max_d2 = d2; mid = i; }
  }
  if (mid == 0 || mid == n - 1) mid = n / 2;   // degenerate: any split works

  Path64 half1(path.begin(), path.begin() + mid + 1);
  Path64 half2(path.begin() + mid, path.end());
  half2.push_back(path[0]);

  Path64 r1 = RamerDouglasPeucker(half1, epsilon);
  Path64 r2 = RamerDouglasPeucker(half2, epsilon);

  Path64 result = std::move(r1);
  // r2 repeats the anchor (its first) and vertex 0 (its last): skip both
  for (size_t i = 1; i + 1 < r2.size(); ++i) result.push_back(r2[i]);
  return result;
}

static Path64 SimplifyContour(const Path64& path, double epsilon,
  double maxAreaRatio, double& epsUsed)
{
  const double area0 = std::fabs(Area(path));
  const double maxAreaErr = maxAreaRatio * area0;

  Path64 best;
  auto feasible = [&](double eps, Path64& out)
  {
    Path64 r = ClosedRDP(path, eps);
    if (r.size() < 3) return false;
    if (std::fabs(std::fabs(Area(r)) - area0) > maxAreaErr) return false;
    if (SelfIntersectionCount(r) != 0) return false;
    out = std::move(r);
    return true;
  };

  // phase 1: shrink epsilon until both constraints hold
  double eps = epsilon;
  int guard = 0;
  while (!feasible(eps, best) && guard++ < 200) eps *= 0.5;
  eps = std::max(eps, 1e-12);

  // phase 2: binary-search refine upward to recover as much thinning as
  // possible while staying within the area / self-intersection limits
  double bestEps = eps;
  if (eps < epsilon)
  {
    double lo = eps, hi = epsilon;
    Path64 cand;
    for (int it = 0; it < 25; ++it)
    {
      double mid = 0.5 * (lo + hi);
      if (feasible(mid, cand)) { best = cand; bestEps = mid; lo = mid; }
      else hi = mid;
    }
  }
  epsUsed = bestEps;   // epsilon that produced 'best' (<= epsilon after correction)
  return best;
}

// ---------------------------------------------------------------------------
// sample contours (small, deterministic, generated in-code)
// ---------------------------------------------------------------------------

// tiny deterministic pseudo-random in [-1, 1]
static double Jitter(unsigned i)
{
  unsigned h = i * 2654435761u;
  h ^= h >> 13; h *= 2246822519u; h ^= h >> 16;
  return (h % 2000) / 1000.0 - 1.0;
}

static Path64 MakeCircle(size_t n, double r, double cx, double cy)
{
  Path64 p;
  for (size_t i = 0; i < n; ++i)
  {
    double a = 2.0 * PI * i / n;
    p.emplace_back(Point64(cx + r * std::cos(a), cy + r * std::sin(a)));
  }
  return p;
}

static Path64 MakeZigzag()
{
  Path64 p;
  const int teeth = 500;
  for (int i = 0; i <= teeth * 2; ++i)          // top: sawtooth, left to right
  {
    double x = i * 20.0;
    double y = 9000 + (i % 2 ? 50 : -50) + Jitter(i) * 4;
    p.emplace_back(Point64(x, y));
  }
  for (int i = teeth * 2; i >= 0; i -= 40)      // bottom: straight, right to left
    p.emplace_back(Point64(i * 20.0, 0.0));
  return p;
}

static Path64 MakeGlyph()                        // block letter "E", x10 scale
{
  const double pts[][2] = {
    {0,0},{300,0},{300,60},{60,60},{60,120},{220,120},
    {220,180},{60,180},{60,240},{300,240},{300,300},{0,300}};
  const size_t n = sizeof(pts) / sizeof(pts[0]);
  Path64 p;
  for (size_t i = 0; i < n; ++i)                 // densify every edge
  {
    double x0 = pts[i][0] * 10, y0 = pts[i][1] * 10;
    double x1 = pts[(i + 1) % n][0] * 10, y1 = pts[(i + 1) % n][1] * 10;
    double len = std::hypot(x1 - x0, y1 - y0);
    int steps = std::max(1, int(len / 5));
    for (int s = 0; s < steps; ++s)
    {
      double t = double(s) / steps;
      p.emplace_back(Point64(x0 + t * (x1 - x0), y0 + t * (y1 - y0)));
    }
  }
  return p;
}

static Path64 MakeNearLine()
{
  Path64 p;
  const int n = 5000;
  for (int i = 0; i <= n; ++i)                   // top: near-straight + noise
  {
    double x = i * 4.0;
    double y = 100 + Jitter(i) * 3 + 2 * std::sin(i * 0.01);
    p.emplace_back(Point64(x, y));
  }
  for (int i = n; i >= 0; i -= 250)              // bottom: straight back
    p.emplace_back(Point64(i * 4.0, 0.0));
  return p;
}

// ---------------------------------------------------------------------------

int main()
{
  struct Case { std::string name; Path64 path; double epsilon; };
  std::vector<Case> cases = {
    // circle epsilon is deliberately oversized (200) so the area-constraint
    // correction is forced to pull it back under the 1% area-error limit
    {"circle (10k)",  MakeCircle(10000, 5000, 10000, 10000), 200.0},
    {"zigzag",        MakeZigzag(),                          60.0},
    {"glyph 'E'",     MakeGlyph(),                           8.0},
    {"near-line",     MakeNearLine(),                        5.0},
  };

  const double kMaxAreaRatio = 0.01;             // area error must stay <= 1%
  bool all_pass = true;

  std::printf("%-14s %9s %9s %8s %8s %10s %10s %8s %9s\n",
    "contour", "verts_in", "verts_out", "eps_in", "eps_used", "max_dev",
    "area_err%", "self_int", "time_ms");
  std::printf("------------------------------------------------------------------------------------\n");

  for (auto& c : cases)
  {
    auto t0 = std::chrono::steady_clock::now();
    double epsUsed = 0;
    Path64 out = SimplifyContour(c.path, c.epsilon, kMaxAreaRatio, epsUsed);
    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    double dev = MaxDeviation(c.path, out);
    double aerr = std::fabs(std::fabs(Area(out)) - std::fabs(Area(c.path)))
      / std::fabs(Area(c.path)) * 100.0;
    size_t si = SelfIntersectionCount(out);

    bool pass = aerr <= 1.0 && si == 0 && dev <= c.epsilon * 1.05 &&
      (c.path.size() < 10000 || ms < 20.0);
    all_pass &= pass;

    std::printf("%-14s %9zu %9zu %8.1f %8.2f %10.3f %10.4f %8zu %9.3f  %s\n",
      c.name.c_str(), c.path.size(), out.size(), c.epsilon, epsUsed,
      dev, aerr, si, ms, pass ? "PASS" : "FAIL");
  }

  std::printf("------------------------------------------------------------------------------------\n");
  std::printf("%s\n", all_pass ? "ALL PASS" : "SOME CASES FAILED");
  return all_pass ? 0 : 1;
}
