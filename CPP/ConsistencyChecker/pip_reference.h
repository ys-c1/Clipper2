// Independent reference point-in-polygon classification.
//
// This file deliberately does NOT include any Clipper2 header: it is the
// second, fully independent implementation that the library output is
// cross-checked against. It is a textbook ray-casting ("even-odd rule")
// algorithm written from first principles, using exact 64-bit integer
// cross products for on-edge detection (no floating point anywhere).
#ifndef CONSISTENCY_CHECKER_PIP_REFERENCE_H
#define CONSISTENCY_CHECKER_PIP_REFERENCE_H

#include <cstdint>
#include <vector>

namespace refcheck {

struct Point {
  int64_t x = 0;
  int64_t y = 0;
  Point() = default;
  Point(int64_t x_, int64_t y_) : x(x_), y(y_) {}
};

using Path = std::vector<Point>;
using Paths = std::vector<Path>;

enum class PipResult { kOutside, kInside, kOn };

// Sign of cross product (b-a) x (p-a), computed with __int128 so that
// products of two int64 deltas cannot overflow.
inline int CrossSign(const Point& a, const Point& b, const Point& p) {
  __int128 dx1 = (__int128)b.x - a.x;
  __int128 dy1 = (__int128)b.y - a.y;
  __int128 dx2 = (__int128)p.x - a.x;
  __int128 dy2 = (__int128)p.y - a.y;
  __int128 cross = dx1 * dy2 - dy1 * dx2;
  return (cross > 0) - (cross < 0);
}

// Ray casting along +x against a single closed path.
// Horizontal edges are skipped (standard convention); a zero cross
// product on a crossing edge means the point lies exactly on the edge.
inline PipResult PointInPath(const Point& p, const Path& poly) {
  const size_t n = poly.size();
  if (n < 3) return PipResult::kOutside;
  bool inside = false;
  for (size_t i = 0, j = n - 1; i < n; j = i++) {
    const Point& a = poly[j];
    const Point& b = poly[i];
    // On-segment test first: bounding box plus zero cross product.
    if (CrossSign(a, b, p) == 0 &&
        p.x >= ((a.x < b.x) ? a.x : b.x) &&
        p.x <= ((a.x < b.x) ? b.x : a.x) &&
        p.y >= ((a.y < b.y) ? a.y : b.y) &&
        p.y <= ((a.y < b.y) ? b.y : a.y)) {
      return PipResult::kOn;
    }
    // Edge straddles the horizontal ray y = p.y (half-open interval).
    if ((a.y > p.y) != (b.y > p.y)) {
      // x coordinate of the edge/ray intersection, compared exactly:
      // p.x < (a.x + (p.y-a.y)*(b.x-a.x)/(b.y-a.y))
      __int128 lhs = (__int128)(p.x - a.x) * (b.y - a.y);
      __int128 rhs = (__int128)(p.y - a.y) * (b.x - a.x);
      // Division-free comparison; sign of (b.y - a.y) decides ordering.
      bool intersects = (b.y > a.y) ? (lhs < rhs) : (lhs > rhs);
      if (intersects) inside = !inside;
    }
  }
  return inside ? PipResult::kInside : PipResult::kOutside;
}

// Even-odd classification over a set of paths. A point on any path is
// reported as kOn, otherwise parity of the winding count decides.
inline PipResult PointInPaths(const Point& p, const Paths& paths) {
  int crossings = 0;
  for (const Path& path : paths) {
    PipResult r = PointInPath(p, path);
    if (r == PipResult::kOn) return PipResult::kOn;
    if (r == PipResult::kInside) ++crossings;
  }
  return (crossings & 1) ? PipResult::kInside : PipResult::kOutside;
}

}  // namespace refcheck

#endif  // CONSISTENCY_CHECKER_PIP_REFERENCE_H
