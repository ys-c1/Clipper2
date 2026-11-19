/*******************************************************************************
* Standalone acceptance test for CleanPolygons (topology cleanup).            *
* Not part of the CMake build — compile & run directly:                       *
*                                                                             *
*   g++ -std=c++17 -O2 -I CPP/Clipper2Lib/include                             *
*       CPP/Tests/CleanTopology.cpp                                           *
*       CPP/Clipper2Lib/src/clipper.engine.cpp                                *
*       CPP/Clipper2Lib/src/clipper.offset.cpp                                *
*       CPP/Clipper2Lib/src/clipper.rectclip.cpp                              *
*       CPP/Clipper2Lib/src/clipper.triangulation.cpp                         *
*       -o /tmp/clean_topology && /tmp/clean_topology                         *
*******************************************************************************/

#include "clipper2/clipper.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>

using namespace Clipper2Lib;

static int failures = 0;

static void Check(bool ok, const std::string& msg)
{
  if (!ok) { ++failures; std::printf("  FAIL: %s\n", msg.c_str()); }
}

// signed area (shoelace), independent of the library's Area()
static double SignedArea(const Path64& p)
{
  double a = 0;
  for (size_t i = 0, j = p.size() - 1; i < p.size(); j = i++)
    a += (double)p[j].x * p[i].y - (double)p[i].x * p[j].y;
  return a * 0.5;
}

static double TotalArea(const Paths64& ps)
{
  double a = 0;
  for (const auto& p : ps) a += SignedArea(p);
  return a;
}

// segments intersect properly (crossing, not just shared endpoints)?
static bool SegsCross(const Point64& a1, const Point64& a2,
  const Point64& b1, const Point64& b2)
{
  auto cross = [](const Point64& o, const Point64& u, const Point64& v) {
    return (double)(u.x - o.x) * (v.y - o.y) - (double)(u.y - o.y) * (v.x - o.x);
  };
  double d1 = cross(b1, b2, a1), d2 = cross(b1, b2, a2);
  double d3 = cross(a1, a2, b1), d4 = cross(a1, a2, b2);
  return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
         ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
}

// a polygon is simple when no two non-adjacent edges cross
static bool IsSimple(const Path64& p)
{
  const size_t n = p.size();
  for (size_t i = 0; i < n; ++i)
    for (size_t j = i + 1; j < n; ++j)
    {
      if (j == i || (j + 1) % n == i || (i + 1) % n == j) continue; // adjacent
      if (SegsCross(p[i], p[(i + 1) % n], p[j], p[(j + 1) % n])) return false;
    }
  return true;
}

static bool HasZeroLengthEdge(const Path64& p)
{
  for (size_t i = 0, j = p.size() - 1; i < p.size(); j = i++)
    if (p[i] == p[j]) return true;
  return false;
}

static bool HasCollinearVertex(const Path64& p)
{
  const size_t n = p.size();
  for (size_t i = 0; i < n; ++i)
  {
    const Point64& a = p[(i + n - 1) % n], & b = p[i], & c = p[(i + 1) % n];
    double cp = (double)(b.x - a.x) * (c.y - b.y) - (double)(b.y - a.y) * (c.x - b.x);
    if (cp == 0) return true;
  }
  return false;
}

static void Report(const char* name, const Paths64& before, const Paths64& after,
  double ref_area)
{
  size_t nv_before = 0, nv_after = 0;
  for (const auto& p : before) nv_before += p.size();
  for (const auto& p : after) nv_after += p.size();
  double area_before = TotalArea(before), area_after = TotalArea(after);
  double err = std::fabs(area_after - ref_area);
  bool winding_ok = (area_before == 0) ||
    ((area_before > 0) == (area_after > 0));
  bool simple = true, zero_edge = false, collinear = false;
  for (const auto& p : after)
  {
    simple = simple && IsSimple(p);
    zero_edge = zero_edge || HasZeroLengthEdge(p);
    collinear = collinear || HasCollinearVertex(p);
  }
  std::printf("%-22s polys %zu->%-3zu verts %4zu->%-4zu area %.1f->%.1f "
    "err=%.2e winding:%s simple:%s\n",
    name, before.size(), after.size(), nv_before, nv_after,
    area_before, area_after, err,
    winding_ok ? "OK" : "FLIPPED", simple ? "yes" : "NO");
  Check(err < 1e-9, std::string(name) + ": area error >= 1e-9");
  Check(winding_ok, std::string(name) + ": winding changed");
  Check(simple, std::string(name) + ": output not simple");
  Check(!zero_edge, std::string(name) + ": zero-length edge remains");
  Check(!collinear, std::string(name) + ": collinear vertex remains");
}

int main()
{
  // -- base shape: concave 8-vertex polygon (CCW); coords are multiples of 4
  // so the collinear mid/quarter points inserted below are exact integers ----
  const Path64 base{ {0,0},{32,0},{32,16},{16,16},{16,32},{0,32},{0,16},{8,16} };

  // group 1: normal (already clean) ----------------------------------------
  Report("1 normal", { base }, CleanPolygons({ base }), SignedArea(base));

  // group 2: degenerate (zero-length edges) ---------------------------------
  Path64 degen;
  for (const auto& pt : base) { degen.push_back(pt); degen.push_back(pt); }
  degen.push_back(degen.front()); // duplicated closing vertex
  Report("2 degenerate edges", { degen }, CleanPolygons({ degen }), SignedArea(base));

  // group 3: collinear vertices ---------------------------------------------
  Path64 coll;
  for (size_t i = 0, j = base.size() - 1; i < base.size(); j = i++)
  {
    coll.push_back(base[j]);
    coll.emplace_back((base[j].x + base[i].x) / 2, (base[j].y + base[i].y) / 2);
    coll.emplace_back((base[j].x + 3 * base[i].x) / 4, (base[j].y + 3 * base[i].y) / 4);
  }
  Report("3 collinear points", { coll }, CleanPolygons({ coll }), SignedArea(base));

  // group 3b: near-collinear with configurable tolerance ---------------------
  Path64 near_coll = base;
  near_coll.insert(near_coll.begin() + 1, { 24, 4 });  // spike 4 units off edge
  near_coll.insert(near_coll.begin() + 3, { 32, 8 });  // exactly on edge
  Paths64 cleaned_tol = CleanPolygons({ near_coll }, 4.5); // tolerance absorbs it
  size_t nv = 0; for (const auto& p : cleaned_tol) nv += p.size();
  std::printf("%-22s polys 1->%-3zu verts %4zu->%-4zu (tolerance 4.5 absorbs "
    "near-collinear spike)\n", "3b tolerance=4.5", cleaned_tol.size(), nv, nv);
  Check(nv == base.size(), "3b: tolerance did not absorb near-collinear vertex");
  Paths64 cleaned_tol0 = CleanPolygons({ near_coll }, 0.0);
  size_t nv0 = 0; for (const auto& p : cleaned_tol0) nv0 += p.size();
  std::printf("%-22s polys 1->%-3zu verts %4zu->%-4zu (tolerance 0 keeps "
    "the spike)\n", "3b tolerance=0", cleaned_tol0.size(), nv0, nv0);
  Check(nv0 == base.size() + 1, "3b: zero tolerance should keep the spike");

  // group 4: self-intersecting bowtie -> two triangles, region area = 4 -----
  const Path64 bowtie{ {0,0},{4,2},{4,0},{0,2} };
  const double bowtie_region_area = 4.0; // two triangles of area 2 each
  Paths64 bt_out = CleanPolygons({ bowtie });
  Report("4 self-intersecting", { bowtie }, bt_out, bowtie_region_area);
  Check(bt_out.size() == 2, "4: bowtie should split into 2 simple polygons");

  // group 4b: clockwise (negative winding) input must stay negative ---------
  Path64 cw_base(base.rbegin(), base.rend());
  Paths64 cw_out = CleanPolygons({ cw_base });
  Report("4b clockwise input", { cw_base }, cw_out, SignedArea(cw_base));
  bool all_neg = true;
  for (const auto& p : cw_out) all_neg = all_neg && SignedArea(p) < 0;
  Check(all_neg, "4b: negative winding not preserved");

  // performance: 10k-vertex gear with degeneracies ---------------------------
  Path64 gear;
  const int teeth = 2500; // 2500 * 4 = 10000 vertices
  for (int i = 0; i < teeth; ++i)
  {
    double a0 = 2 * PI * i / teeth, a1 = 2 * PI * (i + 0.5) / teeth;
    double am = (a0 + a1) / 2;
    int64_t x0 = (int64_t)std::llround(1000000 * std::cos(a0));
    int64_t y0 = (int64_t)std::llround(1000000 * std::sin(a0));
    int64_t x1 = (int64_t)std::llround(1100000 * std::cos(am));
    int64_t y1 = (int64_t)std::llround(1100000 * std::sin(am));
    int64_t x2 = (int64_t)std::llround(1000000 * std::cos(a1));
    int64_t y2 = (int64_t)std::llround(1000000 * std::sin(a1));
    gear.emplace_back(x0, y0);
    gear.emplace_back(x0, y0);          // zero-length edge
    gear.emplace_back(x1, y1);
    gear.emplace_back(x2, y2);
  }
  std::printf("perf input: %zu vertices\n", gear.size());
  auto t0 = std::chrono::steady_clock::now();
  Paths64 gear_out = CleanPolygons({ gear });
  auto t1 = std::chrono::steady_clock::now();
  double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  size_t gv = 0; for (const auto& p : gear_out) gv += p.size();
  std::printf("perf: %.2f ms (%zu output verts, limit 50 ms) %s\n",
    ms, gv, ms < 50 ? "OK" : "TOO SLOW");
  Check(ms < 50, "performance: 10k vertices took >= 50 ms");
  Check(std::fabs(TotalArea(gear_out) - SignedArea(gear)) < 1e-3,
    "performance: gear area changed");

  if (failures == 0) std::printf("\nALL CHECKS PASSED\n");
  else std::printf("\n%d CHECK(S) FAILED\n", failures);
  return failures ? 1 : 0;
}
