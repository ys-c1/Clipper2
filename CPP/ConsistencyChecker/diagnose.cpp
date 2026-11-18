// Diagnostic: for each mismatching grid point, print the distance to the
// nearest solution edge and to the nearest input edge, to establish
// whether mismatches are snap-rounding boundary artifacts.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>

#include "checker.h"
#include "pip_reference.h"

using namespace refcheck;
using Clipper2Lib::Paths64;
using Clipper2Lib::Point64;
using Clipper2Lib::BooleanOp;
using Clipper2Lib::ClipType;
using Clipper2Lib::FillRule;
using Clipper2Lib::PointInPolygon;
using Clipper2Lib::PointInPolygonResult;

static ClipType ToClip(BoolOp o) {
  switch (o) {
    case BoolOp::kIntersection: return ClipType::Intersection;
    case BoolOp::kUnion:        return ClipType::Union;
    case BoolOp::kDifference:   return ClipType::Difference;
    case BoolOp::kXor:          return ClipType::Xor;
  }
  return ClipType::Intersection;
}

static bool RefInside(BoolOp op, PipResult ra, PipResult rb) {
  bool a = (ra == PipResult::kInside), b = (rb == PipResult::kInside);
  switch (op) {
    case BoolOp::kIntersection: return a && b;
    case BoolOp::kUnion:        return a || b;
    case BoolOp::kDifference:   return a && !b;
    case BoolOp::kXor:          return a != b;
  }
  return false;
}

static Paths Convert(const Paths64& in) {
  Paths out;
  for (const auto& p : in) {
    Path rp;
    for (const auto& v : p) rp.emplace_back(v.x, v.y);
    out.push_back(std::move(rp));
  }
  return out;
}

static double DistToPaths(const Point64& p, const Paths64& paths) {
  double best = 1e300;
  for (const auto& path : paths) {
    for (size_t i = 0, j = path.size() - 1; i < path.size(); j = i++) {
      double ax = path[j].x, ay = path[j].y;
      double bx = path[i].x, by = path[i].y;
      double dx = bx - ax, dy = by - ay;
      double len2 = dx * dx + dy * dy;
      double t = len2 ? ((p.x - ax) * dx + (p.y - ay) * dy) / len2 : 0;
      t = std::max(0.0, std::min(1.0, t));
      double ex = ax + t * dx, ey = ay + t * dy;
      best = std::min(best, std::hypot(p.x - ex, p.y - ey));
    }
  }
  return best;
}

int main() {
  CheckerConfig cfg;
  std::mt19937_64 rng(cfg.seed);
  const BoolOp ops[4] = {BoolOp::kIntersection, BoolOp::kUnion,
                         BoolOp::kDifference, BoolOp::kXor};
  for (int ci = 0; ci < cfg.cases; ++ci) {
    Paths64 sub, clip;
    GenerateCase(rng, cfg, sub, clip);
    CaseReport rep = RunCase(ci, sub, clip, cfg);
    if (rep.total_mismatches == 0) continue;

    Paths rs = Convert(sub), rc = Convert(clip);
    Paths64 all = sub;
    all.insert(all.end(), clip.begin(), clip.end());

    for (int k = 0; k < 4; ++k) {
      if (rep.ops[k].mismatches == 0) continue;
      Paths64 sol = BooleanOp(ToClip(ops[k]), FillRule::EvenOdd, sub, clip);
      double gx = double(rep.max_x - rep.min_x) / (cfg.grid_n - 1);
      double gy = double(rep.max_y - rep.min_y) / (cfg.grid_n - 1);
      for (int row = 0; row < cfg.grid_n; ++row)
        for (int col = 0; col < cfg.grid_n; ++col) {
          int64_t ix = std::llround(rep.min_x + col * gx);
          int64_t iy = std::llround(rep.min_y + row * gy);
          Point rp(ix, iy);
          PipResult ra = PointInPaths(rp, rs), rb = PointInPaths(rp, rc);
          if (ra == PipResult::kOn || rb == PipResult::kOn) continue;
          bool ref_in = RefInside(ops[k], ra, rb);
          int w = 0;
          bool api_on = false;
          for (const auto& path : sol) {
            auto r = PointInPolygon(Point64(ix, iy), path);
            if (r == PointInPolygonResult::IsOn) api_on = true;
            if (r == PointInPolygonResult::IsInside) ++w;
          }
          bool api_in = api_on || (w & 1);
          if (ref_in == api_in) continue;
          std::printf(
              "case %d %-12s pt=(%lld,%lld) ref=%d api=%d "
              "dist_solution_edge=%.3f dist_input_edge=%.3f grid=%.1f\n",
              ci, ToString(ops[k]), (long long)ix, (long long)iy,
              (int)ref_in, (int)api_in,
              DistToPaths(Point64(ix, iy), sol),
              DistToPaths(Point64(ix, iy), all), (gx + gy) / 2);
        }
    }
  }
}
