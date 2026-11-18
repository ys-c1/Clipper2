#include "checker.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <unordered_map>

namespace refcheck {

using Clipper2Lib::BooleanOp;
using Clipper2Lib::ClipType;
using Clipper2Lib::FillRule;
using Clipper2Lib::Path64;
using Clipper2Lib::Paths64;
using Clipper2Lib::Point64;
using Clipper2Lib::PointInPolygon;
using Clipper2Lib::PointInPolygonResult;

const char* ToString(BoolOp op) {
  switch (op) {
    case BoolOp::kIntersection: return "intersection";
    case BoolOp::kUnion:        return "union";
    case BoolOp::kDifference:   return "difference";
    case BoolOp::kXor:          return "xor";
  }
  return "?";
}

static ClipType ToClipType(BoolOp op) {
  switch (op) {
    case BoolOp::kIntersection: return ClipType::Intersection;
    case BoolOp::kUnion:        return ClipType::Union;
    case BoolOp::kDifference:   return ClipType::Difference;
    case BoolOp::kXor:          return ClipType::Xor;
  }
  return ClipType::Intersection;
}

// ---------------- random polygon generation -----------------------------

Path64 GenerateSimplePolygon(std::mt19937_64& rng, int vertex_count,
                             int64_t extent) {
  // Radial construction around a random centre: each vertex gets a random
  // angle and radius; sorting by angle guarantees a simple star-shaped
  // polygon (no self-intersections).
  std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
  std::uniform_real_distribution<double> radius_dist(0.25, 1.0);
  std::uniform_int_distribution<int64_t> center_dist(
      extent / 4, extent - extent / 4);
  const double cx = static_cast<double>(center_dist(rng));
  const double cy = static_cast<double>(center_dist(rng));
  const double base_r = extent * 0.18;

  struct V {
    double angle, x, y;
  };
  std::vector<V> verts;
  verts.reserve(vertex_count);
  for (int i = 0; i < vertex_count; ++i) {
    double a = angle_dist(rng);
    double r = base_r * radius_dist(rng);
    verts.push_back({a, cx + r * std::cos(a), cy + r * std::sin(a)});
  }
  std::sort(verts.begin(), verts.end(),
            [](const V& a, const V& b) { return a.angle < b.angle; });

  Path64 out;
  out.reserve(vertex_count);
  for (const V& v : verts)
    out.emplace_back(static_cast<int64_t>(std::llround(v.x)),
                     static_cast<int64_t>(std::llround(v.y)));
  // Remove accidental duplicate consecutive vertices.
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

void GenerateCase(std::mt19937_64& rng, const CheckerConfig& cfg,
                  Paths64& subjects, Paths64& clips) {
  std::uniform_int_distribution<int> vcount_dist(cfg.min_vertices,
                                                 cfg.max_vertices);
  std::uniform_int_distribution<int> subj_dist(1, 3);
  std::uniform_int_distribution<int> clip_dist(1, 2);

  subjects.clear();
  clips.clear();
  int ns = subj_dist(rng);
  for (int i = 0; i < ns; ++i)
    subjects.push_back(
        GenerateSimplePolygon(rng, vcount_dist(rng), cfg.coordinate_extent));
  int nc = clip_dist(rng);
  for (int i = 0; i < nc; ++i)
    clips.push_back(
        GenerateSimplePolygon(rng, vcount_dist(rng), cfg.coordinate_extent));
}

// ---------------- classification paths -----------------------------------

static refcheck::Point ToRef(const Point64& p) {
  return refcheck::Point(p.x, p.y);
}

static Paths ToRefPaths(const Paths64& in) {
  Paths out;
  out.reserve(in.size());
  for (const Path64& p : in) {
    Path rp;
    rp.reserve(p.size());
    for (const Point64& v : p) rp.push_back(ToRef(v));
    out.push_back(std::move(rp));
  }
  return out;
}

// Reference path: independent ray-cast PIP on the two input polygon sets,
// combined with the boolean truth table (even-odd fill rule on both sides).
static bool ReferenceInside(BoolOp op, PipResult ra, PipResult rb) {
  const bool a = (ra == PipResult::kInside);
  const bool b = (rb == PipResult::kInside);
  switch (op) {
    case BoolOp::kIntersection: return a && b;
    case BoolOp::kUnion:        return a || b;
    case BoolOp::kDifference:   return a && !b;
    case BoolOp::kXor:          return a != b;
  }
  return false;
}

// API path: point classification against the boolean solution polygons,
// using Clipper2's own PointInPolygon (even-odd parity across paths).
static bool ApiInside(const Paths64& solution, const Point64& p) {
  int winding = 0;
  for (const Path64& path : solution) {
    PointInPolygonResult r = PointInPolygon(p, path);
    if (r == PointInPolygonResult::IsOn) return true;  // boundary => inside
    if (r == PointInPolygonResult::IsInside) ++winding;
  }
  return (winding & 1) != 0;
}

// Euclidean distance from p to the nearest edge of any input path.
// Computed in double (coordinates ~1e5, squared terms ~1e10, no precision
// issue) and only called at mismatching points, so cost is negligible.
static double DistanceToNearestEdge(const Point64& p, const Paths64& paths) {
  double best = 1e300;
  for (const Path64& path : paths) {
    for (size_t i = 0, j = path.size() - 1; i < path.size(); j = i++) {
      const double ax = static_cast<double>(path[j].x);
      const double ay = static_cast<double>(path[j].y);
      const double bx = static_cast<double>(path[i].x);
      const double by = static_cast<double>(path[i].y);
      const double dx = bx - ax, dy = by - ay;
      const double len2 = dx * dx + dy * dy;
      double t = len2 ? ((static_cast<double>(p.x) - ax) * dx +
                         (static_cast<double>(p.y) - ay) * dy) / len2
                      : 0.0;
      t = std::max(0.0, std::min(1.0, t));
      const double ex = ax + t * dx, ey = ay + t * dy;
      best = std::min(best,
                      std::hypot(static_cast<double>(p.x) - ex,
                                 static_cast<double>(p.y) - ey));
    }
  }
  return best;
}

// ---------------- mismatch region extraction -----------------------------
//
// A mismatch at sample point (row,col) marks one unit cell of an n x n
// lattice. Region boundaries are built from the grid edges that separate
// mismatch cells from non-mismatch cells. Directed edges are then chained
// into closed rings; interior hole rings emerge automatically.

namespace {

struct GridGeometry {
  int64_t min_x, min_y, max_x, max_y;
  int n;
  double dx, dy;
};

GridGeometry MakeGrid(int64_t min_x, int64_t min_y,
                      int64_t max_x, int64_t max_y, int n) {
  GridGeometry g;
  g.min_x = min_x; g.min_y = min_y; g.max_x = max_x; g.max_y = max_y;
  g.n = n;
  g.dx = static_cast<double>(max_x - min_x) / (n - 1);
  g.dy = static_cast<double>(max_y - min_y) / (n - 1);
  return g;
}

// Grid-line vertices (lattice intersections), not cell centres.
Point64 GridVertex(const GridGeometry& g, int col, int row) {
  return Point64(
      static_cast<int64_t>(std::llround(g.min_x + col * g.dx)),
      static_cast<int64_t>(std::llround(g.min_y + row * g.dy)));
}

struct DEdge {  // directed grid edge: from lattice vertex a to b
  int ar, ac, br, bc;
  bool operator==(const DEdge& o) const {
    return ar == o.ar && ac == o.ac && br == o.br && bc == o.bc;
  }
};

struct DEdgeHash {
  size_t operator()(const DEdge& e) const {
    auto h = [](int x) { return static_cast<size_t>(x) * 0x9e3779b1u; };
    return (h(e.ar) ^ (h(e.ac) << 6) ^ (h(e.br) << 12) ^ (h(e.bc) << 18));
  }
};

// Enumerates boundary edges of mismatch cells and chains them into rings
// (each ring is a closed polygon of lattice vertices).
std::vector<Path> ExtractRings(const std::vector<uint8_t>& mask, int n) {
  auto is_mm = [&](int r, int c) {
    if (r < 0 || c < 0 || r >= n || c >= n) return false;
    return mask[r * n + c] != 0;
  };

  // For each mismatch cell, keep only edges facing a non-mismatch
  // neighbour. Edges are directed so the mismatch cell is on the left.
  std::vector<DEdge> edges;
  for (int r = 0; r < n; ++r) {
    for (int c = 0; c < n; ++c) {
      if (!is_mm(r, c)) continue;
      if (!is_mm(r - 1, c))      edges.push_back({r,     c,     r,     c + 1}); // top,    E
      if (!is_mm(r,     c + 1))  edges.push_back({r,     c + 1, r + 1, c + 1}); // right,  S
      if (!is_mm(r + 1, c))      edges.push_back({r + 1, c + 1, r + 1, c    }); // bottom, W
      if (!is_mm(r,     c - 1))  edges.push_back({r + 1, c,     r,     c    }); // left,   N
    }
  }

  // Map each start vertex to its outgoing edge. Boundary edges of a cell
  // set form disjoint simple loops, so each vertex has one outgoing edge.
  std::unordered_map<uint64_t, DEdge> outgoing;
  outgoing.reserve(edges.size() * 2);
  auto vkey = [](int r, int c) -> uint64_t {
    return (static_cast<uint64_t>(static_cast<uint32_t>(r)) << 32) |
           static_cast<uint32_t>(c);
  };
  for (const DEdge& e : edges)
    outgoing.emplace(vkey(e.ar, e.ac), e);

  std::vector<Path> rings;
  std::unordered_map<DEdge, uint8_t, DEdgeHash> used;
  used.reserve(edges.size() * 2);
  for (const DEdge& start : edges) {
    if (used.count(start)) continue;
    Path ring;
    DEdge cur = start;
    do {
      used[cur] = 1;
      // Append the start vertex of each edge; closure re-adds the first.
      ring.emplace_back(cur.ac, cur.ar);  // store as (x=col, y=row)
      auto it = outgoing.find(vkey(cur.br, cur.bc));
      if (it == outgoing.end()) break;     // defensive; loops always close
      cur = it->second;
    } while (!(cur == start));
    rings.push_back(std::move(ring));
  }
  return rings;
}

// Converts a lattice-vertex ring to world coordinates, fusing collinear
// consecutive points so exported polygons stay compact.
Path RingToWorld(const Path& grid_ring, const GridGeometry& g) {
  Path world;
  world.reserve(grid_ring.size());
  for (const Point& gp : grid_ring) {
    Point64 w = GridVertex(g, static_cast<int>(gp.x), static_cast<int>(gp.y));
    world.push_back(refcheck::Point(w.x, w.y));
  }
  // Drop collinear middle points.
  Path compact;
  const size_t m = world.size();
  for (size_t i = 0; i < m; ++i) {
    const Point& a = world[(i + m - 1) % m];
    const Point& b = world[i];
    const Point& c = world[(i + 1) % m];
    if ((b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x) == 0)
      continue;
    compact.push_back(b);
  }
  return compact.empty() ? world : compact;
}

// Counts mismatch cells belonging to the component that contains a ring,
// via flood fill from its first cell (used for the regions' cell_count).
size_t CountComponent(const std::vector<uint8_t>& mask, int n,
                      int seed_row, int seed_col,
                      std::vector<uint8_t>& seen) {
  std::vector<std::pair<int, int>> stack;
  stack.emplace_back(seed_row, seed_col);
  seen[seed_row * n + seed_col] = 1;
  size_t count = 0;
  while (!stack.empty()) {
    auto [r, c] = stack.back();
    stack.pop_back();
    ++count;
    static const int nr[4] = {-1, 1, 0, 0};
    static const int nc[4] = {0, 0, -1, 1};
    for (int d = 0; d < 4; ++d) {
      int r2 = r + nr[d], c2 = c + nc[d];
      if (r2 < 0 || c2 < 0 || r2 >= n || c2 >= n) continue;
      size_t s2 = r2 * n + c2;
      if (mask[s2] && !seen[s2]) {
        seen[s2] = 1;
        stack.emplace_back(r2, c2);
      }
    }
  }
  return count;
}

std::vector<MismatchRegion> BuildRegions(const std::vector<uint8_t>& mask,
                                         int n, const GridGeometry& g) {
  std::vector<MismatchRegion> regions;
  std::vector<Path> rings = ExtractRings(mask, n);
  if (rings.empty()) return regions;

  std::vector<uint8_t> counted(mask.size(), 0);
  int id = 0;
  for (const Path& ring : rings) {
    // First lattice vertex sits at the top-left corner of the mismatch
    // cell just below/right of it; that cell is (row=y, col=x).
    int seed_r = static_cast<int>(ring[0].y);
    int seed_c = static_cast<int>(ring[0].x);
    if (seed_r >= n) seed_r = n - 1;
    if (seed_c >= n) seed_c = n - 1;
    if (!mask[seed_r * n + seed_c]) {  // hole ring: step into the cell
      seed_c = std::min(n - 1, seed_c);
      seed_r = std::max(0, seed_r - 1);
    }
    size_t cells = CountComponent(mask, n, seed_r, seed_c, counted);
    MismatchRegion region;
    region.component_id = id++;
    region.cell_count = cells;
    region.boundary = RingToWorld(ring, g);
    regions.push_back(std::move(region));
  }
  return regions;
}

struct Bounds {
  int64_t min_x, min_y, max_x, max_y;
};

Bounds ComputeBounds(const Paths64& all) {
  auto b = Clipper2Lib::GetBounds(all);
  const int64_t margin_x = (b.right - b.left) / 10 + 10;
  const int64_t margin_y = (b.bottom - b.top) / 10 + 10;
  return {b.left - margin_x, b.top - margin_y,
          b.right + margin_x, b.bottom + margin_y};
}

}  // namespace

// ---------------- per-case driver ----------------------------------------

CaseReport RunCase(int case_index, const Paths64& subjects,
                   const Paths64& clips, const CheckerConfig& cfg) {
  CaseReport rep;
  rep.case_index = case_index;
  rep.subjects = subjects;
  rep.clips = clips;

  Paths64 all = subjects;
  all.insert(all.end(), clips.begin(), clips.end());
  Bounds bd = ComputeBounds(all);
  rep.min_x = bd.min_x; rep.min_y = bd.min_y;
  rep.max_x = bd.max_x; rep.max_y = bd.max_y;
  rep.grid_n = cfg.grid_n;

  const Paths ref_subj = ToRefPaths(subjects);
  const Paths ref_clip = ToRefPaths(clips);
  const GridGeometry grid =
      MakeGrid(rep.min_x, rep.min_y, rep.max_x, rep.max_y, cfg.grid_n);

  const BoolOp ops[4] = {BoolOp::kIntersection, BoolOp::kUnion,
                         BoolOp::kDifference, BoolOp::kXor};

  Paths64 solutions[4];
  for (int k = 0; k < 4; ++k)
    solutions[k] = BooleanOp(ToClipType(ops[k]), FillRule::EvenOdd,
                             subjects, clips);

  std::vector<std::vector<uint8_t>> masks(4);
  std::vector<std::vector<uint8_t>> deep_masks(4);
  rep.ops.resize(4);
  for (int k = 0; k < 4; ++k) {
    rep.ops[k].op = ops[k];
    masks[k].assign(static_cast<size_t>(cfg.grid_n) * cfg.grid_n, 0);
    deep_masks[k].assign(static_cast<size_t>(cfg.grid_n) * cfg.grid_n, 0);
  }

  for (int row = 0; row < cfg.grid_n; ++row) {
    const int64_t iy =
        std::llround(rep.min_y + row * grid.dy);
    for (int col = 0; col < cfg.grid_n; ++col) {
      const int64_t ix =
          std::llround(rep.min_x + col * grid.dx);
      const refcheck::Point rp(ix, iy);
      const Point64 ap(ix, iy);

      const PipResult ra = PointInPaths(rp, ref_subj);
      const PipResult rb = PointInPaths(rp, ref_clip);
      // Boundary points are excluded: the two specs may legitimately
      // differ on boundary convention; interior/exterior is what matters.
      const bool on_boundary = (ra == PipResult::kOn || rb == PipResult::kOn);

      for (int k = 0; k < 4; ++k) {
        if (on_boundary) {
          ++rep.ops[k].boundary;
          continue;
        }
        ++rep.ops[k].compared;
        const bool ref_inside = ReferenceInside(ops[k], ra, rb);
        const bool api_inside = ApiInside(solutions[k], ap);
        if (ref_inside != api_inside) {
          const size_t cell = row * cfg.grid_n + col;
          masks[k][cell] = 1;
          // Classify the mismatch: points within the boundary band of an
          // input edge are expected integer snap-rounding flips; anything
          // farther away is a deep mismatch (a real logic discrepancy).
          const double d = DistanceToNearestEdge(ap, all);
          if (d <= cfg.boundary_band_tolerance) {
            ++rep.ops[k].boundary_band;
          } else {
            deep_masks[k][cell] = 1;
          }
        }
      }
    }
  }

  for (int k = 0; k < 4; ++k) {
    OpReport& op_rep = rep.ops[k];
    op_rep.mismatches =
        std::count(masks[k].begin(), masks[k].end(), uint8_t{1});
    op_rep.regions = BuildRegions(masks[k], cfg.grid_n, grid);
    op_rep.deep_regions = BuildRegions(deep_masks[k], cfg.grid_n, grid);
  }

  for (const OpReport& o : rep.ops) {
    rep.total_compared += o.compared;
    rep.total_mismatches += o.mismatches;
  }
  rep.mismatch_rate = rep.total_compared
                          ? static_cast<double>(rep.total_mismatches) /
                                rep.total_compared
                          : 0.0;
  return rep;
}

// ---------------- mutation injection --------------------------------------

MutationReport RunMutationInjection(const Paths64& subjects,
                                    const Paths64& clips,
                                    BoolOp op,
                                    const CheckerConfig& cfg) {
  MutationReport mr;
  mr.subjects = subjects;
  mr.clips = clips;
  mr.op = op;

  Paths64 solution =
      BooleanOp(ToClipType(op), FillRule::EvenOdd, subjects, clips);
  if (solution.empty() || solution[0].size() < 3) {
    mr.detected = false;  // nothing to mutate for this op
    return mr;
  }
  // Pick the longest path and move one of its vertices outward by a shift
  // larger than two grid spacings, guaranteeing interior points flip.
  size_t best = 0;
  for (size_t i = 1; i < solution.size(); ++i)
    if (solution[i].size() > solution[best].size()) best = i;
  mr.path_index = best;
  Path64& target = solution[best];
  mr.vertex_index = target.size() / 2;

  Paths64 all = subjects;
  all.insert(all.end(), clips.begin(), clips.end());
  Bounds bd = ComputeBounds(all);
  mr.min_x = bd.min_x; mr.min_y = bd.min_y;
  mr.max_x = bd.max_x; mr.max_y = bd.max_y;
  mr.grid_n = cfg.grid_n;
  const GridGeometry grid =
      MakeGrid(mr.min_x, mr.min_y, mr.max_x, mr.max_y, cfg.grid_n);

  const auto b = Clipper2Lib::GetBounds(target);
  int64_t shift_x = std::max<int64_t>(
      static_cast<int64_t>(grid.dx * 3.0), (b.right - b.left) / 4 + 1);
  int64_t shift_y = std::max<int64_t>(
      static_cast<int64_t>(grid.dy * 3.0), (b.bottom - b.top) / 6 + 1);
  mr.old_vertex = target[mr.vertex_index];
  target[mr.vertex_index].x += shift_x;
  target[mr.vertex_index].y += shift_y;
  mr.new_vertex = target[mr.vertex_index];
  mr.mutated_solution = solution;

  // Reference truth is derived from the ORIGINAL inputs (unchanged); the
  // API-side classification now uses the mutated solution.
  const Paths ref_subj = ToRefPaths(subjects);
  const Paths ref_clip = ToRefPaths(clips);
  std::vector<uint8_t> mask(
      static_cast<size_t>(cfg.grid_n) * cfg.grid_n, 0);
  std::vector<uint8_t> deep_mask(mask.size(), 0);

  for (int row = 0; row < cfg.grid_n; ++row) {
    const int64_t iy = std::llround(mr.min_y + row * grid.dy);
    for (int col = 0; col < cfg.grid_n; ++col) {
      const int64_t ix = std::llround(mr.min_x + col * grid.dx);
      const refcheck::Point rp(ix, iy);
      const Point64 ap(ix, iy);
      const PipResult ra = PointInPaths(rp, ref_subj);
      const PipResult rb = PointInPaths(rp, ref_clip);
      if (ra == PipResult::kOn || rb == PipResult::kOn) {
        ++mr.boundary;
        continue;
      }
      ++mr.compared;
      const bool ref_inside = ReferenceInside(mr.op, ra, rb);
      const bool buggy_inside = ApiInside(solution, ap);
      if (ref_inside != buggy_inside) {
        const size_t cell = row * cfg.grid_n + col;
        mask[cell] = 1;
        if (DistanceToNearestEdge(ap, all) > cfg.boundary_band_tolerance)
          deep_mask[cell] = 1;
      }
    }
  }

  mr.mismatches = std::count(mask.begin(), mask.end(), uint8_t{1});
  mr.boundary_band = mr.mismatches -
                     std::count(deep_mask.begin(), deep_mask.end(),
                                uint8_t{1});
  mr.regions = BuildRegions(mask, cfg.grid_n, grid);
  mr.deep_regions = BuildRegions(deep_mask, cfg.grid_n, grid);
  mr.detected = mr.mismatches > 0;
  return mr;
}

std::vector<MutationReport> RunMutationInjectionAll(
    const Paths64& subjects, const Paths64& clips,
    const CheckerConfig& cfg) {
  const BoolOp ops[4] = {BoolOp::kIntersection, BoolOp::kUnion,
                         BoolOp::kDifference, BoolOp::kXor};
  std::vector<MutationReport> out;
  out.reserve(4);
  for (BoolOp op : ops) {
    MutationReport mr = RunMutationInjection(subjects, clips, op, cfg);
    if (!mr.mutated_solution.empty()) out.push_back(std::move(mr));
  }
  return out;
}

// ---------------- export --------------------------------------------------

bool ExportMismatchGeoJson(const std::vector<CaseReport>& reports,
                           const std::string& path) {
  std::ofstream f(path);
  if (!f) return false;
  f << "{\"type\":\"FeatureCollection\",\"features\":[\n";
  bool first = true;
  auto emit_ring = [&](int case_id, const char* op_name,
                       const MismatchRegion& reg, const char* cls) {
    if (!first) f << ",\n";
    first = false;
    f << " {\"type\":\"Feature\",\"properties\":{"
      << "\"case\":" << case_id
      << ",\"operation\":\"" << op_name << "\""
      << ",\"component\":" << reg.component_id
      << ",\"class\":\"" << cls << "\""
      << ",\"mismatch_cells\":" << reg.cell_count
      << "},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[";
    for (size_t i = 0; i < reg.boundary.size(); ++i) {
      const Point& p = reg.boundary[i];
      if (i) f << ",";
      f << "[" << p.x << "," << p.y << "]";
    }
    const Point& p0 = reg.boundary.front();
    f << ",[" << p0.x << "," << p0.y << "]]]}}";
  };
  for (const CaseReport& cr : reports)
    for (const OpReport& op : cr.ops) {
      // A region is "deep" if it also appears in the deep-region set
      // (matched by its first vertex); otherwise it is boundary-band.
      auto is_deep = [&op](const MismatchRegion& reg) {
        for (const MismatchRegion& d : op.deep_regions)
          if (!d.boundary.empty() && !reg.boundary.empty() &&
              d.boundary.front().x == reg.boundary.front().x &&
              d.boundary.front().y == reg.boundary.front().y)
            return true;
        return false;
      };
      for (const MismatchRegion& reg : op.regions)
        emit_ring(cr.case_index, ToString(op.op), reg,
                  is_deep(reg) ? "deep" : "boundary_band");
    }
  f << "\n]}\n";
  return static_cast<bool>(f);
}

bool ExportCaseSvg(const CaseReport& rep, const std::string& path) {
  std::ofstream f(path);
  if (!f) return false;
  const int64_t w = rep.max_x - rep.min_x;
  const int64_t h = rep.max_y - rep.min_y;
  f << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
    << "viewBox=\"" << rep.min_x << " " << rep.min_y << " " << w << " " << h
    << "\" width=\"900\" height=\""
    << std::llround(900.0 * h / w) << "\">\n";
  f << "<rect x=\"" << rep.min_x << "\" y=\"" << rep.min_y
    << "\" width=\"" << w << "\" height=\"" << h
    << "\" fill=\"white\"/>\n";

  auto emit_path = [&](const Path64& p, const char* stroke,
                       const char* fill, double sw) {
    f << "<polygon points=\"";
    for (const Point64& v : p) f << v.x << "," << v.y << " ";
    f << "\" fill=\"" << fill << "\" fill-opacity=\"0.08\" stroke=\""
      << stroke << "\" stroke-width=\"" << sw << "\"/>\n";
  };
  for (const Path64& p : rep.subjects) emit_path(p, "#1f77b4", "#1f77b4", 60);
  for (const Path64& p : rep.clips) emit_path(p, "#2ca02c", "#2ca02c", 60);

  for (const OpReport& op : rep.ops)
    for (const MismatchRegion& reg : op.regions) {
      f << "<polygon points=\"";
      for (const Point& v : reg.boundary) f << v.x << "," << v.y << " ";
      f << "\" fill=\"red\" fill-opacity=\"0.45\" stroke=\"darkred\" "
           "stroke-width=\"40\"/>\n";
    }
  f << "</svg>\n";
  return static_cast<bool>(f);
}

bool ExportMutationSvg(const MutationReport& mr, const std::string& path) {
  std::ofstream f(path);
  if (!f) return false;
  const int64_t w = mr.max_x - mr.min_x;
  const int64_t h = mr.max_y - mr.min_y;
  f << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
    << "viewBox=\"" << mr.min_x << " " << mr.min_y << " " << w << " " << h
    << "\" width=\"900\" height=\""
    << std::llround(900.0 * h / w) << "\">\n";
  f << "<rect x=\"" << mr.min_x << "\" y=\"" << mr.min_y
    << "\" width=\"" << w << "\" height=\"" << h << "\" fill=\"white\"/>\n";
  auto emit = [&](const Path64& p, const char* color, double opacity) {
    f << "<polygon points=\"";
    for (const Point64& v : p) f << v.x << "," << v.y << " ";
    f << "\" fill=\"" << color << "\" fill-opacity=\"" << opacity
      << "\" stroke=\"" << color << "\" stroke-width=\"60\"/>\n";
  };
  for (const Path64& p : mr.subjects) emit(p, "#1f77b4", 0.06);
  for (const Path64& p : mr.clips) emit(p, "#2ca02c", 0.06);
  emit(mr.mutated_solution[mr.path_index], "orange", 0.10);
  for (const MismatchRegion& reg : mr.regions) {
    f << "<polygon points=\"";
    for (const Point& v : reg.boundary) f << v.x << "," << v.y << " ";
    f << "\" fill=\"red\" fill-opacity=\"0.45\" stroke=\"darkred\" "
         "stroke-width=\"40\"/>\n";
  }
  // Mark the moved vertex.
  f << "<circle cx=\"" << mr.new_vertex.x << "\" cy=\"" << mr.new_vertex.y
    << "\" r=\"" << (w / 120) << "\" fill=\"black\"/>\n";
  f << "</svg>\n";
  return static_cast<bool>(f);
}

bool ExportSummary(const std::vector<CaseReport>& reports,
                   double elapsed_seconds, const MutationReport* mutation,
                   const std::string& path) {
  std::ofstream f(path);
  if (!f) return false;
  size_t total_compared = 0, total_mm = 0, total_boundary = 0;
  size_t total_band = 0, total_deep = 0, total_regions = 0;
  for (const CaseReport& cr : reports) {
    total_compared += cr.total_compared;
    for (const OpReport& op : cr.ops) {
      total_boundary += op.boundary;
      total_mm += op.mismatches;
      total_band += op.boundary_band;
      total_deep += op.mismatches - op.boundary_band;
      total_regions += op.regions.size();
    }
  }
  const double rate = total_compared
                          ? 100.0 * total_mm / total_compared : 0.0;
  f << std::fixed << std::setprecision(6);
  f << "{\n";
  f << "  \"cases\": " << reports.size() << ",\n";
  f << "  \"points_compared\": " << total_compared << ",\n";
  f << "  \"points_on_boundary_skipped\": " << total_boundary << ",\n";
  f << "  \"mismatches\": " << total_mm << ",\n";
  f << "  \"boundary_band_mismatches\": " << total_band << ",\n";
  f << "  \"deep_mismatches\": " << total_deep << ",\n";
  f << "  \"mismatch_rate_percent\": " << rate << ",\n";
  f << "  \"mismatch_regions\": " << total_regions << ",\n";
  f << "  \"elapsed_seconds\": " << elapsed_seconds << ",\n";
  if (mutation) {
    f << "  \"mutation_injection\": {\n";
    f << "    \"operation\": \"" << ToString(mutation->op) << "\",\n";
    f << "    \"path_index\": " << mutation->path_index << ",\n";
    f << "    \"vertex_index\": " << mutation->vertex_index << ",\n";
    f << "    \"old_vertex\": [" << mutation->old_vertex.x << ", "
      << mutation->old_vertex.y << "],\n";
    f << "    \"new_vertex\": [" << mutation->new_vertex.x << ", "
      << mutation->new_vertex.y << "],\n";
    f << "    \"points_compared\": " << mutation->compared << ",\n";
    f << "    \"mismatches\": " << mutation->mismatches << ",\n";
    f << "    \"boundary_band_mismatches\": " << mutation->boundary_band
      << ",\n";
    f << "    \"deep_mismatches\": "
      << (mutation->mismatches - mutation->boundary_band) << ",\n";
    f << "    \"regions\": " << mutation->regions.size() << ",\n";
    f << "    \"deep_regions\": " << mutation->deep_regions.size() << ",\n";
    f << "    \"detected\": " << (mutation->detected ? "true" : "false")
      << "\n";
    f << "  },\n";
  } else {
    f << "  \"mutation_injection\": null,\n";
  }
  f << "  \"per_case\": [\n";
  for (size_t i = 0; i < reports.size(); ++i) {
    const CaseReport& cr = reports[i];
    f << "    {\"case\": " << cr.case_index
      << ", \"compared\": " << cr.total_compared
      << ", \"mismatches\": " << cr.total_mismatches
      << ", \"rate_percent\": " << (cr.mismatch_rate * 100.0)
      << ", \"ops\": {";
    for (size_t k = 0; k < cr.ops.size(); ++k) {
      const OpReport& op = cr.ops[k];
      if (k) f << ", ";
      f << "\"" << ToString(op.op) << "\": {\"compared\": " << op.compared
        << ", \"boundary_skipped\": " << op.boundary
        << ", \"mismatches\": " << op.mismatches
        << ", \"boundary_band\": " << op.boundary_band
        << ", \"deep\": " << (op.mismatches - op.boundary_band)
        << ", \"regions\": " << op.regions.size() << "}";
    }
    f << "}}" << (i + 1 == reports.size() ? "\n" : ",\n");
  }
  f << "  ]\n}\n";
  return static_cast<bool>(f);
}

}  // namespace refcheck
