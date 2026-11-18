// Consistency checker: compares the inside/outside verdict of Clipper2's
// boolean operation against an independent per-point point-in-polygon
// reference (see pip_reference.h), on a regular grid of sample points.
#ifndef CONSISTENCY_CHECKER_CHECKER_H
#define CONSISTENCY_CHECKER_CHECKER_H

#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "clipper2/clipper.h"
#include "pip_reference.h"

namespace refcheck {

// One boolean operation under test. Mirrors Clipper2Lib::ClipType.
enum class BoolOp { kIntersection, kUnion, kDifference, kXor };

const char* ToString(BoolOp op);

// A set of mismatching grid cells, exported as a closed boundary polygon
// in the same coordinate system as the input polygons.
struct MismatchRegion {
  int component_id = 0;
  Path boundary;            // boundary polygon (world coordinates)
  size_t cell_count = 0;    // number of mismatching grid cells inside
};

struct OpReport {
  BoolOp op = BoolOp::kIntersection;
  size_t compared = 0;      // grid points actually compared (boundary skipped)
  size_t boundary = 0;      // grid points lying exactly on an edge (excluded)
  size_t mismatches = 0;    // points where the two paths disagree
  // Subset of 'mismatches' lying within boundary_band_tolerance of an
  // input edge: expected snap-rounding flips of the integer pipeline,
  // not logic errors. mismatches - boundary_band = deep mismatches,
  // which must be zero for a correct boolean implementation.
  size_t boundary_band = 0;
  std::vector<MismatchRegion> regions;       // all mismatch regions
  std::vector<MismatchRegion> deep_regions;  // deep-mismatch regions only
};

struct CaseReport {
  int case_index = 0;
  Clipper2Lib::Paths64 subjects;  // kept for export / mutation display
  Clipper2Lib::Paths64 clips;
  int64_t min_x = 0, min_y = 0, max_x = 0, max_y = 0;  // sampling bounds
  int grid_n = 0;                 // grid_n x grid_n sample points
  std::vector<OpReport> ops;
  size_t total_compared = 0;
  size_t total_mismatches = 0;
  double mismatch_rate = 0.0;     // total_mismatches / total_compared
};

// Result of the mutation-injection test: one vertex of the boolean
// solution is deliberately moved (simulating a buggy boolean result), and
// the independent reference classification must flag the affected area.
struct MutationReport {
  BoolOp op = BoolOp::kUnion;
  size_t path_index = 0;
  size_t vertex_index = 0;
  Clipper2Lib::Point64 old_vertex{0, 0};
  Clipper2Lib::Point64 new_vertex{0, 0};
  size_t compared = 0;
  size_t boundary = 0;
  size_t mismatches = 0;
  size_t boundary_band = 0;  // subset of mismatches near an input edge
  std::vector<MismatchRegion> regions;
  std::vector<MismatchRegion> deep_regions;
  Clipper2Lib::Paths64 subjects;
  Clipper2Lib::Paths64 clips;
  Clipper2Lib::Paths64 mutated_solution;
  int64_t min_x = 0, min_y = 0, max_x = 0, max_y = 0;
  int grid_n = 0;
  bool detected = false;
};

struct CheckerConfig {
  int cases = 50;
  int grid_n = 100;               // 100 x 100 = 1e4 sample points per case
  int min_vertices = 6;
  int max_vertices = 30;
  uint64_t seed = 20261119ULL;
  int coordinate_extent = 100000; // polygons generated in [0, extent]
  // A mismatch closer than this (integer units) to any input edge is
  // classified as a boundary-band snap-rounding flip rather than a deep
  // logic error. Integer intersection rounding is bounded by 0.5 unit.
  double boundary_band_tolerance = 1.0;
  std::string out_dir = "artifacts/consistency";
};

// Generates a simple (non self-intersecting) polygon: vertices are sorted
// by angle around the centroid of radially perturbed points, so every
// generated path is a valid star-shaped simple polygon.
Clipper2Lib::Path64 GenerateSimplePolygon(std::mt19937_64& rng, int vertex_count,
                             int64_t extent);

// Generates one test case: 1-3 subject paths and 1-2 clip paths, each a
// simple polygon (some nested, to exercise multi-path even-odd parity).
void GenerateCase(std::mt19937_64& rng, const CheckerConfig& cfg,
                  Clipper2Lib::Paths64& subjects,
                  Clipper2Lib::Paths64& clips);

// Runs both classification paths for every grid point and every boolean
// op, and builds mismatch regions.
CaseReport RunCase(int case_index,
                   const Clipper2Lib::Paths64& subjects,
                   const Clipper2Lib::Paths64& clips,
                   const CheckerConfig& cfg);

// Moves one vertex of the boolean solution for the given operation (the
// API output), then re-classifies the grid: reference truth stays anchored
// to the original inputs, so every interior point whose verdict changed is
// a mismatch the checker must detect. Returns detection details.
MutationReport RunMutationInjection(
    const Clipper2Lib::Paths64& subjects,
    const Clipper2Lib::Paths64& clips,
    BoolOp op,
    const CheckerConfig& cfg);

// Runs the mutation test for every boolean op and returns one report per
// op (ops with empty solutions are omitted).
std::vector<MutationReport> RunMutationInjectionAll(
    const Clipper2Lib::Paths64& subjects,
    const Clipper2Lib::Paths64& clips,
    const CheckerConfig& cfg);

// ---- export -------------------------------------------------------------

// Writes every mismatch region of every op/case as a GeoJSON FeatureCollection.
bool ExportMismatchGeoJson(const std::vector<CaseReport>& reports,
                           const std::string& path);

// Writes one SVG per case: inputs, boolean result and mismatching cells.
bool ExportCaseSvg(const CaseReport& report, const std::string& path);

// Writes an SVG visualising the mutated solution and detected region.
bool ExportMutationSvg(const MutationReport& mutation,
                       const std::string& path);

// Writes a machine-readable summary (one JSON object).
bool ExportSummary(const std::vector<CaseReport>& reports,
                   double elapsed_seconds, const MutationReport* mutation,
                   const std::string& path);

}  // namespace refcheck

#endif  // CONSISTENCY_CHECKER_CHECKER_H
