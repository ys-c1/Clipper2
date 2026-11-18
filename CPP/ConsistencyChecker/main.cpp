// Consistency checker entry point.
//
// Two independent paths decide inside/outside for the same boolean result:
//   * API path:        Clipper2 BooleanOp() + library PointInPolygon()
//   * reference path:  a self-contained ray-casting PIP (pip_reference.h)
//                      combined with the boolean truth table
//
// The mutation-injection test deliberately moves one vertex of the
// boolean SOLUTION (simulating a buggy API result) and verifies that the
// reference path flags the affected area.
//
// Usage:
//   consistency_checker [--cases N] [--grid N] [--seed S] [--out DIR]
//                       [--no-inject]
//
// Exit code is non-zero if the mismatch threshold is exceeded, or if an
// injected mutation is not detected.
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "checker.h"

using namespace refcheck;
using Clipper2Lib::Paths64;

namespace {

struct Args {
  CheckerConfig cfg;
  bool inject = true;
};

void PrintUsage() {
  std::cout <<
      "usage: consistency_checker [--cases N] [--grid N] [--seed S]\n"
      "                           [--out DIR] [--no-inject]\n"
      "  --cases N      number of random polygon pairs (default 50)\n"
      "  --grid N       N x N sample points per case (default 100 => 1e4)\n"
      "  --seed S       deterministic RNG seed\n"
      "  --out DIR      output directory for reports and regions\n"
      "  --band-tolerance T  mismatch within T units of an edge is a\n"
      "                     boundary-band snap-rounding flip (default 1.0)\n"
      "  --no-inject    skip the mutation-injection test\n";
}

Args ParseArgs(int argc, char** argv) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto next = [&](const char* name) -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "missing value for " << name << "\n";
        std::exit(2);
      }
      return argv[++i];
    };
    if (a == "--cases") args.cfg.cases = std::stoi(next("--cases"));
    else if (a == "--grid") args.cfg.grid_n = std::stoi(next("--grid"));
    else if (a == "--seed")
      args.cfg.seed = std::stoull(next("--seed"));
    else if (a == "--out") args.cfg.out_dir = next("--out");
    else if (a == "--band-tolerance")
      args.cfg.boundary_band_tolerance = std::stod(next("--band-tolerance"));
    else if (a == "--no-inject") args.inject = false;
    else if (a == "-h" || a == "--help") { PrintUsage(); std::exit(0); }
    else { std::cerr << "unknown argument: " << a << "\n"; std::exit(2); }
  }
  return args;
}

}  // namespace

int main(int argc, char** argv) {
  Args args = ParseArgs(argc, argv);
  CheckerConfig& cfg = args.cfg;
  std::filesystem::create_directories(cfg.out_dir);

  std::mt19937_64 rng(cfg.seed);
  const auto t0 = std::chrono::steady_clock::now();

  const long long points_per_case =
      static_cast<long long>(cfg.grid_n) * cfg.grid_n;
  std::cout << "Consistency check: " << cfg.cases
            << " random polygon pair(s), " << cfg.grid_n << "x" << cfg.grid_n
            << " = " << points_per_case
            << " sample points per case, 4 boolean ops each\n";

  std::vector<CaseReport> reports;
  reports.reserve(cfg.cases);
  for (int i = 0; i < cfg.cases; ++i) {
    Paths64 subjects, clips;
    GenerateCase(rng, cfg, subjects, clips);
    reports.push_back(RunCase(i, subjects, clips, cfg));
  }

  // Mutation-injection test: a fresh case whose boolean solution has one
  // vertex moved for every operation. The reference path stays anchored to
  // the original inputs and must detect each discrepancy.
  std::vector<MutationReport> mutations;
  if (args.inject) {
    std::cout << "\nMutation-injection test (one moved vertex per op):\n";
    std::mt19937_64 mrng(cfg.seed ^ 0x9e3779b97f4a7c15ULL);
    Paths64 msub, mclip;
    GenerateCase(mrng, cfg, msub, mclip);
    mutations = RunMutationInjectionAll(msub, mclip, cfg);
    for (const MutationReport& m : mutations) {
      std::cout << "  " << ToString(m.op) << ": moved path " << m.path_index
                << " vertex " << m.vertex_index << " (" << m.old_vertex.x
                << "," << m.old_vertex.y << ")->(" << m.new_vertex.x << ","
                << m.new_vertex.y << "); " << m.mismatches
                << " mismatches (" << (m.mismatches - m.boundary_band)
                << " deep) => " << (m.detected ? "DETECTED" : "NOT DETECTED")
                << "\n";
    }
  }

  const auto t1 = std::chrono::steady_clock::now();
  const double elapsed =
      std::chrono::duration<double>(t1 - t0).count();

  size_t total_compared = 0, total_mm = 0, total_regions = 0;
  size_t total_band = 0, total_deep = 0;
  double worst_rate = 0.0;
  for (const CaseReport& cr : reports) {
    total_compared += cr.total_compared;
    total_mm += cr.total_mismatches;
    worst_rate = std::max(worst_rate, cr.mismatch_rate);
    for (const OpReport& op : cr.ops) {
      total_regions += op.regions.size();
      total_band += op.boundary_band;
      total_deep += op.mismatches - op.boundary_band;
    }
  }
  const double rate_percent =
      total_compared ? 100.0 * total_mm / total_compared : 0.0;

  std::cout << "\nResults (50-case sweep):\n";
  std::cout << "  points compared : " << total_compared << "\n";
  std::cout << "  mismatches      : " << total_mm << " ("
            << rate_percent << " %)\n";
  std::cout << "    boundary-band : " << total_band
            << " (snap-rounding flips within 1 integer unit of an edge)\n";
  std::cout << "    deep          : " << total_deep
            << " (logic discrepancies; must be 0)\n";
  std::cout << "  worst per-case  : " << (worst_rate * 100.0) << " %\n";
  std::cout << "  mismatch regions: " << total_regions << "\n";
  std::cout << "  elapsed         : " << elapsed << " s\n";

  const std::string geo = cfg.out_dir + "/mismatch_regions.geojson";
  const std::string sum = cfg.out_dir + "/summary.json";
  if (!ExportMismatchGeoJson(reports, geo))
    std::cerr << "warning: failed to write " << geo << "\n";
  else
    std::cout << "  regions exported: " << geo << "\n";
  if (!mutations.empty()) {
    // Concatenate every op's mutation regions into one GeoJSON file.
    std::ofstream mg(cfg.out_dir + "/mutation_regions.geojson");
    mg << "{\"type\":\"FeatureCollection\",\"features\":[\n";
    bool first = true;
    for (const MutationReport& m : mutations)
      for (const MismatchRegion& reg : m.deep_regions) {
        if (!first) mg << ",\n";
        first = false;
        mg << " {\"type\":\"Feature\",\"properties\":{\"operation\":\""
           << ToString(m.op) << "\",\"component\":" << reg.component_id
           << ",\"mismatch_cells\":" << reg.cell_count
           << "},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[";
        for (size_t i = 0; i < reg.boundary.size(); ++i) {
          if (i) mg << ",";
          mg << "[" << reg.boundary[i].x << "," << reg.boundary[i].y << "]";
        }
        mg << ",[" << reg.boundary.front().x << ","
           << reg.boundary.front().y << "]]]}}";
      }
    mg << "\n]}\n";
  }
  const MutationReport* first_mutation =
      mutations.empty() ? nullptr : &mutations.front();
  if (!ExportSummary(reports, elapsed, first_mutation, sum))
    std::cerr << "warning: failed to write " << sum << "\n";
  else
    std::cout << "  summary written : " << sum << "\n";

  // One SVG per non-trivial case (cap file count).
  int svg_count = 0;
  for (const CaseReport& cr : reports) {
    if (cr.total_mismatches == 0) continue;
    ExportCaseSvg(cr, cfg.out_dir + "/case_" +
                        std::to_string(cr.case_index) + ".svg");
    if (++svg_count >= 20) break;
  }
  for (const MutationReport& m : mutations)
    if (m.detected)
      ExportMutationSvg(m, cfg.out_dir + "/mutation_" +
                                ToString(m.op) + ".svg");

  // Acceptance criteria:
  //  * overall mismatch rate < 0.01 %
  //  * zero DEEP mismatches (any mismatch beyond the snap-rounding band is
  //    a genuine logic error); boundary-band flips are tolerated
  constexpr double kMaxRatePercent = 0.01;
  int rc = 0;
  if (rate_percent >= kMaxRatePercent) {
    std::cerr << "FAIL: mismatch rate " << rate_percent
              << "% reaches/exceeds threshold " << kMaxRatePercent << "%\n";
    rc = 1;
  }
  if (total_deep != 0) {
    std::cerr << "FAIL: " << total_deep
              << " deep mismatches outside the snap-rounding boundary band\n";
    rc = 1;
  }
  if (args.inject) {
    for (const MutationReport& m : mutations)
      if (!m.detected) {
        std::cerr << "FAIL: injected vertex mutation was not detected for "
                  << ToString(m.op) << "\n";
        rc = 1;
      }
  }
  if (rc == 0) std::cout << "\nPASS\n";
  return rc;
}
