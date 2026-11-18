#include <gtest/gtest.h>
#include "clipper2/clipper.batch.h"

#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>

using namespace Clipper2Lib;

namespace {

// Deterministic small-polygon generator (no wall-clock or thread randomness).
class BatchTestGen {
public:
  explicit BatchTestGen(uint64_t seed) : rng_(seed) {}

  Paths64 MakeRandomSubject(int vertex_count, int64_t extent)
  {
    Path64 path;
    path.reserve(static_cast<size_t>(vertex_count));
    std::uniform_int_distribution<int64_t> dist(-extent, extent);
    for (int i = 0; i < vertex_count; ++i)
      path.push_back(Point64(dist(rng_), dist(rng_)));
    return Paths64{ path };
  }

  // A convex-ish blob: points on a perturbed circle, always well-formed.
  Paths64 MakeBlob(int vertex_count, int64_t cx, int64_t cy, int64_t radius)
  {
    Path64 path;
    path.reserve(static_cast<size_t>(vertex_count));
    std::uniform_int_distribution<int64_t> jitter(-radius / 4, radius / 4);
    for (int i = 0; i < vertex_count; ++i) {
      const double a = 2.0 * 3.14159265358979323846 * i / vertex_count;
      const int64_t r = radius + jitter(rng_);
      path.push_back(Point64(
        cx + static_cast<int64_t>(r * cos(a)),
        cy + static_cast<int64_t>(r * sin(a))));
    }
    return Paths64{ path };
  }

  int64_t Coord(int64_t lo, int64_t hi)
  {
    std::uniform_int_distribution<int64_t> dist(lo, hi);
    return dist(rng_);
  }

private:
  std::mt19937_64 rng_;
};

uint64_t HashPaths64(const Paths64& paths, uint64_t h)
{
  const auto mix = [&h](uint64_t v) {
    h ^= v;
    h *= 1099511628211ULL; // FNV-1a prime
  };
  for (const Path64& path : paths) {
    mix(path.size());
    for (const Point64& pt : path) {
      mix(static_cast<uint64_t>(pt.x));
      mix(static_cast<uint64_t>(pt.y));
    }
  }
  return h;
}

uint64_t HashBatch(const std::vector<Paths64>& results)
{
  uint64_t h = 1469598103934665603ULL; // FNV-1a offset basis
  for (const Paths64& paths : results) h = HashPaths64(paths, h);
  return h;
}

std::string Hex(uint64_t v)
{
  char buf[17];
  snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(v));
  return std::string(buf);
}

std::vector<BooleanJob64> MakeJobs(int count, uint64_t seed)
{
  BatchTestGen gen(seed);
  std::vector<BooleanJob64> jobs;
  jobs.reserve(static_cast<size_t>(count));
  const ClipType cts[4] = {
    ClipType::Intersection, ClipType::Union, ClipType::Difference, ClipType::Xor };
  for (int i = 0; i < count; ++i) {
    BooleanJob64 job;
    job.cliptype = cts[i % 4];
    job.fillrule = (i % 2) ? FillRule::EvenOdd : FillRule::NonZero;
    const int64_t cx = gen.Coord(-5000, 5000);
    const int64_t cy = gen.Coord(-5000, 5000);
    job.subject = gen.MakeBlob(8, cx, cy, 100);
    job.clip = gen.MakeBlob(8, cx + gen.Coord(-80, 80), cy + gen.Coord(-80, 80), 100);
    jobs.push_back(std::move(job));
  }
  return jobs;
}

} // namespace

// Acceptance 1: batch results match per-group sequential calls, point by point.
TEST(BatchBoolean, MatchesSequentialPerGroup)
{
  const std::vector<BooleanJob64> jobs = MakeJobs(5, 12345);
  const std::vector<Paths64> batched = BooleanOpBatch(jobs, 4);
  ASSERT_EQ(batched.size(), jobs.size());
  for (size_t i = 0; i < jobs.size(); ++i) {
    const Paths64 expected = BooleanOp(
      jobs[i].cliptype, jobs[i].fillrule, jobs[i].subject, jobs[i].clip);
    EXPECT_EQ(batched[i], expected) << "group " << i << " differs";
  }
}

// Acceptance 2: parallel and serial results are identical (hash compared).
TEST(BatchBoolean, ParallelHashMatchesSerial)
{
  const std::vector<BooleanJob64> jobs = MakeJobs(500, 777);
  const unsigned thread_counts[] = { 1, 2, 4, 8 };
  uint64_t reference = 0;
  std::vector<Paths64> reference_results;
  for (unsigned threads : thread_counts) {
    const std::vector<Paths64> results = BooleanOpBatch(jobs, threads);
    const uint64_t h = HashBatch(results);
    std::cout << "[BatchHash] threads=" << threads
              << " groups=" << jobs.size()
              << " hash=0x" << Hex(h) << std::endl;
    if (threads == 1) {
      reference = h;
      reference_results = results;
    } else {
      EXPECT_EQ(h, reference) << "hash mismatch at threads=" << threads;
      EXPECT_EQ(results, reference_results)
        << "point-by-point mismatch at threads=" << threads;
    }
  }
  const std::vector<Paths64> auto_results = BooleanOpBatch(jobs, 0);
  EXPECT_EQ(HashBatch(auto_results), reference) << "auto-thread hash mismatch";
}

// Acceptance 3: 10,000 small polygon groups complete in under 500 ms.
TEST(BatchBoolean, Performance10kSmallPolygons)
{
  const std::vector<BooleanJob64> jobs = MakeJobs(10000, 2026);
  const auto start = std::chrono::steady_clock::now();
  const std::vector<Paths64> results = BooleanOpBatch(jobs, 0);
  const auto stop = std::chrono::steady_clock::now();
  const auto ms =
    std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
  std::cout << "[BatchPerf] groups=" << jobs.size()
            << " elapsed_ms=" << ms
            << " hash=0x" << Hex(HashBatch(results)) << std::endl;
  EXPECT_EQ(results.size(), jobs.size());
  EXPECT_LT(ms, 500);
}

// Acceptance 4: empty groups (and an empty job list) are handled correctly.
TEST(BatchBoolean, EmptyGroups)
{
  BatchTestGen gen(99);
  const Paths64 subject = gen.MakeBlob(8, 0, 0, 100);
  const Paths64 clip = gen.MakeBlob(8, 50, 0, 100);
  const Paths64 empty;

  std::vector<BooleanJob64> jobs;
  const ClipType cts[4] = {
    ClipType::Intersection, ClipType::Union, ClipType::Difference, ClipType::Xor };
  for (ClipType ct : cts) {
    jobs.push_back(BooleanJob64{ ct, FillRule::NonZero, empty, empty });
    jobs.push_back(BooleanJob64{ ct, FillRule::NonZero, subject, empty });
    jobs.push_back(BooleanJob64{ ct, FillRule::NonZero, empty, clip });
  }

  const std::vector<Paths64> batched = BooleanOpBatch(jobs, 4);
  ASSERT_EQ(batched.size(), jobs.size());
  for (size_t i = 0; i < jobs.size(); ++i) {
    const Paths64 expected = BooleanOp(
      jobs[i].cliptype, jobs[i].fillrule, jobs[i].subject, jobs[i].clip);
    EXPECT_EQ(batched[i], expected) << "empty-group case " << i << " differs";
  }

  // An empty job list yields an empty result, whatever the thread count.
  EXPECT_TRUE(BooleanOpBatch(std::vector<BooleanJob64>{}, 0).empty());
  EXPECT_TRUE(BooleanOpBatch(std::vector<BooleanJob64>{}, 8).empty());
}
