/*******************************************************************************
* Author    :  Angus Johnson                                                   *
* Date      :  19 November 2026                                                *
* Website   :  https://www.angusj.com                                          *
* Copyright :  Angus Johnson 2010-2026                                         *
* Purpose   :  Batched boolean operations over independent job groups          *
* License   :  https://www.boost.org/LICENSE_1_0.txt                           *
*******************************************************************************/

#ifndef CLIPPER2_CLIPPER_BATCH_H
#define CLIPPER2_CLIPPER_BATCH_H

#include "clipper2/clipper.h"

#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

namespace Clipper2Lib {

  // One independent boolean job: a single (subject, clip) group.
  struct BooleanJob64 {
    ClipType cliptype = ClipType::Intersection;
    FillRule fillrule = FillRule::NonZero;
    Paths64 subject;
    Paths64 clip;
  };

  struct BooleanJobD {
    ClipType cliptype = ClipType::Intersection;
    FillRule fillrule = FillRule::NonZero;
    PathsD subject;
    PathsD clip;
    int decimal_prec = 2;
  };

  namespace detail {

    // Runs `fn` over every job, writing each result to its own preallocated
    // slot, so results stay indexed by input group regardless of scheduling.
    // Work is split into fixed contiguous chunks (static partitioning), which
    // keeps the output fully deterministic for any thread count.
    // max_threads: 0 = hardware concurrency, 1 = serial, N = at most N threads.
    template <typename Job, typename Result, typename Fn>
    inline std::vector<Result> RunBatch(
      const std::vector<Job>& jobs, unsigned max_threads, Fn&& fn)
    {
      const size_t count = jobs.size();
      std::vector<Result> results(count);
      if (count == 0) return results;

      size_t thread_count = max_threads;
      if (thread_count == 0) {
        thread_count = std::thread::hardware_concurrency();
        if (thread_count == 0) thread_count = 1;
      }
      thread_count = std::min(thread_count, count);

      if (thread_count <= 1) {
        for (size_t i = 0; i < count; ++i) results[i] = fn(jobs[i]);
        return results;
      }

      const size_t chunk = (count + thread_count - 1) / thread_count;
      std::vector<std::thread> workers;
      workers.reserve(thread_count - 1);
      for (size_t t = 1; t < thread_count; ++t) {
        const size_t begin = t * chunk;
        const size_t end = std::min(count, begin + chunk);
        if (begin >= end) break;
        workers.emplace_back([&jobs, &results, &fn, begin, end]() {
          for (size_t i = begin; i < end; ++i) results[i] = fn(jobs[i]);
        });
      }
      const size_t main_end = std::min(count, chunk);
      for (size_t i = 0; i < main_end; ++i) results[i] = fn(jobs[i]);
      for (std::thread& w : workers) w.join();
      return results;
    }

  } // namespace detail

  // Batch entry point: executes every job group independently and returns one
  // result per group, in input order (results[i] corresponds to jobs[i]).
  // Empty groups (empty subject and/or clip) are valid and follow the same
  // rules as the single-group BooleanOp. Never throws for an empty job list.
  inline std::vector<Paths64> BooleanOpBatch(
    const std::vector<BooleanJob64>& jobs, unsigned max_threads = 0)
  {
    return detail::RunBatch<BooleanJob64, Paths64>(jobs, max_threads,
      [](const BooleanJob64& j) {
        return BooleanOp(j.cliptype, j.fillrule, j.subject, j.clip);
      });
  }

  inline std::vector<PathsD> BooleanOpBatch(
    const std::vector<BooleanJobD>& jobs, unsigned max_threads = 0)
  {
    return detail::RunBatch<BooleanJobD, PathsD>(jobs, max_threads,
      [](const BooleanJobD& j) {
        return BooleanOp(j.cliptype, j.fillrule, j.subject, j.clip, j.decimal_prec);
      });
  }

} // namespace Clipper2Lib

#endif // CLIPPER2_CLIPPER_BATCH_H
