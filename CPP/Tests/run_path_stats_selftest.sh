#!/usr/bin/env bash
# Builds and runs the standalone PathStatsSelfTest.cpp against the Clipper2
# sources in this repository. No network access, no third-party libraries.
#
# Usage: CPP/Tests/run_path_stats_selftest.sh [--keep]
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cpp_root="$(cd "${here}/.." && pwd)"
out_dir="$(mktemp -d "${TMPDIR:-/tmp}/clipper2-pathstats.XXXXXX")"
trap '[[ "${1:-}" == "--keep" ]] || rm -rf "${out_dir}"' EXIT

cxx="${CXX:-g++}"
# C++17, optimised like a release build; USUALLY only the library sources
# (no gtest) are needed since the self test is a plain main().
"${cxx}" -std=c++17 -O2 -Wall -Wextra -Wpedantic \
  -I"${cpp_root}/Clipper2Lib/include" \
  "${here}/PathStatsSelfTest.cpp" \
  "${cpp_root}/Clipper2Lib/src/clipper.engine.cpp" \
  "${cpp_root}/Clipper2Lib/src/clipper.offset.cpp" \
  "${cpp_root}/Clipper2Lib/src/clipper.rectclip.cpp" \
  "${cpp_root}/Clipper2Lib/src/clipper.triangulation.cpp" \
  -o "${out_dir}/path_stats_selftest"

"${out_dir}/path_stats_selftest"
