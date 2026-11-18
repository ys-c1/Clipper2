# Boolean Consistency Checker

Validates Clipper2's boolean operations by computing the same
inside/outside verdict through **two independent paths** and comparing
them on a regular grid of sample points.

## The two independent paths

| Path | How inside/outside is decided |
| --- | --- |
| **API path** | `Clipper2Lib::BooleanOp()` produces the solution polygons; Clipper2's own `PointInPolygon()` classifies each sample point against them (even-odd parity across solution paths). |
| **Reference path** | A self-contained textbook ray-casting point-in-polygon test (`pip_reference.h`, **no Clipper2 header included**) classifies the point against the *original input* polygon sets, then the boolean truth table combines the two verdicts (and/or/difference/xor). |

`pip_reference.h` is deliberately dependency-free and uses only exact
`__int128` integer arithmetic for cross products, so it cannot share a
bug with the library under test.

Points lying exactly on an input edge are skipped: boundary convention
(which side "owns" the edge) is a spec choice, not an interior/exterior
error.

## Boundary-band vs. deep mismatches

Clipper2's integer pipeline rounds polygon intersections to integer
coordinates (snap rounding, bounded by 0.5 unit). A sample point within
one integer unit of an input edge can therefore flip between inside and
outside even though the boolean logic is correct. Every mismatch is
classified:

- **boundary-band** — within `--band-tolerance` (default 1.0 integer
  unit) of an input edge; expected snap-rounding artifact.
- **deep** — farther from every input edge; a genuine logic discrepancy.
  A correct library must produce **zero** deep mismatches.

## Mutation injection

For each of the four boolean operations, one vertex of the generated
solution polygon is moved by several grid spacings (simulating a buggy
boolean result). The reference verdict stays anchored to the original
inputs, so the affected region must show up as mismatches. All four ops
must be detected; in practice the injected regions are entirely *deep*
mismatches.

## Build and run

```sh
make            # produces ./consistency_checker (no external deps)
./consistency_checker
```

Options: `--cases N`, `--grid N` (N×N sample points, default 100 = 1e4),
`--seed S`, `--out DIR`, `--band-tolerance T`, `--no-inject`.

`make diagnose` builds a helper that prints each mismatch's distance to
the nearest input and solution edges.

## Outputs (in `--out`, default `artifacts/consistency/`)

- `summary.json` — aggregate and per-case/per-op counts, including
  boundary-band and deep mismatch breakdowns and the mutation report.
- `mismatch_regions.geojson` — every mismatch region as a closed polygon
  ring with a `class` of `deep` or `boundary_band`, loadable in any GIS
  tool for locating discrepancies.
- `mutation_regions.geojson` — regions produced by the injected errors.
- `case_*.svg`, `mutation_<op>.svg` — visual overlays.

Exit code is non-zero if the overall mismatch rate reaches 0.01 %, if
any deep mismatch occurs, or if an injected mutation goes undetected.

## Acceptance results (default run)

- 50 random polygon pairs × 10 000 grid points × 4 ops = 2 000 000
  comparisons; mismatch rate ≈ 0.00035 %, all boundary-band, **0 deep**.
- 50-case run (including mutation tests) finishes in ~0.2 s (budget 60 s).
- Mutation injection is detected for all four boolean operations.
