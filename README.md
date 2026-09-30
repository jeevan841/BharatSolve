# BharatSolve — prototype

Indigenous LP / MILP solver core for SIH 2026, problem statement **26119**
(MRPL — Indigenous GPU-Accelerated Optimization Solver). This is the
tomorrow-morning prototype: an LP/MILP engine built from scratch (no
external solver libraries, per the problem statement), a CLI, a small
browser UI, a benchmark harness against HiGHS, and a crude-blending demo.

**Not yet built** (see `docs/SRS_v0.3.md` roadmap): QP/interior-point,
sparse LU (dense LU only, so it doesn't scale past a few thousand rows —
see `25fv47` in the benchmark below), GPU acceleration, cutting planes.

## Build

```
make            # builds bin/bharatsolve
make test       # builds and runs tests/test_core.cpp
```

Requires g++ (C++17) and Python 3 with Flask for the UI. No other
dependencies — deliberately, to match the problem statement's "built from
scratch" requirement.

## Run on the command line

```
./bin/bharatsolve data/afiro.mps --json solution.json
./bin/bharatsolve demo/crudeblend.mps --json blend.json --log 2
```

Exit codes: 0 optimal, 1 feasible (limit hit), 2 infeasible, 3 unbounded,
4 unverified, 5 limit with no solution, 10 input error, 11 unsupported (QP).

## Run the UI

```
python3 ui/server.py          # http://127.0.0.1:5055
```

Upload an `.mps` file or pick a sample, solve, and see the log, the
solution table, and the independent verification report.

## Run in a container / on a cloud VM

```
docker build -t bharatsolve .
docker run -p 5055:5055 bharatsolve
```

This was written to run the same way on a rented Linux VM or GitHub
Codespaces as on this development sandbox; it has **not** been
build-tested in this environment (no Docker available here) — test it
once before relying on it live.

## Benchmark evidence

```
python3 tools/benchmark.py     # writes bench_results.csv
```

Compares BharatSolve against SciPy's HiGHS backend on 8 small/medium
Netlib LP instances plus one MILP (`flugpl`). Current result: **7 of 8
match HiGHS's objective to 4+ significant figures**, with independent
verification passing on every solved instance. The 8th (`25fv47`, 821
rows) hits the time limit — this MVP's dense-LU engine doesn't yet scale
that far; sparse LU is the documented next step.

## Demo

`demo/crudeblend.mps` — a small refinery blending LP (minimize cost of
blending 4 streams into a gasoline pool subject to a volume target, a
sulfur ceiling and an octane floor). Meant to be solved live in the UI
during the presentation.
