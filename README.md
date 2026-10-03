# BharatSolve — prototype

Indigenous LP / MILP solver core for SIH 2026, problem statement **26119**
(MRPL — Indigenous GPU-Accelerated Optimization Solver). This is a
working prototype: an LP/MILP engine built from scratch (no
external solver libraries, per the problem statement), a CLI, a small
browser UI, a benchmark harness against HiGHS, and a crude-blending demo.

**Not yet built** (full list in `docs/BUILT_VS_ROADMAP.md`): QP/interior-point,
cutting planes, C++ API and Python bindings. GPU offload of the pricing step is implemented but
not yet validated on real hardware (off by default; see `docs/GPU.md`).
The basis uses a sparse LU (Markowitz pivoting) with product-form updates. The
`--method dual` option runs a dual simplex with cost shifting and finishes with a short
primal clean-up, so its answers are always verified. It is correct on every sample model but is
not faster than the default primal simplex here, so primal stays the default.

## Build

```
make            # builds bin/bharatsolve
make test       # builds and runs tests/test_core.cpp (about 1 second)
make test-slow  # also solves 25fv47 (about 2 seconds)
make test-cuda-mock  # tests the CUDA backend's host code against a CPU stand-in (no GPU needed)
make CUDA=1     # build with the CUDA backend (needs nvcc); run with --gpu on|auto
python3 tools/fuzz.py 500 1   # random LP/MILP models checked against HiGHS (needs SciPy)
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

Upload an `.mps` file or pick a sample (the crude-blend demo is listed first),
solve, and see the log, the solution table, and the independent verification report.

## Run in a container / on a cloud VM

```
docker build -t bharatsolve .
docker run -p 5055:5055 bharatsolve
```

The Docker image itself has **not** been built (Docker was not available where this
was prepared). The same layout (`src`, `Makefile`, `ui`, `data`, `demo`) was
checked outside Docker: it builds and the UI solves `crudeblend.mps`. Run
`docker build` once before relying on it live.

## Benchmark evidence

```
python3 tools/benchmark.py     # writes bench_results.csv
```

Compares BharatSolve against SciPy's HiGHS backend on 8 small/medium
Netlib LP instances plus one MILP (`flugpl`). Current result: **8 of 8 match
HiGHS's objective to 4+ significant figures**, with independent verification
passing on every instance.

BharatSolve is slower than HiGHS, as expected for a prototype.
The largest instance, `25fv47` (821 rows), takes about 1 s against about 0.2 s
for HiGHS (roughly 5 times slower). The MILP `flugpl` takes about 0.2 s against 0.12 s
(branch-and-bound warm-starts each node from its parent's basis with the dual simplex,
and picks branching variables by pseudocosts). Remaining speed work: better pricing (steepest edge), presolve improvements and cutting planes.

If a solve stops early, the status message says why: time limit, iteration limit
or numerical failure (basis became singular).

## Demo

`demo/crudeblend.mps` — a small refinery blending LP (minimize cost of
blending 4 streams into a gasoline pool subject to a volume target, a
sulfur ceiling and an octane floor). Meant to be solved live in the UI
during the presentation.
