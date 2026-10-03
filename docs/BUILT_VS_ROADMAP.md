# BharatSolve: what is built and what is planned

Status as of this release. Every "Built" item was exercised by the regression suite
(`make test`, 126 checks), the HiGHS benchmark (`tools/benchmark.py`) or the random-model
fuzzer (`tools/fuzz.py`). Nothing is listed as built unless it runs.

| Area | Status | Evidence / notes |
|---|---|---|
| C++17 core, standard library only | **Built** | No external solver or linear-algebra libraries. |
| MPS reader (fixed and free format), RANGES, BOUNDS, integer markers, OBJSENSE | **Built** | Malformed numbers give exit code 10, never a crash. |
| QPS / QPLIB input | Roadmap | QP models are reported as unsupported (exit 11). |
| Scaling | **Built** | Applied before solving, undone before verification. |
| Presolve and postsolve | **Built** (basic) | Fixed variables, empty rows; solution restored to the original space. |
| Primal simplex (default) | **Built** | Phase 1 / phase 2, bounded variables, Harris-style two-pass ratio test, cost perturbation. |
| Dual simplex | **Built** | Cost-shifting start, finishes with a short primal clean-up. Correct on all samples and 2,700+ fuzz models. Not faster than primal on the samples, so primal is the default. |
| Sparse LU of the basis | **Built** | Markowitz pivoting with a stability threshold, product-form updates, refactor every 100 pivots. |
| Branch-and-bound (MILP) | **Built** | Warm start from the parent basis with the dual simplex, pseudocost branching, depth-first. |
| Cutting planes (Gomory / MIR / cover) | Roadmap | Not implemented. |
| Primal heuristics (diving, rounding) | Roadmap | Not implemented. |
| Interior-point method and QP | Roadmap | Not implemented. |
| Bland's rule fallback and stalling detection | Roadmap | Perturbation exists; Bland fallback does not. |
| Independent verification layer | **Built** | Checks primal feasibility, dual feasibility, duality gap, complementarity and integrality on the original model. A result is never "optimal" unless it passes. |
| Infeasible / unbounded detection | **Built** | Exit codes 2 and 3, tested including integer models with an unbounded relaxation. |
| Honest limit reporting | **Built** | Time limit, iteration limit and numerical failure are reported separately. |
| Command-line tool, JSON output, exit codes | **Built** | Codes 0 to 5, 10, 11. JSON also reports the pricing backend and time. |
| Web UI | **Built** | Flask app; the crude-blend demo is listed first. |
| Docker image | Written, not built | Dockerfile exists; `docker build` has not been run yet. |
| Regression tests | **Built** | `make test`, plus `make test-slow` for 25fv47. |
| Random-model fuzzer against HiGHS | **Built** | `tools/fuzz.py`, LP and MILP, zero mismatches on thousands of models. |
| Benchmark harness against HiGHS | **Built** | `tools/benchmark.py`: 8 of 8 instances match HiGHS. |
| Continuous integration | Roadmap | No CI workflow in the repository yet. |
| C++ API / Python bindings | Roadmap | Command line and JSON only for now. |
| Multi-core node solves | Roadmap | Single-threaded. |
| GPU offload of pricing (CUDA) | **Implemented, not validated on hardware** | Off by default. Host code tested against a CPU stand-in for CUDA; never compiled with nvcc or run on a GPU. Pricing is only 15 to 21 % of run time on tested sizes, so the speed-up bound there is about 1.2 to 1.3 times. See `docs/GPU.md`. |
| GPU for interior-point / batched node LPs | Roadmap | Where a GPU is expected to pay off most. |
| MIQP / NLP / MINLP | Roadmap | Not implemented. |

## Measured performance (this release)

| Instance | BharatSolve | HiGHS | Note |
|---|---|---|---|
| 25fv47 (LP, 821 rows) | about 1.0 s | about 0.2 s | about 5 times slower |
| flugpl (MILP) | about 0.2 s | about 0.12 s | about 1.6 times slower |
| Other 6 LP instances | under 0.02 s | under 0.01 s | all objectives match |

Small and medium Netlib LP instances and one MILP only. There is no claim of speed parity with
mature solvers, and no result for large-scale industrial models yet.
