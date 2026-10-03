# Corrected slide text (matches the repository)

Tag convention: **[Built]** works today and is tested. **[Roadmap]** planned, not built.

## Slide 2: BHARATSOLVE, Proposed Solution
- **[Built]** Sovereign LP and MILP solver core in C++17, with a command-line tool and a browser UI
- **[Built]** Continuous engine: primal and dual simplex with sparse LU basis updates
- **[Built]** Robustness: scaling, Harris ratio test, cost perturbation
- **[Built]** Presolve and postsolve to the original model
- **[Built]** MILP engine: branch-and-bound, warm-started dual re-solves, pseudocost branching
- **[Built]** Verification layer: feasibility, duality gap, complementarity and integrality; never called optimal unless it passes
- **[Implemented, awaiting GPU hardware validation]** Optional CUDA offload of the pricing step (off by default; CPU fallback; tested against a CPU stand-in)
- **[Roadmap]** Cutting planes and heuristics, convex QP with interior-point (the strongest GPU use case), C++ API and Python bindings, multi-core

## Slide 3: TECHNICAL APPROACH
Technologies
- **[Built]** C++17 core, standard library only
- **[Built]** Sparse model storage and an in-house sparse LU (Markowitz pivoting)
- **[Built]** Interfaces: command line, JSON output, web UI. **[Roadmap]** C++ API and Python bindings
- **[Built]** Input: MPS (fixed and free). **[Roadmap]** QPLIB / QPS
- **[Built]** Testing: regression suite (126 checks) and a random-model fuzzer cross-checked against HiGHS, used only in the harness
- **[Implemented, unvalidated on hardware]** CUDA kernel for pricing (`--gpu`, off by default). **[Roadmap]** CI, multi-core node solves

Methodology
1. Parse, scale, presolve
2. Root LP: primal or dual simplex with anti-degeneracy safeguards
3. Branch and bound: warm-started children, pseudocost branching
4. Verify on the original model
5. Report: status, objective, gap, log, JSON
6. Validate: Netlib LP instances and a MILP against HiGHS. **[Roadmap]** MIPLIB, QPLIB

## Slide 4: FEASIBILITY AND VIABILITY (changes only)
- Prototype status: correct on 8 of 8 benchmark instances, and on thousands of random LP and MILP models, against HiGHS
- Phased build: LP core and MILP **[Built]**, then cuts, QP and interior-point **[Roadmap]**
- Roadmap line: "tuning, cutting planes, QP, MIQP/NLP/MINLP, GPU where it pays off"
- Risk table: change "Bland fallback" to "perturbation (Bland fallback planned)"; change "Slower than mature solvers" mitigation to "Sparse LU, warm starts, honest small and medium benchmarks"

## Slide 5: IMPACT AND BENEFITS (changes only)
- Performance row: **Commercial: Best | Open-source: Good | BharatSolve: Prototype, about 5 times slower than HiGHS on the largest LP tested**
- Keep: "Honest positioning: no claim of speed parity today"
- Add GPU line: "Measured: pricing is 15 to 21 % of run time on tested models, so GPU offload can give at most about 1.2 to 1.3 times there; GPUs are expected to pay off with interior-point on large models"
- Add one line: "Measured: 8 of 8 benchmark instances match HiGHS; 25fv47 solves in about 1 s against 0.2 s"
