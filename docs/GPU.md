# GPU offload: what exists, what was measured, what is not yet proven

## What is offloaded
The **pricing step** of the simplex method: computing the reduced cost of every column,
`d[j] = cost[j] - a_j . y`. It is a sparse matrix-transpose times vector product, which is
naturally data-parallel (one thread per column).

Everything else in a simplex iteration (LU factorization, triangular solves, ratio test) is
sequential or has tiny parallelism, so it stays on the CPU. GPUs are not a good fit there.

## How it works
| Piece | File | Status |
|---|---|---|
| Shared per-column pricing function (same arithmetic for CPU, emulated device and CUDA) | `src/gpu.hpp` | Built, tested |
| Backend interface and chooser (`--gpu off / auto / on / emulate`) | `src/gpu.hpp`, `src/gpu.cpp` | Built, tested |
| Emulated device (runs the device code path on the CPU, separate buffers) | `src/gpu.cpp` | Built, tested |
| CUDA backend (matrix uploaded once, one thread per column, `--fmad=false`, CPU fallback if a CUDA call fails) | `src/gpu_cuda.cu` | **Written; compiled and run only against a CPU stand-in for the CUDA runtime. Not compiled with `nvcc` and not run on a GPU.** |
| Branch-and-bound shares one uploaded matrix across all nodes | `src/solve.cpp` | Built, tested |

Defaults: the GPU is **off**. `--gpu auto` only uses it above `--gpu-min-nnz` (default 500,000
nonzeros). If a GPU is requested but unavailable, the solver says so and uses the CPU.

Build with CUDA: `make CUDA=1` (needs `nvcc`). Run: `./bin/bharatsolve model.mps --gpu on --log 2`.
The JSON output reports `pricing_backend` and `pricing_seconds`.

## What was verified (no GPU was available)
- `make test`: the emulated-device path is bit-identical to a plain reference on 6 models, and a
  full solve through it matches the CPU step for step (same status, same iteration count),
  including the MILP `flugpl`.
- `make test-cuda-mock`: the real `gpu_cuda.cu` host code and kernel compiled as C++ against a CPU
  stand-in for the CUDA runtime. 30 checks pass: matrix upload, vector copies, launch geometry,
  results identical to the emulated device, full solves identical to the CPU.
- Requesting a GPU where none exists falls back to the CPU and still solves correctly.

## What was NOT verified
- Compilation with `nvcc`, execution on real hardware, and any GPU timing. The mock cannot catch
  mistakes that only real CUDA exposes (device memory rules, launch errors, driver issues).
- So this must be described as "implemented, awaiting hardware validation", not "GPU-accelerated".

## Measured: how much of the run time could a GPU even touch?
Wall time spent in the pricing step on this CPU (single core):

| Model | Size | Pricing share of run time |
|---|---|---|
| 25fv47 (Netlib) | 821 x 1,571, 10,400 nonzeros | about 15 % |
| Synthetic random LP | 3,000 x 9,000, about 75,000 nonzeros | about 21 % |

By Amdahl's law, even an infinitely fast GPU pricing step would speed these runs up by only about
**1.2 to 1.3 times**. At these sizes a GPU transfer round trip (tens of microseconds) is also
comparable to the whole CPU pricing cost, so a real GPU would likely show **no gain** here.

Expected crossover (an estimate from typical GPU latency and bandwidth, not a measurement): the
offload can only help when the matrix has on the order of a million or more nonzeros. That is why
`--gpu auto` defaults to 500,000.

## Findings that matter for the slides
- The risk-table line "GPU gives no benefit: optional, off by default, publish timings either way"
  is exactly the situation measured here. Keep it, and add the 15 to 21 % figure.
- The solver itself still struggles above about 10^4 variables (the synthetic 3,000 x 9,000 model
  did not finish in 30 s), so the practical bottleneck at scale is the number of simplex iterations
  and the basis solves, not pricing.

## Where a GPU would really pay (roadmap)
1. **Interior-point method:** each iteration solves normal equations with a large sparse
   factorization, which maps well to GPU sparse linear algebra. This is the strongest case.
2. **Batched strong branching / node LPs** in branch-and-bound.
3. **Pivot-row products** in the dual simplex (same kernel shape as pricing).

## To validate on a real GPU (about 15 minutes)
```
make clean && make CUDA=1
make CUDA=1 test                          # all tests, with the CUDA backend linked in
./bin/bharatsolve data/25fv47.mps --gpu on --log 2
./bin/bharatsolve BIG.mps --gpu on --json out.json   # compare pricing_seconds with --gpu off
```
Expected: identical status, objective and iteration count with and without the GPU, since the
arithmetic is the same (`--fmad=false`). Report `pricing_seconds` for both. Publish the numbers
whichever way they go.
