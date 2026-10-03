// Optional GPU offload of the pricing step of the simplex method.
//
// Pricing computes the reduced cost of every column of [A | -I]:
//     d[j] = cost[j] - a_j . y        (a sparse matrix-transpose times vector product)
// This is the one step of the simplex iteration that is naturally data-parallel, so it is the
// one that can be offloaded. The basis solves (LU factorization, triangular solves) are
// sequential and stay on the CPU.
//
// Backends: the plain CPU loop (default, no backend object), a CUDA backend (built with
// `make CUDA=1`), and an "emulated device" backend that runs the exact device code path on the
// CPU so the offload logic can be tested on machines without a GPU.
#pragma once
#include "bs.hpp"
#include <memory>
#include <string>
#include <vector>

#ifdef __CUDACC__
#define BS_HD __host__ __device__
#else
#define BS_HD
#endif

namespace bs {

// Shared by the CPU reference, the emulated device and the CUDA kernel, so all three compute
// the same arithmetic in the same order.
BS_HD inline double priceColumn(int j, int n, const int* cs, const int* ri, const double* av,
                                const double* cost, const double* y) {
  if (j < n) {
    double s = 0;
    for (int p = cs[j]; p < cs[j + 1]; p++) s += av[p] * y[ri[p]];
    return cost[j] - s;
  }
  return cost[j] - (-y[j - n]);                      // slack column is -e_i
}

class PricingBackend {
 public:
  virtual ~PricingBackend() {}
  virtual const char* name() const = 0;
  // d must have n + m entries; cost n + m; y m.
  virtual void computeD(const double* cost, const double* y, double* d) = 0;
};

// Runs the device code path (separate buffers, per-column kernel function) on the CPU.
std::shared_ptr<PricingBackend> makeEmulatedPricer(const Model& M);
// CUDA backend; returns nullptr and sets `note` if the build has no CUDA or no device is found.
std::shared_ptr<PricingBackend> makeCudaPricer(const Model& M, std::string& note);
// Chooses a backend from Options::gpu (off | auto | on | emulate). nullptr means "use the CPU loop".
std::shared_ptr<PricingBackend> makePricer(const Model& M, const Options& o, std::string& note);

}  // namespace bs
