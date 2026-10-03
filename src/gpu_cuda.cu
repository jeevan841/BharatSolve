// CUDA backend for the pricing step. Compiled only with `make CUDA=1` (needs nvcc).
// NOTE: this file has been written against the CUDA runtime API but has NOT been compiled or run
// on GPU hardware in the environment where it was developed. The kernel's arithmetic is the shared
// priceColumn() function, which IS tested on the CPU (tests/test_core.cpp, emulated device).
#include "gpu.hpp"
#ifdef BS_MOCK_CUDA            // CPU stand-in used by `make test-cuda-mock`
#include "mock_cuda_runtime.h"
MockDim blockIdx, threadIdx, blockDim;
#define BS_LAUNCH(kernel, blocks, threads, ...) mockLaunch(blocks, threads, [&]() { kernel(__VA_ARGS__); })
#else
#include <cuda_runtime.h>
#define BS_LAUNCH(kernel, blocks, threads, ...) kernel<<<blocks, threads>>>(__VA_ARGS__)
#endif

namespace bs {

__global__ void priceKernel(int n, int N, const int* cs, const int* ri, const double* av,
                            const double* cost, const double* y, double* d) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j < N) d[j] = priceColumn(j, n, cs, ri, av, cost, y);
}

namespace {
class CudaPricer : public PricingBackend {
 public:
  explicit CudaPricer(const Model& M) : n_(M.n), m_(M.m), N_(M.n + M.m) {
    int dev = 0;
    if (cudaGetDeviceCount(&dev) != cudaSuccess || dev == 0) { err_ = "no CUDA device found"; return; }
    if (!alloc((void**)&cs_, sizeof(int) * (n_ + 1)) || !alloc((void**)&ri_, sizeof(int) * M.ri.size()) ||
        !alloc((void**)&av_, sizeof(double) * M.av.size()) || !alloc((void**)&cost_, sizeof(double) * N_) ||
        !alloc((void**)&y_, sizeof(double) * m_) || !alloc((void**)&d_, sizeof(double) * N_)) return;
    std::vector<int> cs(M.cs.begin(), M.cs.end()), ri(M.ri.begin(), M.ri.end());
    std::vector<double> av(M.av.begin(), M.av.end());
    if (!copyUp(cs_, cs.data(), sizeof(int) * (n_ + 1)) || !copyUp(ri_, ri.data(), sizeof(int) * ri.size()) ||
        !copyUp(av_, av.data(), sizeof(double) * av.size())) return;
    ok_ = true;
  }
  ~CudaPricer() override { cudaFree(cs_); cudaFree(ri_); cudaFree(av_); cudaFree(cost_); cudaFree(y_); cudaFree(d_); }
  bool ok() const { return ok_; }
  const std::string& error() const { return err_; }
  const char* name() const override { return "cuda"; }
  void computeD(const double* cost, const double* y, double* d) override {
    // cost and y change between calls, so both are copied every time (N + m doubles).
    if (!copyUp(cost_, cost, sizeof(double) * N_) || !copyUp(y_, y, sizeof(double) * m_)) return fail(d, cost, y);
    const int threads = 128, blocks = (N_ + threads - 1) / threads;
    BS_LAUNCH(priceKernel, blocks, threads, n_, N_, cs_, ri_, av_, cost_, y_, d_);
    if (cudaGetLastError() != cudaSuccess) return fail(d, cost, y);
    if (cudaMemcpy(d, d_, sizeof(double) * N_, cudaMemcpyDeviceToHost) != cudaSuccess) return fail(d, cost, y);
  }
 private:
  bool alloc(void** p, size_t bytes) {
    if (bytes == 0) bytes = 1;
    if (cudaMalloc(p, bytes) != cudaSuccess) { err_ = "cudaMalloc failed"; return false; }
    return true;
  }
  bool copyUp(void* dst, const void* src, size_t bytes) {
    if (bytes == 0) return true;
    if (cudaMemcpy(dst, src, bytes, cudaMemcpyHostToDevice) != cudaSuccess) { err_ = "cudaMemcpy failed"; return false; }
    return true;
  }
  // If a CUDA call fails mid-solve, compute this request on the CPU so results stay correct.
  void fail(double* d, const double* cost, const double* y) {
    for (int j = 0; j < N_; j++) d[j] = priceColumn(j, n_, hcs(), hri(), hav(), cost, y);
  }
  const int* hcs() const { return hostCs_.data(); }
  const int* hri() const { return hostRi_.data(); }
  const double* hav() const { return hostAv_.data(); }
  int n_, m_, N_;
  int *cs_ = nullptr, *ri_ = nullptr;
  double *av_ = nullptr, *cost_ = nullptr, *y_ = nullptr, *d_ = nullptr;
  bool ok_ = false;
  std::string err_;
  std::vector<int> hostCs_, hostRi_;
  std::vector<double> hostAv_;
 public:
  void keepHostCopy(const Model& M) { hostCs_.assign(M.cs.begin(), M.cs.end()); hostRi_.assign(M.ri.begin(), M.ri.end()); hostAv_.assign(M.av.begin(), M.av.end()); }
};
}  // namespace

std::shared_ptr<PricingBackend> makeCudaPricer(const Model& M, std::string& note) {
  std::shared_ptr<CudaPricer> p = std::make_shared<CudaPricer>(M);
  if (!p->ok()) { note = p->error(); return nullptr; }
  p->keepHostCopy(M);
  return p;
}

}  // namespace bs
