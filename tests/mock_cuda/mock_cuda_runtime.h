// A tiny CPU stand-in for the CUDA runtime, used ONLY to test the host-side glue of
// src/gpu_cuda.cu (buffer sizes, copies, launch geometry, fallback) on machines without a GPU.
// It does not model real GPU behaviour (memory spaces, concurrency, timing).
#pragma once
#include <cstdlib>
#include <cstring>
#include <functional>

#define __global__
struct MockDim { unsigned x = 0; };
extern MockDim blockIdx, threadIdx, blockDim;
enum cudaError_t { cudaSuccess = 0, cudaErrorMemoryAllocation = 2, cudaErrorInvalidValue = 1 };
enum cudaMemcpyKind { cudaMemcpyHostToDevice, cudaMemcpyDeviceToHost };
inline const char* cudaGetErrorString(cudaError_t) { return "mock cuda error"; }
inline cudaError_t cudaGetDeviceCount(int* n) { *n = 1; return cudaSuccess; }
inline cudaError_t cudaMalloc(void** p, size_t bytes) { *p = std::malloc(bytes); return *p ? cudaSuccess : cudaErrorMemoryAllocation; }
inline cudaError_t cudaFree(void* p) { std::free(p); return cudaSuccess; }
inline cudaError_t cudaMemcpy(void* dst, const void* src, size_t bytes, cudaMemcpyKind) { std::memcpy(dst, src, bytes); return cudaSuccess; }
inline cudaError_t cudaGetLastError() { return cudaSuccess; }
inline void mockLaunch(unsigned blocks, unsigned threads, const std::function<void()>& kernel) {
  blockDim.x = threads;
  for (unsigned b = 0; b < blocks; b++) for (unsigned t = 0; t < threads; t++) { blockIdx.x = b; threadIdx.x = t; kernel(); }
}
