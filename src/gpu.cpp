#include "gpu.hpp"

namespace bs {

namespace {
// Copies the matrix into private buffers and prices column by column, exactly as the GPU kernel
// would (one "thread" per column). Used to test the offload path without a GPU.
class EmulatedDevicePricer : public PricingBackend {
 public:
  explicit EmulatedDevicePricer(const Model& M)
      : n_(M.n), m_(M.m), cs_(M.cs.begin(), M.cs.end()), ri_(M.ri.begin(), M.ri.end()), av_(M.av.begin(), M.av.end()),
        cost_(M.n + M.m), y_(M.m) {}
  const char* name() const override { return "emulated-device"; }
  void computeD(const double* cost, const double* y, double* d) override {
    std::copy(cost, cost + n_ + m_, cost_.begin());     // host-to-device copies
    std::copy(y, y + m_, y_.begin());
    for (int j = 0; j < n_ + m_; j++)                   // "kernel launch"
      d[j] = priceColumn(j, n_, cs_.data(), ri_.data(), av_.data(), cost_.data(), y_.data());
  }
 private:
  int n_, m_;
  std::vector<int> cs_, ri_;
  std::vector<double> av_, cost_, y_;
};
}  // namespace

std::shared_ptr<PricingBackend> makeEmulatedPricer(const Model& M) {
  return std::make_shared<EmulatedDevicePricer>(M);
}

#ifndef BS_WITH_CUDA
std::shared_ptr<PricingBackend> makeCudaPricer(const Model&, std::string& note) {
  note = "this build has no CUDA support (rebuild with: make CUDA=1)";
  return nullptr;
}
#endif

std::shared_ptr<PricingBackend> makePricer(const Model& M, const Options& o, std::string& note) {
  note = "cpu";
  if (o.gpu == "off") return nullptr;
  if (o.gpu == "emulate") { note = "emulated-device (CPU, for testing)"; return makeEmulatedPricer(M); }
  if (o.gpu == "auto" && (long)M.nnz() < o.gpuMinNnz) { note = "cpu (auto: matrix too small for the GPU to pay off)"; return nullptr; }
  if (o.gpu == "on" || o.gpu == "auto") {
    std::string why;
    std::shared_ptr<PricingBackend> p = makeCudaPricer(M, why);
    if (p) { note = "cuda"; return p; }
    note = "cpu (GPU requested but unavailable: " + why + ")";
    return nullptr;
  }
  note = "cpu (unknown --gpu value '" + o.gpu + "', using the CPU)";
  return nullptr;
}

}  // namespace bs
