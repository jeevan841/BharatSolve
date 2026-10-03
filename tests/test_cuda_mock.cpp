// Exercises src/gpu_cuda.cu (host glue + kernel) against a CPU stand-in for the CUDA runtime.
// Run with: make test-cuda-mock
#include "bs.hpp"
#include "gpu.hpp"
#include <cmath>
#include <iostream>
using namespace bs;
static int fails = 0, passes = 0;
#define CHECK(c, msg) do { if (c) passes++; else { fails++; std::cerr << "  FAIL: " << msg << "\n"; } } while (0)

int main() {
  const char* files[] = {"data/afiro.mps", "data/adlittle.mps", "data/e226.mps", "data/israel.mps", "demo/crudeblend.mps", "data/flugpl.mps"};
  for (const char* f : files) {
    Model M; std::string err;
    CHECK(readMps(f, M, err), std::string("read ") + f);
    if (M.n == 0) continue;
    std::string note;
    std::shared_ptr<PricingBackend> cu = makeCudaPricer(M, note);
    CHECK(cu != nullptr, "mock CUDA backend created: " + note);
    if (!cu) continue;
    std::shared_ptr<PricingBackend> emu = makeEmulatedPricer(M);
    std::vector<double> cost(M.n + M.m), y(M.m), a(M.n + M.m), b(M.n + M.m);
    unsigned long long st = 99;
    auto rnd = [&]() { st = st * 6364136223846793005ULL + 1442695040888963407ULL; return ((st >> 33) % 2001) / 100.0 - 10.0; };
    for (double& v : cost) v = rnd();
    for (double& v : y) v = rnd();
    cu->computeD(cost.data(), y.data(), a.data());
    emu->computeD(cost.data(), y.data(), b.data());
    bool same = true; for (size_t j = 0; j < a.size(); j++) if (a[j] != b[j]) same = false;
    CHECK(same, std::string("CUDA glue + kernel matches the emulated device on ") + f);
    Options o; o.log = 0; o.timeLimit = 60;
    Solution P = solve(M, o);
    o.gpu = "on";
    Solution G = solve(M, o);
    CHECK(G.pricing == "cuda", std::string("backend reported as cuda on ") + f);
    CHECK(G.status == P.status && G.iters == P.iters && std::fabs(G.obj - P.obj) <= 1e-12 * std::max(1.0, std::fabs(P.obj)),
          std::string("full solve through the CUDA glue matches the CPU on ") + f);
  }
  std::cout << passes << " checks passed, " << fails << " failed\n";
  return fails ? 1 : 0;
}
