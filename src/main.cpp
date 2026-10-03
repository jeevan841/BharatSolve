#include "bs.hpp"
#include <cstring>
#include <iostream>

using namespace bs;

static void usage() {
  std::cout << "BharatSolve - indigenous LP/MILP solver core (prototype)\n"
             << "Usage: bharatsolve <model.mps> [options]\n"
             << "  --time <sec>      time limit (default 60)\n"
             << "  --method <m>      auto|primal|dual (default auto)\n"
             << "  --no-presolve\n  --no-scale\n"
             << "  --json <path>     write solution as JSON\n"
             << "  --log <0-3>       verbosity (default 1)\n"
             << "  --gpu <mode>      pricing offload: off (default) | auto | on | emulate\n"
             << "  --gpu-min-nnz <n> 'auto' uses the GPU only above this many nonzeros (default 500000)\n";
}

int main(int argc, char** argv) {
  if (argc < 2) { usage(); return 10; }
  std::string path = argv[1];
  Options o;
  std::string jsonOut;
  for (int i = 2; i < argc; i++) {
    std::string a = argv[i];
    if (a == "--time" && i + 1 < argc) o.timeLimit = std::stod(argv[++i]);
    else if (a == "--method" && i + 1 < argc) o.method = argv[++i];
    else if (a == "--gpu" && i + 1 < argc) o.gpu = argv[++i];
    else if (a == "--gpu-min-nnz" && i + 1 < argc) o.gpuMinNnz = std::stol(argv[++i]);
    else if (a == "--no-presolve") o.presolve = false;
    else if (a == "--no-scale") o.scale = false;
    else if (a == "--json" && i + 1 < argc) jsonOut = argv[++i];
    else if (a == "--log" && i + 1 < argc) o.log = std::stoi(argv[++i]);
    else if (a == "--gap" && i + 1 < argc) o.gap = std::stod(argv[++i]);
    else if (a == "-h" || a == "--help") { usage(); return 0; }
  }
  Model M; std::string err;
  if (!readMps(path, M, err)) { std::cerr << "input error: " << err << "\n"; return exitCode(Status::InputError); }
  if (o.log >= 1) std::cerr << "BharatSolve: " << M.name << "  rows=" << M.m << " cols=" << M.n
                             << " nnz=" << M.nnz() << (M.isMip() ? " (MILP)" : (M.hasQuad() ? " (QP)" : " (LP)")) << "\n";
  Solution S = solve(M, o);
  if (o.log >= 1) {
    std::cerr << "status=" << statusName(S.status) << " obj=" << S.obj << " time=" << S.seconds << "s"
               << " iters=" << S.iters << " nodes=" << S.nodes
               << " verify=" << (S.ver.ran ? (S.ver.ok ? "ok" : "FAILED") : "n/a") << "\n";
    if (!S.message.empty() && S.status != Status::Optimal) std::cerr << "  note: " << S.message << "\n";
    if (!S.ver.ok) for (auto& f : S.ver.failures) std::cerr << "  verify: " << f << "\n";
  }
  if (!jsonOut.empty()) writeJson(M, S, jsonOut);
  else if (o.log >= 2) writeJson(M, S, "/dev/stdout");
  return exitCode(S.status);
}
