// BharatSolve core regression tests. Plain asserts, no framework.
// Run from the repository root:  make test
#include "bs.hpp"
#include "gpu.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

using namespace bs;

static int g_fail = 0, g_pass = 0;
#define CHECK(cond, msg) do { if (cond) { g_pass++; } else { g_fail++; std::cerr << "  FAIL: " << msg << "  [" #cond "]\n"; } } while (0)

static Options quiet() { Options o; o.log = 0; o.timeLimit = 60; return o; }

static bool relClose(double a, double b, double tol = 1e-4) {
  return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b));
}

static bool loadString(const std::string& text, Model& M) {
  std::istringstream in(text);
  std::string err;
  return readMpsStream(in, M, err);
}

static Solution solveString(const std::string& text, Model& M) {
  Solution S;
  if (!loadString(text, M)) { S.status = Status::InputError; return S; }
  return solve(M, quiet());
}

static void testFile(const std::string& path, double expectedObj) {
  std::cout << "[file] " << path << "\n";
  Model M; std::string err;
  CHECK(readMps(path, M, err), "cannot read " + path + ": " + err);
  if (M.n == 0) return;
  Solution S = solve(M, quiet());
  CHECK(S.status == Status::Optimal, path + " status is not optimal");
  CHECK(exitCode(S.status) == 0, path + " exit code should be 0");
  CHECK(S.ver.ran && S.ver.ok, path + " verification should pass");
  CHECK(relClose(S.obj, expectedObj), path + " objective " + std::to_string(S.obj) + " expected " + std::to_string(expectedObj));
}

// ---- inline models --------------------------------------------------------
static const char* INFEASIBLE =
  "NAME INF\nROWS\n N obj\n G r1\n L r2\nCOLUMNS\n x obj 1 r1 1\n x r2 1\nRHS\n rhs r1 5 r2 2\nENDATA\n";
static const char* UNB_GE =       // min -x  s.t. x >= 1
  "NAME UNB1\nROWS\n N obj\n G r1\nCOLUMNS\n x obj -1 r1 1\nRHS\n rhs r1 1\nENDATA\n";
static const char* UNB_LE =       // min -x-y  s.t. -x + y <= 4
  "NAME UNB2\nROWS\n N obj\n L r1\nCOLUMNS\n x obj -1 r1 -1\n y obj -1 r1 1\nRHS\n rhs r1 4\nENDATA\n";
static const char* UNB_MIX =      // min -x+y  s.t. x + y >= 1
  "NAME UNB3\nROWS\n N obj\n G r1\nCOLUMNS\n x obj -1 r1 1\n y obj 1 r1 1\nRHS\n rhs r1 1\nENDATA\n";
static const char* MAXIMIZE =     // max x  s.t. x <= 4
  "NAME MAXT\nOBJSENSE\n    MAX\nROWS\n N obj\n L r1\nCOLUMNS\n x obj 1 r1 1\nRHS\n rhs r1 4\nENDATA\n";

static void testStatuses() {
  std::cout << "[inline] status and exit codes\n";
  Model M;
  Solution S = solveString(INFEASIBLE, M);
  CHECK(S.status == Status::Infeasible, "infeasible model must be reported infeasible");
  CHECK(exitCode(S.status) == 2, "infeasible exit code is 2");

  const char* unb[] = {UNB_GE, UNB_LE, UNB_MIX};
  const char* names[] = {"unbounded (x>=1, min -x)", "unbounded (-x+y<=4)", "unbounded (min -x+y, x+y>=1)"};
  for (int i = 0; i < 3; i++) {
    Model U; Solution T = solveString(unb[i], U);
    CHECK(T.status == Status::Unbounded, std::string(names[i]) + " must be unbounded");
    CHECK(exitCode(T.status) == 3, std::string(names[i]) + " exit code is 3");
  }

  Model X; Solution Z = solveString(MAXIMIZE, X);
  CHECK(Z.status == Status::Optimal, "maximize model should be optimal");
  CHECK(relClose(Z.obj, 4.0), "maximize objective should be 4, got " + std::to_string(Z.obj));
}

static void testHonestLimits() {
  std::cout << "[limits] a time limit must be reported as a time limit\n";
  Model M; std::string err;
  CHECK(readMps("data/e226.mps", M, err), "read e226");
  Options o = quiet(); o.timeLimit = 1e-9;
  Solution S = solve(M, o);
  CHECK(S.status == Status::LimitNoSol, "tiny time limit gives limit_no_solution");
  CHECK(exitCode(S.status) == 5, "limit with no solution exits 5");
  CHECK(S.message.find("time limit") != std::string::npos, "message should say time limit, got: " + S.message);
}

static void testRobustness() {
  std::cout << "[robustness] malformed numbers and integer problems\n";
  Model M;
  CHECK(!loadString("NAME X\nROWS\n N o\nCOLUMNS\n x o abc\nENDATA\n", M), "non-numeric coefficient must be an input error, not a crash");
  // integer problem whose LP relaxation is unbounded: min -x, x integer, x >= 1
  Model U; Solution S = solveString(
    "NAME IU\nROWS\n N o\n G r1\nCOLUMNS\n MARKER 'MARKER' 'INTORG'\n x o -1 r1 1\n MARKER 'MARKER' 'INTEND'\nRHS\n rhs r1 1\nENDATA\n", U);
  CHECK(S.status == Status::Unbounded, "MILP with unbounded relaxation must not be reported infeasible");
  // small bounded MILP: max 3x+2y, x+y<=4.5, x,y integer in [0,10] -> 12 at x=4,y=0
  Model K; Solution T = solveString(
    "NAME KN\nOBJSENSE\n    MAX\nROWS\n N o\n L r1\nCOLUMNS\n MARKER 'MARKER' 'INTORG'\n x o 3 r1 1\n y o 2 r1 1\n MARKER 'MARKER' 'INTEND'\nRHS\n rhs r1 4.5\nBOUNDS\n UP b x 10\n UP b y 10\nENDATA\n", K);
  CHECK(T.status == Status::Optimal && relClose(T.obj, 12.0), "small MILP optimum should be 12, got " + std::to_string(T.obj));
}

static void testDualMethod() {
  std::cout << "[dual] the dual simplex must agree with the primal simplex\n";
  Options o = quiet(); o.method = "dual";
  const char* files[] = {"data/afiro.mps", "data/avgas.mps", "data/chip.mps", "data/adlittle.mps",
                         "data/e226.mps", "data/israel.mps", "demo/crudeblend.mps"};
  for (const char* f : files) {
    Model M; std::string err;
    CHECK(readMps(f, M, err), std::string("read ") + f);
    if (M.n == 0) continue;
    Solution P = solve(M, quiet()), D = solve(M, o);
    CHECK(D.status == Status::Optimal && D.ver.ok, std::string("dual optimal and verified on ") + f);
    CHECK(relClose(D.obj, P.obj, 1e-6), std::string("dual objective matches primal on ") + f);
  }
  Model I; CHECK(loadString(INFEASIBLE, I), "parse infeasible");
  CHECK(solve(I, o).status == Status::Infeasible, "dual reports infeasible model as infeasible");
  Model U; CHECK(loadString(UNB_MIX, U), "parse unbounded");
  CHECK(solve(U, o).status == Status::Unbounded, "dual reports unbounded model as unbounded");
}

static void testGpuOffloadPath() {
  std::cout << "[gpu] the offload path (emulated device) must reproduce the CPU exactly\n";
  const char* files[] = {"data/afiro.mps", "data/adlittle.mps", "data/e226.mps", "data/israel.mps", "demo/crudeblend.mps", "data/flugpl.mps"};
  for (const char* f : files) {
    Model M; std::string err;
    CHECK(readMps(f, M, err), std::string("read ") + f);
    if (M.n == 0) continue;
    // 1. the pricing kernel alone, against a straightforward reference, on pseudo-random vectors
    std::shared_ptr<PricingBackend> dev = makeEmulatedPricer(M);
    std::vector<double> cost(M.n + M.m), y(M.m), d(M.n + M.m), ref(M.n + M.m);
    unsigned long long st = 12345;
    auto rnd = [&]() { st = st * 6364136223846793005ULL + 1442695040888963407ULL; return ((st >> 33) % 2001) / 100.0 - 10.0; };
    for (double& v : cost) v = rnd();
    for (double& v : y) v = rnd();
    dev->computeD(cost.data(), y.data(), d.data());
    for (int j = 0; j < M.n; j++) { double sum = 0; for (int p = M.cs[j]; p < M.cs[j + 1]; p++) sum += M.av[p] * y[M.ri[p]]; ref[j] = cost[j] - sum; }
    for (int i = 0; i < M.m; i++) ref[M.n + i] = cost[M.n + i] + y[i];
    bool same = true; for (int j = 0; j < M.n + M.m; j++) if (d[j] != ref[j]) same = false;
    CHECK(same, std::string("emulated device pricing is bit-identical to the reference on ") + f);
    // 2. a full solve through the offload path gives the same answer, step for step
    Options a = quiet(), b = quiet(); b.gpu = "emulate";
    Solution P = solve(M, a), G = solve(M, b);
    CHECK(G.status == P.status && G.iters == P.iters && relClose(G.obj, P.obj, 1e-12),
          std::string("solve with the emulated device matches the CPU on ") + f);
    CHECK(G.pricing.find("emulated") != std::string::npos, "emulated backend reported");
  }
  // 3. asking for a GPU that is not there must fall back to the CPU and still solve correctly
  Model M; std::string err; readMps("data/afiro.mps", M, err);
  Options o = quiet(); o.gpu = "on";
  Solution S = solve(M, o);
  CHECK(S.status == Status::Optimal && S.ver.ok, "--gpu on without a device still solves correctly");
  std::string note; Options au = quiet(); au.gpu = "auto";
  CHECK(makePricer(M, au, note) == nullptr && note.find("too small") != std::string::npos, "auto mode keeps small models on the CPU");
}

static void testBadInput() {
  std::cout << "[inline] bad input\n";
  Model M;
  CHECK(!loadString("this is not an MPS file\n", M), "garbage must be rejected");
  Model E; std::string err;
  CHECK(!readMps("tests/does_not_exist.mps", E, err), "missing file must be rejected");
  CHECK(exitCode(Status::InputError) == 10, "input error exit code is 10");
}

static void testExitCodeMap() {
  std::cout << "[exit codes]\n";
  CHECK(exitCode(Status::Optimal) == 0, "optimal -> 0");
  CHECK(exitCode(Status::Feasible) == 1, "feasible -> 1");
  CHECK(exitCode(Status::Infeasible) == 2, "infeasible -> 2");
  CHECK(exitCode(Status::Unbounded) == 3, "unbounded -> 3");
  CHECK(exitCode(Status::Unverified) == 4, "unverified -> 4");
  CHECK(exitCode(Status::LimitNoSol) == 5, "limit no solution -> 5");
  CHECK(exitCode(Status::InputError) == 10, "input error -> 10");
  CHECK(exitCode(Status::Unsupported) == 11, "unsupported -> 11");
}

static void testTamper() {
  std::cout << "[verify] tampered solutions must fail\n";
  Model M; std::string err;
  CHECK(readMps("data/afiro.mps", M, err), "read afiro");
  Options o = quiet();
  Solution S = solve(M, o);
  CHECK(S.ver.ok, "afiro baseline verifies");
  Model Mm = M; Mm.toMinimize();
  CHECK(verify(Mm, S.x, S.y, S.d, true, o).ok, "untouched solution verifies independently");

  std::vector<double> bad = S.x;
  int j = 0;
  for (; j < Mm.n; j++) if (!isInf(Mm.lb[j])) break;
  if (j < Mm.n) bad[j] = Mm.lb[j] - 1.0;
  Verification v = verify(Mm, bad, S.y, S.d, true, o);
  CHECK(!v.ok, "bound-violating solution must fail verification");

  Model F; CHECK(readMps("data/flugpl.mps", F, err), "read flugpl");
  Solution T = solve(F, o);
  CHECK(T.status == Status::Optimal && T.ver.ok, "flugpl solves and verifies");
  Model Fm = F; Fm.toMinimize();
  std::vector<double> frac = T.x;
  for (int k = 0; k < Fm.n; k++) if (Fm.isint[k]) { frac[k] += 0.5; break; }
  Verification w = verify(Fm, frac, std::vector<double>(), std::vector<double>(), false, o);
  CHECK(!w.ok, "fractional integer variable must fail verification");
}

int main() {
  testFile("data/afiro.mps", -464.753143);
  testFile("data/avgas.mps", -7.75);
  testFile("data/chip.mps", -900.0);
  testFile("data/adlittle.mps", 225494.963162);
  testFile("data/e226.mps", -11.638929);
  testFile("data/israel.mps", -896644.821810);
  testFile("demo/crudeblend.mps", 587360.0);
  testFile("data/flugpl.mps", 1201500.0);   // MILP
  if (std::getenv("BS_SLOW")) testFile("data/25fv47.mps", 5501.845888);   // ~15 s; run with: make test-slow
  testStatuses();
  testHonestLimits();
  testRobustness();
  testDualMethod();
  testGpuOffloadPath();
  testBadInput();
  testExitCodeMap();
  testTamper();
  std::cout << "\n" << g_pass << " checks passed, " << g_fail << " failed\n";
  return g_fail == 0 ? 0 : 1;
}
