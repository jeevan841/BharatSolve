// BharatSolve core declarations. C++17, standard library only (NFR-01).
#pragma once
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace bs {

constexpr double INF = 1e30;
inline bool isInf(double v) { return std::fabs(v) >= 1e20; }

// ---------------------------------------------------------------- model
struct Model {
  std::string name;
  int m = 0, n = 0;
  std::vector<int> cs, ri;          // CSC: column starts (n+1), row indices
  std::vector<double> av;           // CSC values
  std::vector<double> c, lb, ub;    // objective, column bounds
  std::vector<double> rl, ru;       // row bounds  rl <= Ax <= ru
  std::vector<char> isint;
  std::vector<std::string> rname, cname;
  double objoff = 0;
  bool maximize = false;
  // convex quadratic term 0.5 x'Qx (full symmetric triplets), parsed only
  std::vector<int> qi, qj;
  std::vector<double> qv;
  int nnz() const { return (int)ri.size(); }
  bool isMip() const { for (char v : isint) if (v) return true; return false; }
  bool hasQuad() const { return !qv.empty(); }
  void toMinimize();                // negate c if maximize (keeps flag)
};

enum class Status { Optimal, Feasible, Infeasible, Unbounded, Unverified, LimitNoSol, InputError, Unsupported };
const char* statusName(Status s);
int exitCode(Status s);

struct Options {
  double timeLimit = 60;
  double gap = 1e-4;              // relative MIP gap
  long nodeLimit = 10000000;
  std::string method = "auto";    // auto | dual | primal | ipm
  bool presolve = true;
  bool scale = true;
  int seed = 1;
  int threads = 1;                // accepted; single-threaded in this release
  int log = 1;                    // 0 silent, 1 summary, 2 progress, 3 iterations
  double feasTol = 1e-7, optTol = 1e-7, pivTol = 1e-9, intTol = 1e-6;
  double verifyTol = 1e-6;
  long maxIter = 0;               // 0 = automatic
  bool branchPseudo = true;
};

struct Verification {
  bool ran = false, ok = false;
  double primalInf = 0, intInf = 0, dualInf = 0, compl_ = 0, gapRel = 0;
  bool hasDual = false;
  std::vector<std::string> failures;
};

struct Solution {
  Status status = Status::InputError;
  std::string message;
  double obj = 0;                  // in original sense (max or min)
  double bound = 0;                // MIP best bound (original sense)
  double gap = 0;
  std::vector<double> x, y, d;     // primal, row duals, reduced costs (original sense)
  std::vector<double> ray;         // primal (unbounded) or dual (infeasible) ray, if any
  std::string rayKind;
  long iters = 0, nodes = 0;
  double seconds = 0, firstIncumbentSec = -1;
  // presolve / scaling statistics
  int preRowsRemoved = 0, preColsRemoved = 0;
  double scaleMin = 1, scaleMax = 1;
  double condEstimate = 0;
  int numericalEvents = 0;
  Verification ver;
  std::string method;
};

// ---------------------------------------------------------------- io
bool readMps(const std::string& path, Model& out, std::string& err);
bool readMpsStream(std::istream& in, Model& out, std::string& err);

// ---------------------------------------------------------------- lp engine
enum class LPStatus { Optimal, Infeasible, Unbounded, Limit, DualInfeasible, Numerical };

struct Basis {
  std::vector<int> basic;   // m entries
  std::vector<char> stat;   // n+m entries: 0 basic 1 lower 2 upper 3 free
  bool valid() const { return !basic.empty(); }
};

class Simplex {
 public:
  Simplex(const Model& mdl, const Options& o);
  void setColBounds(int j, double l, double u);      // structural j
  void resetBounds();
  void setBasis(const Basis& b);
  Basis getBasis() const;
  void setCosts(const std::vector<double>& c);
  LPStatus solve(const std::string& method, double deadline);
  double objective() const;
  const std::vector<double>& primal() const { return x_; }   // size n+m
  void duals(std::vector<double>& y, std::vector<double>& d) const;
  const std::vector<double>& ray() const { return ray_; }
  const std::string& rayKind() const { return rayKind_; }
  long iterations() const { return iters_; }
  int numericalEvents() const { return numEvents_; }
  double lbOf(int j) const { return lb_[j]; }
  double ubOf(int j) const { return ub_[j]; }
  int n() const { return n_; }
  int m() const { return m_; }
  std::string lastNote;

 private:
  int m_, n_, N_;
  const Model& M_;
  Options opt_;
  std::vector<double> lb_, ub_, lb0_, ub0_, cost_, cost0_, x_, d_;
  std::vector<int> basic_, pos_;
  std::vector<char> stat_;
  unsigned rng_ = 1;
  // dense LU (basis matrix) + product-form eta updates between refactorizations
  std::vector<double> lu_;
  std::vector<int> piv_;
  struct Eta { int r; std::vector<double> a; };
  std::vector<Eta> etas_;
  long iters_ = 0;
  int numEvents_ = 0;
  std::vector<double> ray_;
  std::string rayKind_;
  bool perturbed_ = false;

  double rnd();
  void colOf(int j, std::vector<double>& out) const;
  void placeNonbasic(int j);
  void slackBasis();
  bool factor();
  void ftran(std::vector<double>& v) const;
  void btran(std::vector<double>& v) const;
  double dotCol(int j, const std::vector<double>& y) const;
  void computePrimal();
  void computeDuals(const std::vector<double>& cst);
  void perturbCosts();
  void restoreCosts();
  LPStatus primalPhase(double deadline, bool phase1);
  LPStatus dual(double deadline);
  void finalize();
};

// ---------------------------------------------------------------- presolve/scale
struct PresolveResult {
  Model reduced;
  bool infeasible = false, unbounded = false;
  std::string message;
  std::vector<int> colMap, rowMap;       // reduced -> original
  std::vector<double> fixedVal;          // original col value if removed
  std::vector<char> colRemoved, rowRemoved;
  struct Singleton { int row, col; double a, oldLB, oldUB, newLB, newUB; };
  std::vector<Singleton> singles;        // in elimination order; undo in reverse
  int rowsRemoved = 0, colsRemoved = 0;
};
PresolveResult presolve(const Model& orig, const Options& o);
void postsolve(const Model& orig, const PresolveResult& pr, const std::vector<double>& xr,
               const std::vector<double>& yr, const std::vector<double>& dr, std::vector<double>& x,
               std::vector<double>& y, std::vector<double>& d);

struct Scaling { std::vector<double> R, C; double minS = 1, maxS = 1; };
Scaling computeScaling(const Model& m);
Model applyScaling(const Model& m, const Scaling& s);

// ---------------------------------------------------------------- verification
Verification verify(const Model& minModel, const std::vector<double>& x, const std::vector<double>& y,
                    const std::vector<double>& d, bool lpDuals, const Options& o);

// ---------------------------------------------------------------- drivers
Solution solve(const Model& model, const Options& o);
Solution solveLPCore(const Model& minModel, const Options& o, double deadline, bool wantMip);
void writeJson(const Model& model, const Solution& s, const std::string& path);
double now();

}  // namespace bs
