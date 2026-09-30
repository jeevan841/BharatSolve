#include "bs.hpp"
#include <algorithm>
#include <fstream>
#include <queue>
#include <sstream>

namespace bs {

static void unscaleSolution(const Scaling& sc, std::vector<double>& x, std::vector<double>& y, std::vector<double>& d) {
  for (size_t j = 0; j < x.size(); j++) x[j] *= sc.C[j];
  for (size_t i = 0; i < y.size(); i++) y[i] *= sc.R[i];
  for (size_t j = 0; j < d.size(); j++) d[j] /= sc.C[j];
}

// Solve the (already minimized, possibly scaled/presolved) LP relaxation once.
static LPStatus solveOnceLP(const Model& mdl, const Options& o, double deadline, Simplex& sx, std::string method) {
  return sx.solve(method, deadline);
}

Solution solveLPCore(const Model& minModel, const Options& o, double deadline, bool /*wantMip*/) {
  Solution S;
  PresolveResult pr = presolve(minModel, o);
  if (pr.infeasible) { S.status = Status::Infeasible; S.message = pr.message; return S; }
  S.preRowsRemoved = pr.rowsRemoved; S.preColsRemoved = pr.colsRemoved;
  Model& R = pr.reduced;
  Scaling sc;
  Model scaled = R;
  if (o.scale && R.n > 0 && R.m > 0) { sc = computeScaling(R); scaled = applyScaling(R, sc); S.scaleMin = sc.minS; S.scaleMax = sc.maxS; }
  else { sc.R.assign(R.m, 1.0); sc.C.assign(R.n, 1.0); }

  if (R.n == 0) {
    // everything fixed by presolve; trivial solution
    std::vector<double> xr, yr, dr;
    postsolve(minModel, pr, xr, yr, dr, S.x, S.y, S.d);
    S.status = Status::Optimal;
    double obj = minModel.objoff; for (int j = 0; j < minModel.n; j++) obj += minModel.c[j] * S.x[j];
    S.obj = obj; return S;
  }

  Simplex sx(scaled, o);
  LPStatus st = sx.solve(o.method, deadline);
  S.iters = sx.iterations(); S.numericalEvents = sx.numericalEvents();
  if (st == LPStatus::Unbounded) { S.status = Status::Unbounded; return S; }
  if (st == LPStatus::Infeasible) { S.status = Status::Infeasible; return S; }
  if (st == LPStatus::Limit || st == LPStatus::Numerical) { S.status = Status::LimitNoSol; S.message = "iteration or time limit before feasibility"; return S; }

  std::vector<double> xr(R.n), yr(R.m), dr(R.n);
  const auto& xf = sx.primal();
  for (int j = 0; j < R.n; j++) xr[j] = xf[j];
  sx.duals(yr, dr);
  unscaleSolution(sc, xr, yr, dr);
  postsolve(minModel, pr, xr, yr, dr, S.x, S.y, S.d);
  double obj = minModel.objoff; for (int j = 0; j < minModel.n; j++) obj += minModel.c[j] * S.x[j];
  S.obj = obj; S.status = Status::Optimal;
  return S;
}

// ---------------- Branch and bound (depth-first, most-fractional branching) ----------------
struct Node { std::vector<double> lb, ub; double bound; int depth; };

Solution solve(const Model& modelIn, const Options& o) {
  double t0 = now(); double deadline = t0 + o.timeLimit;
  Model M = modelIn; M.toMinimize();
  Solution S;
  if (M.n == 0) { S.status = Status::InputError; S.message = "empty model"; return S; }

  if (!M.isMip() && !M.hasQuad()) {
    S = solveLPCore(M, o, deadline, false);
    S.method = (o.method == "dual") ? "dual-simplex" : "primal-simplex(phase1+phase2)";
  } else if (M.hasQuad()) {
    S.status = Status::Unsupported;
    S.message = "QP solve (interior-point) is not implemented in this release; see roadmap. Model was parsed and its convexity/structure can be reported.";
    return S;
  } else {
    // MILP: depth-first branch and bound on top of the LP relaxation.
    Options relO = o; relO.presolve = false;  // node LPs reuse the (light) presolved root only
    PresolveResult pr = presolve(M, o);
    if (pr.infeasible) { S.status = Status::Infeasible; S.message = pr.message; return S; }
    Model R = pr.reduced;
    Scaling sc; Model scaled = R;
    if (o.scale && R.n > 0) { sc = computeScaling(R); scaled = applyScaling(R, sc); } else { sc.R.assign(R.m,1.0); sc.C.assign(R.n,1.0); }
    std::vector<char> isint = R.isint;

    double bestObj = INF; std::vector<double> bestX;
    long nodes = 0; bool anyFeasible = false;
    std::vector<Node> stack;
    Node root{R.lb, R.ub, -INF, 0}; stack.push_back(root);
    while (!stack.empty()) {
      if (now() > deadline) break;
      nodes++;
      Node cur = stack.back(); stack.pop_back();
      // Fresh Simplex per node: warm-starting the same object across nodes is
      // unsafe here because tightening a bound via setColBounds does not snap
      // an already-nonbasic variable's stored value to the new bound, and
      // Phase 1 only repairs infeasible BASIC variables, not a stale nonbasic
      // one — that combination produced silent false "infeasible" nodes.
      // Rebuilding from a clean slack basis each node is correct by
      // construction; for the model sizes in this release the extra
      // factorizations are inexpensive.
      Simplex sx(scaled, o);
      for (int j = 0; j < R.n; j++) { double l = cur.lb[j], u = cur.ub[j]; double c = sc.C.empty()?1.0:sc.C[j]; sx.setColBounds(j, isInf(l)?l:l/c, isInf(u)?u:u/c); }
      LPStatus st = sx.solve("primal", std::min(deadline, now() + std::max(1.0, o.timeLimit / 20)));
      if (st != LPStatus::Optimal) continue;
      double relObj = sx.objective();
      if (relObj > bestObj + 1e-9) continue;  // pruned by bound
      const auto& xf = sx.primal();
      int fracJ = -1; double fracDist = 1e-6;
      for (int j = 0; j < R.n; j++) if (isint[j]) {
        double v = xf[j] * (sc.C.empty()?1.0:sc.C[j]);
        double f = v - std::floor(v);
        double dist = std::min(f, 1 - f);
        if (dist > fracDist) { fracDist = dist; fracJ = j; }
      }
      if (fracJ < 0) {
        anyFeasible = true;
        if (relObj < bestObj) { bestObj = relObj; bestX.assign(xf.begin(), xf.begin() + R.n); S.iters += sx.iterations(); }
        continue;
      }
      double v = xf[fracJ] * (sc.C.empty()?1.0:sc.C[fracJ]);
      Node up = cur, down = cur;
      down.ub[fracJ] = std::min(down.ub[fracJ], std::floor(v));
      up.lb[fracJ] = std::max(up.lb[fracJ], std::ceil(v));
      stack.push_back(up); stack.push_back(down);
      if (nodes > o.nodeLimit) break;
    }
    S.nodes = nodes;
    if (!anyFeasible) { S.status = (now() > deadline) ? Status::LimitNoSol : Status::Infeasible; S.message = "no integer-feasible solution found"; return S; }
    std::vector<double> xr(bestX), yr(R.m, 0.0), dr(R.n, 0.0);
    unscaleSolution(sc, xr, yr, dr);
    postsolve(M, pr, xr, yr, dr, S.x, S.y, S.d);
    double obj = M.objoff; for (int j = 0; j < M.n; j++) obj += M.c[j] * S.x[j];
    S.obj = obj; S.bound = bestObj + M.objoff; S.gap = std::fabs(S.obj - S.bound) / std::max(1.0, std::fabs(S.obj));
    S.status = (now() > deadline) ? Status::Feasible : Status::Optimal;
    S.method = "branch-and-bound(primal-simplex)";
  }

  S.seconds = now() - t0;
  if (S.status == Status::Optimal || S.status == Status::Feasible) {
    bool lpDuals = !modelIn.isMip();
    S.ver = verify(M, S.x, S.y, S.d, lpDuals, o);
    if (!S.ver.ok && S.status == Status::Optimal) S.status = Status::Unverified;
  }
  if (modelIn.maximize) { S.obj = -S.obj; S.bound = -S.bound; }
  return S;
}

static std::string esc(const std::string& s) { std::string o; for (char c : s) { if (c=='"'||c=='\\') o+='\\'; o+=c; } return o; }

void writeJson(const Model& model, const Solution& s, const std::string& path) {
  std::ofstream f(path);
  f << "{\n";
  f << "  \"solver\": \"BharatSolve\",\n";
  f << "  \"model\": \"" << esc(model.name) << "\",\n";
  f << "  \"status\": \"" << statusName(s.status) << "\",\n";
  f << "  \"exit_code\": " << exitCode(s.status) << ",\n";
  f << "  \"method\": \"" << esc(s.method) << "\",\n";
  f << "  \"objective\": " << s.obj << ",\n";
  f << "  \"bound\": " << s.bound << ",\n";
  f << "  \"gap\": " << s.gap << ",\n";
  f << "  \"iterations\": " << s.iters << ",\n";
  f << "  \"nodes\": " << s.nodes << ",\n";
  f << "  \"seconds\": " << s.seconds << ",\n";
  f << "  \"presolve\": {\"rows_removed\": " << s.preRowsRemoved << ", \"cols_removed\": " << s.preColsRemoved << "},\n";
  f << "  \"scaling\": {\"min\": " << s.scaleMin << ", \"max\": " << s.scaleMax << "},\n";
  f << "  \"numerical_events\": " << s.numericalEvents << ",\n";
  f << "  \"verification\": {\"ran\": " << (s.ver.ran ? "true" : "false") << ", \"ok\": " << (s.ver.ok ? "true" : "false")
    << ", \"primal_infeasibility\": " << s.ver.primalInf << ", \"integer_infeasibility\": " << s.ver.intInf
    << ", \"dual_infeasibility\": " << s.ver.dualInf << "},\n";
  f << "  \"variables\": [";
  for (int j = 0; j < (int)s.x.size(); j++) { if (j) f << ", "; f << "{\"name\": \"" << esc(model.cname[j]) << "\", \"value\": " << s.x[j] << "}"; }
  f << "],\n";
  f << "  \"message\": \"" << esc(s.message) << "\"\n";
  f << "}\n";
}

}  // namespace bs
