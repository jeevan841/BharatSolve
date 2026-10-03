#include <cstdlib>
#include <iostream>
#include <memory>
#include "bs.hpp"
#include "gpu.hpp"
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
  std::string pricingNote;
  std::shared_ptr<PricingBackend> pricer = makePricer(scaled, o, pricingNote);
  sx.setPricer(pricer);
  S.pricing = pricingNote;
  if (o.log >= 2) std::cerr << "pricing backend: " << pricingNote << "\n";
  LPStatus st = sx.solve(o.method, deadline);
  S.iters = sx.iterations(); S.numericalEvents = sx.numericalEvents(); S.pricingSeconds = sx.pricingSeconds();
  if (st == LPStatus::Unbounded) { S.status = Status::Unbounded; return S; }
  if (st == LPStatus::Infeasible) { S.status = Status::Infeasible; return S; }
  if (st == LPStatus::Limit) {
    S.status = Status::LimitNoSol;
    S.message = (now() > deadline) ? "time limit reached before an optimal basis was found"
                                   : "iteration limit reached before an optimal basis was found";
    return S;
  }
  if (st == LPStatus::Numerical) {
    S.status = Status::LimitNoSol;
    S.message = "numerical failure: the basis became singular, so no reliable solution was produced";
    return S;
  }

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
struct Node { std::vector<double> lb, ub; double bound; int depth; std::shared_ptr<Basis> basis;
            int bvar = -1; int bdir = 0; double bfrac = 0, pobj = 0; };  // basis: parent's optimal basis (warm start)

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
    long nodes = 0, warmNodes = 0, coldRetries = 0; bool anyFeasible = false;
    std::string bbNote;
    std::shared_ptr<PricingBackend> bbPricer = makePricer(scaled, o, bbNote);   // matrix uploaded once, shared by every node
    S.pricing = bbNote;
    // Pseudocosts: average objective increase per unit of rounding, learned from solved children.
    std::vector<double> pcSum[2] = {std::vector<double>(R.n, 0.0), std::vector<double>(R.n, 0.0)};
    std::vector<long> pcCnt[2] = {std::vector<long>(R.n, 0), std::vector<long>(R.n, 0)};
    double pcTotSum[2] = {0, 0}; long pcTotCnt[2] = {0, 0};
    bool incomplete = false;   // true if any part of the search tree was skipped, so optimality cannot be claimed
    std::vector<Node> stack;
    Node root{R.lb, R.ub, -INF, 0, nullptr, -1, 0, 0, 0}; stack.push_back(root);
    while (!stack.empty()) {
      if (now() > deadline) { incomplete = true; break; }
      nodes++;
      Node cur = stack.back(); stack.pop_back();
      // Each node re-solves the LP relaxation. The root is solved from scratch with the primal
      // simplex. Children start from their parent's optimal basis: tightening one bound keeps that
      // basis dual feasible, so the dual simplex repairs the single violated variable in a few
      // pivots. A fresh Simplex object is used per node, and the dual start snaps every nonbasic
      // variable to its (possibly new) bound, so no stale values can leak between nodes.
      // If a warm start fails numerically, the node is re-solved cold with the primal simplex.
      double nodeDeadline = std::min(deadline, now() + std::max(1.0, o.timeLimit / 20));
      auto makeSimplex = [&]() {
        std::unique_ptr<Simplex> p(new Simplex(scaled, o));
        p->setPricer(bbPricer);
        for (int j = 0; j < R.n; j++) { double l = cur.lb[j], u = cur.ub[j]; double c = sc.C.empty()?1.0:sc.C[j]; p->setColBounds(j, isInf(l)?l:l/c, isInf(u)?u:u/c); }
        return p;
      };
      std::unique_ptr<Simplex> sxp = makeSimplex();
      LPStatus st;
      if (cur.basis) {
        sxp->setBasis(*cur.basis);
        st = sxp->solve("dual", nodeDeadline);
        S.iters += sxp->iterations(); warmNodes++;
        if (st != LPStatus::Optimal && st != LPStatus::Infeasible) {   // warm start failed: cold retry
          sxp = makeSimplex();
          st = sxp->solve("primal", nodeDeadline);
          S.iters += sxp->iterations(); coldRetries++;
        }
      } else {
        st = sxp->solve("primal", nodeDeadline);
        S.iters += sxp->iterations();
      }
      Simplex& sx = *sxp;
      if (st == LPStatus::Infeasible) continue;                       // node proven infeasible: safe to prune
      if (st == LPStatus::Unbounded) {
        if (nodes == 1) { S.status = Status::Unbounded; S.nodes = nodes; S.message = "LP relaxation is unbounded, so the integer problem is unbounded or has no integer-feasible point"; S.seconds = now() - t0; if (modelIn.maximize) { S.obj = -S.obj; S.bound = -S.bound; } return S; }
        incomplete = true; continue;                                  // cannot happen for a bounded root; treat as untrusted
      }
      if (st != LPStatus::Optimal) { incomplete = true; continue; }   // limit or numerical failure: subtree NOT explored
      double relObj = sx.objective();
      if (relObj > bestObj + 1e-9) continue;  // pruned by bound
      if (cur.bvar >= 0 && cur.bfrac > 1e-9) {                        // learn from this child's objective gain
        int d = cur.bdir > 0 ? 1 : 0;
        double g = std::max(0.0, relObj - cur.pobj) / cur.bfrac;
        pcSum[d][cur.bvar] += g; pcCnt[d][cur.bvar]++; pcTotSum[d] += g; pcTotCnt[d]++;
      }
      const auto& xf = sx.primal();
      int fracJ = -1; double bestScore = -1;
      double avg[2] = {pcTotCnt[0] ? pcTotSum[0] / pcTotCnt[0] : 1.0, pcTotCnt[1] ? pcTotSum[1] / pcTotCnt[1] : 1.0};
      for (int j = 0; j < R.n; j++) if (isint[j]) {
        double v = xf[j] * (sc.C.empty()?1.0:sc.C[j]);
        double f = v - std::floor(v);
        double dist = std::min(f, 1 - f);
        if (dist <= 1e-6) continue;
        double pd = pcCnt[0][j] ? pcSum[0][j] / pcCnt[0][j] : avg[0];
        double pu = pcCnt[1][j] ? pcSum[1][j] / pcCnt[1][j] : avg[1];
        double score = std::max(pd * f, 1e-6) * std::max(pu * (1 - f), 1e-6);   // product rule
        if (score > bestScore) { bestScore = score; fracJ = j; }
      }
      if (fracJ < 0) {
        anyFeasible = true;
        if (relObj < bestObj) { bestObj = relObj; bestX.assign(xf.begin(), xf.begin() + R.n); }
        continue;
      }
      double v = xf[fracJ] * (sc.C.empty()?1.0:sc.C[fracJ]);
      Node up = cur, down = cur;
      down.ub[fracJ] = std::min(down.ub[fracJ], std::floor(v));
      up.lb[fracJ] = std::max(up.lb[fracJ], std::ceil(v));
      double fdn = v - std::floor(v);
      down.bvar = up.bvar = fracJ; down.bdir = -1; up.bdir = 1; down.bfrac = fdn; up.bfrac = 1 - fdn;
      down.pobj = up.pobj = relObj;
      std::shared_ptr<Basis> bp = std::make_shared<Basis>(sx.getBasis());
      up.basis = bp; down.basis = bp;
      if (fdn < 0.5) { stack.push_back(up); stack.push_back(down); } else { stack.push_back(down); stack.push_back(up); }
      if (nodes > o.nodeLimit) { incomplete = true; break; }
    }
    S.nodes = nodes;
    if (!anyFeasible) {
      S.status = incomplete ? Status::LimitNoSol : Status::Infeasible;
      S.message = incomplete ? "search stopped (time, node limit or numerical trouble) before any integer-feasible solution was found"
                             : "no integer-feasible solution exists";
      S.seconds = now() - t0;
      return S;
    }
    std::vector<double> xr(bestX), yr(R.m, 0.0), dr(R.n, 0.0);
    unscaleSolution(sc, xr, yr, dr);
    postsolve(M, pr, xr, yr, dr, S.x, S.y, S.d);
    double obj = M.objoff; for (int j = 0; j < M.n; j++) obj += M.c[j] * S.x[j];
    S.obj = obj; S.bound = bestObj + M.objoff; S.gap = std::fabs(S.obj - S.bound) / std::max(1.0, std::fabs(S.obj));
    S.status = incomplete ? Status::Feasible : Status::Optimal;
    if (incomplete) S.message = "incumbent found but the search was not completed, so optimality is not proven";
    S.method = "branch-and-bound(dual-simplex warm start)";
    if (getenv("BS_DEBUG")) std::cerr << "[bb] nodes=" << nodes << " warm=" << warmNodes << " coldRetries=" << coldRetries << "\n";
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
  f << "  \"pricing_backend\": \"" << esc(s.pricing) << "\",\n";
  f << "  \"pricing_seconds\": " << s.pricingSeconds << ",\n";
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
