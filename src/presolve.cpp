// Presolve (first release, "Must" reductions only): remove fixed columns,
// drop empty rows, and eliminate singleton rows by tightening the column
// bound they imply. Every reduction records what postsolve needs to undo it.
// Roadmap: dominated columns, duplicate rows, forcing rows with multiple vars.
#include "bs.hpp"
#include <cmath>

namespace bs {

PresolveResult presolve(const Model& orig, const Options& o) {
  PresolveResult pr;
  int n = orig.n, m = orig.m;
  std::vector<char> colDead(n, 0), rowDead(m, 0);
  std::vector<double> lb = orig.lb, ub = orig.ub, rl = orig.rl, ru = orig.ru;
  std::vector<double> fixedVal(n, 0.0);

  if (!o.presolve) { pr.reduced = orig; pr.colMap.resize(n); pr.rowMap.resize(m); for (int j=0;j<n;j++) pr.colMap[j]=j; for (int i=0;i<m;i++) pr.rowMap[i]=i; pr.colRemoved.assign(n,0); pr.rowRemoved.assign(m,0); return pr; }

  // Pass 1: fixed columns (lb==ub) -> remove, folding value into row bounds & objective offset.
  double objAdj = 0;
  bool changed = true;
  int rounds = 0;
  // column-wise adjacency built once (structure doesn't change; we only mask rows/cols)
  while (changed && rounds < 20) {
    changed = false; rounds++;
    for (int j = 0; j < n; j++) {
      if (colDead[j]) continue;
      if (lb[j] == ub[j]) {
        double v = lb[j];
        fixedVal[j] = v;
        colDead[j] = 1;
        objAdj += orig.c[j] * v;
        for (int p = orig.cs[j]; p < orig.cs[j + 1]; p++) {
          int i = orig.ri[p];
          if (rowDead[i]) continue;
          double a = orig.av[p];
          if (!isInf(rl[i])) rl[i] -= a * v;
          if (!isInf(ru[i])) ru[i] -= a * v;
        }
        changed = true;
      }
    }
    // Pass 2: singleton rows -> tighten the one remaining column's bound, then drop row.
    std::vector<int> singCount(m, 0), singCol(m, -1); std::vector<double> singA(m, 0.0);
    for (int j = 0; j < n; j++) {
      if (colDead[j]) continue;
      for (int p = orig.cs[j]; p < orig.cs[j + 1]; p++) {
        int i = orig.ri[p];
        if (rowDead[i]) continue;
        singCount[i]++; singCol[i] = j; singA[i] = orig.av[p];
      }
    }
    for (int i = 0; i < m; i++) {
      if (rowDead[i]) continue;
      if (singCount[i] == 0) {
        // empty row: check feasibility of 0 in [rl,ru], then drop.
        if ((!isInf(rl[i]) && rl[i] > o.feasTol) || (!isInf(ru[i]) && ru[i] < -o.feasTol)) { pr.infeasible = true; pr.message = "empty row infeasible after fixing: " + orig.rname[i]; return pr; }
        rowDead[i] = 1; pr.rowsRemoved++; changed = true;
      } else if (singCount[i] == 1) {
        // Singleton-row bound tightening is DISABLED in this release: its
        // postsolve dual reconstruction (recovering y_i for the eliminated
        // row) has an unresolved correctness issue that showed up as a
        // dual-infeasibility on ADLITTLE during verification. Re-enable once
        // that recovery is fixed and re-verified; for now these rows are
        // simply left in the model (correct, just slightly less reduction).
        continue;
      }
    }
  }
  for (int j = 0; j < n; j++) if (lb[j] > ub[j] + 1e-7) { pr.infeasible = true; pr.message = "column " + orig.cname[j] + " has lb>ub"; return pr; }

  pr.colMap.clear(); pr.rowMap.clear();
  std::vector<int> newColOf(n, -1), newRowOf(m, -1);
  Model& R = pr.reduced;
  for (int j = 0; j < n; j++) if (!colDead[j]) { newColOf[j] = (int)pr.colMap.size(); pr.colMap.push_back(j); }
  for (int i = 0; i < m; i++) if (!rowDead[i]) { newRowOf[i] = (int)pr.rowMap.size(); pr.rowMap.push_back(i); }
  R.n = (int)pr.colMap.size(); R.m = (int)pr.rowMap.size();
  R.name = orig.name; R.maximize = orig.maximize; R.objoff = orig.objoff + objAdj;
  R.c.resize(R.n); R.lb.resize(R.n); R.ub.resize(R.n); R.isint.resize(R.n); R.cname.resize(R.n);
  for (int j = 0; j < n; j++) if (!colDead[j]) { int nj = newColOf[j]; R.c[nj] = orig.c[j]; R.lb[nj] = lb[j]; R.ub[nj] = ub[j]; R.isint[nj] = orig.isint[j]; R.cname[nj] = orig.cname[j]; }
  R.rl.resize(R.m); R.ru.resize(R.m); R.rname.resize(R.m);
  for (int i = 0; i < m; i++) if (!rowDead[i]) { int ni = newRowOf[i]; R.rl[ni] = rl[i]; R.ru[ni] = ru[i]; R.rname[ni] = orig.rname[i]; }
  R.cs.assign(R.n + 1, 0);
  std::vector<std::vector<std::pair<int, double>>> cols(R.n);
  for (int j = 0; j < n; j++) {
    if (colDead[j]) continue;
    int nj = newColOf[j];
    for (int p = orig.cs[j]; p < orig.cs[j + 1]; p++) { int i = orig.ri[p]; if (rowDead[i]) continue; cols[nj].emplace_back(newRowOf[i], orig.av[p]); }
  }
  for (int j = 0; j < R.n; j++) R.cs[j + 1] = R.cs[j] + (int)cols[j].size();
  R.ri.assign(R.cs[R.n], 0); R.av.assign(R.cs[R.n], 0.0);
  for (int j = 0; j < R.n; j++) { int p = R.cs[j]; for (auto& e : cols[j]) { R.ri[p] = e.first; R.av[p] = e.second; p++; } }

  pr.colRemoved = colDead; pr.rowRemoved = rowDead; pr.fixedVal = fixedVal;
  pr.colsRemoved = (int)(n - R.n);
  return pr;
}

void postsolve(const Model& orig, const PresolveResult& pr, const std::vector<double>& xr,
               const std::vector<double>& yr, const std::vector<double>& dr, std::vector<double>& x,
               std::vector<double>& y, std::vector<double>& d) {
  int n = orig.n, m = orig.m;
  x.assign(n, 0.0); y.assign(m, 0.0); d.assign(n, 0.0);
  for (int j = 0; j < n; j++) if (pr.colRemoved.size() == (size_t)n && pr.colRemoved[j]) x[j] = pr.fixedVal[j];
  for (size_t k = 0; k < pr.colMap.size(); k++) x[pr.colMap[k]] = xr[k];
  for (size_t k = 0; k < pr.rowMap.size(); k++) y[pr.rowMap[k]] = yr[k];
  if (!dr.empty()) for (size_t k = 0; k < pr.colMap.size(); k++) d[pr.colMap[k]] = dr[k];
  if (!yr.empty()) {
    // Column reduced costs in the reduced problem's own numbering, needed to
    // attribute duals when a singleton row's implied bound is the active one.
    // Recomputed here rather than threaded through, since it is cheap: for
    // each singleton, if x[col] sits at the bound THIS row introduced (and
    // that bound differs from what existed before it), the row absorbs the
    // reduced cost that would otherwise sit unexplained on the column.
    for (auto it = pr.singles.rbegin(); it != pr.singles.rend(); ++it) {
      int j = it->col, i = it->row;
      bool tightenedLower = it->newLB > it->oldLB + 1e-12;
      bool tightenedUpper = it->newUB < it->oldUB - 1e-12;
      bool atLower = tightenedLower && std::fabs(x[j] - it->newLB) < 1e-6;
      bool atUpper = tightenedUpper && std::fabs(x[j] - it->newUB) < 1e-6;
      if (atLower || atUpper) {
        double rc = d[j];      // reduced cost currently attributed to column j
        y[i] = rc / it->a;
        d[j] = 0.0;
      } else {
        y[i] = 0.0;
      }
    }
  }
}

Scaling computeScaling(const Model& m) {
  Scaling s; s.R.assign(m.m, 1.0); s.C.assign(m.n, 1.0);
  for (int pass = 0; pass < 20; pass++) {
    std::vector<double> rmin(m.m, INF), rmax(m.m, 0), cmin(m.n, INF), cmax(m.n, 0);
    for (int j = 0; j < m.n; j++) for (int p = m.cs[j]; p < m.cs[j + 1]; p++) {
      double v = std::fabs(m.av[p]) * s.R[m.ri[p]] * s.C[j];
      if (v <= 0) continue;
      rmin[m.ri[p]] = std::min(rmin[m.ri[p]], v); rmax[m.ri[p]] = std::max(rmax[m.ri[p]], v);
      cmin[j] = std::min(cmin[j], v); cmax[j] = std::max(cmax[j], v);
    }
    double worst = 1.0;
    for (int i = 0; i < m.m; i++) if (rmax[i] > 0) { double g = std::sqrt(rmin[i] * rmax[i]); s.R[i] /= g; worst = std::max(worst, rmax[i] / std::max(rmin[i], 1e-30)); }
    for (int j = 0; j < m.n; j++) if (cmax[j] > 0) { double g = std::sqrt(cmin[j] * cmax[j]); s.C[j] /= g; }
    if (worst < 1e3) break;
  }
  s.minS = INF; s.maxS = 0;
  for (double v : s.R) { s.minS = std::min(s.minS, v); s.maxS = std::max(s.maxS, v); }
  for (double v : s.C) { s.minS = std::min(s.minS, v); s.maxS = std::max(s.maxS, v); }
  return s;
}

Model applyScaling(const Model& m, const Scaling& s) {
  Model r = m;
  for (int j = 0; j < r.n; j++) for (int p = r.cs[j]; p < r.cs[j + 1]; p++) r.av[p] *= s.R[r.ri[p]] * s.C[j];
  for (int j = 0; j < r.n; j++) {
    double c = s.C[j];
    r.c[j] *= c;
    if (!isInf(r.lb[j])) r.lb[j] /= c;
    if (!isInf(r.ub[j])) r.ub[j] /= c;
  }
  for (int i = 0; i < r.m; i++) {
    double rr = s.R[i];
    if (!isInf(r.rl[i])) r.rl[i] *= rr;
    if (!isInf(r.ru[i])) r.ru[i] *= rr;
  }
  return r;
}

}  // namespace bs
