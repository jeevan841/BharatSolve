#include <climits>
#include <algorithm>
// Bounded revised simplex over the model's own m rows and n structural columns,
// using slack variables indexed n..n+m-1 for row i (slack bounds -ru<=s<=-rl-ish
// handled via the general form l <= Ax <= u through slack = Ax with bounds [rl,ru]).
// First release uses dense LU refactorization each `refactorEvery` iterations
// (VR/SX roadmap: replace with sparse LU + Forrest-Tomlin updates in a later pass).
#include "bs.hpp"
#include "gpu.hpp"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <iostream>

namespace bs {

double now() {
  using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count();
}

Simplex::Simplex(const Model& mdl, const Options& o) : M_(mdl), opt_(o) {
  n_ = mdl.n; m_ = mdl.m; N_ = n_ + m_;
  lb_.assign(N_, 0); ub_.assign(N_, 0);
  for (int j = 0; j < n_; j++) { lb_[j] = mdl.lb[j]; ub_[j] = mdl.ub[j]; }
  for (int i = 0; i < m_; i++) { lb_[n_ + i] = mdl.rl[i]; ub_[n_ + i] = mdl.ru[i]; }
  lb0_ = lb_; ub0_ = ub_;
  cost_.assign(N_, 0);
  for (int j = 0; j < n_; j++) cost_[j] = mdl.c[j];
  cost0_ = cost_;
  x_.assign(N_, 0); d_.assign(N_, 0);
  basic_.assign(m_, 0); pos_.assign(N_, -1);
  stat_.assign(N_, 1);
  rng_ = 12345u + (unsigned)o.seed;
}

double Simplex::rnd() { rng_ = rng_ * 1103515245u + 12345u; return ((rng_ >> 8) & 0xFFFFFF) / double(0xFFFFFF); }

void Simplex::setColBounds(int j, double l, double u) { lb_[j] = l; ub_[j] = u; }
void Simplex::resetBounds() { lb_ = lb0_; ub_ = ub0_; }
void Simplex::setCosts(const std::vector<double>& c) { for (int j = 0; j < n_; j++) cost_[j] = c[j]; cost0_ = cost_; }

void Simplex::colOf(int j, std::vector<double>& out) const {
  out.assign(m_, 0.0);
  if (j < n_) {
    for (int p = M_.cs[j]; p < M_.cs[j + 1]; p++) out[M_.ri[p]] = M_.av[p];
  } else {
    out[j - n_] = -1.0;  // slack column: Ax - s = 0, s in [rl,ru]
  }
}

void Simplex::placeNonbasic(int j) {
  // Place each nonbasic structural variable at the bound consistent with the
  // sign of its cost, so the slack basis is dual-feasible at zero duals and
  // the dual simplex can be used as the phase-1/phase-2 method throughout
  // (rather than requiring a separately feasible starting basis).
  bool hasL = !isInf(lb_[j]), hasU = !isInf(ub_[j]);
  if (hasL && hasU) {
    if (cost_[j] >= 0) stat_[j] = 1; else stat_[j] = 2;
  } else if (hasL) stat_[j] = 1;
  else if (hasU) stat_[j] = 2;
  else stat_[j] = 3;
  x_[j] = stat_[j] == 1 ? lb_[j] : (stat_[j] == 2 ? ub_[j] : 0.0);
}

void Simplex::slackBasis() {
  std::fill(pos_.begin(), pos_.end(), -1);
  for (int j = 0; j < n_; j++) { placeNonbasic(j); }
  for (int i = 0; i < m_; i++) {
    int s = n_ + i;
    basic_[i] = s; pos_[s] = i; stat_[s] = 0;
  }
}

// Sparse LU of the current basis B (m x m); column k of B is column basic_[k] of [A | -I].
// Right-looking elimination with Markowitz pivot choice and a column-wise stability threshold.
bool Simplex::factor() {
  const int m = m_;
  etas_.clear();
  typedef std::pair<int, double> Ent;
  std::vector<std::vector<Ent>> rows(m);          // active rows: (column id, value)
  std::vector<std::vector<int>> colRows(m);       // rows that may hold the column (may be stale)
  std::vector<int> colCnt(m, 0);
  for (int k = 0; k < m; k++) {
    int j = basic_[k];
    if (j < n_) {
      for (int p = M_.cs[j]; p < M_.cs[j + 1]; p++) {
        if (M_.av[p] == 0.0) continue;
        rows[M_.ri[p]].push_back({k, M_.av[p]}); colRows[k].push_back(M_.ri[p]); colCnt[k]++;
      }
    } else {
      rows[j - n_].push_back({k, -1.0}); colRows[k].push_back(j - n_); colCnt[k]++;
    }
  }
  std::vector<char> rowDone(m, 0), colDone(m, 0);
  std::vector<int> pos(m, -1);
  prow_.assign(m, 0); pcol_.assign(m, 0); pivVal_.assign(m, 0.0);
  lPtr_.assign(m + 1, 0); lIdx_.clear(); lVal_.clear();
  uRowPtr_.assign(m + 1, 0); uRowIdx_.clear(); uRowVal_.clear();
  const double thr = 0.1;                         // threshold pivoting parameter

  auto findEntry = [&](int i, int c) -> int {
    const auto& r = rows[i];
    for (int q = 0; q < (int)r.size(); q++) if (r[q].first == c) return q;
    return -1;
  };

  // Examines column c and returns the best Markowitz pivot row for it (or -1).
  auto bestInColumn = [&](int c, long& score, double& absv) -> int {
    std::vector<int>& cr = colRows[c];
    std::vector<int> valid; valid.reserve(cr.size());
    std::vector<double> vals;
    double cmax = 0;
    for (int i : cr) {
      if (rowDone[i]) continue;
      int q = findEntry(i, c); if (q < 0) continue;
      bool dup = false; for (int v : valid) if (v == i) { dup = true; break; }
      if (dup) continue;
      valid.push_back(i); vals.push_back(std::fabs(rows[i][q].second)); cmax = std::max(cmax, vals.back());
    }
    cr = valid;                                   // prune stale entries
    int bestI = -1; score = LONG_MAX; absv = 0;
    for (size_t q = 0; q < valid.size(); q++) {
      double a = vals[q];
      if (a < thr * cmax || a < 1e-11) continue;
      long sc = (long)((int)rows[valid[q]].size() - 1) * (long)(colCnt[c] - 1);
      if (sc < score || (sc == score && a > absv)) { score = sc; absv = a; bestI = valid[q]; }
    }
    return bestI;
  };

  std::vector<std::pair<int, int>> cand;
  for (int t = 0; t < m; t++) {
    int bestR = -1, bestC = -1; long bestScore = LONG_MAX; double bestAbs = 0;
    cand.clear();
    for (int c = 0; c < m; c++) if (!colDone[c]) {
      if (colCnt[c] == 0) { numEvents_++; return false; }   // structurally singular
      cand.push_back({colCnt[c], c});
    }
    size_t take = std::min<size_t>(4, cand.size());
    std::partial_sort(cand.begin(), cand.begin() + take, cand.end());
    for (size_t q = 0; q < take; q++) {
      long sc; double av; int r = bestInColumn(cand[q].second, sc, av);
      if (r >= 0 && (sc < bestScore || (sc == bestScore && av > bestAbs))) { bestScore = sc; bestAbs = av; bestR = r; bestC = cand[q].second; }
      if (bestScore == 0) break;
    }
    if (bestR < 0) {                               // fall back to every remaining column
      for (size_t q = take; q < cand.size(); q++) {
        long sc; double av; int r = bestInColumn(cand[q].second, sc, av);
        if (r >= 0 && (sc < bestScore || (sc == bestScore && av > bestAbs))) { bestScore = sc; bestAbs = av; bestR = r; bestC = cand[q].second; }
      }
    }
    if (bestR < 0) { numEvents_++; return false; }          // numerically singular
    int r = bestR, c = bestC;
    double p = rows[r][findEntry(r, c)].second;
    prow_[t] = r; pcol_[t] = c; pivVal_[t] = p;
    // U row: the rest of the pivot row
    for (const Ent& e : rows[r]) { colCnt[e.first]--; if (e.first != c) { uRowIdx_.push_back(e.first); uRowVal_.push_back(e.second); } }
    uRowPtr_[t + 1] = (int)uRowIdx_.size();
    rowDone[r] = 1; colDone[c] = 1;
    // eliminate column c from the other active rows
    std::vector<int> targets;
    for (int i : colRows[c]) if (!rowDone[i] && findEntry(i, c) >= 0) {
      bool dup = false; for (int v : targets) if (v == i) { dup = true; break; }
      if (!dup) targets.push_back(i);
    }
    const std::vector<Ent> prowEnts = rows[r];
    for (int i : targets) {
      auto& ri = rows[i];
      int q = findEntry(i, c);
      double l = ri[q].second / p;
      ri[q] = ri.back(); ri.pop_back(); colCnt[c]--;
      lIdx_.push_back(i); lVal_.push_back(l);
      for (int u = 0; u < (int)ri.size(); u++) pos[ri[u].first] = u;
      for (const Ent& e : prowEnts) {
        if (e.first == c) continue;
        int u = pos[e.first];
        if (u >= 0) ri[u].second -= l * e.second;
        else { ri.push_back({e.first, -l * e.second}); pos[e.first] = (int)ri.size() - 1; colCnt[e.first]++; colRows[e.first].push_back(i); }
      }
      for (const Ent& e : ri) pos[e.first] = -1;
      for (int u = 0; u < (int)ri.size();) {                // drop entries that cancelled to zero
        if (std::fabs(ri[u].second) < 1e-14) { colCnt[ri[u].first]--; ri[u] = ri.back(); ri.pop_back(); } else u++;
      }
    }
    lPtr_[t + 1] = (int)lIdx_.size();
    rows[r].clear(); rows[r].shrink_to_fit();
  }
  // U stored by column (indexed by column id) with row ids, for column-oriented back substitution.
  uColPtr_.assign(m + 1, 0);
  for (int c2 : uRowIdx_) uColPtr_[c2 + 1]++;
  for (int c2 = 0; c2 < m; c2++) uColPtr_[c2 + 1] += uColPtr_[c2];
  uColIdx_.assign(uRowIdx_.size(), 0); uColVal_.assign(uRowIdx_.size(), 0.0);
  std::vector<int> fillp(uColPtr_.begin(), uColPtr_.end() - 1);
  for (int t = 0; t < m; t++)
    for (int q = uRowPtr_[t]; q < uRowPtr_[t + 1]; q++) {
      int c2 = uRowIdx_[q]; int w = fillp[c2]++;
      uColIdx_[w] = prow_[t]; uColVal_[w] = uRowVal_[q];
    }
  return true;
}

// Solve B x = v. Input is indexed by row, output by basis position.
void Simplex::ftran(std::vector<double>& v) const {
  for (int t = 0; t < m_; t++) {                       // forward: apply the L multipliers
    double vr = v[prow_[t]];
    if (vr == 0.0) continue;
    for (int q = lPtr_[t]; q < lPtr_[t + 1]; q++) v[lIdx_[q]] -= lVal_[q] * vr;
  }
  std::vector<double> x(m_, 0.0);
  for (int t = m_ - 1; t >= 0; t--) {                  // backward: column-oriented U solve
    double rv = v[prow_[t]];
    if (rv == 0.0) continue;
    int c = pcol_[t];
    double xv = rv / pivVal_[t];
    x[c] = xv;
    for (int q = uColPtr_[c]; q < uColPtr_[c + 1]; q++) v[uColIdx_[q]] -= uColVal_[q] * xv;
  }
  v.swap(x);
  for (const Eta& e : etas_) {                         // product-form updates, oldest first
    double t = v[e.r];
    if (t == 0.0) continue;
    for (size_t q = 0; q < e.idx.size(); q++) v[e.idx[q]] += e.val[q] * t;
    v[e.r] = t * e.ar;
  }
}

// Solve B^T y = v. Input is indexed by basis position, output by row.
void Simplex::btran(std::vector<double>& v) const {
  for (auto it = etas_.rbegin(); it != etas_.rend(); ++it) {
    double s = it->ar * v[it->r];
    for (size_t q = 0; q < it->idx.size(); q++) s += it->val[q] * v[it->idx[q]];
    v[it->r] = s;
  }
  std::vector<double> z(m_, 0.0);
  for (int t = 0; t < m_; t++) {                       // U^T solve, column-oriented
    int c = pcol_[t];
    double wc = v[c];
    if (wc == 0.0) continue;
    double zr = wc / pivVal_[t];
    z[prow_[t]] = zr;
    for (int q = uRowPtr_[t]; q < uRowPtr_[t + 1]; q++) v[uRowIdx_[q]] -= uRowVal_[q] * zr;
  }
  for (int t = m_ - 1; t >= 0; t--) {                  // apply L^T
    double s = 0.0;
    for (int q = lPtr_[t]; q < lPtr_[t + 1]; q++) s += lVal_[q] * z[lIdx_[q]];
    z[prow_[t]] -= s;
  }
  v.swap(z);
}

double Simplex::dotCol(int j, const std::vector<double>& y) const {
  if (j < n_) { double s = 0; for (int p = M_.cs[j]; p < M_.cs[j + 1]; p++) s += M_.av[p] * y[M_.ri[p]]; return s; }
  return -y[j - n_];
}

void Simplex::computePrimal() {
  std::vector<double> rhs(m_, 0.0);
  for (int j = 0; j < N_; j++) if (pos_[j] < 0 || basic_[pos_[j]] != j) {
    if (stat_[j] == 0) continue;
    double xv = x_[j];
    if (xv == 0) continue;
    if (j < n_) for (int p = M_.cs[j]; p < M_.cs[j + 1]; p++) rhs[M_.ri[p]] -= M_.av[p] * xv;
    else rhs[j - n_] -= -xv;
  }
  ftran(rhs);
  for (int i = 0; i < m_; i++) x_[basic_[i]] = rhs[i];
}

// d[j] = cst[j] - a_j . y for every column. Runs on the GPU backend when one is attached.
void Simplex::priceAll(const std::vector<double>& cst, const std::vector<double>& y) {
  double t0 = now();
  if (pricer_) pricer_->computeD(cst.data(), y.data(), d_.data());
  else for (int j = 0; j < N_; j++) d_[j] = cst[j] - dotCol(j, y);
  pricingSeconds_ += now() - t0;
}

void Simplex::computeDuals(const std::vector<double>& cst) {
  std::vector<double> cb(m_);
  for (int i = 0; i < m_; i++) cb[i] = cst[basic_[i]];
  std::vector<double> y = cb;
  btran(y);
  priceAll(cst, y);
}

void Simplex::perturbCosts() {
  if (perturbed_) return;
  for (int j = 0; j < N_; j++) cost_[j] = cost0_[j] + 1e-7 * (rnd() - 0.5) * (1.0 + std::fabs(cost0_[j]));
  perturbed_ = true;
}
void Simplex::restoreCosts() { cost_ = cost0_; perturbed_ = false; }

// ---------- Primal simplex, phase 1 (composite objective) and phase 2 ----------
// Phase 1 minimizes total bound infeasibility of the current basic variables.
// Each iteration, every basic variable gets a fresh phase-1 cost and an
// effective bound based on which side (if any) it currently violates, so a
// variable that becomes feasible mid-phase is immediately treated normally
// again (Wolfe's composite-objective method). Phase 2 is the same pivot code
// with phase1=false, using the model's real costs and bounds throughout.
LPStatus Simplex::primalPhase(double deadline, bool phase1) {
  std::vector<double> col;
  std::vector<double> c1(N_), elb(N_), eub(N_);
  for (;;) {
    if (now() > deadline) return LPStatus::Limit;
    if (iters_ % 200 == 0) { if (!factor()) return LPStatus::Numerical; computePrimal(); }
    const std::vector<double>* costPtr;
    if (phase1) {
      for (int j = 0; j < N_; j++) { c1[j] = 0; elb[j] = lb_[j]; eub[j] = ub_[j]; }
      double infeas = 0;
      for (int i = 0; i < m_; i++) {
        int bj = basic_[i];
        if (!isInf(lb_[bj]) && x_[bj] < lb_[bj] - opt_.feasTol) { c1[bj] = -1; elb[bj] = -INF; infeas += lb_[bj] - x_[bj]; }
        else if (!isInf(ub_[bj]) && x_[bj] > ub_[bj] + opt_.feasTol) { c1[bj] = 1; eub[bj] = INF; infeas += x_[bj] - ub_[bj]; }
      }
      if (infeas < opt_.feasTol) return LPStatus::Optimal;  // phase 1 done, feasible
      costPtr = &c1;
    } else {
      elb = lb_; eub = ub_;
      costPtr = &cost_;
    }
    const std::vector<double>& cst = *costPtr;
    // duals from current basis costs
    std::vector<double> cb(m_); for (int i = 0; i < m_; i++) cb[i] = cst[basic_[i]];
    std::vector<double> y = cb; btran(y);
    priceAll(cst, y);

    int enter = -1; double best = opt_.optTol; int dir = 0;
    for (int j = 0; j < N_; j++) {
      if (pos_[j] >= 0 && basic_[pos_[j]] == j) continue;
      double dj = d_[j];
      bool hasL = !isInf(elb[j]), hasU = !isInf(eub[j]);
      if (stat_[j] == 1 && hasL && dj < -best) { best = -dj; enter = j; dir = 1; }
      else if (stat_[j] == 2 && hasU && dj > best) { best = dj; enter = j; dir = -1; }
      else if (stat_[j] == 3 && std::fabs(dj) > best) { best = std::fabs(dj); enter = j; dir = dj < 0 ? 1 : -1; }
    }
    if (enter < 0) { if (phase1) return LPStatus::Infeasible; finalize(); return LPStatus::Optimal; }
    colOf(enter, col); ftran(col);
    double step = INF; int leave = -1; int leaveTo = 1; bool boundFlip = false;
    double flipLimit = (isInf(elb[enter]) || isInf(eub[enter])) ? INF : (eub[enter] - elb[enter]);
    // Two-pass (Harris-style) ratio test: pass 1 finds the largest step that keeps every
    // basic variable within a small tolerance of its bound; pass 2 picks, among the rows
    // that block no later than that step, the one with the largest pivot magnitude.
    // This avoids tiny pivots that make the basis near-singular.
    struct Cand { int i; double t; double a; int to; };
    std::vector<Cand> cands;
    double harris = INF;
    const double ratioTol = 1e-9;
    for (int i = 0; i < m_; i++) {
      double a = col[i] * dir;
      int bj = basic_[i];
      double t = INF; int to = 0;
      if (phase1 && c1[bj] == -1 && a < -opt_.pivTol) { t = (lb_[bj] - x_[bj]) / (-a); to = 1; }
      else if (phase1 && c1[bj] == 1 && a > opt_.pivTol) { t = (x_[bj] - ub_[bj]) / a; to = 2; }
      else if (a > opt_.pivTol) { if (isInf(elb[bj])) continue; t = (x_[bj] - elb[bj]) / a; to = 1; }
      else if (a < -opt_.pivTol) { if (isInf(eub[bj])) continue; t = (eub[bj] - x_[bj]) / (-a); to = 2; }
      else continue;
      if (t < 0) t = 0;
      cands.push_back({i, t, std::fabs(a), to});
      harris = std::min(harris, t + ratioTol / std::fabs(a));
    }
    double bestA = 0;
    for (auto& c : cands) if (c.t <= harris && c.a > bestA) { bestA = c.a; leave = c.i; leaveTo = c.to; step = c.t; }
    if (step > flipLimit) { step = flipLimit; boundFlip = true; }
    if (step >= INF) {
      if (phase1) return LPStatus::Numerical;  // phase 1 objective cannot be unbounded on a well-posed model
      ray_.assign(N_, 0); ray_[enter] = dir; for (int i = 0; i < m_; i++) ray_[basic_[i]] = -col[i] * dir; rayKind_ = "primal_unbounded"; return LPStatus::Unbounded;
    }
    x_[enter] += step * dir;
    for (int i = 0; i < m_; i++) x_[basic_[i]] -= step * dir * col[i];
    if (boundFlip) { stat_[enter] = stat_[enter] == 1 ? 2 : 1; }
    else {
      int leaveVar = basic_[leave];
      stat_[leaveVar] = leaveTo; x_[leaveVar] = leaveTo == 1 ? (isInf(elb[leaveVar]) ? lb_[leaveVar] : elb[leaveVar]) : (isInf(eub[leaveVar]) ? ub_[leaveVar] : eub[leaveVar]);
      pos_[leaveVar] = -1; basic_[leave] = enter; pos_[enter] = leave; stat_[enter] = 0;
      {
        Eta e; e.r = leave; e.ar = 1.0 / col[leave];
        for (int i = 0; i < m_; i++) if (i != leave && col[i] != 0.0) { e.idx.push_back(i); e.val.push_back(-col[i] / col[leave]); }
        etas_.push_back(std::move(e));
        if (etas_.size() > 100) { if (!factor()) return LPStatus::Numerical; }
      }
    }
    iters_++;
    if (iters_ > (opt_.maxIter > 0 ? opt_.maxIter : 20000 + 50 * N_)) return LPStatus::Limit;
  }
}

// ---------- Dual simplex ----------
// Dual simplex with cost shifting.
// 1. Make the start dual feasible: boxed variables sit at the bound matching the sign of their
//    reduced cost; any remaining dual infeasibility is removed by shifting that variable's cost.
// 2. Iterate: pick the most infeasible basic variable, choose the entering variable by the dual
//    ratio test (Harris two-pass, largest pivot among near-ties), pivot.
// 3. At primal feasibility the caller restores the true costs and finishes with primal phase 2,
//    so the shifts can never produce a wrong answer; they only cost a few extra primal pivots.
LPStatus Simplex::dual(double deadline) {
  const double dualTol = 1e-9;
  auto isBasic = [&](int j) { return pos_[j] >= 0 && basic_[pos_[j]] == j; };
  auto fixedVar = [&](int j) { return !isInf(lb_[j]) && !isInf(ub_[j]) && lb_[j] == ub_[j]; };

  cost_ = cost0_;
  computeDuals(cost_);
  for (int j = 0; j < N_; j++) {
    if (isBasic(j)) continue;
    bool hasL = !isInf(lb_[j]), hasU = !isInf(ub_[j]);
    double dj = d_[j];
    if (hasL && hasU) stat_[j] = (dj >= 0 || lb_[j] == ub_[j]) ? 1 : 2;
    else if (hasL) { stat_[j] = 1; if (dj < -dualTol) { cost_[j] -= dj; d_[j] = 0; } }
    else if (hasU) { stat_[j] = 2; if (dj > dualTol) { cost_[j] -= dj; d_[j] = 0; } }
    else { stat_[j] = 3; if (std::fabs(dj) > dualTol) { cost_[j] -= dj; d_[j] = 0; } }
    x_[j] = stat_[j] == 1 ? lb_[j] : (stat_[j] == 2 ? ub_[j] : 0.0);
  }
  computePrimal();

  std::vector<double> col, rho(m_), alphaRow(N_);
  for (;;) {
    if (now() > deadline) return LPStatus::Limit;
    if (iters_ % 100 == 0 || etas_.size() > 100) {
      if (!factor()) return LPStatus::Numerical;
      computePrimal(); computeDuals(cost_);
    }
    // leaving row: largest bound violation among basic variables
    int leave = -1; double worst = opt_.feasTol; int leaveDir = 0;
    for (int i = 0; i < m_; i++) {
      int bj = basic_[i];
      if (!isInf(lb_[bj]) && x_[bj] < lb_[bj] - worst) { worst = lb_[bj] - x_[bj]; leave = i; leaveDir = 1; }
      else if (!isInf(ub_[bj]) && x_[bj] > ub_[bj] + worst) { worst = x_[bj] - ub_[bj]; leave = i; leaveDir = -1; }
    }
    if (leave < 0) { return LPStatus::Optimal; }   // primal feasible for the (shifted) costs
    std::fill(rho.begin(), rho.end(), 0.0); rho[leave] = 1.0; btran(rho);
    // pivot row alpha_j = rho . a_j for every nonbasic j; normalise so that alpha > 0 means
    // "increasing x_j moves the leaving variable toward its violated bound"
    const double sgn = (leaveDir == 1) ? -1.0 : 1.0;
    double harris = INF;
    for (int j = 0; j < N_; j++) {
      alphaRow[j] = 0.0;
      if (isBasic(j) || fixedVar(j)) continue;
      double a = dotCol(j, rho); alphaRow[j] = a;
      double an = a * sgn;
      bool up = (stat_[j] == 1 || stat_[j] == 3), down = (stat_[j] == 2 || stat_[j] == 3);
      if ((an > opt_.pivTol && up) || (an < -opt_.pivTol && down))
        harris = std::min(harris, (std::fabs(d_[j]) + dualTol) / std::fabs(an));
    }
    int enter = -1; double bestA = 0;
    for (int j = 0; j < N_; j++) {
      if (isBasic(j) || fixedVar(j)) continue;
      double an = alphaRow[j] * sgn;
      bool up = (stat_[j] == 1 || stat_[j] == 3), down = (stat_[j] == 2 || stat_[j] == 3);
      if (!((an > opt_.pivTol && up) || (an < -opt_.pivTol && down))) continue;
      if (std::fabs(d_[j]) / std::fabs(an) <= harris && std::fabs(an) > bestA) { bestA = std::fabs(an); enter = j; }
    }
    if (enter < 0) {
      // No variable can repair this row: the model is primal infeasible, whatever the costs are.
      ray_.assign(m_, 0); for (int i = 0; i < m_; i++) ray_[i] = rho[i] * (leaveDir == 1 ? 1 : -1);
      rayKind_ = "dual_infeasible_primal";
      return LPStatus::Infeasible;
    }
    colOf(enter, col); ftran(col);
    double alphaRq = col[leave];
    if (std::fabs(alphaRq) < 1e-9) { if (!factor()) return LPStatus::Numerical; computePrimal(); computeDuals(cost_); continue; }
    int leaveVar = basic_[leave];
    double target = leaveDir == 1 ? lb_[leaveVar] : ub_[leaveVar];
    double step = (x_[leaveVar] - target) / alphaRq;
    // dual update (reduced costs)
    double theta = d_[enter] / alphaRq;
    for (int j = 0; j < N_; j++) if (alphaRow[j] != 0.0) d_[j] -= theta * alphaRow[j];
    d_[enter] = 0.0; d_[leaveVar] = -theta;
    // primal update
    x_[enter] += step;
    for (int i = 0; i < m_; i++) if (i != leave) x_[basic_[i]] -= step * col[i];
    x_[leaveVar] = target;
    stat_[leaveVar] = leaveDir == 1 ? 1 : 2;
    pos_[leaveVar] = -1; basic_[leave] = enter; pos_[enter] = leave; stat_[enter] = 0;
    {
      Eta e; e.r = leave; e.ar = 1.0 / alphaRq;
      for (int i = 0; i < m_; i++) if (i != leave && col[i] != 0.0) { e.idx.push_back(i); e.val.push_back(-col[i] / alphaRq); }
      etas_.push_back(std::move(e));
    }
    iters_++;
    if (iters_ > (opt_.maxIter > 0 ? opt_.maxIter : 20000 + 50 * N_)) return LPStatus::Limit;
  }
}

void Simplex::finalize() { computePrimal(); computeDuals(cost_); }

LPStatus Simplex::solve(const std::string& method, double deadline) {
  bool haveBasis = false;
  for (int i = 0; i < m_; i++) if (pos_[basic_[i]] == i) { haveBasis = true; break; }
  if (!haveBasis) slackBasis();
  if (!factor()) { slackBasis(); factor(); }
  computePrimal();
  // "dual" is kept available for development/comparison; it requires a
  // dual-feasible start that most industrial models don't naturally have
  // (see placeNonbasic), so the default path below always uses the primal
  // engine with an explicit phase 1 for feasibility, which is correct
  // regardless of cost/bound sign combinations.
  if (method == "dual") {
    long it0 = iters_;
    LPStatus st = dual(deadline);
    if (getenv("BS_DEBUG")) std::cerr << "[dual] dual-simplex iterations: " << (iters_ - it0) << "\n";
    restoreCosts();                                       // drop the cost shifts
    if (st == LPStatus::Optimal) {
      st = primalPhase(deadline, false);                  // finish with the true costs
      if (st == LPStatus::Numerical) { numEvents_++; slackBasis(); factor(); computePrimal(); return solve("primal", deadline); }
      return st;
    }
    if (st == LPStatus::Numerical) { numEvents_++; slackBasis(); factor(); computePrimal(); return solve("primal", deadline); }
    return st;                                            // infeasible or limit
  }

  bool primalFeas = true;
  for (int i = 0; i < m_; i++) { int bj = basic_[i]; if ((!isInf(lb_[bj]) && x_[bj] < lb_[bj] - opt_.feasTol) || (!isInf(ub_[bj]) && x_[bj] > ub_[bj] + opt_.feasTol)) { primalFeas = false; break; } }
  LPStatus st = LPStatus::Optimal;
  if (!primalFeas) {
    st = primalPhase(deadline, true);
    if (st != LPStatus::Optimal) {
      if (st == LPStatus::Numerical) { numEvents_++; }
      return (st == LPStatus::Infeasible) ? LPStatus::Infeasible : st;
    }
  }
  st = primalPhase(deadline, false);
  if (st == LPStatus::Numerical) {
    numEvents_++;
    perturbCosts();
    slackBasis(); factor(); computePrimal();
    if (!primalFeas) { st = primalPhase(deadline, true); if (st == LPStatus::Optimal) st = primalPhase(deadline, false); }
    else st = primalPhase(deadline, false);
    restoreCosts();
    if (st == LPStatus::Optimal) computeDuals(cost_);
  }
  return st;
}

double Simplex::objective() const { double s = 0; for (int j = 0; j < n_; j++) s += cost0_[j] * x_[j]; return s; }

void Simplex::duals(std::vector<double>& y, std::vector<double>& d) const {
  y.assign(m_, 0); d.assign(n_, 0);
  for (int i = 0; i < m_; i++) y[i] = d_[n_ + i];
  for (int j = 0; j < n_; j++) d[j] = d_[j];
}

void Simplex::setBasis(const Basis& b) { basic_ = b.basic; stat_ = b.stat; std::fill(pos_.begin(), pos_.end(), -1); for (int i = 0; i < m_; i++) pos_[basic_[i]] = i; }
Basis Simplex::getBasis() const { Basis b; b.basic = basic_; b.stat = stat_; return b; }

}  // namespace bs
