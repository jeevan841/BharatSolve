// Bounded revised simplex over the model's own m rows and n structural columns,
// using slack variables indexed n..n+m-1 for row i (slack bounds -ru<=s<=-rl-ish
// handled via the general form l <= Ax <= u through slack = Ax with bounds [rl,ru]).
// First release uses dense LU refactorization each `refactorEvery` iterations
// (VR/SX roadmap: replace with sparse LU + Forrest-Tomlin updates in a later pass).
#include "bs.hpp"
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

// Dense LU with partial pivoting on the current basis matrix B (m x m), columns = colOf(basic_[k]).
bool Simplex::factor() {
  lu_.assign((size_t)m_ * m_, 0.0);
  piv_.assign(m_, 0);
  std::vector<double> col;
  for (int k = 0; k < m_; k++) {
    colOf(basic_[k], col);
    for (int i = 0; i < m_; i++) lu_[(size_t)i * m_ + k] = col[i];
  }
  for (int i = 0; i < m_; i++) piv_[i] = i;
  for (int k = 0; k < m_; k++) {
    int piv = k; double best = std::fabs(lu_[(size_t)k * m_ + k]);
    for (int i = k + 1; i < m_; i++) {
      double v = std::fabs(lu_[(size_t)i * m_ + k]);
      if (v > best) { best = v; piv = i; }
    }
    if (best < 1e-11) { numEvents_++; return false; }
    if (piv != k) {
      for (int j = 0; j < m_; j++) std::swap(lu_[(size_t)k * m_ + j], lu_[(size_t)piv * m_ + j]);
      std::swap(piv_[k], piv_[piv]);
    }
    double diag = lu_[(size_t)k * m_ + k];
    for (int i = k + 1; i < m_; i++) {
      double f = lu_[(size_t)i * m_ + k] / diag;
      lu_[(size_t)i * m_ + k] = f;
      if (f != 0.0) for (int j = k + 1; j < m_; j++) lu_[(size_t)i * m_ + j] -= f * lu_[(size_t)k * m_ + j];
    }
  }
  etas_.clear();
  return true;
}

void Simplex::ftran(std::vector<double>& v) const {
  std::vector<double> pv(m_);
  for (int i = 0; i < m_; i++) pv[i] = v[piv_[i]];
  for (int i = 0; i < m_; i++) for (int j = 0; j < i; j++) pv[i] -= lu_[(size_t)i * m_ + j] * pv[j];
  for (int i = m_ - 1; i >= 0; i--) {
    for (int j = i + 1; j < m_; j++) pv[i] -= lu_[(size_t)i * m_ + j] * pv[j];
    pv[i] /= lu_[(size_t)i * m_ + i];
  }
  v = pv;
  for (auto& e : etas_) {
    double t = v[e.r];
    for (int i = 0; i < m_; i++) if (i != e.r) v[i] += e.a[i] * t;
    v[e.r] = t * e.a[e.r];
  }
}

void Simplex::btran(std::vector<double>& v) const {
  for (auto it = etas_.rbegin(); it != etas_.rend(); ++it) {
    double s = 0; for (int i = 0; i < m_; i++) s += it->a[i] * v[i];
    v[it->r] = s;
  }
  std::vector<double> y(m_);
  for (int i = 0; i < m_; i++) {
    double s = v[i];
    for (int j = 0; j < i; j++) s -= lu_[(size_t)j * m_ + i] * y[j];
    y[i] = s / lu_[(size_t)i * m_ + i];
  }
  std::vector<double> z(m_);
  for (int i = m_ - 1; i >= 0; i--) {
    double s = y[i];
    for (int j = i + 1; j < m_; j++) s -= lu_[(size_t)j * m_ + i] * z[j];
    z[i] = s;
  }
  std::vector<double> out(m_);
  for (int i = 0; i < m_; i++) out[piv_[i]] = z[i];
  v = out;
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

void Simplex::computeDuals(const std::vector<double>& cst) {
  std::vector<double> cb(m_);
  for (int i = 0; i < m_; i++) cb[i] = cst[basic_[i]];
  std::vector<double> y = cb;
  btran(y);
  for (int j = 0; j < N_; j++) d_[j] = cst[j] - dotCol(j, y);
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
    for (int j = 0; j < N_; j++) d_[j] = cst[j] - dotCol(j, y);

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
    for (int i = 0; i < m_; i++) {
      double a = col[i] * dir;
      int bj = basic_[i];
      if (a > opt_.pivTol) {
        if (isInf(elb[bj])) continue;
        double t = (x_[bj] - elb[bj]) / a;
        if (t < step - 1e-12) { step = t; leave = i; leaveTo = 1; }
      } else if (a < -opt_.pivTol) {
        if (isInf(eub[bj])) continue;
        double t = (eub[bj] - x_[bj]) / (-a);
        if (t < step - 1e-12) { step = t; leave = i; leaveTo = 2; }
      }
    }
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
      stat_[leaveVar] = leaveTo; x_[leaveVar] = leaveTo == 1 ? elb[leaveVar] : eub[leaveVar];
      pos_[leaveVar] = -1; basic_[leave] = enter; pos_[enter] = leave; stat_[enter] = 0;
      if (m_ <= 400) { std::vector<double> a(m_); for (int i = 0; i < m_; i++) a[i] = (i == leave) ? 1.0 / col[i] : -col[i] / col[leave]; etas_.push_back({leave, a}); if (etas_.size() > 60) { if (!factor()) return LPStatus::Numerical; } }
      else { if (!factor()) return LPStatus::Numerical; }
    }
    iters_++;
    if (iters_ > (opt_.maxIter > 0 ? opt_.maxIter : 20000 + 50 * N_)) return LPStatus::Limit;
  }
}

// ---------- Dual simplex ----------
LPStatus Simplex::dual(double deadline) {
  std::vector<double> row(m_), col;
  if (getenv("BS_DEBUG")) {
    computeDuals(cost_);
    double worstDF = 0; int worstJ = -1;
    for (int j = 0; j < N_; j++) {
      if (pos_[j] >= 0 && basic_[pos_[j]] == j) continue;
      double dj = d_[j], v = 0;
      if (stat_[j] == 1) v = std::max(0.0, -dj);
      else if (stat_[j] == 2) v = std::max(0.0, dj);
      else v = std::fabs(dj);
      if (v > worstDF) { worstDF = v; worstJ = j; }
    }
    std::cerr << "[dual] initial dual infeasibility=" << worstDF << " at j=" << worstJ
               << " stat=" << (worstJ>=0?(int)stat_[worstJ]:-1) << "\n";
  }
  for (;;) {
    if (now() > deadline) return LPStatus::Limit;
    if (iters_ % 200 == 0) { if (!factor()) return LPStatus::Numerical; }
    computePrimal();
    int leave = -1; double worst = opt_.feasTol; int leaveDir = 0;
    for (int i = 0; i < m_; i++) {
      int bj = basic_[i];
      if (!isInf(lb_[bj]) && x_[bj] < lb_[bj] - worst) { worst = lb_[bj] - x_[bj]; leave = i; leaveDir = 1; }
      else if (!isInf(ub_[bj]) && x_[bj] > ub_[bj] + worst) { worst = x_[bj] - ub_[bj]; leave = i; leaveDir = -1; }
    }
    if (leave < 0) { finalize(); return LPStatus::Optimal; }
    // row `leave` of B^-1 A
    std::vector<double> e(m_, 0.0); e[leave] = 1.0; btran(e);
    computeDuals(cost_);
    int enter = -1; double bestRatio = INF;
    for (int j = 0; j < N_; j++) {
      if (pos_[j] >= 0 && basic_[pos_[j]] == j) continue;
      double alpha = dotCol(j, e) * (leaveDir == 1 ? -1.0 : 1.0);
      bool canIncrease = (stat_[j] == 1) || (stat_[j] == 3);
      bool canDecrease = (stat_[j] == 2) || (stat_[j] == 3);
      if (alpha > opt_.pivTol && canIncrease) { double r = d_[j] / alpha; if (r < bestRatio - 1e-12) { bestRatio = r; enter = j; } }
      else if (alpha < -opt_.pivTol && canDecrease) { double r = d_[j] / alpha; if (r < bestRatio - 1e-12) { bestRatio = r; enter = j; } }
    }
    if (enter < 0) {
      if (getenv("BS_DEBUG")) std::cerr << "[dual] no entering var found; leave=" << leave << " leaveDir=" << leaveDir << " worst=" << worst << " iter=" << iters_ << "\n";
      ray_.assign(m_, 0); for (int i = 0; i < m_; i++) ray_[i] = e[i] * (leaveDir == 1 ? 1 : -1); rayKind_ = "dual_infeasible_primal"; return LPStatus::Infeasible;
    }
    if (getenv("BS_DEBUG")) std::cerr << "[dual] iter=" << iters_ << " leave=" << leave << "(" << leaveDir << ") enter=" << enter << " ratio=" << bestRatio << " obj=" << objective() << "\n";
    colOf(enter, col); ftran(col);
    int leaveVar = basic_[leave];
    double target = leaveDir == 1 ? lb_[leaveVar] : ub_[leaveVar];
    double step = (x_[leaveVar] - target) / col[leave];
    x_[enter] += step;
    for (int i = 0; i < m_; i++) if (i != leave) x_[basic_[i]] -= step * col[i];
    x_[leaveVar] = target;
    stat_[leaveVar] = leaveDir == 1 ? 1 : 2;
    pos_[leaveVar] = -1; basic_[leave] = enter; pos_[enter] = leave; stat_[enter] = 0;
    if (!factor()) return LPStatus::Numerical;
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
  if (method == "dual") return dual(deadline);

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
