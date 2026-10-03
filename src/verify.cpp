// Independent verification (VR requirements): re-derives feasibility and
// optimality residuals directly from the ORIGINAL model and the returned
// solution, without reusing any solver-internal state (basis, factorization).
#include "bs.hpp"
#include <cmath>

namespace bs {

Verification verify(const Model& M, const std::vector<double>& x, const std::vector<double>& y,
                    const std::vector<double>& /*d*/, bool lpDuals, const Options& o) {
  Verification v; v.ran = true;
  int n = M.n, m = M.m;
  // Primal feasibility: bounds + row activity, from scratch.
  std::vector<double> Ax(m, 0.0);
  for (int j = 0; j < n; j++) for (int p = M.cs[j]; p < M.cs[j + 1]; p++) Ax[M.ri[p]] += M.av[p] * x[j];
  double maxInf = 0;
  for (int i = 0; i < m; i++) {
    if (!isInf(M.rl[i]) && Ax[i] < M.rl[i] - 1e-9) maxInf = std::max(maxInf, M.rl[i] - Ax[i]);
    if (!isInf(M.ru[i]) && Ax[i] > M.ru[i] + 1e-9) maxInf = std::max(maxInf, Ax[i] - M.ru[i]);
  }
  for (int j = 0; j < n; j++) {
    if (!isInf(M.lb[j]) && x[j] < M.lb[j] - 1e-9) maxInf = std::max(maxInf, M.lb[j] - x[j]);
    if (!isInf(M.ub[j]) && x[j] > M.ub[j] + 1e-9) maxInf = std::max(maxInf, x[j] - M.ub[j]);
  }
  v.primalInf = maxInf;
  // Integrality
  double intInf = 0;
  for (int j = 0; j < n; j++) if (M.isint[j]) intInf = std::max(intInf, std::fabs(x[j] - std::round(x[j])));
  v.intInf = intInf;
  if (v.primalInf > o.verifyTol) v.failures.push_back("primal infeasibility " + std::to_string(v.primalInf) + " exceeds tolerance");
  if (intInf > o.verifyTol) v.failures.push_back("integrality violation " + std::to_string(intInf));

  if (lpDuals && !y.empty()) {
    v.hasDual = true;
    // Dual feasibility & complementary slackness, recomputed independently:
    // reduced cost r_j = c_j - sum_i A_ij y_i must match sign of active bound.
    std::vector<double> r(n);
    for (int j = 0; j < n; j++) { double s = M.c[j]; for (int p = M.cs[j]; p < M.cs[j + 1]; p++) s -= M.av[p] * y[M.ri[p]]; r[j] = s; }
    double dualInf = 0, comp = 0;
    for (int j = 0; j < n; j++) {
      bool atLower = !isInf(M.lb[j]) && std::fabs(x[j] - M.lb[j]) < 1e-6;
      bool atUpper = !isInf(M.ub[j]) && std::fabs(x[j] - M.ub[j]) < 1e-6;
      if (atLower && !atUpper) { if (r[j] < -1e-6) dualInf = std::max(dualInf, -r[j]); }
      else if (atUpper && !atLower) { if (r[j] > 1e-6) dualInf = std::max(dualInf, r[j]); }
      else if (!atLower && !atUpper) { dualInf = std::max(dualInf, std::fabs(r[j])); }
      comp = std::max(comp, std::fabs(r[j] * (x[j] - (atLower ? M.lb[j] : (atUpper ? M.ub[j] : x[j])))));
    }
    v.dualInf = dualInf; v.compl_ = comp;
    if (dualInf > o.verifyTol) v.failures.push_back("dual infeasibility " + std::to_string(dualInf));
  }
  v.ok = v.failures.empty();
  return v;
}

}  // namespace bs
