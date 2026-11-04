// time_walk.cpp
//
// Log-likelihood of the event times under an exponential Hawkes kernel with
// jumps G_j and an optional fast component with unit jumps,
//   lambda(t) = mu(t) + K * sum_{t_j < t} G_j * exp(-beta * (t - t_j))
//                     + K2 * sum_{t_j < t} exp(-beta2 * (t - t_j)),
// excluding the baseline integral, which the caller adds from its quadrature.
// With excitation fixed at birth G_j = g_j = exp(z_j' a). With path-dependent
// excitation G_j also carries the changes to earlier events' g at t_j,
// discounted to t_j (see path_jumps_cpp()), so it depends on beta too.
// The caller passes dG (n x p, dG_j / da_m) and dGb (dG_j / dbeta, or length
// 0 for none). The recursions are, with S+ etc. the values just after the
// jump at t_{i-1},
//   S_i   = exp(-beta * dt) * S+_{i-1}                 (S+ = S + G)
//   R_i   = exp(-beta * dt) * (R_{i-1} + dt * S+_{i-1})
//   U_m,i = exp(-beta * dt) * U+_m,{i-1}               (U+ = U + dG_m)
//   V_i   = exp(-beta * dt) * V+_{i-1}                 (V+ = V + dGb)
// so dS_i / dbeta = V_i - R_i and dS_i / da_m = U_m,i; S2 and R2 are the
// same for the fast component with G = 1. K2 = 0 drops it.

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Compensator of one exponential component, (1 / beta) * sum_i g_i (1 - exp(-beta tau_i)),
// and the derivative of the sum with respect to beta. 1 - exp(-x) is taken as
// -expm1(-x): for small beta * tau the naive form rounds to 0, which removes
// the penalty on K and makes the likelihood unbounded as beta -> 0.
static void exp_compensator(NumericVector t, NumericVector g, double beta, double T_end,
                            double& comp, double& dcomp_dbeta) {
  comp = 0.0;
  dcomp_dbeta = 0.0;
  for (int i = 0; i < t.size(); ++i) {
    const double tau = T_end - t[i];
    const double x = beta * tau;
    const double one_minus_e = -std::expm1(-x);
    comp += g[i] * one_minus_e;
    // d/dbeta of (1 - exp(-beta tau)) / beta is tau^2 (x e^-x - (1 - e^-x)) / x^2,
    // whose two terms cancel for small x; use its series there.
    const double r = x < 1e-4 ? -0.5 + x / 3.0 - x * x / 8.0
                              : (x * std::exp(-x) - one_minus_e) / (x * x);
    dcomp_dbeta += g[i] * tau * tau * r;
  }
}

// [[Rcpp::export]]
List time_ll_cpp(NumericVector t, NumericVector mu_ev, NumericVector G,
                 NumericMatrix dG, NumericVector dGb, double K, double beta, double T_end,
                 double K2 = 0.0, double beta2 = 1.0) {
  const int n = t.size();
  const int p = dG.ncol();
  const bool has_b = dGb.size() == n;
  NumericVector inv_lambda(n), log_lambda(n), da(p);
  std::vector<double> U(p, 0.0), dcomp_da(p, 0.0);

  double S = 0.0, R = 0.0, V = 0.0, sum_log = 0.0, dK = 0.0, dbeta = 0.0;
  double S2 = 0.0, R2 = 0.0, dK2 = 0.0, dbeta2 = 0.0;
  const bool fast = K2 > 0.0;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(G[i])) return List::create(_["ll"] = R_NegInf);
    if (i > 0) {
      const double dt = t[i] - t[i - 1];
      if (dt < 0) stop("times must be increasing");
      const double e = std::exp(-beta * dt);
      R = e * (R + dt * S);
      S = e * S;
      V = e * V;
      for (int m = 0; m < p; ++m) U[m] *= e;
      if (fast) {
        const double e2 = std::exp(-beta2 * dt);
        R2 = e2 * (R2 + dt * S2);
        S2 = e2 * S2;
      }
    }
    const double lam = mu_ev[i] + K * S + K2 * S2;
    if (!std::isfinite(lam) || lam <= 0.0) {
      return List::create(_["ll"] = R_NegInf);
    }
    log_lambda[i] = std::log(lam);
    inv_lambda[i] = 1.0 / lam;
    sum_log += log_lambda[i];
    dK += S / lam;
    dbeta += K * (V - R) / lam;
    for (int m = 0; m < p; ++m) da[m] += K * U[m] / lam;
    dK2 += S2 / lam;
    dbeta2 -= K2 * R2 / lam;

    // The jump at t_i adds to S, U and V but not to R (zero lag).
    S += G[i];
    if (has_b) V += dGb[i];
    for (int m = 0; m < p; ++m) U[m] += dG(i, m);
    if (fast) S2 += 1.0;
  }

  // Excitation parts of the compensator.
  double comp, dcomp_dbeta, comp2 = 0.0, dcomp2_dbeta2 = 0.0;
  exp_compensator(t, G, beta, T_end, comp, dcomp_dbeta);
  if (fast) {
    exp_compensator(t, NumericVector(n, 1.0), beta2, T_end, comp2, dcomp2_dbeta2);
  }
  double dcomp_jump_b = 0.0;
  for (int i = 0; i < n; ++i) {
    const double one_minus_e = -std::expm1(-beta * (T_end - t[i]));
    for (int m = 0; m < p; ++m) dcomp_da[m] += dG(i, m) * one_minus_e;
    if (has_b) dcomp_jump_b += dGb[i] * one_minus_e;
  }
  for (int m = 0; m < p; ++m) da[m] -= (K / beta) * dcomp_da[m];

  return List::create(
    _["ll"] = sum_log - (K / beta) * comp - (fast ? (K2 / beta2) * comp2 : 0.0),
    _["dK"] = dK - comp / beta,
    _["dbeta"] = dbeta - K * dcomp_dbeta - (K / beta) * dcomp_jump_b,
    _["da"] = da,
    _["dK2"] = dK2 - (fast ? comp2 / beta2 : 0.0),
    _["dbeta2"] = dbeta2 - K2 * dcomp2_dbeta2,
    _["inv_lambda"] = inv_lambda,
    _["log_lambda"] = log_lambda
  );
}

// Effective jumps for path-dependent excitation. Event j's excitation at
// time t is g_j(c) = exp(alpha * logn_j + gamma * h(c_j(t))), where c_j(t)
// starts at c0_j and gains one at each later linked event (updates from
// connectedness_updates_cpp()). Its contribution K g_j(c_j(t)) e^{-beta (t - t_j)}
// equals that of a jump g_j(c0_j) at t_j plus, at each update time t_k, a jump
// (g_j(c + 1) - g_j(c)) e^{-beta (t_k - t_j)}. G_k collects the jumps at t_k.
// Returns G, dG (columns alpha if size, gamma) and dGb = dG / dbeta.
// [[Rcpp::export]]
List path_jumps_cpp(NumericVector t, NumericVector logn, IntegerVector c0,
                    IntegerVector upd_ptr, IntegerVector upd_j, IntegerVector upd_c,
                    double alpha, double gamma, bool size, bool log1p_scale, double beta) {
  const int n = t.size();
  const int p = size ? 2 : 1;
  const int gcol = size ? 1 : 0;
  NumericVector G(n), dGb(n);
  NumericMatrix dG(n, p);
  auto h = [log1p_scale](double c) { return log1p_scale ? std::log1p(c) : c; };
  for (int k = 0; k < n; ++k) {
    const double a = size ? alpha * logn[k] : 0.0;
    const double g0 = std::exp(a + gamma * h(c0[k]));
    G[k] = g0;
    if (size) dG(k, 0) = logn[k] * g0;
    dG(k, gcol) = h(c0[k]) * g0;
    for (int u = upd_ptr[k]; u < upd_ptr[k + 1]; ++u) {
      const int j = upd_j[u];
      const double c = upd_c[u];
      const double aj = size ? alpha * logn[j] : 0.0;
      const double lag = t[k] - t[j];
      const double disc = std::exp(-beta * lag);
      const double g_new = std::exp(aj + gamma * h(c + 1.0));
      const double g_old = std::exp(aj + gamma * h(c));
      const double d = (g_new - g_old) * disc;
      G[k] += d;
      if (size) dG(k, 0) += logn[j] * d;
      dG(k, gcol) += (h(c + 1.0) * g_new - h(c) * g_old) * disc;
      dGb[k] -= lag * d;
    }
  }
  return List::create(_["G"] = G, _["dG"] = dG, _["dGb"] = dGb);
}
