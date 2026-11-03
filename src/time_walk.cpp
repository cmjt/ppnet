// time_walk.cpp
//
// Log-likelihood of the event times under an exponential Hawkes kernel with
// mark-dependent jumps and an optional fast component with unit jumps,
//   lambda(t) = mu(t) + K * sum_{t_j < t} g_j * exp(-beta * (t - t_j))
//                     + K2 * sum_{t_j < t} exp(-beta2 * (t - t_j)),
//   g_j = exp(z_j' a),
// excluding the baseline integral, which the caller adds from its quadrature.
// The recursions are, with S+ and U+ the values just after the jump at t_{i-1},
//   S_i   = exp(-beta * dt) * S+_{i-1}                 (S+ = S + g)
//   R_i   = exp(-beta * dt) * (R_{i-1} + dt * S+_{i-1})
//   U_m,i = exp(-beta * dt) * U+_m,{i-1}               (U+ = U + z_m * g)
// where R_i = -dS_i / dbeta and U_m,i = dS_i / da_m; S2 and R2 are the same
// for the fast component with g = 1. K2 = 0 drops it.

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
List time_ll_cpp(NumericVector t, NumericVector mu_ev, NumericVector g,
                 NumericMatrix Z, double K, double beta, double T_end,
                 double K2 = 0.0, double beta2 = 1.0) {
  const int n = t.size();
  const int p = Z.ncol();
  NumericVector inv_lambda(n), log_lambda(n), da(p);
  std::vector<double> U(p, 0.0), dcomp_da(p, 0.0);

  double S = 0.0, R = 0.0, sum_log = 0.0, dK = 0.0, dbeta = 0.0;
  double S2 = 0.0, R2 = 0.0, dK2 = 0.0, dbeta2 = 0.0;
  const bool fast = K2 > 0.0;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(g[i])) return List::create(_["ll"] = R_NegInf);
    if (i > 0) {
      const double dt = t[i] - t[i - 1];
      if (dt < 0) stop("times must be increasing");
      const double e = std::exp(-beta * dt);
      R = e * (R + dt * S);
      S = e * S;
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
    dbeta -= K * R / lam;
    for (int m = 0; m < p; ++m) da[m] += K * U[m] / lam;
    dK2 += S2 / lam;
    dbeta2 -= K2 * R2 / lam;

    // The jump at t_i adds to S and U but not to R (zero lag).
    S += g[i];
    for (int m = 0; m < p; ++m) U[m] += Z(i, m) * g[i];
    if (fast) S2 += 1.0;
  }

  // Excitation parts of the compensator.
  double comp, dcomp_dbeta, comp2 = 0.0, dcomp2_dbeta2 = 0.0;
  exp_compensator(t, g, beta, T_end, comp, dcomp_dbeta);
  if (fast) {
    exp_compensator(t, NumericVector(n, 1.0), beta2, T_end, comp2, dcomp2_dbeta2);
  }
  for (int i = 0; i < n; ++i) {
    const double one_minus_e = -std::expm1(-beta * (T_end - t[i]));
    for (int m = 0; m < p; ++m) dcomp_da[m] += Z(i, m) * g[i] * one_minus_e;
  }
  for (int m = 0; m < p; ++m) da[m] -= (K / beta) * dcomp_da[m];

  return List::create(
    _["ll"] = sum_log - (K / beta) * comp - (fast ? (K2 / beta2) * comp2 : 0.0),
    _["dK"] = dK - comp / beta,
    _["dbeta"] = dbeta - K * dcomp_dbeta,
    _["da"] = da,
    _["dK2"] = dK2 - (fast ? comp2 / beta2 : 0.0),
    _["dbeta2"] = dbeta2 - K2 * dcomp2_dbeta2,
    _["inv_lambda"] = inv_lambda,
    _["log_lambda"] = log_lambda
  );
}
