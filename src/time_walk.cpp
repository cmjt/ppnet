// time_walk.cpp
//
// Log-likelihood of the event times under an exponential Hawkes kernel with
// mark-dependent jumps,
//   lambda(t) = mu(t) + K * sum_{t_j < t} g_j * exp(-beta * (t - t_j)),
//   g_j = exp(z_j' a),
// excluding the baseline integral, which the caller adds from its quadrature.
// The recursions are, with S+ and U+ the values just after the jump at t_{i-1},
//   S_i   = exp(-beta * dt) * S+_{i-1}                 (S+ = S + g)
//   R_i   = exp(-beta * dt) * (R_{i-1} + dt * S+_{i-1})
//   U_m,i = exp(-beta * dt) * U+_m,{i-1}               (U+ = U + z_m * g)
// where R_i = -dS_i / dbeta and U_m,i = dS_i / da_m.

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List time_ll_cpp(NumericVector t, NumericVector mu_ev, NumericVector g,
                 NumericMatrix Z, double K, double beta, double T_end) {
  const int n = t.size();
  const int p = Z.ncol();
  NumericVector inv_lambda(n), log_lambda(n), da(p);
  std::vector<double> U(p, 0.0), dcomp_da(p, 0.0);

  double S = 0.0, R = 0.0, sum_log = 0.0, dK = 0.0, dbeta = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(g[i])) return List::create(_["ll"] = R_NegInf);
    if (i > 0) {
      const double dt = t[i] - t[i - 1];
      if (dt < 0) stop("times must be increasing");
      const double e = std::exp(-beta * dt);
      R = e * (R + dt * S);
      S = e * S;
      for (int m = 0; m < p; ++m) U[m] *= e;
    }
    const double lam = mu_ev[i] + K * S;
    if (!std::isfinite(lam) || lam <= 0.0) {
      return List::create(_["ll"] = R_NegInf);
    }
    log_lambda[i] = std::log(lam);
    inv_lambda[i] = 1.0 / lam;
    sum_log += log_lambda[i];
    dK += S / lam;
    dbeta -= K * R / lam;
    for (int m = 0; m < p; ++m) da[m] += K * U[m] / lam;

    // The jump at t_i adds to S and U but not to R (zero lag).
    S += g[i];
    for (int m = 0; m < p; ++m) U[m] += Z(i, m) * g[i];
  }

  // Excitation part of the compensator: (K / beta) * sum_i g_i (1 - exp(-beta * tau_i)).
  double comp = 0.0, dcomp_dbeta = 0.0;
  for (int i = 0; i < n; ++i) {
    const double tau = T_end - t[i];
    const double e = std::exp(-beta * tau);
    comp += g[i] * (1.0 - e);
    dcomp_dbeta += g[i] * (-(1.0 - e) / (beta * beta) + tau * e / beta);
    for (int m = 0; m < p; ++m) dcomp_da[m] += Z(i, m) * g[i] * (1.0 - e);
  }
  for (int m = 0; m < p; ++m) da[m] -= (K / beta) * dcomp_da[m];

  return List::create(
    _["ll"] = sum_log - (K / beta) * comp,
    _["dK"] = dK - comp / beta,
    _["dbeta"] = dbeta - K * dcomp_dbeta,
    _["da"] = da,
    _["inv_lambda"] = inv_lambda,
    _["log_lambda"] = log_lambda
  );
}
