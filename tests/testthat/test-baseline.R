# The baseline integral is checked against stats::integrate(), and event-time
# log-intensities against a direct sum, without using the cached evaluators.

brute_baseline_time_ll <- function(t, p, spec, T0, T_end) {
  b <- unlist(p[baseline_coef_names(spec)])
  mu_fun <- function(s) p$mu * exp(drop(.baseline_design(spec, s) %*% b))
  lam <- vapply(seq_along(t), function(i) {
    mu_fun(t[i]) + p$K * sum(exp(-p$beta * (t[i] - t[seq_len(i - 1L)])))
  }, numeric(1))
  base_int <- stats::integrate(mu_fun, T0, T_end, subdivisions = 10000L,
                               rel.tol = 1e-12)$value
  sum(log(lam)) - base_int - (p$K / p$beta) * sum(1 - exp(-p$beta * (T_end - t)))
}

spec <- baseline_spec(trend = TRUE, periods = c(daily = 1, weekly = 7),
                      harmonics = c(2, 1), trend_scale = 10)
p_base <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7,
               b_trend = 0.4, b_daily_cos1 = 0.6, b_daily_sin1 = -0.3,
               b_daily_cos2 = 0.2, b_daily_sin2 = 0.1,
               b_weekly_cos1 = -0.5, b_weekly_sin1 = 0.25)

test_that("baseline coefficient names follow the design columns", {
  expect_equal(baseline_coef_names(spec), names(p_base)[-(1:5)])
  expect_equal(colnames(.baseline_design(spec, c(0, 1.3))), baseline_coef_names(spec))
  expect_length(baseline_coef_names(NULL), 0L)
})

test_that("baseline loglik matches a brute-force calculation", {
  set.seed(21)
  ev <- sim_hawkesNet(p_base, T_end = 15, mark_type = "ba_bip", baseline = spec)$ev
  skip_if(nrow(ev$times) < 5L)

  # Mark part is unchanged by the baseline, so take it from the constant fit.
  p0 <- p_base
  p0[baseline_coef_names(spec)] <- 0
  mark <- loglik(ev, p_base[1:5], mark_type = "ba_bip", T_end = 15) -
    brute_baseline_time_ll(ev$times$t, p0, spec, 0, 15)

  expect_equal(
    loglik(ev, p_base, mark_type = "ba_bip", T_end = 15, baseline = spec),
    brute_baseline_time_ll(ev$times$t, p_base, spec, 0, 15) + mark,
    tolerance = 1e-8
  )
})

test_that("a zero-coefficient baseline equals the constant baseline", {
  set.seed(22)
  ev <- sim_hawkesNet(p_base[1:5], T_end = 12, mark_type = "ba_bip")$ev
  p0 <- p_base
  p0[baseline_coef_names(spec)] <- 0
  expect_equal(loglik(ev, p0, mark_type = "ba_bip", T_end = 12, baseline = spec),
               loglik(ev, p_base[1:5], mark_type = "ba_bip", T_end = 12),
               tolerance = 1e-10)
})

test_that("missing baseline coefficients are an error", {
  set.seed(23)
  ev <- sim_hawkesNet(p_base[1:5], T_end = 8, mark_type = "ba_bip")$ev
  expect_error(loglik(ev, p_base[1:6], mark_type = "ba_bip", baseline = spec),
               "Missing baseline coefficients")
  expect_error(fit_hawkesNet(ev, p_base[1:6], mark_type = "ba_bip", baseline = spec),
               "Missing baseline coefficients")
})

test_that("the thinning bound holds over the window", {
  s <- seq(0, 15, by = 0.001)
  b <- unlist(p_base[baseline_coef_names(spec)])
  mu_s <- p_base$mu * exp(drop(.baseline_design(spec, s) %*% b))
  expect_true(all(mu_s <= .baseline_upper(spec, p_base, 0, 15)))
})

test_that("simulated event counts follow the baseline", {
  # With K = 0 the process is Poisson, so counts on the high and low halves of
  # the daily cycle should match their integrated rates.
  set.seed(24)
  sp <- baseline_spec(periods = c(daily = 1))
  p <- list(mu = 20, K = 0, beta = 1, beta_edges = 0.3, lambda_new = 0.7,
            b_daily_cos1 = 1, b_daily_sin1 = 0)
  t <- sim_hawkesNet(p, T_end = 30, mark_type = "ba_bip", baseline = sp)$ev$times$t
  high <- (t %% 1) < 0.25 | (t %% 1) > 0.75
  mu_fun <- function(s) 20 * exp(cos(2 * pi * s))
  expected_high <- 30 * 2 * stats::integrate(mu_fun, 0, 0.25)$value
  expected_total <- 30 * stats::integrate(mu_fun, 0, 1)$value
  expect_lt(abs(sum(high) - expected_high), 4 * sqrt(expected_high))
  expect_lt(abs(length(t) - expected_total), 4 * sqrt(expected_total))
})

test_that("fit_hawkesNet recovers baseline coefficients", {
  set.seed(25)
  sp <- baseline_spec(trend = TRUE, periods = c(weekly = 7), trend_scale = 100)
  p <- list(mu = 3, K = 0.4, beta = 2, beta_edges = 0.3, lambda_new = 0.7,
            b_trend = 0.5, b_weekly_cos1 = 0.6, b_weekly_sin1 = -0.4)
  ev <- sim_hawkesNet(p, T_end = 100, mark_type = "ba_bip", baseline = sp)$ev
  init <- modifyList(p, list(mu = 1, K = 0.2, beta = 1,
                             b_trend = 0, b_weekly_cos1 = 0, b_weekly_sin1 = 0))

  fit <- suppressWarnings(fit_hawkesNet(
    ev, init, fixed_params = c("beta_edges", "lambda_new"),
    mark_type = "ba_bip", T_end = 100, baseline = sp, hessian = TRUE
  ))
  est <- unlist(fit$par[c("b_trend", "b_weekly_cos1", "b_weekly_sin1")])
  se <- sqrt(diag(solve(fit$hessian)))[names(est)]
  expect_true(all(abs(est - unlist(p[names(est)])) < 4 * se))
})
