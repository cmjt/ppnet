# diag_temporal() is checked against a directly coded intensity, integrated
# numerically over each gap.

test_that("diag_temporal matches a direct calculation with all temporal options", {
  local_min_one_perp()
  set.seed(111)
  sp <- baseline_spec(trend = TRUE, periods = c(daily = 1), trend_scale = 10)
  exc <- excitation_spec(size = TRUE, connectedness = TRUE)
  p <- list(mu = 1, K = 0.3, beta = 0.5, K_fast = 2, beta_fast = 6, beta_edges = 0.3,
            lambda_new = 0.9, alpha = 0.4, gamma = 0.3,
            b_trend = 0.3, b_daily_cos1 = 0.5, b_daily_sin1 = -0.2)
  ev <- sim_hawkesNet(p, T_end = 12, mark_type = "ba_bip", baseline = sp,
                      excitation = exc, kernel = "exp2")$ev
  t <- ev$times$t
  skip_if(length(t) < 10L)

  s <- event_summary(ev)
  g <- s$n_perps^p$alpha * (1 + s$connectedness)^p$gamma
  lam <- function(u) vapply(u, function(x) {
    j <- t < x
    p$mu * exp(p$b_trend * x / 10 + p$b_daily_cos1 * cos(2 * pi * x) + p$b_daily_sin1 * sin(2 * pi * x)) +
      p$K * sum(g[j] * exp(-p$beta * (x - t[j]))) + p$K_fast * sum(exp(-p$beta_fast * (x - t[j])))
  }, numeric(1))

  d <- diag_temporal(ev, p, T_end = 13, baseline = sp, excitation = exc, kernel = "exp2")
  expect_equal(d$lambda, lam(t), tolerance = 1e-10)
  gaps <- mapply(function(a, b) stats::integrate(lam, a, b, rel.tol = 1e-10)$value,
                 c(0, t[-length(t)]), t)
  expect_equal(d$comp_inc, gaps, tolerance = 1e-7)
  expect_equal(d$meta$comp_tail, stats::integrate(lam, max(t), 13, rel.tol = 1e-10)$value,
               tolerance = 1e-7)

  # Same total as the likelihood's time part
  ll_fn <- .build_loglik_fn(ev, T_end = 13, mark_type = "ba_bip", baseline = sp,
                            excitation = exc, kernel = "exp2")
  mark <- ll_fn(modifyList(p, list(K = 0, K_fast = 1e-300, alpha = 0, gamma = 0))) -
    (sum(log(p$mu * exp(drop(.baseline_design(sp, t) %*% c(0.3, 0.5, -0.2))))) -
       p$mu * .baseline_integral(sp, c(0.3, 0.5, -0.2), 0, 13))
  expect_equal(as.numeric(ll_fn(p)), sum(log(d$lambda)) - sum(d$comp_inc) - d$meta$comp_tail + mark,
               tolerance = 1e-7)
})

test_that("rescaled times are Exp(1) under the true model", {
  set.seed(112)
  sp <- baseline_spec(periods = c(daily = 1))
  p <- list(mu = 3, K = 0.3, beta = 0.5, K_fast = 2, beta_fast = 6, m = 0.5, beta_edges = 0.3,
            lambda_new = 0.6, lambda_new0 = 0.5, b_daily_cos1 = 0.7, b_daily_sin1 = 0)
  att <- attachment_spec(degree = FALSE, min_perps = 1)
  ev <- sim_hawkesNet(p, T_end = 150, mark_type = "ba_bip", baseline = sp, attachment = att,
                      kernel = "exp2", return_net = FALSE)$ev
  d <- diag_temporal(ev, p, baseline = sp, kernel = "exp2")
  expect_gt(suppressWarnings(ks.test(d$tau, "pexp"))$p.value, 0.001)

  # A constant-rate, one-kernel model misses the daily cycle and fast clustering
  d_bad <- diag_temporal(ev, p)
  expect_lt(suppressWarnings(ks.test(d_bad$tau, "pexp"))$p.value, 0.001)
})

test_that("plot(fit) uses the fit's temporal settings", {
  set.seed(113)
  sp <- baseline_spec(periods = c(daily = 1))
  p <- list(mu = 2, K = 0.3, beta = 1, beta_edges = 0.3, lambda_new = 0.7,
            b_daily_cos1 = 0.5, b_daily_sin1 = 0)
  ev <- sim_hawkesNet(p, T_end = 20, mark_type = "ba_bip", baseline = sp)$ev
  fit <- suppressWarnings(fit_hawkesNet(ev, p, fixed_params = c("beta_edges", "lambda_new"),
                                        mark_type = "ba_bip", T_end = 25, baseline = sp))
  expect_identical(fit$model$baseline, sp)
  expect_equal(fit$model$T_end, 25)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  d <- plot(fit, type = "tau_qq")
  expect_equal(d$tau, diag_temporal(ev, fit$par, T_end = 25, baseline = sp)$tau)
})
