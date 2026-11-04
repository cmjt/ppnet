# Analytic gradients are checked against central differences of the
# likelihood itself, on the natural scale (num_grad() and check_grad() are in
# helper-brute.R).

test_that("BA gradient matches central differences", {
  set.seed(51)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3)
  ev <- sim_hawkesNet(p, T_end = 15, mark_type = "ba")$ev
  check_grad(.build_loglik_fn(ev, T_end = 15, mark_type = "ba"),
             list(mu = 0.8, K = 0.6, beta = 1.4, beta_edges = 0.5))
})

test_that("BA-bip gradient matches central differences", {
  set.seed(52)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
  ev <- sim_hawkesNet(p, T_end = 15, mark_type = "ba_bip")$ev
  check_grad(.build_loglik_fn(ev, T_end = 15, mark_type = "ba_bip"),
             list(mu = 1.3, K = 0.3, beta = 0.7, beta_edges = 0.1, lambda_new = 1.2))
})

test_that("baseline gradient matches central differences", {
  set.seed(53)
  sp <- baseline_spec(trend = TRUE, periods = c(daily = 1, weekly = 7),
                      harmonics = c(2, 1), trend_scale = 10)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7,
            b_trend = 0.4, b_daily_cos1 = 0.6, b_daily_sin1 = -0.3,
            b_daily_cos2 = 0.2, b_daily_sin2 = 0.1,
            b_weekly_cos1 = -0.5, b_weekly_sin1 = 0.25)
  ev <- sim_hawkesNet(p, T_end = 15, mark_type = "ba_bip", baseline = sp)$ev
  check_grad(.build_loglik_fn(ev, T_end = 16, mark_type = "ba_bip", baseline = sp),
             modifyList(p, list(b_trend = -0.2, b_weekly_sin1 = 0.5)))
})

test_that("a window starting after 0 has the right gradient", {
  set.seed(54)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
  ev <- sim_hawkesNet(p, T0 = 3, T_end = 15, mark_type = "ba_bip")$ev
  check_grad(.build_loglik_fn(ev, T0 = 3, T_end = 15, mark_type = "ba_bip"), p)
})

test_that("models without an analytic mark gradient return none", {
  set.seed(55)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, node_lambda = 1.2, CS_edges = 0.4)
  ev <- sim_hawkesNet(p, T_end = 6, mark_type = "cs", formula_rhs = "edges")$ev
  ll_fn <- .build_loglik_fn(ev, T_end = 6, mark_type = "cs", formula_rhs = "edges")
  expect_null(attr(ll_fn(p, gradient = TRUE), "gradient"))
  expect_equal(as.numeric(ll_fn(p, gradient = TRUE)), ll_fn(p))
})

test_that("the gradient is near zero at the fitted optimum", {
  set.seed(56)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
  ev <- sim_hawkesNet(p, T_end = 30, mark_type = "ba_bip")$ev
  fit <- suppressWarnings(fit_hawkesNet(ev, p, mark_type = "ba_bip", T_end = 30,
                                        control = list(reltol = 1e-12)))
  g <- attr(.build_loglik_fn(ev, T_end = 30, mark_type = "ba_bip")(fit$par, gradient = TRUE),
            "gradient")
  expect_lt(max(abs(g * unlist(fit$par[names(g)]))), 1e-3)
  expect_lt(fit$fit$counts[["function"]], 200)
})

test_that("excitation gradient matches central differences", {
  local_min_one_perp()
  set.seed(57)
  sp <- excitation_spec(size = TRUE, connectedness = TRUE, connectedness_window = 4)
  bl <- baseline_spec(periods = c(weekly = 7))
  p <- list(mu = 1, K = 0.4, beta = 1, beta_edges = 0.3, lambda_new = 0.9,
            alpha = 0.5, gamma = 0.3, b_weekly_cos1 = 0.3, b_weekly_sin1 = -0.2)
  ev <- sim_hawkesNet(p, T_end = 15, mark_type = "ba_bip", excitation = sp, baseline = bl)$ev
  check_grad(.build_loglik_fn(ev, T_end = 16, mark_type = "ba_bip", excitation = sp, baseline = bl),
             modifyList(p, list(alpha = 0.2, gamma = -0.1)))
})

test_that("attachment gradient matches central differences", {
  set.seed(58)
  sp <- attachment_spec(degree = TRUE, coffenders = TRUE, coffender_window = 3)
  p <- list(mu = 1, K = 0.5, beta = 1, m = 0.8, beta_edges = 0.4, delta = 0.5,
            theta = 0.7, lambda_new = 0.9)
  ev <- sim_hawkesNet(p, T_end = 20, mark_type = "ba_bip", attachment = sp)$ev
  check_grad(.build_loglik_fn(ev, T_end = 20, mark_type = "ba_bip", attachment = sp),
             modifyList(p, list(m = 1.3, beta_edges = 0.2, delta = 2, theta = -0.3)))
  check_grad(.build_loglik_fn(ev, T_end = 20, mark_type = "ba_bip",
                              attachment = attachment_spec(degree = FALSE)),
             p[c("mu", "K", "beta", "m", "beta_edges", "lambda_new")])
})

test_that("attachment gradient with 3-paths matches central differences", {
  set.seed(60)
  sp <- attachment_spec(degree = TRUE, coffenders = TRUE, paths = TRUE, coffender_window = 3)
  p <- list(mu = 1, K = 0.5, beta = 1, m = 0.8, beta_edges = 0.4, delta = 0.5,
            theta = 0.7, theta_paths = -0.3, lambda_new = 0.9)
  ev <- sim_hawkesNet(p, T_end = 20, mark_type = "ba_bip", attachment = sp)$ev
  check_grad(.build_loglik_fn(ev, T_end = 20, mark_type = "ba_bip", attachment = sp),
             modifyList(p, list(theta = -0.2, theta_paths = 0.4)))
  sp_re <- attachment_spec(degree = TRUE, coffenders = TRUE, paths = TRUE, coffender_window = 3,
                           random_effect = TRUE)
  check_grad(.build_loglik_fn(ev, T_end = 20, mark_type = "ba_bip", attachment = sp_re),
             c(p, random_effect_var = 0.7))
})

test_that("the time likelihood stays bounded as beta goes to 0", {
  # With 1 - exp(-beta * tau) computed naively, a tiny beta drops the
  # compensator's K term and a huge K then looks ever better.
  set.seed(59)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
  ev <- sim_hawkesNet(p, T_end = 15, mark_type = "ba_bip")$ev
  ll_fn <- .build_loglik_fn(ev, T_end = 15, mark_type = "ba_bip")
  t <- ev$times$t

  # As beta -> 0 the kernel is flat: lambda_i = mu + K (i - 1), compensator
  # mu T + K sum(T - t_i).
  limit <- function(K) sum(log(1 + K * (seq_along(t) - 1))) - 15 - K * sum(15 - t)
  mark <- ll_fn(p) - brute_time_ll(t, 1, 0.5, 1, 15)
  for (beta in c(1e-6, 1e-12, 1e-200)) {
    q <- modifyList(p, list(beta = beta, K = 2 * beta^0))
    expect_equal(as.numeric(ll_fn(q)), limit(2) + mark, tolerance = 1e-5)
  }
  # beta * tau below 1e-4 uses the series for the beta derivative; step
  # relative to beta so the difference stays inside (0, Inf).
  q <- modifyList(p, list(beta = 1e-6))
  g <- attr(ll_fn(q, gradient = TRUE), "gradient")
  expect_equal(g[["beta"]], num_grad(ll_fn, q, "beta", h = 1e-9)[[1]], tolerance = 1e-5)
})

test_that("exp2 gradient matches central differences", {
  local_min_one_perp()
  set.seed(60)
  exc <- excitation_spec(size = TRUE, connectedness = TRUE)
  p <- list(mu = 1, K = 0.3, beta = 0.5, K_fast = 2, beta_fast = 8,
            beta_edges = 0.3, lambda_new = 0.9, alpha = 0.4, gamma = 0.2)
  ev <- sim_hawkesNet(p, T_end = 15, mark_type = "ba_bip", kernel = "exp2", excitation = exc)$ev
  check_grad(.build_loglik_fn(ev, T_end = 16, mark_type = "ba_bip", kernel = "exp2",
                              excitation = exc),
             modifyList(p, list(K_fast = 1.2, beta_fast = 5)))
})

test_that("count-model gradient matches central differences", {
  set.seed(61)
  sp <- attachment_spec(degree = TRUE, coffenders = TRUE, new_perps = "negbin", min_perps = 1)
  p <- list(mu = 1, K = 0.5, beta = 1, m = 0.8, beta_edges = 0.4, delta = 0.5, theta = 0.3,
            lambda_new = 0.6, lambda_new0 = 0.8, size_new = 1.5)
  ev <- sim_hawkesNet(p, T_end = 20, mark_type = "ba_bip", attachment = sp)$ev
  check_grad(.build_loglik_fn(ev, T_end = 20, mark_type = "ba_bip", attachment = sp), p)
})
