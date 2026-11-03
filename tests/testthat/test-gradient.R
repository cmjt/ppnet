# Analytic gradients are checked against central differences of the
# likelihood itself, on the natural scale.

num_grad <- function(f, p, nms, h = 1e-5) {
  vapply(nms, function(nm) {
    up <- p; dn <- p
    up[[nm]] <- p[[nm]] + h
    dn[[nm]] <- p[[nm]] - h
    (f(up) - f(dn)) / (2 * h)
  }, numeric(1))
}

check_grad <- function(ll_fn, p, tol = 1e-5) {
  g <- attr(ll_fn(p, gradient = TRUE), "gradient")
  expect_false(is.null(g))
  expect_setequal(names(g), names(p))
  expect_equal(g[names(p)], num_grad(ll_fn, p, names(p)), tolerance = tol)
}

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
