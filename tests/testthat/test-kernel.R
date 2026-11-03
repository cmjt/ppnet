# Two-component kernel: the slow component carries the mark weights g_i, the
# fast one has unit jumps.

brute_exp2_time_ll <- function(t, p, T_end, g = rep(1, length(t))) {
  lam <- vapply(seq_along(t), function(i) {
    j <- seq_len(i - 1L)
    p$mu + p$K * sum(g[j] * exp(-p$beta * (t[i] - t[j]))) +
      p$K_fast * sum(exp(-p$beta_fast * (t[i] - t[j])))
  }, numeric(1))
  sum(log(lam)) - p$mu * T_end -
    (p$K / p$beta) * sum(g * (1 - exp(-p$beta * (T_end - t)))) -
    (p$K_fast / p$beta_fast) * sum(1 - exp(-p$beta_fast * (T_end - t)))
}

p2 <- list(mu = 1, K = 0.3, beta = 0.5, K_fast = 2, beta_fast = 8,
           beta_edges = 0.3, lambda_new = 0.9)

test_that("exp2 loglik matches a brute-force calculation", {
  set.seed(81)
  ev <- sim_hawkesNet(p2, T_end = 15, mark_type = "ba_bip", kernel = "exp2")$ev
  skip_if(nrow(ev$times) < 10L)
  t <- ev$times$t
  mark <- loglik(ev, p2, mark_type = "ba_bip", T_end = 15) -
    brute_exp2_time_ll(t, modifyList(p2, list(K_fast = 0)), 15)

  expect_equal(loglik(ev, p2, mark_type = "ba_bip", T_end = 15, kernel = "exp2"),
               brute_exp2_time_ll(t, p2, 15) + mark, tolerance = 1e-8)
})

test_that("excitation scales only the slow component", {
  local_min_one_perp()
  set.seed(82)
  exc <- excitation_spec(size = TRUE, connectedness = TRUE)
  p <- c(p2, alpha = 0.5, gamma = 0.3)
  ev <- sim_hawkesNet(p, T_end = 15, mark_type = "ba_bip", kernel = "exp2", excitation = exc)$ev
  skip_if(nrow(ev$times) < 10L)
  t <- ev$times$t
  s <- event_summary(ev)
  g <- s$n_perps^0.5 * (1 + s$connectedness)^0.3
  mark <- loglik(ev, p2, mark_type = "ba_bip", T_end = 15) -
    brute_exp2_time_ll(t, modifyList(p2, list(K_fast = 0)), 15)

  expect_equal(loglik(ev, p, mark_type = "ba_bip", T_end = 15, kernel = "exp2", excitation = exc),
               brute_exp2_time_ll(t, p, 15, g = g) + mark, tolerance = 1e-8)
})

test_that("a vanishing fast component gives the one-component model", {
  set.seed(83)
  ev <- sim_hawkesNet(p2, T_end = 12, mark_type = "ba_bip")$ev
  expect_equal(loglik(ev, modifyList(p2, list(K_fast = 1e-300)), mark_type = "ba_bip",
                      T_end = 12, kernel = "exp2"),
               loglik(ev, p2, mark_type = "ba_bip", T_end = 12), tolerance = 1e-10)
})

test_that("exp2 needs its parameters and beta_fast > beta", {
  set.seed(84)
  ev <- sim_hawkesNet(p2, T_end = 8, mark_type = "ba_bip")$ev
  expect_equal(loglik(ev, modifyList(p2, list(beta_fast = 0.4)), mark_type = "ba_bip",
                      kernel = "exp2"), -Inf)
  expect_error(loglik(ev, p2[-4], mark_type = "ba_bip", kernel = "exp2"), "K_fast and beta_fast")
  expect_error(fit_hawkesNet(ev, modifyList(p2, list(beta_fast = 0.1)), mark_type = "ba_bip",
                             kernel = "exp2"), "beta_fast > beta")
  expect_error(sim_hawkesNet(modifyList(p2, list(beta_fast = 0.1)), T_end = 2,
                             mark_type = "ba_bip", kernel = "exp2"), "beta_fast>beta")
})

test_that("fit_hawkesNet recovers both components", {
  set.seed(85)
  p <- list(mu = 2, K = 0.15, beta = 0.4, K_fast = 3, beta_fast = 12,
            m = 1, beta_edges = 0.3, lambda_new = 1)
  att <- attachment_spec(degree = FALSE)
  ev <- sim_hawkesNet(p, T_end = 150, mark_type = "ba_bip", kernel = "exp2", attachment = att)$ev
  init <- modifyList(p, list(mu = 1, K = 0.1, beta = 1, K_fast = 1, beta_fast = 5))

  fit <- suppressWarnings(fit_hawkesNet(
    ev, init, fixed_params = c("m", "beta_edges", "lambda_new"), mark_type = "ba_bip",
    T_end = 150, attachment = att, kernel = "exp2", hessian = TRUE
  ))
  nm <- c("mu", "K", "beta", "K_fast", "beta_fast")
  se <- sqrt(diag(solve(fit$hessian)))[nm]
  est <- log(unlist(fit$par[nm]))
  expect_true(all(abs(est - log(unlist(p[nm]))) < 4 * se))
})
