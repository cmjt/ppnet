test_that("compare_fits tabulates fits and tests nested models", {
  set.seed(101)
  p <- list(mu = 2, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
  ev <- sim_hawkesNet(p, T_end = 30, mark_type = "ba_bip")$ev
  sp <- baseline_spec(periods = c(daily = 1))
  b0 <- list(b_daily_cos1 = 0, b_daily_sin1 = 0)
  fixed <- c("beta_edges", "lambda_new")

  f0 <- suppressWarnings(fit_hawkesNet(ev, p, fixed_params = fixed, mark_type = "ba_bip",
                                       T_end = 30, hessian = TRUE))
  f1 <- suppressWarnings(fit_hawkesNet(ev, c(p, b0), fixed_params = fixed, mark_type = "ba_bip",
                                       T_end = 30, baseline = sp, hessian = TRUE))
  tab <- compare_fits(constant = f0, daily = f1)

  expect_equal(tab$model, c("constant", "daily"))
  expect_equal(tab$n_par, c(3, 5))
  expect_equal(tab$loglik, c(f0$loglik, f1$loglik))
  expect_equal(tab$AIC, -2 * tab$loglik + 2 * tab$n_par)
  expect_equal(tab$df, c(NA, 2))
  expect_equal(tab$LR[2], 2 * (f1$loglik - f0$loglik))
  expect_equal(tab$p[2], pchisq(tab$LR[2], 2, lower.tail = FALSE))
  expect_true(all(tab$hessian_ok))

  expect_warning(compare_fits(daily = f1, constant = f0), "increasing order")
  expect_error(compare_fits(f0, list()), "fit_hawkesNet")
})
