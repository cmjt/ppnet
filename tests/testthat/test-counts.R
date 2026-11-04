# New-perp counts given old-perp counts, checked against direct calls to the
# count densities.

brute_count_ll <- function(n_new, n_old, p, family, min_perps) {
  d <- function(y, mu) {
    if (family == "poisson") dpois(y, mu, log = TRUE) else dnbinom(y, size = p$size_new, mu = mu, log = TRUE)
  }
  shifted <- min_perps == 1 & n_old == 0
  sum(ifelse(shifted, d(n_new - 1, p$lambda_new0), d(n_new, p$lambda_new)))
}

p_cnt <- list(mu = 1, K = 0.5, beta = 1, m = 0.8, beta_edges = 0.4, delta = 0.5,
              lambda_new = 0.6, lambda_new0 = 0.8, size_new = 1.5)

test_that("count models match the densities directly", {
  set.seed(91)
  sp <- attachment_spec(degree = TRUE, new_perps = "negbin", min_perps = 1)
  ev <- sim_hawkesNet(p_cnt, T_end = 20, mark_type = "ba_bip", attachment = sp)$ev
  s <- event_summary(ev)
  expect_true(all(s$n_perps >= 1))
  expect_true(any(s$n_old == 0) && any(s$n_old > 0))

  base <- attachment_spec(degree = TRUE)
  attach_only <- compile_logpmf_attach(ev, base)(p_cnt) -
    sum(dpois(s$n_new, p_cnt$lambda_new, log = TRUE))
  for (family in c("poisson", "negbin")) {
    spec <- attachment_spec(degree = TRUE, new_perps = family, min_perps = 1)
    expect_equal(compile_logpmf_attach(ev, spec)(p_cnt),
                 attach_only + brute_count_ll(s$n_new, s$n_old, p_cnt, family, 1),
                 tolerance = 1e-10)
  }
})

test_that("min_perps = 1 rejects events with no perps", {
  set.seed(92)
  sp0 <- attachment_spec(degree = TRUE)
  ev <- sim_hawkesNet(modifyList(p_cnt, list(lambda_new = 0.2)), T_end = 20,
                      mark_type = "ba_bip", attachment = sp0)$ev
  skip_if(all(event_summary(ev)$n_perps > 0))
  expect_equal(compile_logpmf_attach(ev, attachment_spec(min_perps = 1))(p_cnt), -Inf)
})

test_that("negative binomial counts approach Poisson as size grows", {
  set.seed(93)
  sp <- attachment_spec(degree = TRUE, min_perps = 1)
  ev <- sim_hawkesNet(p_cnt, T_end = 15, mark_type = "ba_bip", attachment = sp)$ev
  nb <- attachment_spec(degree = TRUE, new_perps = "negbin", min_perps = 1)
  expect_equal(compile_logpmf_attach(ev, nb)(modifyList(p_cnt, list(size_new = 1e8))),
               compile_logpmf_attach(ev, sp)(p_cnt), tolerance = 1e-6)
})

test_that("fit_hawkesNet recovers the count parameters", {
  set.seed(94)
  sp <- attachment_spec(degree = FALSE, new_perps = "negbin", min_perps = 1)
  p <- list(mu = 3, K = 0.3, beta = 1, m = 0.8, beta_edges = 0.3,
            lambda_new = 0.5, lambda_new0 = 1.2, size_new = 2)
  ev <- sim_hawkesNet(p, T_end = 200, mark_type = "ba_bip", attachment = sp)$ev
  init <- modifyList(p, list(lambda_new = 1, lambda_new0 = 1, size_new = 1))
  fit <- suppressWarnings(fit_hawkesNet(ev, init, fixed_params = c("mu", "K", "beta"),
                                        mark_type = "ba_bip", T_end = 200,
                                        attachment = sp, hessian = TRUE))
  nm <- c("lambda_new", "lambda_new0", "size_new")
  se <- sqrt(diag(solve(fit$hessian)))[nm]
  expect_true(all(abs(log(unlist(fit$par[nm])) - log(unlist(p[nm]))) < 4 * se))
})

test_that("a random effect makes perps differ in how often they return", {
  # Same expected weights, but a random effect spreads attendance counts out.
  sp <- attachment_spec(degree = FALSE)
  p <- list(mu = 4, K = 0.2, beta = 1, m = 1, beta_edges = 0.05, lambda_new = 0.7)
  attendance <- function(shape) {
    set.seed(95)
    ev <- sim_hawkesNet(p, T_end = 100, mark_type = "ba_bip", attachment = sp,
                        random_effect_shape = shape)$ev
    as.numeric(table(ev$edges$j))
  }
  a0 <- attendance(Inf)
  a1 <- attendance(0.3)
  expect_gt(var(a1) / mean(a1), 1.5 * var(a0) / mean(a0))
})
