test_that("sim_perp_events gives a valid bipartite history", {
  set.seed(11)
  s <- sim_perp_events(T_end = 104, random_effect_shape = 1, K_self = 0.3, excite_on = "all",
                       extra_mean = 1, size_link = 1)
  ev <- s$ev
  es <- event_summary(ev)
  expect_true(nrow(es) > 50)
  expect_true(all(es$n_perps >= 1))
  expect_true(all(ev$times$t >= 0 & ev$times$t <= 104))
  # No perp twice at one event
  expect_false(any(duplicated(ev$edges[, c("i", "j")])))
  # Perp table matches the edges
  seen <- s$perps[!is.na(s$perps$id), ]
  att <- table(ev$edges$j)
  expect_equal(as.integer(att[seen$id]), seen$n_events)
  # The fitted models accept it
  spec <- attachment_spec(degree = TRUE, new_perps = "negbin", min_perps = 1)
  p <- list(mu = 1, K = 0.2, beta = 0.5, alpha = 0, gamma = 0, m = 0.5, beta_edges = 0.05,
            delta = 0.5, lambda_new = 0.6, lambda_new0 = 0.5, size_new = 2)
  ll <- loglik(ev, p, mark_type = "ba_bip", T_end = 104, attachment = spec,
               excitation = excitation_spec(size = TRUE, connectedness = TRUE))
  expect_true(is.finite(ll))
})

test_that("sim_perp_events without co-offenders gives one perp per event", {
  set.seed(12)
  s <- sim_perp_events(T_end = 52, extra_mean = 0)
  expect_true(all(event_summary(s$ev)$n_perps == 1))
})

test_that("sim_perp_events matches the expected event count", {
  # No random effect, self-excitation or burn-in, careers that never end: the
  # active count grows as entry_rate * t, so E[N] = rate * entry_rate * T^2 / 2.
  set.seed(13)
  n <- replicate(40, nrow(sim_perp_events(T_end = 100, entry_rate = 1, career_mean = Inf,
                                          rate = 0.05, burn_in = 0, extra_mean = 0)$ev$times))
  expect_equal(mean(n), 250, tolerance = 0.06)
})

test_that("sim_perp_events self-excitation scales the rate by 1 / (1 - K_self)", {
  # Stationary active set (long burn-in): rate per active perp is
  # rate / (1 - K_self) on average.
  set.seed(14)
  args <- list(T_end = 200, entry_rate = 1, career_mean = 50, rate = 0.05, extra_mean = 0)
  n0 <- replicate(20, nrow(do.call(sim_perp_events, args)$ev$times))
  n1 <- replicate(20, nrow(do.call(sim_perp_events, c(args, K_self = 0.5))$ev$times))
  expect_equal(mean(n0), 0.05 * 50 * 200, tolerance = 0.06)
  expect_equal(mean(n1) / mean(n0), 2, tolerance = 0.1)
})
