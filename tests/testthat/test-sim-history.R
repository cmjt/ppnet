sp_h <- attachment_spec(coffenders = TRUE, paths = TRUE, coffender_window = 6, episode_gap = 0.1,
                        recent = TRUE, new_perps = "negbin", min_perps = 1)
p_h <- list(mu = 2, K = 0.5, beta = 1.5, m = 0.6, beta_edges = 0.3, delta = 0.5, theta = 0.3,
            theta_paths = 0.1, theta_recent = 0.8, lambda_new = 0.6, lambda_new0 = 0.4, size_new = 2)

set.seed(301)
hist_h <- sim_hawkesNet(p_h, T_end = 25, mark_type = "ba_bip", attachment = sp_h, return_net = FALSE)$ev

test_that("runaway simulations stop at max_events or max_seconds and say so", {
  expect_warning(s <- sim_hawkesNet(p_h, T_end = 1e6, mark_type = "ba_bip", attachment = sp_h,
                                    max_events = 50, return_net = FALSE), "max_events")
  expect_true(s$exploded)
  expect_identical(s$stop, "max_events")
  expect_equal(s$events$n, 50L)
  expect_equal(nrow(s$ev$times), 50L)
  expect_lt(s$t_stop, 1e6)
  expect_warning(s <- sim_hawkesNet(p_h, T_end = 1e7, mark_type = "ba_bip", attachment = sp_h,
                                    max_seconds = 0.5, return_net = FALSE), "max_seconds")
  expect_identical(s$stop, "max_seconds")
  ok <- sim_hawkesNet(p_h, T_end = 5, mark_type = "ba_bip", attachment = sp_h, return_net = FALSE)
  expect_false(ok$exploded)
  expect_identical(ok$stop, "T_end")
})

test_that("a history is kept as the start of the output", {
  s <- sim_hawkesNet(p_h, T_end = 35, mark_type = "ba_bip", attachment = sp_h, history = hist_h,
                     return_net = FALSE)
  n_h <- nrow(hist_h$times)
  expect_equal(s$n_history, n_h)
  expect_equal(s$history_end, max(hist_h$times$t))
  expect_equal(s$ev$times$t[seq_len(n_h)], hist_h$times$t)
  expect_true(all(s$events$t > s$history_end))
  expect_equal(s$events$n, nrow(s$ev$times) - n_h)
  old <- s$ev$nodes[s$ev$nodes$event_id <= n_h, c("id", "role")]
  expect_equal(old$id, as.character(hist_h$nodes$id))
  expect_true(all(startsWith(s$ev$nodes$id[s$ev$nodes$event_id > n_h], "s")))
  # Later events link back to history perps.
  late <- s$ev$edges[s$ev$edges$event_id > n_h, ]
  expect_true(any(c(late$i, late$j) %in% hist_h$nodes$id))
  expect_true(is.finite(loglik(s$ev, p_h, mark_type = "ba_bip", attachment = sp_h)))
  expect_error(sim_hawkesNet(p_h, T_end = 35, T0 = 1, mark_type = "ba_bip", attachment = sp_h,
                             history = hist_h), "after the last history event")
  expect_error(sim_hawkesNet(p_h, T_end = 35, mark_type = "ba_bip", attachment = sp_h,
                             history = hist_h, burn_in = 5), "not both")
})

test_that("a burn-in becomes the history", {
  s <- sim_hawkesNet(p_h, T_end = 10, mark_type = "ba_bip", attachment = sp_h, burn_in = 15,
                     return_net = FALSE)
  expect_equal(s$history_end, 0)
  expect_equal(s$n_history, sum(s$ev$times$t < 0))
  expect_true(s$n_history > 0)
  expect_true(all(s$events$t >= 0))
  expect_lt(min(s$ev$times$t), 0)
})

test_that("continuing a history keeps the time part calibrated", {
  # Under the model, N - Lambda over the continuation has mean 0, with the
  # history's events (and path-dependent jumps) already in the intensity.
  exc <- excitation_spec(size = TRUE, connectedness = TRUE, connectedness_window = 5,
                         path_dependent = TRUE)
  p <- c(p_h, alpha = 0.3, gamma = 0.4, K_fast = 1, beta_fast = 6)
  T1 <- 30
  set.seed(302)
  d <- vapply(1:200, function(r) {
    s <- sim_hawkesNet(p, T_end = T1, mark_type = "ba_bip", attachment = sp_h, excitation = exc,
                       kernel = "exp2", history = hist_h, return_net = FALSE)
    dg <- diag_temporal(s$ev, p, T_end = T1, excitation = exc, kernel = "exp2")
    post <- seq_along(dg$t) > s$n_history
    s$events$n - (sum(dg$comp_inc[post]) + dg$meta$comp_tail)
  }, numeric(1))
  expect_lt(abs(mean(d) / (sd(d) / sqrt(length(d)))), 4)
})

test_that("continuing a history keeps the mark part calibrated", {
  # Per-event scores of the mark likelihood have mean 0 given the past, so
  # summed over the simulated events they average 0 across runs only if the
  # history's perps were loaded into the simulator correctly.
  mark <- function(ev) compile_logpmf_attach(ev, sp_h)
  n_h <- nrow(hist_h$times)
  set.seed(303)
  sc <- t(vapply(1:200, function(r) {
    s <- sim_hawkesNet(p_h, T_end = 30, mark_type = "ba_bip", attachment = sp_h,
                       history = hist_h, return_net = FALSE)
    S <- attr(mark(s$ev)(p_h, scores = TRUE), "scores")
    colSums(S[seq_len(nrow(S)) > n_h, c("m", "beta_edges", "delta", "theta", "theta_paths",
                                        "theta_recent"), drop = FALSE])
  }, numeric(6)))
  z <- colMeans(sc) / (apply(sc, 2, sd) / sqrt(nrow(sc)))
  expect_true(all(abs(z) < 4), info = paste(names(z), round(z, 2), collapse = ", "))
})

test_that("simulate(from = \"end\") continues the data", {
  fit <- fit_hawkesNet(hist_h, p_h, mark_type = "ba_bip", attachment = sp_h)
  expect_error(simulate(fit, from = "end"), "give `T_end`")
  sims <- simulate(fit, nsim = 2, seed = 1, from = "end", T_end = 30)
  expect_equal(attr(sims, "history_end"), max(hist_h$times$t))
  expect_equal(attr(sims, "exploded"), c(FALSE, FALSE))
  expect_equal(sims[[1]]$times$t[seq_len(nrow(hist_h$times))], hist_h$times$t)
})
