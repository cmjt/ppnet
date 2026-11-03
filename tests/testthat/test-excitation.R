# The excitation weights are recomputed from the raw edge lists, without
# event_summary(), so the likelihood and its design are checked independently.

brute_excitation_ll <- function(ev, p, T_end, window = Inf, scale = "log1p") {
  arrivals <- nodes_by_event(ev)
  edges <- edges_by_event(ev)
  t <- ev$times$t
  perps <- lapply(seq_along(t), function(k) {
    neighbours_of(as.character(arrivals[[k]]$id[arrivals[[k]]$role == "event"]), edges[[k]])
  })
  conn <- vapply(seq_along(t), function(k) {
    sum(vapply(seq_len(k - 1L), function(m) {
      any(perps[[m]] %in% perps[[k]]) && t[k] - t[m] <= window
    }, logical(1)))
  }, numeric(1))
  h <- if (scale == "log1p") log1p(conn) else conn
  g <- lengths(perps)^p$alpha * exp(p$gamma * h)

  lam <- vapply(seq_along(t), function(i) {
    j <- seq_len(i - 1L)
    p$mu + p$K * sum(g[j] * exp(-p$beta * (t[i] - t[j])))
  }, numeric(1))
  sum(log(lam)) - p$mu * T_end - (p$K / p$beta) * sum(g * (1 - exp(-p$beta * (T_end - t))))
}

p_exc <- list(mu = 1, K = 0.4, beta = 1, beta_edges = 0.3, lambda_new = 0.9,
              alpha = 0.5, gamma = 0.3)
exc <- excitation_spec(size = TRUE, connectedness = TRUE)

# The mark part is unchanged by the excitation, so take it from the plain model.
mark_part <- function(ev, p, T_end) {
  p0 <- modifyList(p, list(alpha = 0, gamma = 0))
  loglik(ev, p0[c("mu", "K", "beta", "beta_edges", "lambda_new")], mark_type = "ba_bip", T_end = T_end) -
    brute_excitation_ll(ev, p0, T_end)
}

test_that("excitation loglik matches a brute-force calculation", {
  local_min_one_perp()
  set.seed(61)
  ev <- sim_hawkesNet(p_exc, T_end = 15, mark_type = "ba_bip", excitation = exc)$ev
  skip_if(nrow(ev$times) < 5L)
  mark <- mark_part(ev, p_exc, 15)

  expect_equal(loglik(ev, p_exc, mark_type = "ba_bip", T_end = 15, excitation = exc),
               brute_excitation_ll(ev, p_exc, 15) + mark, tolerance = 1e-8)

  lin <- excitation_spec(size = FALSE, connectedness = TRUE,
                         connectedness_scale = "linear", connectedness_window = 2)
  p_lin <- modifyList(p_exc, list(alpha = 0, gamma = 0.1))
  expect_equal(loglik(ev, p_lin[-6], mark_type = "ba_bip", T_end = 15, excitation = lin),
               brute_excitation_ll(ev, p_lin, 15, window = 2, scale = "linear") + mark,
               tolerance = 1e-8)
})

test_that("zero excitation coefficients give the plain Hawkes model", {
  local_min_one_perp()
  set.seed(62)
  ev <- sim_hawkesNet(p_exc[1:5], T_end = 12, mark_type = "ba_bip")$ev
  p0 <- modifyList(p_exc, list(alpha = 0, gamma = 0))
  expect_equal(loglik(ev, p0, mark_type = "ba_bip", T_end = 12, excitation = exc),
               loglik(ev, p_exc[1:5], mark_type = "ba_bip", T_end = 12), tolerance = 1e-10)
})

test_that("event_summary applies the connectedness window", {
  # Events at times 1..5; event 3 links to 1 and 2, event 4 to 3.
  ev <- make_events(
    edges = data.frame(i = c("e1", "e1", "e2", "e2", "e3", "e3", "e3", "e4", "e5"),
                       j = c("a", "b", "b", "c", "a", "c", "d", "d", "f"),
                       time = c(1, 1, 2, 2, 3, 3, 3, 4, 5)),
    nodes = data.frame(id = c("e1", "a", "b", "e2", "c", "e3", "d", "e4", "e5", "f"),
                       time = c(1, 1, 1, 2, 2, 3, 3, 4, 5, 5),
                       role = c("event", "perp", "perp", "event", "perp", "event", "perp",
                                "event", "event", "perp")),
    allow_implicit_birth = FALSE
  )
  expect_equal(event_summary(ev, connectedness_window = 1)$connectedness, c(0, 1, 1, 1, 0))
  expect_equal(event_summary(ev, connectedness_window = 0.5)$connectedness, rep(0, 5))
})

test_that("excitation needs a bipartite history and its coefficients", {
  set.seed(63)
  ev <- sim_hawkesNet(p_exc[1:4], T_end = 6, mark_type = "ba")$ev
  expect_error(loglik(ev, c(p_exc[1:4], alpha = 0.1), mark_type = "ba",
                      excitation = excitation_spec()), "roles")
  expect_error(sim_hawkesNet(c(p_exc[1:4], alpha = 0.1), T_end = 2, mark_type = "ba",
                             excitation = excitation_spec()), "bipartite")

  ev_bip <- sim_hawkesNet(p_exc[1:5], T_end = 6, mark_type = "ba_bip")$ev
  expect_error(fit_hawkesNet(ev_bip, p_exc[1:6], mark_type = "ba_bip", excitation = exc),
               "Missing excitation coefficients")
})

test_that("fit_hawkesNet recovers excitation coefficients", {
  local_min_one_perp()
  set.seed(64)
  sp <- excitation_spec(size = TRUE, connectedness = TRUE, connectedness_window = 5)
  p <- list(mu = 1.5, K = 0.4, beta = 1.5, beta_edges = 0.3, lambda_new = 1,
            alpha = 0.6, gamma = 0.4)
  ev <- sim_hawkesNet(p, T_end = 200, mark_type = "ba_bip", excitation = sp)$ev
  init <- modifyList(p, list(mu = 1, K = 0.2, beta = 1, alpha = 0, gamma = 0))

  fit <- suppressWarnings(fit_hawkesNet(
    ev, init, fixed_params = c("beta_edges", "lambda_new"),
    mark_type = "ba_bip", T_end = 200, excitation = sp, hessian = TRUE
  ))
  est <- unlist(fit$par[c("alpha", "gamma")])
  se <- sqrt(diag(solve(fit$hessian)))[names(est)]
  expect_true(all(abs(est - c(0.6, 0.4)) < 4 * se))
})
