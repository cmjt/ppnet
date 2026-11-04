# The excitation weights are recomputed from the raw edge lists, without
# event_summary(), so the likelihood and its design are checked independently.

brute_excitation_ll <- function(ev, p, T_end, window = Inf, scale = "log1p", min_gap = 0) {
  arrivals <- nodes_by_event(ev)
  edges <- edges_by_event(ev)
  t <- ev$times$t
  perps <- lapply(seq_along(t), function(k) {
    neighbours_of(as.character(arrivals[[k]]$id[arrivals[[k]]$role == "event"]), edges[[k]])
  })
  conn <- vapply(seq_along(t), function(k) {
    sum(vapply(seq_len(k - 1L), function(m) {
      any(perps[[m]] %in% perps[[k]]) && t[k] - t[m] <= window && t[k] - t[m] > min_gap
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

test_that("connectedness_min_gap leaves out recent events in the likelihood", {
  local_min_one_perp()
  set.seed(63)
  ev <- sim_hawkesNet(p_exc, T_end = 15, mark_type = "ba_bip", excitation = exc)$ev
  skip_if(nrow(ev$times) < 5L)
  mark <- mark_part(ev, p_exc, 15)
  gapped <- excitation_spec(size = TRUE, connectedness = TRUE, connectedness_window = 4,
                            connectedness_min_gap = 0.5)
  expect_equal(loglik(ev, p_exc, mark_type = "ba_bip", T_end = 15, excitation = gapped),
               brute_excitation_ll(ev, p_exc, 15, window = 4, min_gap = 0.5) + mark,
               tolerance = 1e-8)
  expect_error(excitation_spec(connectedness = TRUE, connectedness_window = 1,
                               connectedness_min_gap = 1))
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
  # A minimum gap of 1 drops links to the event just before, keeping 3 -> 1.
  expect_equal(event_summary(ev, connectedness_window = 2, connectedness_min_gap = 1)$connectedness,
               c(0, 0, 1, 0, 0))
  expect_equal(event_summary(ev, connectedness_min_gap = 1.5)$connectedness,
               event_summary(ev, connectedness_window = Inf, connectedness_min_gap = 1.5)$connectedness)
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

# Path-dependent excitation, from its definition: c_i(t) counts the events
# before t, other than i, sharing a perp with i, at a distance in
# (min_gap, window] either side of t_i. The compensator integrates each
# event's term piece by piece between the times its c_i changes.
brute_path_ll <- function(ev, p, T_end, window = Inf, scale = "log1p", min_gap = 0) {
  arrivals <- nodes_by_event(ev)
  edges <- edges_by_event(ev)
  t <- ev$times$t
  n <- length(t)
  perps <- lapply(seq_len(n), function(k) {
    neighbours_of(as.character(arrivals[[k]]$id[arrivals[[k]]$role == "event"]), edges[[k]])
  })
  linked <- outer(seq_len(n), seq_len(n), Vectorize(function(i, m) {
    i != m && any(perps[[i]] %in% perps[[m]]) && abs(t[i] - t[m]) <= window &&
      abs(t[i] - t[m]) > min_gap
  }))
  h <- function(c) if (scale == "log1p") log1p(c) else c
  g <- function(i, c) lengths(perps)[i]^p$alpha * exp(p$gamma * h(c))
  c_at <- function(i, tt) sum(linked[i, ] & t < tt)

  lam <- vapply(seq_len(n), function(k) {
    i <- seq_len(k - 1L)
    p$mu + p$K * sum(vapply(i, function(j) g(j, c_at(j, t[k])), numeric(1)) *
                       exp(-p$beta * (t[k] - t[i])))
  }, numeric(1))
  comp <- sum(vapply(seq_len(n), function(i) {
    br <- sort(unique(c(t[i], t[linked[i, ] & t > t[i]], T_end)))
    sum(vapply(seq_len(length(br) - 1L), function(s) {
      a <- br[s]; b <- br[s + 1L]
      g(i, sum(linked[i, ] & t <= a)) *
        (exp(-p$beta * (a - t[i])) - exp(-p$beta * (b - t[i]))) / p$beta
    }, numeric(1)))
  }, numeric(1)))
  sum(log(lam)) - p$mu * T_end - p$K * comp
}

test_that("path-dependent excitation loglik matches a brute-force calculation", {
  local_min_one_perp()
  set.seed(65)
  pd <- excitation_spec(size = TRUE, connectedness = TRUE, path_dependent = TRUE)
  ev <- sim_hawkesNet(p_exc, T_end = 15, mark_type = "ba_bip", excitation = pd)$ev
  skip_if(nrow(ev$times) < 5L)
  mark <- mark_part(ev, p_exc, 15)
  expect_equal(loglik(ev, p_exc, mark_type = "ba_bip", T_end = 15, excitation = pd),
               brute_path_ll(ev, p_exc, 15) + mark, tolerance = 1e-8)
  # It differs from excitation fixed at birth once gamma != 0 ...
  expect_gt(abs(loglik(ev, p_exc, mark_type = "ba_bip", T_end = 15, excitation = pd) -
                  loglik(ev, p_exc, mark_type = "ba_bip", T_end = 15, excitation = exc)), 1e-3)
  # ... and equals it at gamma = 0.
  p0 <- modifyList(p_exc, list(gamma = 0))
  expect_equal(loglik(ev, p0, mark_type = "ba_bip", T_end = 15, excitation = pd),
               loglik(ev, p0, mark_type = "ba_bip", T_end = 15, excitation = exc), tolerance = 1e-10)

  win <- excitation_spec(size = FALSE, connectedness = TRUE, connectedness_scale = "linear",
                         connectedness_window = 3, connectedness_min_gap = 0.3,
                         path_dependent = TRUE)
  p_lin <- modifyList(p_exc, list(alpha = 0, gamma = -0.2))
  expect_equal(loglik(ev, p_lin[-6], mark_type = "ba_bip", T_end = 15, excitation = win),
               brute_path_ll(ev, p_lin, 15, window = 3, scale = "linear", min_gap = 0.3) + mark,
               tolerance = 1e-8)
  expect_error(excitation_spec(size = TRUE, path_dependent = TRUE), "connectedness")
})

test_that("the path-dependent simulator matches its compensator", {
  # N(T) - Lambda(T) is a martingale with variance E N(T), so pooled over
  # replicates it is about N(0, sum N) under the model that simulated. With
  # gamma < 0 later links cut an event's excitation, so the compensator with
  # excitation fixed at birth is clearly off (z near -7 here) and a simulator
  # that skipped the updates would fail; gamma < 0 also cannot explode.
  set.seed(66)
  pd <- excitation_spec(size = FALSE, connectedness = TRUE, connectedness_scale = "linear",
                        connectedness_window = 5, path_dependent = TRUE)
  fx <- excitation_spec(size = FALSE, connectedness = TRUE, connectedness_scale = "linear",
                        connectedness_window = 5)
  p <- list(mu = 1, K = 0.8, beta = 1, beta_edges = 0.3, lambda_new = 0.5, gamma = -1)
  comp <- function(ev, sp) {
    d <- diag_temporal(ev, p, T_end = 20, excitation = sp)
    sum(d$comp_inc) + d$meta$comp_tail
  }
  tot <- vapply(1:400, function(r) {
    ev <- suppressMessages(sim_hawkesNet(p, T_end = 20, mark_type = "ba_bip", excitation = pd,
                                         return_net = FALSE))$ev
    c(n = nrow(ev$times), pd = comp(ev, pd), fx = comp(ev, fx))
  }, numeric(3))
  z <- (sum(tot["n", ]) - rowSums(tot[c("pd", "fx"), ])) / sqrt(sum(tot["n", ]))
  expect_lt(abs(z[["pd"]]), 4)
  expect_lt(z[["fx"]], -4)
})

test_that("fit_hawkesNet recovers path-dependent excitation coefficients", {
  local_min_one_perp()
  set.seed(67)
  sp <- excitation_spec(size = TRUE, connectedness = TRUE, connectedness_window = 5,
                        path_dependent = TRUE)
  p <- list(mu = 1.5, K = 0.4, beta = 1.5, beta_edges = 0.3, lambda_new = 1,
            alpha = 0.6, gamma = 0.4)
  ev <- suppressMessages(sim_hawkesNet(p, T_end = 200, mark_type = "ba_bip", excitation = sp,
                                       return_net = FALSE))$ev
  init <- modifyList(p, list(mu = 1, K = 0.2, beta = 1, alpha = 0, gamma = 0))
  fit <- suppressWarnings(fit_hawkesNet(
    ev, init, fixed_params = c("beta_edges", "lambda_new"),
    mark_type = "ba_bip", T_end = 200, excitation = sp, hessian = TRUE
  ))
  est <- unlist(fit$par[c("alpha", "gamma")])
  se <- sqrt(diag(solve(fit$hessian)))[names(est)]
  expect_true(all(abs(est - c(0.6, 0.4)) < 4 * se))
})
