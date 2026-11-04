# Every old perp is scored directly with log p or log(1 - p), so this checks
# the running-sum evaluator and its sum(q) = m shortcut independently.

brute_attach_mark_ll <- function(ev, p, spec) {
  arrivals <- nodes_by_event(ev)
  edges <- edges_by_event(ev)
  t <- ev$times$t
  h <- function(n) if (spec$coffender_scale == "log1p") log1p(n) else n
  deg <- numeric(0); last <- numeric(0); ties <- list(); past <- list()

  ll <- 0
  for (k in seq_along(t)) {
    nodes_k <- arrivals[[k]]
    event_id <- as.character(nodes_k$id[nodes_k$role == "event"])
    new_perps <- as.character(nodes_k$id[nodes_k$role == "perp"])
    perps <- neighbours_of(event_id, edges[[k]])
    ll <- ll + dpois(length(new_perps), p$lambda_new, log = TRUE)

    old <- names(deg)
    if (length(old) > 0L) {
      x <- vapply(old, function(v) sum(t[k] - (ties[[v]] %||% numeric(0)) <= spec$coffender_window),
                  numeric(1))
      w <- exp(-p$beta_edges * (t[k] - last[old]))
      if (spec$degree) w <- w * (deg[old] + p$delta)
      if (spec$coffenders) w <- w * exp(p$theta * h(x))
      if (isTRUE(spec$paths)) {
        # 3-paths k - v - j - u: sum of (n_j - 1) over v's events j in the window
        y <- vapply(old, function(v) {
          pj <- past[[v]]
          if (is.null(pj)) 0 else sum(pj$n1[t[k] - pj$t <= spec$coffender_window])
        }, numeric(1))
        w <- w * exp(p$theta_paths * h(y))
      }
      pr <- 1 - exp(-p$m * w / sum(w))
      ll <- ll + bern_ll(pr, old %in% perps)
    }

    deg[new_perps] <- 0
    deg[perps] <- deg[perps] + 1
    last[perps] <- t[k]
    for (v in perps) {
      past[[v]] <- rbind(past[[v]], data.frame(t = t[k], n1 = length(perps) - 1))
      tv <- ties[[v]] %||% numeric(0)
      tv[setdiff(perps, v)] <- t[k]
      ties[[v]] <- tv
    }
  }
  ll
}

att <- attachment_spec(degree = TRUE, coffenders = TRUE, coffender_window = 3)
p_att <- list(mu = 1, K = 0.5, beta = 1, m = 0.8, beta_edges = 0.4, delta = 0.5,
              theta = 0.7, theta_paths = -0.4, lambda_new = 0.9)

test_that("attachment loglik matches a brute-force calculation", {
  set.seed(71)
  ev <- sim_hawkesNet(p_att, T_end = 20, mark_type = "ba_bip", attachment = att)$ev
  skip_if(nrow(ev$times) < 10L)
  expect_gt(sum(event_summary(ev)$n_old), 5)

  time_ll <- brute_time_ll(ev$times$t, 1, 0.5, 1, 20)
  for (spec in list(att,
                    attachment_spec(degree = FALSE, coffenders = FALSE),
                    attachment_spec(degree = TRUE, coffenders = TRUE, coffender_scale = "linear"),
                    attachment_spec(degree = TRUE, paths = TRUE, coffender_window = 3),
                    attachment_spec(degree = TRUE, coffenders = TRUE, paths = TRUE,
                                    coffender_window = 3, coffender_scale = "linear"))) {
    q <- p_att[c("mu", "K", "beta", attachment_param_names(spec))]
    expect_equal(loglik(ev, q, mark_type = "ba_bip", T_end = 20, attachment = spec),
                 time_ll + brute_attach_mark_ll(ev, q, spec), tolerance = 1e-8)
  }
})

test_that("delta = 0 is allowed: joins proportional to past attendance", {
  set.seed(72)
  ev <- sim_hawkesNet(p_att, T_end = 20, mark_type = "ba_bip", attachment = att)$ev
  q <- modifyList(p_att, list(delta = 0))
  time_ll <- brute_time_ll(ev$times$t, 1, 0.5, 1, 20)
  expect_equal(loglik(ev, q, mark_type = "ba_bip", T_end = 20, attachment = att),
               time_ll + brute_attach_mark_ll(ev, q, att), tolerance = 1e-8)
  expect_equal(loglik(ev, modifyList(q, list(delta = -0.1)), mark_type = "ba_bip",
                      T_end = 20, attachment = att), -Inf)
  re <- attachment_spec(degree = TRUE, coffenders = TRUE, coffender_window = 3,
                        random_effect = TRUE)
  ll0 <- compile_logpmf_attach(ev, re)(c(q, random_effect_var = 0.5))
  ll_eps <- compile_logpmf_attach(ev, re)(c(modifyList(q, list(delta = 1e-9)),
                                            random_effect_var = 0.5))
  expect_equal(ll0, ll_eps, tolerance = 1e-6)
})

test_that("running sums stay accurate over long histories", {
  # More than 512 events exercises the periodic recompute; strong recency
  # decay makes the sums span many orders of magnitude.
  set.seed(72)
  p <- modifyList(p_att, list(mu = 3, beta_edges = 2))
  ev <- sim_hawkesNet(p, T_end = 250, mark_type = "ba_bip", attachment = att)$ev
  expect_gt(nrow(ev$times), 600)
  mark <- compile_logpmf_attach(ev, att)
  expect_equal(mark(p), brute_attach_mark_ll(ev, p, att), tolerance = 1e-8)
})

test_that("attachment needs BA-bip and all its parameters", {
  set.seed(73)
  ev <- sim_hawkesNet(p_att, T_end = 5, mark_type = "ba_bip", attachment = att)$ev
  expect_error(loglik(ev, p_att, mark_type = "cs_bip", attachment = att), "ba_bip")
  expect_error(fit_hawkesNet(ev, p_att[-7], mark_type = "ba_bip", attachment = att),
               "Missing attachment parameters: theta")
  expect_equal(loglik(ev, modifyList(p_att, list(m = -1)), mark_type = "ba_bip", attachment = att), -Inf)
})

test_that("fit_hawkesNet recovers attachment parameters", {
  set.seed(74)
  p <- list(mu = 3, K = 0.4, beta = 1.5, m = 1.2, beta_edges = 0.3, delta = 1,
            theta = 0.6, lambda_new = 1)
  sp <- attachment_spec(degree = TRUE, coffenders = TRUE, coffender_window = 10)
  ev <- sim_hawkesNet(p, T_end = 300, mark_type = "ba_bip", attachment = sp)$ev
  init <- modifyList(p, list(m = 0.5, beta_edges = 1, delta = 0.3, theta = 0))

  fit <- suppressWarnings(fit_hawkesNet(
    ev, init, fixed_params = c("mu", "K", "beta"), mark_type = "ba_bip",
    T_end = 300, attachment = sp, hessian = TRUE
  ))
  nm <- c("m", "beta_edges", "delta", "theta", "lambda_new")
  se_work <- sqrt(diag(solve(fit$hessian)))[nm]
  est_work <- c(log(unlist(fit$par[c("m", "beta_edges", "delta")])), fit$par$theta,
                log(fit$par$lambda_new))
  true_work <- c(log(c(1.2, 0.3, 1)), 0.6, log(1))
  expect_true(all(abs(est_work - true_work) < 4 * se_work))
})

test_that("the C++ simulator joins perps with the model's probabilities", {
  # Given the simulated history, recompute every old perp's join probability
  # with the brute-force walk and compare observed joins with their sum,
  # overall and for the most and least active perps.
  set.seed(75)
  p <- list(mu = 3, K = 0.3, beta = 1, m = 1.5, beta_edges = 0.3, delta = 0.5,
            theta = 0.8, theta_paths = -0.3, lambda_new = 0.6, lambda_new0 = 0.7)
  sp <- attachment_spec(degree = TRUE, coffenders = TRUE, paths = TRUE, coffender_window = 4,
                        min_perps = 1)
  ev <- sim_hawkesNet(p, T_end = 120, mark_type = "ba_bip", attachment = sp, return_net = FALSE)$ev

  arrivals <- nodes_by_event(ev); edges <- edges_by_event(ev); t <- ev$times$t
  deg <- numeric(0); last <- numeric(0); ties <- list(); past <- list()
  obs <- pr <- dg <- numeric(0)
  for (k in seq_along(t)) {
    nodes_k <- arrivals[[k]]
    perps <- neighbours_of(as.character(nodes_k$id[nodes_k$role == "event"]), edges[[k]])
    old <- names(deg)
    if (length(old) > 0L) {
      x <- vapply(old, function(v) sum(t[k] - (ties[[v]] %||% numeric(0)) <= 4), numeric(1))
      y <- vapply(old, function(v) {
        pj <- past[[v]]
        if (is.null(pj)) 0 else sum(pj$n1[t[k] - pj$t <= 4])
      }, numeric(1))
      w <- (deg[old] + p$delta) * exp(p$theta * log1p(x) + p$theta_paths * log1p(y)) *
        exp(-p$beta_edges * (t[k] - last[old]))
      pr <- c(pr, 1 - exp(-p$m * w / sum(w)))
      obs <- c(obs, old %in% perps)
      dg <- c(dg, deg[old])
    }
    new <- as.character(nodes_k$id[nodes_k$role == "perp"])
    deg[new] <- 0
    deg[perps] <- deg[perps] + 1
    last[perps] <- t[k]
    for (v in perps) past[[v]] <- rbind(past[[v]], data.frame(t = t[k], n1 = length(perps) - 1))
    for (v in perps) { tv <- ties[[v]] %||% numeric(0); tv[setdiff(perps, v)] <- t[k]; ties[[v]] <- tv }
  }
  z <- function(keep) sum(obs[keep] - pr[keep]) / sqrt(sum(pr[keep] * (1 - pr[keep])))
  expect_gt(sum(obs), 200)
  expect_lt(abs(z(TRUE)), 4)
  expect_lt(abs(z(dg >= 3)), 4)
  expect_lt(abs(z(dg <= 1)), 4)
})

test_that("fit_hawkesNet recovers the 3-path coefficient", {
  set.seed(76)
  p <- list(mu = 3, K = 0.4, beta = 1.5, m = 1.2, beta_edges = 0.3, delta = 1,
            theta_paths = 0.5, lambda_new = 1)
  sp <- attachment_spec(degree = TRUE, paths = TRUE, coffender_window = 10)
  ev <- sim_hawkesNet(p, T_end = 300, mark_type = "ba_bip", attachment = sp, return_net = FALSE)$ev
  fit <- suppressWarnings(fit_hawkesNet(
    ev, modifyList(p, list(m = 0.5, beta_edges = 1, delta = 0.3, theta_paths = 0)),
    fixed_params = c("mu", "K", "beta"), mark_type = "ba_bip", T_end = 300,
    attachment = sp, hessian = TRUE
  ))
  se <- sqrt(diag(solve(fit$hessian)))[["theta_paths"]]
  # 3-paths are roughly attendance times (group size - 1), so theta_paths is
  # strongly correlated with delta and its standard error is wide.
  expect_lt(abs(fit$par$theta_paths - 0.5), 4 * se)
})
