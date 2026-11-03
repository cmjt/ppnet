# Perp frailty in the attachment model. The brute force scores every old perp
# at every event, then integrates each perp's frailty out numerically.

brute_frailty_mark_ll <- function(ev, p, spec) {
  arrivals <- nodes_by_event(ev)
  edges <- edges_by_event(ev)
  t <- ev$times$t
  h <- function(n) if (spec$coffender_scale == "log1p") log1p(n) else n
  deg <- numeric(0); last <- numeric(0); ties <- list()
  Q <- numeric(0); qj <- list()

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
      q <- p$m * w / sum(w)
      joined <- old %in% perps
      Q[old] <- Q[old] + q
      for (v in old[joined]) qj[[v]] <- c(qj[[v]], q[[v]])
    }

    deg[new_perps] <- 0
    Q[new_perps] <- 0
    deg[perps] <- deg[perps] + 1
    last[perps] <- t[k]
    for (v in perps) {
      tv <- ties[[v]] %||% numeric(0)
      tv[setdiff(perps, v)] <- t[k]
      ties[[v]] <- tv
    }
  }

  # Each perp: log E[prod_joins (1 - exp(-u q)) exp(-u (Q - sum q_joins))]
  # for u ~ Gamma(a, a), integrated numerically in log u after scaling by
  # the integrand's maximum. Perps that were never candidates contribute 0.
  a <- 1 / p$frailty_var
  for (v in names(Q)) {
    qv <- qj[[v]] %||% numeric(0)
    lf <- function(s) {
      vapply(s, function(si) sum(log(-expm1(-exp(si) * qv))) - exp(si) * (Q[[v]] - sum(qv)) +
               dgamma(exp(si), a, a, log = TRUE) + si, numeric(1))
    }
    top <- optimize(lf, c(-30, 10), maximum = TRUE, tol = 1e-10)
    lo <- top$maximum - max(60 / (a + length(qv)), 3)
    hi <- top$maximum + 5
    ll <- ll + top$objective +
      log(integrate(function(s) exp(lf(s) - top$objective), lo, hi, rel.tol = 1e-12,
                    subdivisions = 1000L)$value)
  }
  ll
}

fr <- attachment_spec(degree = TRUE, coffenders = TRUE, coffender_window = 3, frailty = TRUE)
p_fr <- list(mu = 1, K = 0.5, beta = 1, m = 0.8, beta_edges = 0.4, delta = 0.5,
             theta = 0.7, lambda_new = 0.9, frailty_var = 0.6)

test_that("frailty loglik matches a brute-force calculation", {
  set.seed(81)
  ev <- sim_hawkesNet(p_fr, T_end = 20, mark_type = "ba_bip", attachment = fr)$ev
  skip_if(nrow(ev$times) < 10L)
  expect_gt(sum(event_summary(ev)$n_old), 5)

  for (spec in list(fr, attachment_spec(degree = FALSE, frailty = TRUE),
                    attachment_spec(degree = TRUE, coffenders = TRUE, coffender_scale = "linear",
                                    frailty = TRUE))) {
    q <- p_fr[c("m", "beta_edges", attachment_param_names(spec))]
    for (fv in c(0.6, 2.5)) {
      q$frailty_var <- fv
      expect_equal(compile_logpmf_attach(ev, spec)(q), brute_frailty_mark_ll(ev, q, spec),
                   tolerance = 1e-8)
    }
  }
})

test_that("frailty sums stay accurate across reference-time epochs", {
  # beta_edges * T_end is well over 300, so the cumulative sums restart
  # several times and perps' stretches of unchanged state span epochs.
  set.seed(82)
  p <- modifyList(p_fr, list(mu = 3, beta_edges = 4, m = 0.5))
  ev <- sim_hawkesNet(p, T_end = 250, mark_type = "ba_bip", attachment = fr)$ev
  expect_gt(nrow(ev$times), 600)
  expect_equal(compile_logpmf_attach(ev, fr)(p), brute_frailty_mark_ll(ev, p, fr), tolerance = 1e-8)
})

test_that("frailty_var = 0 is the model without frailty", {
  set.seed(83)
  ev <- sim_hawkesNet(p_fr, T_end = 20, mark_type = "ba_bip", attachment = fr)$ev
  no_fr <- attachment_spec(degree = TRUE, coffenders = TRUE, coffender_window = 3)
  ll0 <- compile_logpmf_attach(ev, no_fr)(p_fr)
  expect_equal(compile_logpmf_attach(ev, fr)(modifyList(p_fr, list(frailty_var = 0))), ll0)
  # Near zero the large terms of the frailty integral must cancel stably.
  for (fv in c(1e-3, 1e-5, 1e-7)) {
    ll_fv <- compile_logpmf_attach(ev, fr)(modifyList(p_fr, list(frailty_var = fv)))
    expect_lt(abs(ll_fv - ll0), 1e3 * fv)
  }
  expect_equal(compile_logpmf_attach(ev, fr)(modifyList(p_fr, list(frailty_var = 1e-300))), ll0)
  expect_equal(compile_logpmf_attach(ev, fr)(modifyList(p_fr, list(frailty_var = -1))), -Inf)
})

test_that("frailty gradient matches central differences", {
  set.seed(84)
  ev <- sim_hawkesNet(p_fr, T_end = 20, mark_type = "ba_bip", attachment = fr)$ev
  check_grad(.build_loglik_fn(ev, T_end = 20, mark_type = "ba_bip", attachment = fr),
             modifyList(p_fr, list(m = 1.1, beta_edges = 0.2, delta = 0.8, theta = 0.3,
                                   frailty_var = 1.3)))

  set.seed(85)
  p <- modifyList(p_fr, list(mu = 3, beta_edges = 4, m = 0.5))
  ev <- sim_hawkesNet(p, T_end = 150, mark_type = "ba_bip", attachment = fr)$ev
  mark <- compile_logpmf_attach(ev, fr)
  q <- p[attachment_param_names(fr)]
  g <- attr(mark(q, gradient = TRUE), "gradient")
  expect_equal(g[names(q)], num_grad(mark, q, names(q), h = 1e-6), tolerance = 1e-5)
})

test_that("simulator frailty spreads returns and matches frailty_shape", {
  p <- list(mu = 4, K = 0.3, beta = 1, m = 1.5, beta_edges = 0.02, lambda_new = 0.8)
  disp <- function(spec, params, ...) {
    set.seed(86)
    ev <- sim_hawkesNet(params, T_end = 150, mark_type = "ba_bip", attachment = spec,
                        return_net = FALSE, ...)$ev
    a <- as.numeric(table(ev$edges$j))
    var(a) / mean(a)
  }
  base <- attachment_spec(degree = FALSE)
  with_fr <- attachment_spec(degree = FALSE, frailty = TRUE)
  d0 <- disp(base, p)
  d1 <- disp(with_fr, c(p, frailty_var = 0.8))
  expect_equal(d1, disp(base, p, frailty_shape = 1 / 0.8))
  expect_gt(d1, 1.2 * d0)
  expect_error(sim_hawkesNet(c(p, frailty_var = 1), T_end = 5, mark_type = "ba_bip",
                             attachment = with_fr, frailty_shape = 2),
               "only when")
})

test_that("frailty_var can be fixed at 0 in a fit", {
  set.seed(87)
  ev <- sim_hawkesNet(p_fr, T_end = 20, mark_type = "ba_bip", attachment = fr)$ev
  no_fr <- attachment_spec(degree = TRUE, coffenders = TRUE, coffender_window = 3)
  f0 <- suppressWarnings(fit_hawkesNet(
    ev, modifyList(p_fr, list(frailty_var = 0)), fixed_params = c("mu", "K", "beta", "frailty_var"),
    mark_type = "ba_bip", T_end = 20, attachment = fr))
  f1 <- suppressWarnings(fit_hawkesNet(
    ev, p_fr[names(p_fr) != "frailty_var"], fixed_params = c("mu", "K", "beta"),
    mark_type = "ba_bip", T_end = 20, attachment = no_fr))
  expect_equal(f0$par$frailty_var, 0)
  expect_equal(f0$loglik, f1$loglik, tolerance = 1e-6)
})
