# 4-cycles in the attachment model, by pseudo-likelihood. The brute force
# scores every old perp at every event by its conditional given the others.

# Conditional join probabilities and outcomes for every (event, old perp),
# with z the 4-cycle statistic (before the scale h).
brute_cycle_terms <- function(ev, p, spec) {
  arrivals <- nodes_by_event(ev)
  edges <- edges_by_event(ev)
  t <- ev$times$t
  win <- spec$coffender_window
  gap <- spec$episode_gap %||% 0
  h <- function(n) if (spec$coffender_scale == "log1p") log1p(n) else n
  deg <- numeric(0); last <- numeric(0); ties <- list(); shared <- list(); own <- list()
  out <- list()
  for (k in seq_along(t)) {
    nodes_k <- arrivals[[k]]
    event_id <- as.character(nodes_k$id[nodes_k$role == "event"])
    perps <- neighbours_of(event_id, edges[[k]])
    old <- names(deg)
    if (length(old) > 0L) {
      w <- exp(-p$beta_edges * (t[k] - last[old]))
      if (spec$degree) w <- w * (deg[old] + p$delta)
      if (spec$coffenders) {
        x <- vapply(old, function(v) sum(t[k] - (ties[[v]] %||% numeric(0)) <= win), numeric(1))
        w <- w * exp(p$theta * h(x))
      }
      if (isTRUE(spec$recent)) {
        r <- vapply(old, function(v) sum(t[k] - own[[v]] <= gap), numeric(1))
        w <- w * exp(p$theta_recent * h(r))
      }
      q <- p$m * w / sum(w)
      y <- old %in% perps
      joiners <- old[y]
      # shared events with the other old joiners: within the episode gap
      # (short) and beyond it, up to the window (long)
      zz <- vapply(old, function(v) {
        rowSums(vapply(setdiff(joiners, v), function(u) {
          tt <- shared[[paste(sort(c(u, v)), collapse = "|")]] %||% numeric(0)
          d <- t[k] - tt
          c(sum(d <= win & d > gap), sum(d <= gap))
        }, numeric(2)))
      }, numeric(2))
      z <- zz[1, ]; zs <- zz[2, ]
      eta <- log(expm1(q)) + (p$theta_cycles %||% 0) * h(z) +
        (p$theta_cycles_short %||% 0) * h(zs)
      out[[k]] <- data.frame(k = k, v = old, y = y, z = z, zs = zs, q = q, p = plogis(eta))
    }
    new <- as.character(nodes_k$id[nodes_k$role == "perp"])
    deg[new] <- 0
    deg[perps] <- deg[perps] + 1
    last[perps] <- t[k]
    for (v in perps) own[[v]] <- c(own[[v]], t[k])
    for (v in perps) {
      tv <- ties[[v]] %||% numeric(0)
      tv[setdiff(perps, v)] <- t[k]
      ties[[v]] <- tv
    }
    if (length(perps) > 1L) {
      for (pr in utils::combn(sort(perps), 2, simplify = FALSE)) {
        key <- paste(pr, collapse = "|")
        shared[[key]] <- c(shared[[key]], t[k])
      }
    }
  }
  do.call(rbind, out)
}

brute_cycle_mark_ll <- function(ev, p, spec) {
  d <- brute_cycle_terms(ev, p, spec)
  n_new <- vapply(nodes_by_event(ev), function(n) sum(n$role == "perp"), numeric(1))
  sum(dpois(n_new, p$lambda_new, log = TRUE)) +
    sum(ifelse(d$y, log(d$p), log1p(-d$p)))
}

cyc <- attachment_spec(degree = TRUE, coffenders = TRUE, coffender_window = 4, cycles = TRUE)
p_cyc <- list(mu = 2, K = 0.5, beta = 1, m = 1.5, beta_edges = 0.3, delta = 0.5,
              theta = 0.4, lambda_new = 0.8, theta_cycles = 0.6)

test_that("4-cycle pseudo-likelihood matches a brute-force calculation", {
  set.seed(91)
  ev <- sim_hawkesNet(p_cyc, T_end = 25, mark_type = "ba_bip", attachment = cyc,
                      return_net = FALSE)$ev
  d <- brute_cycle_terms(ev, p_cyc, cyc)
  expect_gt(sum(d$z > 0 & d$y), 5)
  for (spec in list(cyc,
                    attachment_spec(degree = FALSE, cycles = TRUE, coffender_window = 4,
                                    coffender_scale = "linear"),
                    attachment_spec(degree = TRUE, cycles = TRUE, coffender_window = Inf))) {
    q <- p_cyc[c("m", "beta_edges", attachment_param_names(spec))]
    expect_equal(compile_logpmf_attach(ev, spec)(q), brute_cycle_mark_ll(ev, q, spec),
                 tolerance = 1e-8)
  }
})

test_that("theta_cycles = 0 is the model without 4-cycles", {
  set.seed(92)
  ev <- sim_hawkesNet(p_cyc, T_end = 25, mark_type = "ba_bip", attachment = cyc,
                      return_net = FALSE)$ev
  no_cyc <- attachment_spec(degree = TRUE, coffenders = TRUE, coffender_window = 4)
  expect_equal(compile_logpmf_attach(ev, cyc)(modifyList(p_cyc, list(theta_cycles = 0))),
               compile_logpmf_attach(ev, no_cyc)(p_cyc))
  expect_error(attachment_spec(cycles = TRUE, random_effect = TRUE), "cannot be combined")
})

test_that("4-cycle gradient and per-event scores are right", {
  set.seed(93)
  ev <- sim_hawkesNet(p_cyc, T_end = 25, mark_type = "ba_bip", attachment = cyc,
                      return_net = FALSE)$ev
  check_grad(.build_loglik_fn(ev, T_end = 25, mark_type = "ba_bip", attachment = cyc),
             modifyList(p_cyc, list(theta_cycles = 0.5, theta = -0.2)))
  sp <- attachment_spec(degree = TRUE, coffenders = TRUE, paths = TRUE, cycles = TRUE,
                        coffender_window = 4, new_perps = "negbin", min_perps = 1)
  q <- c(p_cyc, theta_paths = 0.2, lambda_new0 = 0.9, size_new = 1.5)
  ev <- sim_hawkesNet(q, T_end = 25, mark_type = "ba_bip", attachment = sp, return_net = FALSE)$ev
  check_grad(.build_loglik_fn(ev, T_end = 25, mark_type = "ba_bip", attachment = sp), q)
  ll <- compile_logpmf_attach(ev, sp)(q, gradient = TRUE, scores = TRUE)
  expect_equal(colSums(attr(ll, "scores")), attr(ll, "gradient"), tolerance = 1e-8)
})

# Positive theta_cycles feeds on itself (past co-offenders pull each other in,
# making bigger groups and more shared events), so like ERGMs with closure
# terms the model runs away to huge events beyond a point; theta_cycles = 0.9
# with m = 1.5 already does. The values here stay below that.

test_that("the simulator draws joins with the model's conditional probabilities", {
  # Gibbs sampling targets the per-event autologistic model, whose full
  # conditionals are the pseudo-likelihood terms: observed joins should
  # match their sum, overall and among perps with past co-offenders present.
  set.seed(94)
  p <- modifyList(p_cyc, list(mu = 3))
  ev <- sim_hawkesNet(p, T_end = 80, mark_type = "ba_bip", attachment = cyc,
                      return_net = FALSE)$ev
  d <- brute_cycle_terms(ev, p, cyc)
  z <- function(keep) sum(d$y[keep] - d$p[keep]) / sqrt(sum(d$p[keep] * (1 - d$p[keep])))
  expect_gt(sum(d$y & d$z > 0), 50)
  expect_lt(abs(z(TRUE)), 4)
  expect_lt(abs(z(d$z > 0)), 4)
  expect_lt(abs(z(d$z == 0)), 4)
})

test_that("fit_hawkesNet recovers theta_cycles, with sandwich standard errors", {
  set.seed(95)
  p <- list(mu = 3, K = 0.4, beta = 1.5, m = 1.2, beta_edges = 0.3, delta = 1,
            lambda_new = 1, theta_cycles = 0.6)
  sp <- attachment_spec(degree = TRUE, cycles = TRUE, coffender_window = 10)
  ev <- sim_hawkesNet(p, T_end = 300, mark_type = "ba_bip", attachment = sp, return_net = FALSE)$ev
  fit <- suppressWarnings(fit_hawkesNet(
    ev, modifyList(p, list(m = 0.5, beta_edges = 1, delta = 0.3, theta_cycles = 0)),
    fixed_params = c("mu", "K", "beta"), mark_type = "ba_bip", T_end = 300,
    attachment = sp, hessian = TRUE
  ))
  V <- attachment_vcov(fit)
  se <- attr(V, "se")
  expect_equal(names(se), c("m", "beta_edges", "delta", "lambda_new", "theta_cycles"))
  expect_lt(abs(fit$par$theta_cycles - 0.6), 4 * se[["theta_cycles"]])
  expect_lt(se[["theta_cycles"]], 0.5)
  # Natural scale: estimate * se for the log-transformed parameters, unchanged
  # for theta_cycles (untransformed).
  se_nat <- attr(V, "se_natural")
  expect_equal(names(se_nat), names(se))
  logp <- c("m", "beta_edges", "delta", "lambda_new")
  expect_equal(se_nat[logp], unlist(fit$par[logp]) * se[logp])
  expect_equal(se_nat[["theta_cycles"]], se[["theta_cycles"]])
  # Without 4-cycles it is a likelihood, so sandwich and Hessian agree roughly.
  sp0 <- attachment_spec(degree = TRUE, coffender_window = 10)
  fit0 <- suppressWarnings(fit_hawkesNet(
    ev, p[names(p) != "theta_cycles"], fixed_params = c("mu", "K", "beta"),
    mark_type = "ba_bip", T_end = 300, attachment = sp0, hessian = TRUE
  ))
  ratio <- attr(attachment_vcov(fit0), "se") / sqrt(diag(solve(fit0$hessian)))[c("m", "beta_edges", "delta", "lambda_new")]
  expect_true(all(ratio > 0.7 & ratio < 1.4))
})

# Episode gap: within-episode vs separate-occasion 4-cycles, and recent
# attendance. Events here are about 0.3 apart, so a gap of 0.5 separates
# bursts from returns.
epi <- attachment_spec(degree = TRUE, coffenders = TRUE, coffender_window = 4, cycles = TRUE,
                       episode_gap = 0.5, recent = TRUE)
p_epi <- c(p_cyc, theta_cycles_short = 0.4, theta_recent = 0.5)

test_that("episode-gap terms match a brute-force calculation", {
  set.seed(96)
  ev <- sim_hawkesNet(p_epi, T_end = 25, mark_type = "ba_bip", attachment = epi,
                      return_net = FALSE)$ev
  d <- brute_cycle_terms(ev, p_epi, epi)
  expect_gt(sum(d$zs > 0 & d$y), 3)
  expect_gt(sum(d$z > 0 & d$y), 3)
  rec <- attachment_spec(degree = TRUE, coffender_window = 4, episode_gap = 0.5, recent = TRUE)
  for (spec in list(epi, rec,
                    attachment_spec(degree = FALSE, cycles = TRUE, coffender_window = 4,
                                    episode_gap = 0.5, coffender_scale = "linear"))) {
    q <- p_epi[c("m", "beta_edges", attachment_param_names(spec))]
    expect_equal(compile_logpmf_attach(ev, spec)(q), brute_cycle_mark_ll(ev, q, spec),
                 tolerance = 1e-8)
  }
  expect_equal(attachment_param_names(epi),
               c("m", "beta_edges", "delta", "theta", "theta_recent", "lambda_new",
                 "theta_cycles", "theta_cycles_short"))
  expect_error(attachment_spec(recent = TRUE), "episode_gap")
  expect_error(attachment_spec(coffender_window = 1, episode_gap = 1), "less than")
})

test_that("episode-gap gradient and scores are right", {
  set.seed(97)
  ev <- sim_hawkesNet(p_epi, T_end = 25, mark_type = "ba_bip", attachment = epi,
                      return_net = FALSE)$ev
  check_grad(.build_loglik_fn(ev, T_end = 25, mark_type = "ba_bip", attachment = epi),
             modifyList(p_epi, list(theta_recent = -0.3, theta_cycles_short = 0.8)))
  ll <- compile_logpmf_attach(ev, epi)(p_epi, gradient = TRUE, scores = TRUE)
  expect_equal(colSums(attr(ll, "scores")), attr(ll, "gradient"), tolerance = 1e-8)
})

test_that("the simulator draws episode-gap joins with the model's conditionals", {
  set.seed(98)
  p <- modifyList(p_epi, list(mu = 3, theta_cycles_short = 1, theta_cycles = 0.3))
  ev <- sim_hawkesNet(p, T_end = 80, mark_type = "ba_bip", attachment = epi,
                      return_net = FALSE)$ev
  d <- brute_cycle_terms(ev, p, epi)
  z <- function(keep) sum(d$y[keep] - d$p[keep]) / sqrt(sum(d$p[keep] * (1 - d$p[keep])))
  expect_gt(sum(d$y & d$zs > 0), 30)
  expect_gt(sum(d$y & d$z > 0), 30)
  expect_lt(abs(z(TRUE)), 4)
  expect_lt(abs(z(d$zs > 0)), 4)
  expect_lt(abs(z(d$z > 0 & d$zs == 0)), 4)
  expect_lt(abs(z(d$z == 0 & d$zs == 0)), 4)
})
