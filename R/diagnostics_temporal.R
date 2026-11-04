# diagnostics_temporal.R

#' Temporal diagnostics for HawkesNet
#'
#' Computes time-series diagnostics for the temporal part of the model at the
#' observed event times: inter-arrival times, the intensity just before each
#' event, and the compensator over each gap. Under a correctly specified
#' model the compensator increments `tau` are independent Exp(1) (the
#' time-rescaling theorem), which `plot_temporal_diag()` checks.
#'
#' @param events An events object created by `make_events()`.
#' @param params Named list of model parameters on the natural scale, as for
#'   [loglik()]: `mu`, `K`, `beta`, plus any baseline coefficients,
#'   excitation coefficients and `K_fast`, `beta_fast`.
#' @param T_end Optional end time for diagnostics window. Defaults to last event.
#' @param T0 Start time for diagnostics window.
#' @param baseline,excitation,kernel As for [loglik()].
#'
#' @return A list with components:
#' \describe{
#'   \item{t}{Event times.}
#'   \item{dt}{Inter-arrival times (first is `t[1] - T0`).}
#'   \item{lambda}{Intensity just before each event time.}
#'   \item{comp_inc}{Compensator over each gap ending at an event.}
#'   \item{tau}{Time-rescaled increments, equal to `comp_inc`; Exp(1) under
#'     a correct temporal model.}
#'   \item{meta}{List with `T0`, `T_end`, `comp_tail` (compensator from the
#'     last event to `T_end`) and `ok`.}
#' }
#' @export
diag_temporal <- function(events, params, T_end = NULL, T0 = 0,
                          baseline = NULL, excitation = NULL, kernel = c("exp", "exp2")) {
  kernel <- match.arg(kernel)
  if (is.null(events) || is.null(events$times) || !"t" %in% names(events$times)) {
    stop("diag_temporal(): `events$times$t` is required.", call. = FALSE)
  }

  t <- as.numeric(events$times$t)
  if (length(t) == 0L) {
    return(list(
      t = numeric(0),
      dt = numeric(0),
      lambda = numeric(0),
      comp_inc = numeric(0),
      tau = numeric(0),
      meta = list(T0 = T0, T_end = T_end, ok = TRUE, note = "No events.")
    ))
  }
  if (any(!is.finite(t))) stop("diag_temporal(): non-finite event times.", call. = FALSE)
  if (any(diff(t) < 0)) stop("diag_temporal(): event times must be nondecreasing.", call. = FALSE)
  if (is.null(T_end)) T_end <- max(t)
  if (!is.finite(T_end) || T_end < max(t)) stop("diag_temporal(): T_end must be >= max(t).", call. = FALSE)

  n <- length(t)
  dt <- c(t[1] - T0, diff(t))
  invalid <- function(note) {
    list(t = t, dt = dt, lambda = rep(NA_real_, n), comp_inc = rep(NA_real_, n),
         tau = rep(NA_real_, n), meta = list(T0 = T0, T_end = T_end, ok = FALSE, note = note))
  }

  mu <- params$mu
  K <- params$K
  beta <- params$beta
  if (!is.finite(mu) || !is.finite(K) || !is.finite(beta) || mu <= 0 || K < 0 || beta <= 0) {
    return(invalid("Invalid mu/K/beta."))
  }

  # Baseline at the events, and its integral over each gap (last: to T_end)
  lo <- c(T0, t)
  hi <- c(t, T_end)
  if (is.null(baseline)) {
    mu_ev <- rep.int(mu, n)
    base_inc <- mu * (hi - lo)
  } else {
    b <- .baseline_coefs(baseline, params)
    mu_ev <- mu * exp(drop(.baseline_design(baseline, t) %*% b))
    base_inc <- mu * .baseline_integral(baseline, b, lo, hi)
  }

  # Kernel components: (K, beta, jump sizes); excitation scales the slow one
  # With path-dependent excitation these are the effective jumps, which
  # carry the later changes in each event's excitation.
  g <- if (is.null(excitation)) rep.int(1, n) else {
    .excitation_jumps(excitation, .excitation_design(excitation, events),
                      .excitation_coefs(excitation, params), beta)$G
  }
  comps <- list(list(K = K, beta = beta, g = g))
  if (kernel == "exp2") {
    if (is.null(params$K_fast) || is.null(params$beta_fast)) {
      stop("diag_temporal(): kernel = 'exp2' needs K_fast and beta_fast.", call. = FALSE)
    }
    comps[[2L]] <- list(K = params$K_fast, beta = params$beta_fast, g = rep.int(1, n))
  }

  exc_inc <- numeric(n + 1L)
  lambda_exc <- numeric(n)
  for (cp in comps) {
    S <- 0  # excitation sum just after the previous event
    for (k in seq_len(n + 1L)) {
      d <- hi[k] - lo[k]
      exc_inc[k] <- exc_inc[k] + (cp$K / cp$beta) * S * (-expm1(-cp$beta * d))
      if (k <= n) {
        S <- S * exp(-cp$beta * d)
        lambda_exc[k] <- lambda_exc[k] + cp$K * S
        S <- S + cp$g[k]
      }
    }
  }

  comp <- base_inc + exc_inc
  comp_inc <- comp[seq_len(n)]
  list(
    t = t,
    dt = dt,
    lambda = mu_ev + lambda_exc,
    comp_inc = comp_inc,
    tau = comp_inc,
    meta = list(T0 = T0, T_end = T_end, comp_tail = comp[n + 1L], ok = TRUE)
  )
}


#' Plot temporal diagnostics from `diag_temporal()`
#'
#' @param diag Output of `diag_temporal()`.
#' @param which Character vector of plots to show. Any of:
#'   `"count"`, `"dt"`, `"lambda"`, `"tau_hist"`, `"tau_qq"`.
#' @param ask Logical; if `TRUE`, prompt between plots.
#' @export
plot_temporal_diag <- function(diag,
                               which = c("count", "dt", "lambda", "tau_hist", "tau_qq"),
                               ask = FALSE) {
  if (is.null(diag) || is.null(diag$t)) stop("plot_temporal_diag(): invalid `diag`.", call. = FALSE)

  t <- diag$t
  if (length(t) == 0L) {
    plot.new()
    text(0.5, 0.5, "No events to plot.")
    return(invisible(NULL))
  }

  old_ask <- par(ask = ask)
  on.exit(par(old_ask), add = TRUE)

  which <- unique(which)

  if ("count" %in% which) {
    plot(t, seq_along(t), type = "s", xlab = "time", ylab = "N(t)",
         main = "Cumulative event count")
  }

  if ("dt" %in% which) {
    dt <- diag$dt
    hist(dt, breaks = "FD", xlab = "Δt", main = "Inter-arrival times", col = "grey90", border = "grey40")
  }

  if ("lambda" %in% which) {
    lam <- diag$lambda
    plot(t, lam, type = "l", xlab = "time", ylab = expression(lambda(t[k]^-1)),
         main = "Intensity at event times", lty = 1)
    # lines(t, lam, lty = 1)
  }

  if ("tau_hist" %in% which) {
    tau <- diag$tau
    hist(tau, breaks = "FD", probability = TRUE, xlab = expression(tau[k]),
         main = "Time-rescaled increments", col = "grey90", border = "grey40")
    xx <- seq(0, max(tau, na.rm = TRUE), length.out = 200)
    lines(xx, dexp(xx), lty = 1)
    legend("topright", legend = c("hist", "Exp(1) density"), bty = "n", lty = c(NA, 1), pch = c(15, NA))
  }

  if ("tau_qq" %in% which) {
    tau <- diag$tau
    tau <- tau[is.finite(tau) & tau >= 0]
    if (length(tau) < 5) {
      plot.new()
      text(0.5, 0.5, "Not enough tau values for QQ plot.")
    } else {
      n <- length(tau)
      theo <- stats::qexp(ppoints(n))
      samp <- sort(tau)
      plot(theo, samp, xlab = "Exp(1) theoretical quantiles", ylab = "Sample quantiles",
           main = "QQ plot: tau vs Exp(1)")
      abline(0, 1, lty = 2)
    }
  }

  invisible(NULL)
}
