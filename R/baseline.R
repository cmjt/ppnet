# baseline.R

#' Time-varying baseline specification
#'
#' Describes a baseline rate of the form
#' \deqn{\mu(t) = \mu \exp\{x(t)^\top b\}}
#' where \eqn{x(t)} holds an optional log-linear trend and Fourier terms for
#' each period. With all coefficients at zero the baseline is the constant
#' `mu`, so a constant-baseline fit is nested in every specification.
#'
#' Coefficients are passed as ordinary scalar parameters named by
#' [baseline_coef_names()]: `b_trend`, then `b_<period>_cos<h>` and
#' `b_<period>_sin<h>` for each period and harmonic `h`.
#'
#' @param trend Logical; include the trend term `(t - trend_origin) / trend_scale`.
#'   `mu` is then the rate at `trend_origin` (with periodic terms at zero).
#' @param periods Named numeric vector of cycle lengths in the time units of
#'   the events, e.g. `c(daily = 1, weekly = 7, yearly = 365.25)`. Phases are
#'   measured from time 0.
#' @param harmonics Number of Fourier harmonics, either one value for all
#'   periods or one per period.
#' @param trend_origin,trend_scale Origin and unit of the trend covariate. They
#'   are fixed by the user, not the data window, so the simulator and the
#'   likelihood agree on what `b_trend` means.
#'
#' @return An object of class `baseline_spec`.
#' @export
#'
#' @examples
#' spec <- baseline_spec(trend = TRUE, periods = c(weekly = 7), harmonics = 2)
#' baseline_coef_names(spec)
baseline_spec <- function(trend = FALSE, periods = NULL, harmonics = 1L,
                          trend_origin = 0, trend_scale = 1) {
  periods <- periods %||% numeric(0)
  if (length(periods) > 0L) {
    if (is.null(names(periods)) || any(names(periods) == "") || anyDuplicated(names(periods))) {
      stop("`periods` must be a named numeric vector with unique names.", call. = FALSE)
    }
    if (any(!is.finite(periods) | periods <= 0)) stop("`periods` must be > 0.", call. = FALSE)
    harmonics <- rep_len(as.integer(harmonics), length(periods))
    if (any(harmonics < 1L)) stop("`harmonics` must be >= 1.", call. = FALSE)
    names(harmonics) <- names(periods)
  } else {
    harmonics <- integer(0)
  }
  if (!is.finite(trend_scale) || trend_scale <= 0) stop("`trend_scale` must be > 0.", call. = FALSE)

  structure(
    list(trend = isTRUE(trend), periods = periods, harmonics = harmonics,
         trend_origin = trend_origin, trend_scale = trend_scale),
    class = "baseline_spec"
  )
}

#' Names of the baseline coefficients
#'
#' @param spec A `baseline_spec`, or `NULL` for a constant baseline.
#' @return Character vector of parameter names, in design-matrix column order.
#' @export
baseline_coef_names <- function(spec) {
  if (is.null(spec)) return(character(0))
  nms <- if (spec$trend) "b_trend" else character(0)
  for (p in names(spec$periods)) {
    h <- seq_len(spec$harmonics[[p]])
    nms <- c(nms, as.vector(rbind(paste0("b_", p, "_cos", h), paste0("b_", p, "_sin", h))))
  }
  nms
}

# Internal: design matrix x(t), one row per time, columns as baseline_coef_names().
.baseline_design <- function(spec, t) {
  cols <- list()
  if (spec$trend) cols[[1L]] <- (t - spec$trend_origin) / spec$trend_scale
  for (p in names(spec$periods)) {
    for (h in seq_len(spec$harmonics[[p]])) {
      a <- 2 * pi * h * t / spec$periods[[p]]
      cols <- c(cols, list(cos(a), sin(a)))
    }
  }
  X <- matrix(unlist(cols), nrow = length(t))
  colnames(X) <- baseline_coef_names(spec)
  X
}

# Internal: pull the coefficient vector out of a params list, in column order.
.baseline_coefs <- function(spec, params) {
  nms <- baseline_coef_names(spec)
  missing <- setdiff(nms, names(params))
  if (length(missing) > 0L) {
    stop("Missing baseline coefficients: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  vapply(nms, function(nm) as.numeric(params[[nm]]), numeric(1))
}

# Internal: Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch).
.gauss_legendre <- function(n) {
  i <- seq_len(n - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1L, ]^2)
}

# Internal: quadrature for the integral of mu(t) over [T0, T_end]. Panels are
# at most half the shortest cycle, so 8 nodes per panel are accurate to
# roughly 1e-10 relative for moderate seasonal amplitudes.
.baseline_quadrature <- function(spec, T0, T_end, nodes = 8L) {
  span <- T_end - T0
  if (span <= 0) return(list(X = .baseline_design(spec, numeric(0)), w = numeric(0)))

  h <- span
  if (length(spec$periods) > 0L) h <- min(h, min(spec$periods / spec$harmonics) / 2)
  n_panels <- ceiling(span / h)
  h <- span / n_panels

  gl <- .gauss_legendre(nodes)
  left <- T0 + (seq_len(n_panels) - 1L) * h
  s <- as.vector(outer((gl$x + 1) * h / 2, left, `+`))
  w <- rep(gl$w * h / 2, n_panels)
  list(X = .baseline_design(spec, s), w = w)
}

# Internal: an upper bound on mu(t) over [T0, T_end], for thinning.
.baseline_upper <- function(spec, params, T0, T_end) {
  b <- .baseline_coefs(spec, params)
  bound <- 0
  if (spec$trend) {
    ends <- (c(T0, T_end) - spec$trend_origin) / spec$trend_scale
    bound <- max(b[["b_trend"]] * ends)
  }
  for (p in names(spec$periods)) {
    for (h in seq_len(spec$harmonics[[p]])) {
      bound <- bound + sqrt(b[[paste0("b_", p, "_cos", h)]]^2 + b[[paste0("b_", p, "_sin", h)]]^2)
    }
  }
  params$mu * exp(bound)
}

# Internal: integral of exp(x(t)' b) over each interval [lo_i, hi_i], with
# the same panel width and Gauss-Legendre rule as .baseline_quadrature().
.baseline_integral <- function(spec, b, lo, hi, nodes = 8L) {
  len <- hi - lo
  h <- max(len, 0)
  if (length(spec$periods) > 0L) h <- min(h, min(spec$periods / spec$harmonics) / 2)
  if (h <= 0) return(numeric(length(lo)))
  n_panels <- pmax(1L, ceiling(len / h))
  width <- len / n_panels

  idx <- rep(seq_along(lo), n_panels)
  left <- lo[idx] + (sequence(n_panels) - 1L) * width[idx]
  gl <- .gauss_legendre(nodes)
  s <- as.vector(outer((gl$x + 1) / 2, width[idx]) + rep(left, each = nodes))
  w <- as.vector(outer(gl$w / 2, width[idx]))
  vals <- w * exp(drop(.baseline_design(spec, s) %*% b))
  out <- numeric(length(lo))
  sums <- rowsum(vals, rep(idx, each = nodes))
  out[as.integer(rownames(sums))] <- sums[, 1]
  out
}
