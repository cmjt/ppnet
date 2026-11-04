# excitation.R

#' Mark-dependent excitation specification
#'
#' Lets each event's contribution to the intensity depend on its marks:
#' \deqn{\lambda(t) = \mu(t) + K \sum_{t_i < t} g_i e^{-\beta (t - t_i)},
#'   \qquad g_i = n_i^{\alpha} \exp\{\gamma\, h(c_i)\}}
#' where \eqn{n_i} is the number of perps at event \eqn{i} and \eqn{c_i} its
#' connectedness, the number of distinct earlier events sharing a perp with it
#' (see [event_summary()]). With `alpha = gamma = 0` this is the constant
#' excitation `K`, so the plain Hawkes model is nested in every specification.
#'
#' Only bipartite event-perp histories (`mark_type = "ba_bip"` or `"cs_bip"`)
#' have these marks. The marks are observed history for the event times, so
#' the time and mark parts of the likelihood still separate.
#'
#' @param size Logical; include \eqn{n_i^\alpha} (parameter `alpha`).
#' @param connectedness Logical; include \eqn{\exp\{\gamma h(c_i)\}} (parameter `gamma`).
#' @param connectedness_scale `"log1p"` (default) for \eqn{h(c) = \log(1 + c)},
#'   giving \eqn{(1 + c_i)^\gamma}, or `"linear"` for \eqn{h(c) = c}. The log
#'   scale keeps heavily connected events from dominating the fit.
#' @param connectedness_window Only earlier events at most this long before
#'   event \eqn{i} count towards \eqn{c_i}. With the default `Inf`,
#'   connectedness keeps growing for prolific perps as the history lengthens,
#'   which can make a long process explosive when `gamma > 0`.
#' @param connectedness_min_gap Only earlier events more than this long before
#'   event \eqn{i} count towards \eqn{c_i}. Records of the same episode (e.g.
#'   minutes apart) are then left to the kernel, typically the fast component
#'   of `kernel = "exp2"`, instead of making an event look connected. `0` (the
#'   default) counts all.
#' @param path_dependent Logical. `FALSE` (the default) fixes \eqn{c_i} at
#'   event \eqn{i}, from the events before it. `TRUE` makes the excitation
#'   mark-path dependent: \eqn{c_i(t)} also counts the later events, before
#'   \eqn{t}, that share a perp with event \eqn{i} (with the same window and
#'   minimum gap, now after \eqn{t_i}), so an event's effect on the intensity
#'   grows as its perps go on to offend again,
#'   \deqn{\lambda(t) = \mu(t) + K \sum_{t_i < t} n_i^{\alpha}
#'     \exp\{\gamma\, h(c_i(t))\}\, e^{-\beta (t - t_i)}.}
#'   Needs `connectedness = TRUE`. The likelihood stays exact: each change in
#'   \eqn{c_i} acts as an extra jump at the time of the later event. `gamma`
#'   then measures how much an event's effect depends on its perps' later
#'   offending together, which is partly the outcome itself, so it is harder
#'   to read as an effect of connectedness on later offending. Stability is
#'   weaker too: with `gamma > 0` and an infinite window an old event's effect
#'   can keep growing; a finite `connectedness_window` stops each event's
#'   \eqn{c_i} changing once it is that old.
#'
#' @return An object of class `excitation_spec`.
#' @export
#'
#' @examples
#' excitation_coef_names(excitation_spec(size = TRUE, connectedness = TRUE))
excitation_spec <- function(size = TRUE, connectedness = FALSE,
                            connectedness_scale = c("log1p", "linear"),
                            connectedness_window = Inf, connectedness_min_gap = 0,
                            path_dependent = FALSE) {
  connectedness_scale <- match.arg(connectedness_scale)
  stopifnot(is.numeric(connectedness_window), length(connectedness_window) == 1L,
            connectedness_window >= 0, is.numeric(connectedness_min_gap),
            length(connectedness_min_gap) == 1L, connectedness_min_gap >= 0,
            connectedness_min_gap < connectedness_window)
  if (!isTRUE(size) && !isTRUE(connectedness)) {
    stop("excitation_spec() needs `size`, `connectedness` or both; use NULL for constant excitation.",
         call. = FALSE)
  }
  if (isTRUE(path_dependent) && !isTRUE(connectedness)) {
    stop("`path_dependent` needs `connectedness = TRUE`.", call. = FALSE)
  }
  structure(
    list(size = isTRUE(size), connectedness = isTRUE(connectedness),
         connectedness_scale = connectedness_scale,
         connectedness_window = connectedness_window,
         connectedness_min_gap = connectedness_min_gap,
         path_dependent = isTRUE(path_dependent)),
    class = "excitation_spec"
  )
}

#' Names of the excitation coefficients
#'
#' @param spec An `excitation_spec`, or `NULL` for constant excitation.
#' @return Character vector: `"alpha"` and/or `"gamma"`, in design-matrix order.
#' @export
excitation_coef_names <- function(spec) {
  if (is.null(spec)) return(character(0))
  c(if (spec$size) "alpha", if (spec$connectedness) "gamma")
}

# Internal: excitation covariates z_i from per-event perp counts and
# connectedness, one row per event, columns as excitation_coef_names().
.excitation_z <- function(spec, n_perps, connectedness) {
  cols <- list()
  if (spec$size) {
    if (any(n_perps < 1L)) {
      stop("Excitation by size needs at least one perp at every event.", call. = FALSE)
    }
    cols$alpha <- log(n_perps)
  }
  if (spec$connectedness) {
    cols$gamma <- if (spec$connectedness_scale == "log1p") log1p(connectedness) else connectedness
  }
  Z <- matrix(unlist(cols), nrow = length(n_perps))
  colnames(Z) <- names(cols)
  Z
}

# Internal: the excitation design for an events object: Z, the covariates at
# each event's birth, and with `path_dependent` the connectedness updates
# from connectedness_updates_cpp().
.excitation_design <- function(spec, events) {
  window <- spec$connectedness_window
  min_gap <- spec$connectedness_min_gap %||% 0
  s <- event_summary(events, connectedness_window = window, connectedness_min_gap = min_gap)
  out <- list(Z = .excitation_z(spec, s$n_perps, s$connectedness))
  if (isTRUE(spec$path_dependent)) {
    a <- .event_arrays(events)
    f <- .bip_focal(a)
    out$path <- c(list(t = a$t, logn = log(pmax(s$n_perps, 1L))),
                  connectedness_updates_cpp(a$t, a$node_ptr, a$edge_ptr, a$edge_i, a$edge_j,
                                            f$focal, length(a$ids), window, min_gap))
  }
  out
}

# Internal: the jumps of the slow kernel component at each event, G, with
# dG = dG / da (columns as excitation_coef_names()) and dGb = dG / dbeta
# (length 0 when G does not depend on beta). Without path dependence G is the
# excitation at birth, exp(Z a).
.excitation_jumps <- function(spec, design, a, beta) {
  if (is.null(design$path)) {
    g <- exp(drop(design$Z %*% a))
    return(list(G = g, dG = design$Z * g, dGb = numeric(0)))
  }
  pd <- design$path
  path_jumps_cpp(pd$t, pd$logn, pd$c0, pd$upd_ptr, pd$upd_j, pd$upd_c,
                 if (spec$size) a[["alpha"]] else 0, a[["gamma"]], spec$size,
                 spec$connectedness_scale == "log1p", beta)
}

# Internal: pull the coefficient vector out of a params list, in column order.
.excitation_coefs <- function(spec, params) {
  nms <- excitation_coef_names(spec)
  missing <- setdiff(nms, names(params))
  if (length(missing) > 0L) {
    stop("Missing excitation coefficients: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  vapply(nms, function(nm) as.numeric(params[[nm]]), numeric(1))
}
