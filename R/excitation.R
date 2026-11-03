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
#'
#' @return An object of class `excitation_spec`.
#' @export
#'
#' @examples
#' excitation_coef_names(excitation_spec(size = TRUE, connectedness = TRUE))
excitation_spec <- function(size = TRUE, connectedness = FALSE,
                            connectedness_scale = c("log1p", "linear"),
                            connectedness_window = Inf) {
  connectedness_scale <- match.arg(connectedness_scale)
  stopifnot(is.numeric(connectedness_window), length(connectedness_window) == 1L,
            connectedness_window >= 0)
  if (!isTRUE(size) && !isTRUE(connectedness)) {
    stop("excitation_spec() needs `size`, `connectedness` or both; use NULL for constant excitation.",
         call. = FALSE)
  }
  structure(
    list(size = isTRUE(size), connectedness = isTRUE(connectedness),
         connectedness_scale = connectedness_scale,
         connectedness_window = connectedness_window),
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

# Internal: the excitation design for an events object.
.excitation_design <- function(spec, events) {
  s <- event_summary(events, connectedness_window = spec$connectedness_window)
  .excitation_z(spec, s$n_perps, s$connectedness)
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
