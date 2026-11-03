# plot_methods.R

.temporal_types <- c("count", "dt", "lambda", "tau_hist", "tau_qq")
.network_types  <- c("size", "edges_vs_nodes", "degree_hist", "degree_ccdf")

.is_events <- function(x) {
  inherits(x, "events") || (!is.null(x$times) && is.data.frame(x$times) && "t" %in% names(x$times))
}

#' Plot diagnostics for an events object
#'
#' @param x An events object created by `make_events()`.
#' @param type Which diagnostic plot to produce. Temporal: `count`, `dt`, `lambda`,
#'   `tau_hist`, `tau_qq`. Network: `size`, `edges_vs_nodes`, `degree_hist`,
#'   `degree_ccdf`.
#' @param params Optional named list of parameters. Required for `lambda`,
#'   `tau_hist`, and `tau_qq`.
#' @param T_end Optional end time for temporal diagnostics.
#' @param T0 Start time for temporal diagnostics.
#' @param net_init_fun Function returning an initial `network` for network diagnostics.
#' @param net0 Optional initial network override.
#' @param net_step_fun Network update function for network diagnostics.
#' @param loglog For `degree_ccdf`, plot on log-log axes.
#' @param baseline,excitation,kernel Temporal model settings, as for
#'   [loglik()]; used by `lambda`, `tau_hist` and `tau_qq`.
#' @param ... Passed through to diagnostic helpers.
#'
#' @return Invisibly returns the computed diagnostics object.
#' @export
#' @method plot events
plot.events <- function(x,
                        type = c(.temporal_types, .network_types),
                        params = NULL,
                        T_end = NULL,
                        T0 = 0,
                        net_init_fun = net_init,
                        net0 = NULL,
                        net_step_fun = net_add_event,
                        loglog = TRUE,
                        baseline = NULL,
                        excitation = NULL,
                        kernel = "exp",
                        ...) {
  type <- match.arg(type)

  if (!.is_events(x)) {
    stop("plot.events(): `x` does not look like an events object (needs `x$times$t`).", call. = FALSE)
  }

  # ---- Temporal diagnostics ----
  if (type %in% .temporal_types) {
    if (type %in% c("lambda", "tau_hist", "tau_qq")) {
      if (is.null(params)) {
        stop("plot(events, type = '", type, "'): please supply `params` (needs mu, K, beta).",
             call. = FALSE)
      }
      d <- diag_temporal(events = x, params = params, T_end = T_end, T0 = T0,
                         baseline = baseline, excitation = excitation, kernel = kernel)
    } else {
      # count/dt don't need params; give a light-weight diag object
      t <- as.numeric(x$times$t)
      d <- list(t = t, dt = c(t[1] - T0, diff(t)))
    }

    plot_temporal_diag(d, which = type, ask = FALSE)
    return(invisible(d))
  }

  # ---- Network diagnostics ----
  d <- diag_network(events = x, net_init_fun = net_init_fun, net0 = net0, net_step_fun = net_step_fun)
  if (type == "degree_ccdf") {
    plot_network_diag(d, which = type, ask = FALSE, loglog = loglog)
  } else {
    plot_network_diag(d, which = type, ask = FALSE)
  }
  invisible(d)
}


#' Plot diagnostics for a fitted HawkesNet model
#'
#' The temporal model settings (`baseline`, `excitation`, `kernel`, `T0`,
#' `T_end`) are taken from the fit unless given.
#'
#' @param x A fitted model object. Recommended structure is the list returned by
#'   `fit_hawkesNet()` with `class(x) = 'hawkesNet_fit'` and the original events
#'   stored as `x$events`.
#' @inheritParams plot.events
#'
#' @return Invisibly returns the computed diagnostics object.
#' @export
#' @method plot hawkesNet_fit
plot.hawkesNet_fit <- function(x,
                               type = c(.temporal_types, .network_types),
                               params = NULL,
                               T_end = NULL,
                               T0 = 0,
                               net_init_fun = net_init,
                               net0 = NULL,
                               net_step_fun = net_add_event,
                               loglog = TRUE,
                               baseline = x$model$baseline,
                               excitation = x$model$excitation,
                               kernel = x$model$kernel %||% "exp",
                               ...) {
  type <- match.arg(type)
  if (missing(T_end) && !is.null(x$model$T_end)) T_end <- x$model$T_end
  if (missing(T0) && !is.null(x$model$T0)) T0 <- x$model$T0

  ev <- x$events %||% x$ev
  if (is.null(ev) || !.is_events(ev)) {
    stop(
      "plot(fit, ...): couldn't find events on `fit`.\n",
      "Recommended: store the events used for fitting as `fit$events` and set class to 'hawkesNet_fit'.\n",
      "Workaround: call `plot(events, type = ..., params = fit$par)`.",
      call. = FALSE
    )
  }

  if (is.null(params)) {
    params <- x$par %||% x$params
  }

  # Temporal plots that need params
  if (type %in% c("lambda", "tau_hist", "tau_qq") && is.null(params)) {
    stop("plot(fit, type = '", type, "'): couldn't find parameters on `fit` (expected `fit$par`).",
         call. = FALSE)
  }

  plot(ev,
       type = type,
       params = params,
       T_end = T_end,
       T0 = T0,
       net_init_fun = net_init_fun,
       net0 = net0,
       net_step_fun = net_step_fun,
       loglog = loglog,
       baseline = baseline,
       excitation = excitation,
       kernel = kernel,
       ...)
}
