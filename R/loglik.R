# loglik.R

# Internal: initialize recursion cache for an exponential Hawkes kernel.
time_cache_init <- function(T0) {
  list(last_t = T0, S = 0.0, n = 0L)
}

time_cache_step <- function(cache, t, beta) {
  dt <- t - cache$last_t
  if (dt < 0) stop("times must be increasing")
  cache$S <- exp(-beta * dt) * cache$S
  cache$last_t <- t
  cache
}

temporal_intensity_from_cache <- function(cache, mu, K, beta, kernel_density = FALSE) {
  scale <- if (kernel_density) beta else 1.0
  mu + K * scale * cache$S
}

#' Log-likelihood for a HawkesNet model
#'
#' Computes the log-likelihood for the observed event sequence under an
#' exponential temporal Hawkes kernel and a chosen mark model (`mark_type`).
#'
#' @inheritParams hawkesnet-lik-args
#' @param params Named list of model parameters on the natural scale. Must
#'   contain `mu`, `K`, and `beta`, plus any mark parameters required by
#'   `mark_type` (e.g. `beta_edges` for BA).
#' @param net_init_fun Function to initialize the pre-event network if `net0` is
#'   `NULL`. Defaults to `net_init`.
#' @param net0 Optional initial network object to start recursion from.
#' @param net_step_fun Function used to update the network after each event.
#'   Defaults to `net_add_event`.
#' @param formula_rhs String, the CS kernel formula
#'
#' @return A numeric scalar log-likelihood. Returns `-Inf` for invalid parameter
#'   values or if the likelihood is undefined.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' params_true <- list(mu = 0.5, K = 0.5, beta = 0.5, beta_edges = 2)
#' sim <- sim_hawkesNet(params = params_true, T_end = 10, mark_type = "ba")
#' loglik(sim$ev, params_true, mark_type = "ba")
# Internal: build a compiled log-likelihood function that reuses mark caches
# (critical for optimizers that evaluate the objective many times).
.build_loglik_fn <- function(
    events, T_end = NULL, T0 = 0,
    mark_type = "ba",
    debug = FALSE,
    net_init_fun = net_init, net0 = NULL,
    net_step_fun = net_add_event,
    formula_rhs = NULL,
    truncation = Inf,
    baseline = NULL,
    excitation = NULL,
    attachment = NULL,
    kernel = c("exp", "exp2"),
    ...
) {
  stopifnot(inherits(events, "events"))
  event_times <- events$times$t
  kernel <- match.arg(kernel)

  # In the case of no events
  if (length(event_times) == 0L) {
    return(function(params) 0.0)
  }

  # If the user does not specify T_end, assume the event window lasts exactly
  # up until the last event
  if (is.null(T_end)) T_end <- max(event_times)
  if (T_end < max(event_times)) stop("T_end must be >= max event time", call. = FALSE)

  if (!is.null(baseline) && !inherits(baseline, "baseline_spec")) {
    stop("`baseline` must be NULL or a baseline_spec().", call. = FALSE)
  }
  # The baseline design is fixed by the data, so build it once: rows at the
  # event times for the log-intensity, and quadrature rows for its integral.
  if (!is.null(baseline)) {
    X_ev <- .baseline_design(baseline, event_times)
    quad <- .baseline_quadrature(baseline, T0, T_end)
  }

  if (!is.null(excitation) && !inherits(excitation, "excitation_spec")) {
    stop("`excitation` must be NULL or an excitation_spec().", call. = FALSE)
  }
  # Excitation covariates are marks of the observed history, fixed by the data.
  exc_design <- if (!is.null(excitation)) .excitation_design(excitation, events)

  # Grab the lists of node and edge arrivals by event
  arrivals <- nodes_by_event(events)
  edges <- edges_by_event(events)

  # Compile the mark kernel function.
  # This allow us to calculate the mark-likelihood contribution as a simple
  # function of the parameters.
  # The reason this works is that the network changes are fixed regardless of
  # the parameter values.
  # The default history (empty start, net_add_event steps) is walked with
  # arrays in C++; a custom start or step function needs the network objects.
  use_arrays <- is.null(net0) && identical(net_init_fun, net_init) &&
    identical(net_step_fun, net_add_event)

  if (!is.null(attachment)) {
    if (!inherits(attachment, "attachment_spec")) {
      stop("`attachment` must be NULL or an attachment_spec().", call. = FALSE)
    }
    if (mark_type != "ba_bip") stop("`attachment` needs mark_type = 'ba_bip'.", call. = FALSE)
    mark_type <- "attach"
  }

  mark_logpmf <- switch(
    mark_type,
    "attach" = compile_logpmf_attach(events, attachment),
    "ba" = if (use_arrays) compile_logpmf_ba_array(events, delta = 0.001, ...) else compile_logpmf_ba(
      events = events,
      net_init_fun = net_init_fun,
      net0 = net0,
      net_step_fun = net_step_fun,
      delta = 0.001,
      ...
    ),
    "ba_bip" = if (use_arrays) compile_logpmf_ba_bip_array(events, delta = 0.001, ...) else compile_logpmf_ba_bip(
      events = events,
      net_init_fun = net_init_fun,
      net0 = net0,
      net_step_fun = net_step_fun,
      delta = 0.001,
      ...
    ),
    "cs" = compile_logpmf_cs(
      events = events,
      formula_rhs = formula_rhs,
      truncation = truncation,
      net_init_fun = net_init_fun,
      net0 = net0,
      net_step_fun = net_step_fun,
      model_cache = NULL,
      ...
    ),
    "cs_bip" = compile_logpmf_cs_bip(
      events = events,
      formula_rhs = formula_rhs,
      truncation = truncation,
      net_init_fun = net_init_fun,
      net0 = net0,
      net_step_fun = net_step_fun,
      model_cache = NULL,
      ...
    ),
    stop("mark_type must be one of 'ba', 'cs', 'ba_bip', 'cs_bip'", call. = FALSE)
  )



  # Debug: precompute layout from full network derived from events
  debug_layout <- NULL
  if (isTRUE(debug)) {
    full_net <- create_net_from_events(events)
    debug_layout <- .debug_prepare_full_layout(full_net)
  }

  # Mark evaluators that can give analytic gradients take a `gradient` argument
  mark_has_grad <- "gradient" %in% names(formals(mark_logpmf))

  # Return the likelihood function over params (natural scale). With
  # `gradient = TRUE` the value carries attribute "gradient", the derivatives
  # with respect to mu, K, beta, the baseline and excitation coefficients and the mark
  # parameters; it is NULL if the mark model has no analytic gradient.
  function(params, gradient = FALSE) {
    mu <- params$mu
    K <- params$K
    beta <- params$beta

    if (!is.finite(mu) || !is.finite(K) || !is.finite(beta) || mu <= 0 || K < 0 || beta <= 0)
      return(-Inf)

    # Fast component: must decay faster than the slow one, so they can't swap.
    K_fast <- 0
    beta_fast <- 1
    if (kernel == "exp2") {
      K_fast <- params$K_fast
      beta_fast <- params$beta_fast
      if (is.null(K_fast) || is.null(beta_fast)) {
        stop("kernel = 'exp2' needs parameters K_fast and beta_fast.", call. = FALSE)
      }
      if (!is.finite(K_fast) || !is.finite(beta_fast) || K_fast <= 0 || beta_fast <= beta)
        return(-Inf)
    }

    if (is.null(baseline)) {
      mu_ev <- rep.int(mu, length(event_times))
      baseline_int <- mu * (T_end - T0)
    } else {
      b <- .baseline_coefs(baseline, params)
      if (any(!is.finite(b))) return(-Inf)
      mu_ev <- mu * exp(drop(X_ev %*% b))
      mu_quad <- mu * quad$w * exp(drop(quad$X %*% b))
      baseline_int <- sum(mu_quad)
    }

    # mark loglik ONCE (it already sums over events)
    mark_ll_total <- if (mark_has_grad) mark_logpmf(params, gradient = gradient) else mark_logpmf(params)
    if (!is.finite(mark_ll_total)) return(-Inf)

    if (is.null(excitation)) {
      jumps <- list(G = rep.int(1, length(event_times)),
                    dG = matrix(0, nrow = length(event_times), ncol = 0L), dGb = numeric(0))
    } else {
      a <- .excitation_coefs(excitation, params)
      if (any(!is.finite(a))) return(-Inf)
      jumps <- .excitation_jumps(excitation, exc_design, a, beta)
    }

    tm <- time_ll_cpp(event_times, mu_ev, jumps$G, jumps$dG, jumps$dGb, K, beta, T_end,
                      K_fast, beta_fast)
    if (!is.finite(tm$ll)) return(-Inf)

    # Debug plot (two panels), replaying the history one event at a time
    if (isTRUE(debug)) {
      net <- if (is.null(net0)) net_init_fun() else net0
      ll_time <- cumsum(tm$log_lambda)
      for (idx in seq_along(event_times)) {
        net <- net_step_fun(net, arrivals[[idx]], edges[[idx]], event_times[idx])
        .debug_plot_step(net, debug_layout, event_times[idx], 0, ll_time[idx],
                         ll_time[idx] + mark_ll_total, arrivals[[idx]], edges[[idx]],
                         event_id = events$times$event_id[idx], edge_probs = NULL)
        readline("Press Enter / Return to continue to next event...")
      }
    }

    ll <- tm$ll - baseline_int + as.numeric(mark_ll_total)

    if (gradient) {
      mark_grad <- attr(mark_ll_total, "gradient")
      if (mark_has_grad && !is.null(mark_grad)) {
        # d/d log mu_i of sum(log lambda) is mu_i / lambda_i
        g_ev <- mu_ev * tm$inv_lambda
        g <- c(mu = (sum(g_ev) - baseline_int) / mu, K = tm$dK, beta = tm$dbeta)
        if (!is.null(baseline)) {
          g_b <- drop(crossprod(X_ev, g_ev) - crossprod(quad$X, mu_quad))
          g <- c(g, g_b)
        }
        if (!is.null(excitation)) g <- c(g, stats::setNames(tm$da, excitation_coef_names(excitation)))
        if (kernel == "exp2") g <- c(g, K_fast = tm$dK2, beta_fast = tm$dbeta2)
        attr(ll, "gradient") <- c(g, mark_grad)
      }
    }
    ll
  }
}

#' Log-likelihood for a HawkesNet model
#'
#' Computes the log-likelihood for the observed event sequence under an
#' exponential temporal Hawkes kernel and a chosen mark model (`mark_type`).
#'
#' @inheritParams hawkesnet-lik-args
#' @param params Named list of model parameters on the natural scale. Must
#'   contain `mu`, `K`, and `beta`, plus any mark parameters required by
#'   `mark_type` (e.g. `beta_edges` for BA).
#' @param net_init_fun Function to initialize the pre-event network if `net0` is
#'   `NULL`. Defaults to `net_init`.
#' @param net0 Optional initial network object to start recursion from.
#' @param net_step_fun Function used to update the network after each event.
#'   Defaults to `net_add_event`.
#' @param formula_rhs String, the CS kernel formula
#' @param truncation CS only: number of most recent old nodes each new node
#'   can attach to.
#' @param baseline Optional [baseline_spec()] for a time-varying baseline
#'   \eqn{\mu(t) = \mu \exp\{x(t)^\top b\}}. `params` must then also
#'   contain the coefficients named by [baseline_coef_names()]. `NULL` (the
#'   default) gives the constant baseline `mu`.
#' @param excitation Optional [excitation_spec()] letting each event's
#'   excitation depend on its size and connectedness (bipartite mark types
#'   only). `params` must then also contain `alpha` and/or `gamma`. `NULL`
#'   (the default) gives constant excitation `K`.
#' @param attachment Optional [attachment_spec()] replacing the BA-bip rule
#'   for which old perps join each event (`mark_type = "ba_bip"` only).
#'   `params` must then contain [attachment_param_names()].
#' @param kernel `"exp"` (default) for one exponential component, or `"exp2"`
#'   to add a fast component \eqn{K_{fast} \sum_i e^{-\beta_{fast}(t - t_i)}}
#'   with unit jumps; `K`, `beta` and any `excitation` stay with the slow
#'   component. `params` must then contain `K_fast` and `beta_fast`, and the
#'   likelihood is `-Inf` unless `beta_fast > beta`.
#'
#' @return A numeric scalar log-likelihood. Returns `-Inf` for invalid parameter
#'   values or if the likelihood is undefined.
#'
#' @export
loglik <- function(
    events, params, T_end = NULL, T0 = 0,
    mark_type = "ba",
    debug = FALSE,
    net_init_fun = net_init, net0 = NULL,
    net_step_fun = net_add_event,
    formula_rhs = NULL,
    truncation = Inf,
    baseline = NULL,
    excitation = NULL,
    attachment = NULL,
    kernel = c("exp", "exp2"),
    ...
) {
  ll_fn <- .build_loglik_fn(
    events = events, T_end = T_end, T0 = T0,
    mark_type = mark_type, debug = debug,
    net_init_fun = net_init_fun, net0 = net0,
    net_step_fun = net_step_fun,
    formula_rhs = formula_rhs,
    truncation = truncation,
    baseline = baseline,
    excitation = excitation,
    attachment = attachment,
    kernel = kernel,
    ...
  )
  ll_fn(params)
}
