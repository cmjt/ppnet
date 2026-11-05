#' Simulate a HawkesNet model
#'
#' Simulates event times by thinning and, at each event, its marks (new nodes
#' and edges) from the chosen mark model.
#'
#' @param params Named list of parameters on the natural scale, as for
#'   [loglik()].
#' @param T_end End of the simulation window.
#' @param T0 Start of the simulation window.
#' @param delta Degree offset for the original BA mark models.
#' @param mark_type One of `"ba"`, `"cs"`, `"ba_bip"`, `"cs_bip"`.
#' @param return_ev Logical; return the simulated history as an events object.
#' @param debug Logical; enable debug plotting.
#' @param net_init_fun,net_step_fun Functions that create the starting network
#'   and add one event to it.
#' @param baseline Optional [baseline_spec()] for a time-varying baseline; its
#'   coefficients go in `params`.
#' @param excitation Optional [excitation_spec()] for size- and
#'   connectedness-dependent excitation (bipartite mark types only); its
#'   coefficients go in `params`.
#' @param attachment Optional [attachment_spec()] for which old perps join
#'   each event (`mark_type = "ba_bip"` only); its parameters go in `params`.
#' @param kernel `"exp"` or `"exp2"`; see [loglik()].
#' @param random_effect_shape With an `attachment` spec without `random_effect`
#'   only:
#'   finite values give each perp a fixed propensity
#'   \eqn{u \sim \text{Gamma}(\text{shape}, \text{shape})} (mean 1, variance
#'   1 / shape) multiplying its join rate, as `random_effect_var` does in a spec with
#'   `random_effect = TRUE`. The fitted model then has no such term, so this
#'   simulates unmodelled differences between perps, e.g. to check what the
#'   co-offender effect picks up.
#' @param return_net Logical; also return the final `network` object. With
#'   `attachment` it is built once at the end, which is slow for long
#'   histories, so set `FALSE` when only `ev` is needed.
#' @param history Optional events object (bipartite, as from [make_events()]
#'   or a fit's `events`) to simulate on from: the simulation starts at `T0`
#'   (by default the last history time) with the history's events in the
#'   intensity and its perps, attendance and co-offending in the mark state,
#'   so it continues the observed process (for forecasts and goodness-of-fit
#'   checks). The returned `ev` holds the history followed by the simulated
#'   events, whose node ids start with `"s"` (attachment model). With a
#'   random effect, the history perps' effects are drawn from their
#'   distribution, not conditioned on their history.
#' @param burn_in Length of a burn-in before `T0`: the simulation starts at
#'   `T0 - burn_in` with nothing before it, and the events before `T0` are
#'   returned as the history, so the process at `T0` is not starting from
#'   empty. Not with `history`.
#' @param max_events,max_seconds Stop the simulation when it has simulated
#'   this many events (after any history) or run this long, e.g. for a model
#'   that explodes. The output then has `exploded = TRUE` and `stop` saying
#'   which limit was reached, and holds the events up to that point.
#' @param verbose Logical; if `TRUE`, report the simulation time with
#'   [message()].
#' @param ... Passed to the mark simulator (e.g. `formula_rhs` for the CS
#'   mark types).
#'
#' @return A list with `events` (the number `n` and times `t` of the
#'   simulated events after any history or burn-in), `accept_probs` (thinning
#'   acceptance probabilities), `ev` (an events object with any history and the
#'   simulated events, if `return_ev`), `net` (the final network, if
#'   `return_net`), `history_end` and `n_history` (the time the history or
#'   burn-in ends and its number of events; `NULL` and 0 without one),
#'   `stop` (`"T_end"`, `"max_events"` or `"max_seconds"`), `exploded`
#'   (`TRUE` unless `T_end` was reached) and `t_stop` (the time reached).
#'   [simulate()] on a fit calls this with the fitted parameters.
#' @export
#'
#' @examples
#' sim <- sim_hawkesNet(list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3,
#'                           lambda_new = 0.7), T_end = 10, mark_type = "ba_bip")
#' sim$events$n
sim_hawkesNet <- function(
    params, T_end, T0 = 0,
    delta = 0.001,
    mark_type = "ba",
    return_ev = TRUE,
    debug = FALSE,
    net_init_fun = net_init,
    net_step_fun = net_add_event,
    baseline = NULL,
    excitation = NULL,
    attachment = NULL,
    kernel = c("exp", "exp2"),
    random_effect_shape = Inf,
    return_net = TRUE,
    history = NULL,
    burn_in = 0,
    max_events = 1e6,
    max_seconds = Inf,
    verbose = FALSE,
    ...
) {
  start_time <- proc.time()
  stopifnot(is.numeric(burn_in), length(burn_in) == 1L, burn_in >= 0,
            is.numeric(max_events), max_events >= 1, is.numeric(max_seconds), max_seconds > 0)

  # A history (observed events) or a burn-in comes before T0; the simulation
  # proper starts at T0 and its own events are those after the history.
  history_end <- NULL
  if (!is.null(history)) {
    if (burn_in > 0) stop("Give either `history` or `burn_in`, not both.", call. = FALSE)
    if (!inherits(history, "events")) stop("`history` must be an events object.", call. = FALSE)
    if (!mark_type %in% c("ba_bip", "cs_bip")) {
      stop("`history` needs a bipartite mark_type ('ba_bip' or 'cs_bip').", call. = FALSE)
    }
    t_hist_max <- max(history$times$t)
    if (missing(T0)) T0 <- t_hist_max
    if (T0 < t_hist_max) stop("`T0` must be at or after the last history event.", call. = FALSE)
    history_end <- T0
  } else if (burn_in > 0) {
    history_end <- T0
    T0 <- T0 - burn_in
  }

  mu <- params$mu
  K <- params$K
  beta <- params$beta
  if (!is.finite(mu) || mu <= 0) stop("mu>0 required")
  if (!is.finite(K)  || K  < 0) stop("K>=0 required")
  if (!is.finite(beta) || beta <= 0) stop("beta>0 required")
  kernel <- match.arg(kernel)
  K_fast <- 0
  beta_fast <- 1
  if (kernel == "exp2") {
    K_fast <- params$K_fast
    beta_fast <- params$beta_fast
    if (is.null(K_fast) || !is.finite(K_fast) || K_fast < 0) stop("K_fast>=0 required")
    if (is.null(beta_fast) || !is.finite(beta_fast) || beta_fast <= beta) stop("beta_fast>beta required")
  }
  # Slow component (cache) plus fast component (cache_f, unit jumps).
  intensity <- function(cache, cache_f, mu_t) {
    temporal_intensity_from_cache(cache, mu = mu_t, K = K, beta = beta) + K_fast * cache_f$S
  }

  # Baseline rate at time t, and a bound on it over the window for thinning.
  # The excitation only decays between events, so mu_max + K * S bounds the
  # intensity from any time until the next accepted event.
  if (is.null(baseline)) {
    mu_at <- function(t) mu
    mu_max <- mu
  } else {
    b <- .baseline_coefs(baseline, params)
    mu_at <- function(t) mu * exp(sum(.baseline_design(baseline, t) * b))
    mu_max <- .baseline_upper(baseline, params, T0, T_end)
  }

  # Mark-dependent jump sizes need each event's perps and the earlier events
  # they were at, kept as perp id -> event indices.
  if (!is.null(excitation)) {
    if (!mark_type %in% c("ba_bip", "cs_bip")) {
      stop("`excitation` needs a bipartite mark_type ('ba_bip' or 'cs_bip').", call. = FALSE)
    }
    a_exc <- .excitation_coefs(excitation, params)
    perp_history <- new.env(hash = TRUE)
    # Path dependence: each event's perp count and current connectedness, so
    # its excitation can be updated when later events link to it.
    path_dep <- isTRUE(excitation$path_dependent)
    n_perps_ev <- numeric(64)
    conn_ev <- numeric(64)
    g_of <- function(n, cn) exp(sum(.excitation_z(excitation, n, cn) * a_exc))
  }

  # Mark kernel dispatch
  if (mark_type == "ba") {
    sim_mark_fun <- sim_mark_ba
  } else if (mark_type == "ba_bip") {
    sim_mark_fun <- sim_mark_ba_bip
  } else if (mark_type == "cs") {
    sim_mark_fun <- sim_mark_cs
  } else if (mark_type == "cs_bip") {
    sim_mark_fun <- sim_mark_cs_bip
  } else {
    stop("mark_type must be one of 'ba', 'cs', 'ba_bip', 'cs_bip'")
  }
  if (!is.null(attachment)) {
    if (mark_type != "ba_bip") stop("`attachment` needs mark_type = 'ba_bip'.", call. = FALSE)
    sim_mark_fun <- .make_sim_mark_attach(attachment, random_effect_shape = random_effect_shape)
  }


  # CS: model cache for ernm
  cs_model_cache <- NULL
  if (mark_type %in% c("cs", "cs_bip")) {
    cs_model_cache <- list(rhs = NULL, model = NULL)
  }

  t <- T0

  # Grown by doubling: appending with c() each iteration is quadratic.
  events_t <- numeric(64)
  n_events <- 0L
  accept_probs <- numeric(64)
  n_props <- 0L
  grow <- function(x, n) if (n > length(x)) c(x, numeric(length(x))) else x
  nodes_rows <- list()
  edges_rows <- list()
  cache <- time_cache_init(T0 = T0)
  cache_f <- time_cache_init(T0 = T0)
  net <- NULL
  n_history <- 0L

  if (!is.null(history)) {
    th <- history$times$t
    n_history <- length(th)
    events_t <- c(th, numeric(64))
    n_events <- n_history
    # The history's events in the intensity at T0 (its excitation jumps
    # include, with path dependence, the later changes up to T0).
    G <- if (is.null(excitation)) rep(1, n_history) else {
      .excitation_jumps(excitation, .excitation_design(excitation, history), a_exc, beta)$G
    }
    cache$S <- sum(G * exp(-beta * (T0 - th)))
    cache_f$S <- sum(exp(-beta_fast * (T0 - th)))
    # Its nodes and edges, with times, start the output.
    at <- function(event_id) th[match(event_id, history$times$event_id)]
    nodes_rows[[1L]] <- data.frame(id = as.character(history$nodes$id),
                                   role = as.character(history$nodes$role),
                                   time = at(history$nodes$event_id), stringsAsFactors = FALSE)
    edges_rows[[1L]] <- data.frame(i = as.character(history$edges$i),
                                   j = as.character(history$edges$j),
                                   time = at(history$edges$event_id), stringsAsFactors = FALSE)
    if (!is.null(excitation)) {
      hp <- .history_perps(history)
      for (k in seq_len(n_history)) {
        for (v in c(hp$old[[k]], hp$new[[k]])) {
          assign(v, c(perp_history[[v]], k), envir = perp_history)
        }
      }
      if (path_dep) {
        d <- .excitation_design(excitation, history)$path
        n_perps_ev <- c(lengths(hp$old) + lengths(hp$new), numeric(64))
        conn_ev <- c(d$c0 + tabulate(d$upd_j + 1L, nbins = n_history), numeric(64))
      }
    }
    if (!is.null(attachment)) {
      sim_mark_fun(NULL, T0, params, history = history)
    } else {
      net <- create_net_from_events(history)
    }
  }
  if (is.null(net)) net <- net_init_fun()

  # compute initial bound M = lambda(t)
  M <- intensity(cache, cache_f, mu_max)
  stop_reason <- "T_end"

  while (t < T_end) {
    if (!is.finite(M) || M <= 0) break

    if (n_events - n_history >= max_events) {
      stop_reason <- "max_events"
      break
    }
    if (is.finite(max_seconds) && n_props %% 256L == 0L &&
        (proc.time() - start_time)[3] > max_seconds) {
      stop_reason <- "max_seconds"
      break
    }

    # propose next candidate time
    dt <- stats::rexp(1, rate = M)
    t_prop <- t + dt
    if (t_prop > T_end) break

    # move cache to candidate time and compute true intensity there
    cache_prop <- time_cache_step(cache, t_prop, beta = beta)
    cache_f_prop <- time_cache_step(cache_f, t_prop, beta = beta_fast)
    lambda_prop <- intensity(cache_prop, cache_f_prop, mu_at(t_prop))

    accepted <- FALSE
    a <- lambda_prop / M
    n_props <- n_props + 1L
    accept_probs <- grow(accept_probs, n_props)
    accept_probs[n_props] <- a
    if (debug) message("t=", t_prop, " M=", M, " lambda=", lambda_prop, " a=", a)

    if (stats::runif(1) < a) {
      accepted <- TRUE
      # accepted: sample mark given *pre-event* state at t_prop
      if (mark_type %in% c("cs", "cs_bip")) {
        mk <- sim_mark_fun(net = net, t_k = t_prop, params = params, model_cache = cs_model_cache, ...)
      } else {
        mk <- sim_mark_fun(net = net, t_k = t_prop, params = params, delta = delta, ...)
      }

      # CS: carry forward ERNM model cache (to avoid rebuilding each event)
      if (mark_type %in% c("cs", "cs_bip") && !is.null(mk$model_cache)) {
        cs_model_cache <- mk$model_cache
      }

      arrivals_k <- mk$arrivals
      edges_k <- mk$edges

      # Normalise arrivals to a data.frame with at least `id`
      if (is.null(arrivals_k)) {
        arrivals_k <- data.frame(id = character(0), stringsAsFactors = FALSE)
      } else if (is.character(arrivals_k)) {
        arrivals_k <- data.frame(id = as.character(arrivals_k), stringsAsFactors = FALSE)
      }
      if (!is.data.frame(arrivals_k) || !("id" %in% names(arrivals_k))) {
        stop("sim_mark_fun must return `arrivals` as a data.frame with column `id`.")
      }

      # BA-bip: ensure roles exist (event first, then perps) if not provided by sim_mark_fun
      if (mark_type == "ba_bip" && !("role" %in% names(arrivals_k))) {
        if (nrow(arrivals_k) < 1L) stop("BA-bip strict: sim_mark_fun produced no arrivals.")
        arrivals_k$role <- c("event", rep("perp", nrow(arrivals_k) - 1L))
      }

      if (nrow(arrivals_k) > 0) {
        df_nodes <- arrivals_k
        df_nodes$time <- rep(t_prop, nrow(arrivals_k))
        nodes_rows[[length(nodes_rows) + 1L]] <- df_nodes
      }

      if (!is.null(edges_k) && nrow(edges_k) > 0) {
        df_edges <- data.frame(
          i = as.character(edges_k$i),
          j = as.character(edges_k$j),
          time = rep(t_prop, nrow(edges_k)),
          stringsAsFactors = FALSE
        )
        edges_rows[[length(edges_rows) + 1L]] <- df_edges
      }

      # The attachment simulator keeps its own state, so the network is only
      # built at the end, if asked for.
      if (is.null(attachment)) net <- net_step_fun(net, arrivals_k, edges_k, t_prop)

      # Hawkes jump at an accepted event (event-driven ground process)
      if (is.null(excitation)) {
        g_k <- 1
      } else {
        event_id <- as.character(arrivals_k$id[arrivals_k$role == "event"])
        ei <- as.character(edges_k$i)
        ej <- as.character(edges_k$j)
        perps <- unique(c(ej[ei == event_id], ei[ej == event_id]))
        past <- unique(unlist(mget(perps, envir = perp_history, ifnotfound = list(integer(0)))))
        gap <- t_prop - events_t[past]
        past <- past[gap <= excitation$connectedness_window &
                       gap > (excitation$connectedness_min_gap %||% 0)]
        g_k <- g_of(length(perps), length(past))
        k <- n_events + 1L
        for (v in perps) assign(v, c(perp_history[[v]], k), envir = perp_history)
        if (path_dep) {
          # Earlier linked events gain one connectedness now; their change in
          # excitation enters the slow sum discounted to t_prop.
          for (j in past) {
            dg <- g_of(n_perps_ev[j], conn_ev[j] + 1) - g_of(n_perps_ev[j], conn_ev[j])
            cache_prop$S <- cache_prop$S + dg * exp(-beta * (t_prop - events_t[j]))
            conn_ev[j] <- conn_ev[j] + 1
          }
          n_perps_ev <- grow(n_perps_ev, k)
          conn_ev <- grow(conn_ev, k)
          n_perps_ev[k] <- length(perps)
          conn_ev[k] <- length(past)
        }
      }
      cache_prop$S <- cache_prop$S + g_k
      cache_f_prop$S <- cache_f_prop$S + 1
      n_events <- n_events + 1L
      events_t <- grow(events_t, n_events)
      events_t[n_events] <- t_prop

      M_after <- intensity(cache_prop, cache_f_prop, mu_max)
    }

    # advance time regardless; update cache + bound for next iteration
    t <- t_prop
    cache <- cache_prop
    cache_f <- cache_f_prop
    M <- if (accepted) M_after else intensity(cache_prop, cache_f_prop, mu_max)
  }

  events_t <- events_t[seq_len(n_events)]
  # With a burn-in, its events form the history.
  if (is.null(history) && !is.null(history_end)) n_history <- sum(events_t < history_end)
  sim_t <- events_t[seq_len(n_events) > n_history]
  out <- list(
    events = list(n = length(sim_t), t = sim_t),
    accept_probs = accept_probs[seq_len(n_props)],
    history_end = history_end,
    n_history = n_history,
    stop = stop_reason,
    exploded = stop_reason != "T_end",
    t_stop = if (stop_reason == "T_end") T_end else t
  )
  if (out$exploded) {
    warning("Simulation stopped at t = ", signif(out$t_stop, 4), " (", stop_reason,
            "): the model may be explosive.", call. = FALSE)
  }

  if (return_ev || (!is.null(attachment) && return_net)) {
    nodes <- if (length(nodes_rows)) .bind_rows_fast(nodes_rows) else data.frame(id=character(0), time=numeric(0))
    edges <- if (length(edges_rows)) .bind_rows_fast(edges_rows) else data.frame(i=character(0), j=character(0), time=numeric(0))
    ev <- make_events(
      nodes = nodes,
      edges = edges,
      allow_implicit_birth = (mark_type != "ba_bip")
    )
    if (return_ev) out$ev <- ev
  }

  if (!is.null(attachment)) net <- if (return_net) create_net_from_events(ev) else NULL
  out$net <- net

  if (verbose) message("Simulation took ", round((proc.time() - start_time)[3], 2), " seconds")

  out
}

# Internal: rbind for many small data.frames with the same columns, column by
# column (do.call(rbind, ...) is slow for thousands of rows).
.bind_rows_fast <- function(rows) {
  cols <- names(rows[[1L]])
  out <- lapply(cols, function(cl) unlist(lapply(rows, `[[`, cl), use.names = FALSE))
  names(out) <- cols
  as.data.frame(out, stringsAsFactors = FALSE)
}
