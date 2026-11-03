# event_arrays.R

# Internal: encode an events object as 0-based integer arrays grouped by event,
# the layout used by the C++ walks in src/event_walk.cpp.
.event_arrays <- function(events) {
  stopifnot(inherits(events, "events"))
  t <- events$times$t
  n_ev <- length(t)

  nd <- events$nodes[order(events$nodes$event_id), , drop = FALSE]
  ids <- as.character(nd$id)
  ed <- events$edges[order(events$edges$event_id), , drop = FALSE]
  ei <- match(as.character(ed$i), ids)
  ej <- match(as.character(ed$j), ids)
  if (anyNA(ei) || anyNA(ej)) stop("Edge endpoint without a node birth.", call. = FALSE)

  ptr <- function(event_id) c(0L, cumsum(tabulate(event_id, nbins = n_ev)))

  list(
    t = t,
    ids = ids,
    role = if ("role" %in% names(nd)) as.character(nd$role) else rep(NA_character_, length(ids)),
    node_event = as.integer(nd$event_id),
    node_ptr = as.integer(ptr(nd$event_id)),
    node_birth = t[nd$event_id],
    edge_event = as.integer(ed$event_id),
    edge_ptr = as.integer(ptr(ed$event_id)),
    edge_i = as.integer(ei - 1L),
    edge_j = as.integer(ej - 1L)
  )
}

# Internal: for bipartite histories, the 0-based index of each event's event
# node (-1 if there is not exactly one), and whether each event is valid
# (one event node, joined to all of its new perps).
.bip_focal <- function(a) {
  n_ev <- length(a$t)
  is_event <- !is.na(a$role) & a$role == "event"
  is_perp <- !is.na(a$role) & a$role == "perp"

  n_event_nodes <- tabulate(a$node_event[is_event], nbins = n_ev)
  focal <- rep(-1L, n_ev)
  one <- which(is_event & n_event_nodes[a$node_event] == 1L)
  focal[a$node_event[one]] <- as.integer(one - 1L)

  # A new perp is joined to its event node if an edge of its arrival event
  # links them; an edge at a later event does not count.
  f <- focal[a$edge_event]
  at_i <- a$edge_i == f
  at_j <- a$edge_j == f
  linked <- paste(c(a$edge_j[at_i], a$edge_i[at_j]) + 1L,
                  c(a$edge_event[at_i], a$edge_event[at_j]))
  new_perp <- which(is_perp)
  missing <- new_perp[!(paste(new_perp, a$node_event[new_perp]) %in% linked)]

  valid <- n_event_nodes == 1L
  valid[a$node_event[missing]] <- FALSE

  list(focal = focal, valid = valid,
       n_new = tabulate(a$node_event[is_perp], nbins = n_ev))
}

# Array-based BA-bip mark evaluator; same likelihood as compile_logpmf_ba_bip().
compile_logpmf_ba_bip_array <- function(events, delta = 0.001, eps = 1e-12) {
  stopifnot(is.numeric(delta), length(delta) == 1L, is.finite(delta), delta >= 0)
  stopifnot(is.numeric(eps), length(eps) == 1L, is.finite(eps), eps > 0, eps < 0.5)

  a <- .event_arrays(events)
  if (all(is.na(a$role))) stop("BA-bip compile: expected `role` column in new_nodes.", call. = FALSE)
  f <- .bip_focal(a)
  cand <- !is.na(a$role) & a$role == "perp"
  invalid <- !all(f$valid)

  # With `gradient = TRUE` the value carries attribute "gradient", the
  # derivatives with respect to beta_edges and lambda_new.
  function(params, gradient = FALSE) {
    beta_edges <- params[["beta_edges"]]
    lambda_new <- params[["lambda_new"]]

    if (!(is.numeric(beta_edges) && length(beta_edges) == 1L && is.finite(beta_edges) && beta_edges >= 0)) return(-Inf)
    if (!(is.numeric(lambda_new) && length(lambda_new) == 1L && is.finite(lambda_new) && lambda_new >= 0)) return(-Inf)
    if (invalid) return(-Inf)

    attach <- mark_ll_attach_cpp(a$t, a$node_ptr, cand, a$node_birth, a$edge_ptr,
                                 a$edge_i, a$edge_j, f$focal, beta_edges, delta, eps)
    ll <- sum(stats::dpois(f$n_new, lambda_new, log = TRUE)) + attach[1L]
    if (gradient && is.finite(ll)) {
      attr(ll, "gradient") <- c(beta_edges = attach[2L],
                                lambda_new = sum(f$n_new) / lambda_new - length(f$n_new))
    }
    ll
  }
}

# Array-based BA mark evaluator; same likelihood as compile_logpmf_ba().
compile_logpmf_ba_array <- function(events, delta = 0.001, eps = 1e-12) {
  stopifnot(is.numeric(delta), length(delta) == 1L, is.finite(delta), delta >= 0)
  stopifnot(is.numeric(eps), length(eps) == 1L, is.finite(eps), eps > 0)

  a <- .event_arrays(events)
  n_ev <- length(a$t)
  if (any(tabulate(a$node_event, nbins = n_ev) != 1L)) {
    stop("BA compile: expected exactly 1 arriving node per event.", call. = FALSE)
  }
  if (n_ev > 0L && a$edge_ptr[2L] > 0L) {
    stop("BA event has no eligible existing nodes, but edges were observed.", call. = FALSE)
  }
  cand <- rep(TRUE, length(a$ids))
  focal <- rep(-1L, n_ev)

  function(params, gradient = FALSE) {
    beta_edges <- params[["beta_edges"]]
    if (!(is.numeric(beta_edges) && length(beta_edges) == 1L && is.finite(beta_edges) && beta_edges >= 0)) {
      return(-Inf)
    }
    attach <- mark_ll_attach_cpp(a$t, a$node_ptr, cand, a$node_birth, a$edge_ptr,
                                 a$edge_i, a$edge_j, focal, beta_edges, delta, eps)
    ll <- attach[1L]
    if (gradient && is.finite(ll)) attr(ll, "gradient") <- c(beta_edges = attach[2L])
    ll
  }
}

#' Per-event summaries of a bipartite event-perp history
#'
#' Walks the events in time order and records, for each event, quantities that
#' depend only on the history up to and including that event.
#'
#' @param events An events object with node roles `"event"` and `"perp"`, as
#'   produced for `mark_type = "ba_bip"` or `"cs_bip"`.
#' @param connectedness_window Only earlier events at most this long before
#'   an event count towards its connectedness. `Inf` (the default) counts all.
#'
#' @return A data.frame with one row per event: `event_id`, `t`, `n_perps`
#'   (perps joined to the event), `n_new` and `n_old` (perps first seen at this
#'   event, or earlier), and `connectedness`, the number of distinct earlier
#'   events (within `connectedness_window`) that share at least one perp with
#'   this one.
#' @export
event_summary <- function(events, connectedness_window = Inf) {
  stopifnot(is.numeric(connectedness_window), length(connectedness_window) == 1L,
            connectedness_window >= 0)
  a <- .event_arrays(events)
  if (all(is.na(a$role))) stop("event_summary() needs node roles 'event' and 'perp'.", call. = FALSE)
  f <- .bip_focal(a)
  if (any(f$focal < 0L)) stop("Every event needs exactly one event node.", call. = FALSE)

  s <- event_summary_cpp(a$t, a$node_ptr, a$edge_ptr, a$edge_i, a$edge_j, f$focal,
                         length(a$ids), connectedness_window)
  data.frame(event_id = events$times$event_id, t = a$t, s)
}
