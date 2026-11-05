# R/sim_ba.R

# Internal: generate a new node id from the given state
next_node_id <- function(net) {
  existing <- network::network.vertex.names(net)
  if (is.null(existing) || length(existing) == 0L) return("1")

  suppressWarnings(ints <- as.integer(existing))
  ints <- ints[is.finite(ints)]
  if (length(ints) == 0L) {
    # fallback for non-integer IDs
    return(paste0("n", length(existing) + 1L))
  }
  as.character(max(ints) + 1L)
}

#' Simulate BA mark (arrivals + edges) for one event
#'
#' Paper BA: exactly one arriving node at each accepted event time.
#' Given pre-event state, connect the new node to each old node v
#' independently with probability p_v from edge_probs_ba().
#'
#' @param net Pre-event state (must include named deg, born).
#' @param t_k Numeric scalar event time.
#' @param params List containing at least beta_edges and delta.
#' @param delta Degree offset: weights are proportional to degree + `delta`.
#' @param new_node_id Optional character scalar to force the arriving node id
#'   (useful for tests). If NULL, generated from state.
#'
#' @return list(arrivals = <chr>, edges = data.frame(i,j))
#' @export
sim_mark_ba <- function(net, t_k, params, delta = 0.001, new_node_id = NULL) {
  stopifnot(length(t_k) == 1L, is.numeric(t_k), is.finite(t_k))
  stopifnot(is.list(params))

  stopifnot("beta_edges" %in% names(params))
  beta_edges <- params$beta_edges

  stopifnot(length(beta_edges) == 1L, is.numeric(beta_edges), is.finite(beta_edges), beta_edges >= 0)
  stopifnot(length(delta) == 1L, is.numeric(delta), is.finite(delta), delta >= 0)

  # 1) Choose/generate arriving node id
  if (is.null(new_node_id)) {
    new_node_id <- next_node_id(net)
  }
  stopifnot(length(new_node_id) == 1L, nzchar(new_node_id))
  new_node_id <- as.character(new_node_id)
  new_nodes <- data.frame(id = new_node_id)

  # Compute attachment probabilities over old nodes
  p_old <- edge_probs_ba(net, t_k = t_k, beta_edges = beta_edges, delta = delta)
  old_nodes <- names(p_old)

  # No old nodes => no edges
  if (length(old_nodes) == 0L) {
    return(list(
      arrivals = new_nodes,
      edges = data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
    ))
  }

  add <- stats::runif(length(p_old)) < p_old
  if (!any(add)) {
    return(list(
      arrivals = new_nodes,
      edges = data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
    ))
  }

  edges <- data.frame(
    i = rep.int(new_node_id, sum(add)),
    j = old_nodes[add],
    stringsAsFactors = FALSE
  )

  list(arrivals = new_nodes, edges = edges)
}

