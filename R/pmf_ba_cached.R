# Internal: compute BA attachment probabilities over existing (old) nodes.
# Returns a named numeric vector p_old with names = old node IDs.
edge_probs_ba <- function(net, t_k, beta_edges = 0.5, delta = 0.001) {
  if (net %n% 'n' == 0) {
    return(numeric(0))
  }

  # Birth-time age of old nodes
  age <- t_k - network::get.vertex.attribute(net, "time")
  # Degrees of old nodes
  deg <- sna::degree(net)

  # Attachment attractiveness score (must be >= 0)
  #   delta=0 results in pure degree-based preference
  w <- (deg + delta) * exp(-beta_edges * age)

  # Normalizing term
  W <- sum(w)
  if (!is.finite(W) || W <= 0) {
    stop("Invalid BA weights: sum(weights) must be > 0.", call. = FALSE)
  }

  # Normalise probabilities
  #   IMPORTANT: if only 1 old node, this implies probability of attachment = 1)
  p_old <- w / W
  names(p_old) <- network::network.vertex.names(net)

  # Numerical safety
  p_old <- pmin(pmax(p_old, 0), 1)

  p_old
}

# Internal: validate that edges_k are BA-valid and return oriented endpoints.
#
# Paper BA (strict):
# - Exactly one arriving node `arrival_k` (character scalar).
# - Every edge must be between that new node and one old node in `old_nodes`.
# - No endpoints outside {arrival_k} U old_nodes.
# - No new-new or old-old edges.
# - No duplicate (new, old) pairs within the event.
#
# Returns:
#   list(new_end = <chr>, old_end = <chr>)
# one element per row of edges_k (same order).
validate_edges_ba <- function(edges_k, arrival_k, old_nodes) {
  stopifnot(is.data.frame(edges_k), all(c("i", "j") %in% names(edges_k)))

  if (nrow(edges_k) == 0L) {
    return(list(new_end = character(0), old_end = character(0)))
  }

  # Exactly one new node for BA
  stopifnot(length(arrival_k) == 1L, nzchar(arrival_k))
  new_node <- as.character(arrival_k[[1]])

  old_nodes <- as.character(old_nodes)
  # if old_nodes is empty, any edge is impossible
  if (length(old_nodes) == 0L) {
    stop("BA event has edges but there are no eligible old nodes.", call. = FALSE)
  }

  i <- as.character(edges_k$i)
  j <- as.character(edges_k$j)

  i_is_new <- (i == new_node)
  j_is_new <- (j == new_node)
  i_is_old <- (i %in% old_nodes)
  j_is_old <- (j %in% old_nodes)

  # Any endpoint must be either the new node or an eligible old node
  endpoint_ok <- (i_is_new | i_is_old) & (j_is_new | j_is_old)
  if (any(!endpoint_ok)) {
    stop("BA event has edge endpoint not in {arrival_k} or old_nodes.", call. = FALSE)
  }

  # Edges must be new-old: exactly one endpoint is the new node, the other is old
  is_new_old <- (i_is_new & j_is_old) | (j_is_new & i_is_old)
  if (any(!is_new_old)) {
    stop("BA event has non-BA edges (must be between the new node and an old node).", call. = FALSE)
  }

  # Orient edges: new_end is always the new node; old_end is the old node
  new_end <- ifelse(i_is_new, i, j)
  old_end <- ifelse(i_is_old, i, j)

  # Sanity: all new_end should equal new_node (given strict checks)
  # (keep as stopifnot for internal consistency)
  stopifnot(all(new_end == new_node))

  # No duplicates: Bernoulli model can't contain repeated (new, old) within one event
  if (anyDuplicated(old_end) > 0L) {
    stop("Duplicate BA edges detected within event (same old node repeated).", call. = FALSE)
  }

  list(new_end = new_end, old_end = old_end)
}

# Compile BA mark log-likelihood into a cheap function(params) evaluator.
#
# After compilation, evaluation is all vectorised arithmetic (no network calls).
compile_logpmf_ba <- function(
    events,
    T0 = 0,
    net_init_fun = net_init,
    net0 = NULL,
    net_step_fun = net_add_event,
    delta = 0.001,
    eps = 1e-12
) {
  stopifnot(is.numeric(delta), length(delta) == 1L, is.finite(delta), delta >= 0)
  stopifnot(is.numeric(eps), length(eps) == 1L, is.finite(eps), eps > 0)

  arrivals <- nodes_by_event(events)
  edges_by <- edges_by_event(events)
  event_times <- events$times$t

  # One pass through the data (network manipulation happens ONCE here)
  net <- if (is.null(net0)) net_init_fun() else net0

  cache <- vector("list", length(event_times))

  for (idx in seq_along(event_times)) {
    t_k <- event_times[idx]
    new_nodes <- arrivals[[idx]]
    new_edges <- edges_by[[idx]]

    # Your BA kernel requires exactly 1 arriving node
    if (!is.null(new_nodes) && nrow(new_nodes) != 1L) {
      stop("BA compile: expected exactly 1 arriving node per event.", call. = FALSE)
    }

    old_ids <- network::network.vertex.names(net)
    if (is.null(old_ids)) old_ids <- character(0)

    if (length(old_ids) == 0L) {
      # Must match your current behaviour:
      # - if no existing nodes and no edges: logp = 0
      # - if edges exist: invalid
      if (!is.null(new_edges) && nrow(new_edges) > 0L) {
        stop("BA event has no eligible existing nodes, but edges were observed.", call. = FALSE)
      }
      cache[[idx]] <- list(type = "empty")  # contributes 0
    } else {
      # Precompute quantities that DO NOT depend on beta_edges:
      # degrees + ages at t_k for all old nodes.
      age <- t_k - network::get.vertex.attribute(net, "time")
      deg <- sna::degree(net)

      # Determine which old nodes were attached-to in this event.
      # Your current BA PMF uses:
      # present <- intersect(existing_nodes, c(new_edges$i, new_edges$j))
      present_ids <- character(0)
      if (!is.null(new_edges) && nrow(new_edges) > 0L) {
        present_ids <- intersect(old_ids, c(new_edges$i, new_edges$j))
      }
      present <- old_ids %in% present_ids

      cache[[idx]] <- list(
        type = "ba",
        ids = old_ids,
        deg = as.numeric(deg),
        age = as.numeric(age),
        present = present
      )
    }

    # Advance net along observed history
    net <- net_step_fun(net, new_nodes, new_edges, t_k)
  }

  # Return a cheap evaluator: function(params) -> total mark loglik
  force(cache); force(delta); force(eps)
  function(params) {
    beta_edges <- params[["beta_edges"]]
    if (!(is.numeric(beta_edges) && length(beta_edges) == 1L && is.finite(beta_edges) && beta_edges >= 0)) {
      return(-Inf)
    }

    ll <- 0.0
    for (idx in seq_along(cache)) {
      c <- cache[[idx]]
      if (identical(c$type, "empty")) next

      # weights w_i = (deg_i + delta) * exp(-beta_edges * age_i)
      w <- (c$deg + delta) * exp(-beta_edges * c$age)
      W <- sum(w)

      if (!is.finite(W) || W <= 0) return(-Inf)

      p <- w / W

      # Numerical safety (avoid log(0) / log(negative))
      p <- pmin(pmax(p, eps), 1 - eps)

      # logp = sum_{present} log(p) + sum_{absent} log(1-p)
      ll <- ll + sum(ifelse(c$present, log(p), log1p(-p)))
    }

    as.numeric(ll)
  }
}
