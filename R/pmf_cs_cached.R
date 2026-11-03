# Compile CS mark log-likelihood into a cheap function(params) evaluator.
#
# Mirrors log_pmf_cs():
#   - M ~ Poisson(node_lambda) (unless t_k > max_node_time)
#   - Candidate edges are deterministic (.cs_candidates_new_old)
#   - Edge probs use ERNM change-stats via mdl$computeChangeStats()
#   - p = plogis(C %*% theta) * exp(-beta_edges * age_old)
#   - Bernoulli-product loglik over all candidates
#
# After compilation, evaluation uses only cached matrices/vectors.

compile_logpmf_cs <- function(
    events,
    formula_rhs,
    truncation = Inf,
    max_node_time = 1e10,
    net_init_fun = net_init,
    net0 = NULL,
    net_step_fun = net_add_event,
    model_cache = NULL,
    eps = 1e-12
) {
  stopifnot(is.list(events))
  stopifnot(is.character(formula_rhs), length(formula_rhs) == 1L, nzchar(formula_rhs))
  stopifnot(is.numeric(max_node_time), length(max_node_time) == 1L, is.finite(max_node_time))
  stopifnot(is.numeric(eps), length(eps) == 1L, is.finite(eps), eps > 0, eps < 0.5)

  arrivals <- nodes_by_event(events)
  edges_by <- edges_by_event(events)
  event_times <- events$times$t

  net <- if (is.null(net0)) net_init_fun() else net0

  cache <- vector("list", length(event_times))
  stat_names_global <- NULL

  for (k in seq_along(event_times)) {
    t_k <- event_times[k]
    new_nodes <- arrivals[[k]]
    new_edges <- edges_by[[k]]

    # arrivals_k (same normalisation as log_pmf_cs)
    if (is.null(new_nodes) || nrow(new_nodes) == 0L) {
      arrivals_k <- character(0)
    } else {
      stopifnot(is.data.frame(new_nodes), "id" %in% names(new_nodes))
      arrivals_k <- unique(as.character(new_nodes$id))
      arrivals_k <- arrivals_k[!is.na(arrivals_k) & nzchar(arrivals_k)]
    }
    M <- length(arrivals_k)

    # old nodes from current net state
    old_nodes <- network::network.vertex.names(net)
    if (is.null(old_nodes)) old_nodes <- character(0)
    old_nodes <- unique(old_nodes[!is.na(old_nodes) & nzchar(old_nodes)])

    # new nodes must not already exist (same check)
    if (length(intersect(arrivals_k, old_nodes)) > 0L) {
      bad <- intersect(arrivals_k, old_nodes)
      stop("compile_logpmf_cs(): new_nodes contains nodes already present in net: ",
           paste(bad, collapse = ", "), call. = FALSE)
    }

    # Poisson term usage flag (same behaviour)
    use_poisson <- isTRUE(t_k <= max_node_time)

    # Handle no old nodes
    if (length(old_nodes) == 0L) {
      # then there are no candidate edges; must observe none
      if (!is.null(new_edges) && nrow(new_edges) > 0L) {
        cache[[k]] <- list(type = "invalid", reason = "edges_with_no_old")
      } else {
        cache[[k]] <- list(type = "no_old", use_poisson = use_poisson, M = M)
      }

      net <- net_step_fun(net, new_nodes, new_edges, t_k)
      next
    }

    # Ensure edges_k is a df(i,j) even if empty/null
    if (is.null(new_edges)) {
      edges_k <- data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
    } else {
      edges_k <- new_edges
      stopifnot(is.data.frame(edges_k), all(c("i", "j") %in% names(edges_k)))
      if (nrow(edges_k) == 0L) {
        edges_k <- data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
      }
    }

    # If no arrivals, then no edges allowed (same behaviour)
    if (M == 0L) {
      if (nrow(edges_k) > 0L) {
        cache[[k]] <- list(type = "invalid", reason = "edges_with_no_arrivals")
      } else {
        cache[[k]] <- list(type = "no_arrivals", use_poisson = use_poisson, M = 0L)
      }

      net <- net_step_fun(net, new_nodes, new_edges, t_k)
      next
    }

    # Augment net with new nodes (no edges) for change stats (same as log_pmf_cs)
    net2 <- .cs_augment_net_with_new_nodes(net, arrivals_k, t_k = t_k)

    # Deterministic candidate set (tails/heads are character ids)
    cand <- .cs_candidates_new_old(arrivals_k, old_nodes, truncation = truncation)
    tails <- cand$tail
    heads <- cand$head

    if (length(tails) == 0L) {
      # no candidates => must observe no edges
      if (nrow(edges_k) > 0L) {
        cache[[k]] <- list(type = "invalid", reason = "edges_with_no_candidates")
      } else {
        cache[[k]] <- list(type = "no_candidates", use_poisson = use_poisson, M = M)
      }

      net <- net_step_fun(net, new_nodes, new_edges, t_k)
      next
    }

    # build idx mapping in net2 (must include arrivals_k)
    idx <- .cs_make_idx(net2)
    tail_idx <- idx[tails]
    head_idx <- idx[heads]
    if (any(is.na(tail_idx)) || any(is.na(head_idx))) {
      stop("compile_logpmf_cs(): missing idx mapping for some candidate nodes.", call. = FALSE)
    }

    # Get ERNM model (cached) and compute change stats once
    res_mdl <- .cs_get_model(net2, formula_rhs = formula_rhs, model_cache = model_cache)
    mdl <- res_mdl$mdl
    model_cache <- res_mdl$model_cache

    C <- mdl$computeChangeStats(as.integer(tail_idx), as.integer(head_idx))
    if (!is.matrix(C)) C <- as.matrix(C)

    cs_stats <- mdl$statistics()
    stat_names <- names(cs_stats)
    if (is.null(stat_names) || any(!nzchar(stat_names))) {
      stop("compile_logpmf_cs(): change-stat matrix has no valid column names.", call. = FALSE)
    }

    # Enforce consistent change-stat columns across events
    if (is.null(stat_names_global)) {
      stat_names_global <- stat_names
    } else if (!identical(stat_names, stat_names_global)) {
      stop("compile_logpmf_cs(): change-stat column names changed across events.", call. = FALSE)
    }

    # Ages of OLD endpoint (heads are old nodes by construction)
    times <- .cs_get_vertex_times(net2)
    names(times) <- network::network.vertex.names(net2)
    age <- t_k - times[heads]
    if (any(!is.finite(age))) stop("compile_logpmf_cs(): non-finite ages.", call. = FALSE)
    if (any(age < 0)) stop("compile_logpmf_cs(): negative ages encountered.", call. = FALSE)

    # Validate observed edges and orient to (new, old), then map to candidate set
    oriented <- validate_edges_cs(edges_k, new_nodes = arrivals_k, old_nodes = old_nodes)
    obs_keys <- paste(oriented$new_end, oriented$old_end, sep = "|")

    cand_keys <- paste(tails, heads, sep = "|")

    # Observed edges must be subset of candidate set (same as log_pmf_cs)
    if (length(setdiff(obs_keys, cand_keys)) > 0L) {
      # With a finite truncation this means the data can't be fit at that
      # truncation, which would otherwise show up as a flat -Inf objective.
      if (is.finite(truncation)) {
        stop("compile_logpmf_cs(): event ", k, " has an edge to an old node outside the ",
             truncation, " most recent; increase `truncation`.", call. = FALSE)
      }
      cache[[k]] <- list(type = "invalid", reason = "obs_not_subset_candidates")
      net <- net_step_fun(net, new_nodes, new_edges, t_k)
      next
    }

    is_obs <- cand_keys %in% obs_keys

    cache[[k]] <- list(
      type = "cs",
      use_poisson = use_poisson,
      M = M,
      C = C,
      age = as.numeric(age),
      is_obs = is_obs
    )

    # advance observed network
    net <- net_step_fun(net, new_nodes, new_edges, t_k)
  }

  force(cache); force(stat_names_global); force(eps)

  function(params) {
    # same parameter requirements as log_pmf_cs
    if (!is.list(params)) return(-Inf)

    beta_edges <- params[["beta_edges"]]
    node_lambda <- params[["node_lambda"]]

    if (!(is.numeric(beta_edges) && length(beta_edges) == 1L && is.finite(beta_edges) && beta_edges >= 0)) return(-Inf)
    if (!(is.numeric(node_lambda) && length(node_lambda) == 1L && is.finite(node_lambda) && node_lambda > 0)) return(-Inf)

    # Extract CS_* thetas once (matching ERNM column names)
    # If we never computed any change-stats during compilation, there are no CS thetas to use.
    if (is.null(stat_names_global) || length(stat_names_global) == 0L) {
      thetas <- numeric(0)
    } else {
      thetas <- .cs_extract_thetas(params, stat_names = stat_names_global)
    }

    ll <- 0.0

    for (k in seq_along(cache)) {
      c <- cache[[k]]

      if (identical(c$type, "invalid")) return(-Inf)

      # Poisson arrivals term
      if (isTRUE(c$use_poisson)) {
        ll <- ll + stats::dpois(c$M, lambda = node_lambda, log = TRUE)
      } else {
        # if past cutoff, must have zero arrivals (as per log_pmf_cs)
        if (c$M != 0L) return(-Inf)
      }

      if (!identical(c$type, "cs")) next

      eta <- as.vector(c$C %*% unname(thetas))
      if (any(!is.finite(eta))) return(-Inf)

      p0 <- stats::plogis(eta)

      p <- p0 * exp(-beta_edges * c$age)
      p <- pmin(pmax(p, eps), 1 - eps)

      ll <- ll + sum(ifelse(c$is_obs, log(p), log1p(-p)))
    }

    as.numeric(ll)
  }
}
