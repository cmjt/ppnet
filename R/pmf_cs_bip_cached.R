# R/pmf_cs_bip_cached.R
#
# Cache-based CS-bip mark log-likelihood (fast evaluator).
#
# Spec (mirrors BA-bip semantics, but CS edge model):
# - Each event must introduce exactly one new node with role == "event"
# - It may introduce K_new >= 0 new nodes with role == "perp"
# - K_new ~ Poisson(lambda_new) (optionally only for t_k <= max_node_time)
# - Deterministic edges (prob 1): new_event -- each new_perp (must be observed)
# - Optional edges (modeled): new_event -- each OLD perp (Bernoulli-product)
#     p = plogis(C %*% theta) * exp(-beta_edges * age_oldperp)
#   where C is computed via ERNM change-stats on the net augmented with new nodes
#
# No edges other than:
#   (new_event, new_perp) and (new_event, old_perp)
# are allowed within a mark. Others are treated as invalid.

compile_logpmf_cs_bip <- function(
    events,
    formula_rhs,
    truncation = Inf,          # kept for API symmetry; currently unused (candidate set is new_event x old_perps)
    max_node_time = Inf,       # if finite: only apply Poisson arrivals for t_k <= max_node_time, else require K_new == 0
    net_init_fun = net_init,
    net0 = NULL,
    net_step_fun = net_add_event,
    model_cache = NULL,
    eps = 1e-12
) {
  stopifnot(is.list(events))
  stopifnot(is.character(formula_rhs), length(formula_rhs) == 1L, nzchar(formula_rhs))
  stopifnot(is.numeric(max_node_time), length(max_node_time) == 1L, is.finite(max_node_time) || is.infinite(max_node_time))
  stopifnot(is.numeric(eps), length(eps) == 1L, is.finite(eps), eps > 0, eps < 0.5)

  arrivals <- nodes_by_event(events)
  edges_by <- edges_by_event(events)
  event_times <- .get_event_times(events)

  net <- if (is.null(net0)) net_init_fun() else net0

  cache <- vector("list", length(event_times))
  stat_names_global <- NULL

  for (k in seq_along(event_times)) {
    t_k <- event_times[k]
    new_nodes <- arrivals[[k]]
    new_edges <- edges_by[[k]]

    if (is.null(new_nodes)) new_nodes <- data.frame(id = character(0), stringsAsFactors = FALSE)
    stopifnot(is.data.frame(new_nodes), "id" %in% names(new_nodes))

    if (!("role" %in% names(new_nodes))) {
      cache[[k]] <- list(type = "invalid", reason = "missing_role_column")
      net <- net_step_fun(net, new_nodes, new_edges, t_k)
      next
    }

    # edges_k normalisation
    if (is.null(new_edges)) {
      edges_k <- data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
    } else {
      edges_k <- new_edges
      stopifnot(is.data.frame(edges_k), all(c("i", "j") %in% names(edges_k)))
      if (nrow(edges_k) == 0L) edges_k <- data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
    }

    new_nodes$role <- as.character(new_nodes$role)
    new_ids <- unique(as.character(new_nodes$id))
    new_ids <- new_ids[!is.na(new_ids) & nzchar(new_ids)]

    is_event <- !is.na(new_nodes$role) & new_nodes$role == "event"
    is_perp  <- !is.na(new_nodes$role) & new_nodes$role == "perp"

    if (sum(is_event) != 1L) {
      cache[[k]] <- list(type = "invalid", reason = "not_exactly_one_event")
      net <- net_step_fun(net, new_nodes, new_edges, t_k)
      next
    }

    new_event <- as.character(new_nodes$id[which(is_event)[1]])
    new_perps <- unique(as.character(new_nodes$id[which(is_perp)]))
    new_perps <- new_perps[!is.na(new_perps) & nzchar(new_perps)]
    K_new <- length(new_perps)

    # Old perps (by role) from current net
    roles <- network::get.vertex.attribute(net, "role")
    ids <- network::network.vertex.names(net)
    if (is.null(ids)) ids <- character(0)
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]

    old_perp_ids <- character(0)
    if (length(ids) > 0L && !is.null(roles)) {
      roles_chr <- as.character(roles)
      old_perp_idx <- which(!is.na(roles_chr) & roles_chr == "perp")
      if (length(old_perp_idx) > 0L) old_perp_ids <- ids[old_perp_idx]
    }

    # Deterministic edges: new_event -- each new_perp must be present
    if (K_new > 0L) {
      edge_keys_all <- .edge_keys_from_df_undirected(edges_k)
      need <- .edge_key_undirected(rep.int(new_event, K_new), new_perps)
      if (!all(need %in% edge_keys_all)) {
        cache[[k]] <- list(type = "invalid", reason = "missing_event_newperp_edges")
        net <- net_step_fun(net, new_nodes, new_edges, t_k)
        next
      }
    }

    # Enforce: no illegal edges in this mark.
    # Allowed: (new_event, new_perp) and (new_event, old_perp)
    if (nrow(edges_k) > 0L) {
      a <- as.character(edges_k$i)
      b <- as.character(edges_k$j)

      ok_newperp <- (a == new_event & b %in% new_perps) | (b == new_event & a %in% new_perps)
      ok_oldperp <- (a == new_event & b %in% old_perp_ids) | (b == new_event & a %in% old_perp_ids)

      if (!all(ok_newperp | ok_oldperp)) {
        cache[[k]] <- list(type = "invalid", reason = "illegal_edge_present")
        net <- net_step_fun(net, new_nodes, new_edges, t_k)
        next
      }
    }

    # Poisson usage (optional cutoff; defaults to always on)
    use_poisson <- isTRUE(t_k <= max_node_time)

    # If no old perps, then there are no modeled edges this event.
    # We still store K_new for the Poisson term.
    if (length(old_perp_ids) == 0L) {
      cache[[k]] <- list(type = "no_old_perps", use_poisson = use_poisson, K_new = K_new)
      net <- net_step_fun(net, new_nodes, new_edges, t_k)
      next
    }

    # Build augmented net with ALL new nodes (event + perps), but no edges yet.
    net2 <- .cs_augment_net_with_new_nodes(net, new_ids, t_k = t_k)

    # Candidate set: new_event -> each old perp
    tails <- rep.int(new_event, length(old_perp_ids))
    heads <- old_perp_ids

    # idx mapping into net2
    idx <- .cs_make_idx(net2)
    tail_idx <- idx[tails]
    head_idx <- idx[heads]
    if (any(is.na(tail_idx)) || any(is.na(head_idx))) {
      stop("compile_logpmf_cs_bip(): missing idx mapping for candidate nodes.", call. = FALSE)
    }

    # ERNM model + change stats
    res_mdl <- .cs_get_model(net2, formula_rhs = formula_rhs, model_cache = model_cache)
    mdl <- res_mdl$mdl
    model_cache <- res_mdl$model_cache

    C <- mdl$computeChangeStats(as.integer(tail_idx), as.integer(head_idx))
    if (!is.matrix(C)) C <- as.matrix(C)

    cs_stats <- mdl$statistics()
    stat_names <- names(cs_stats)
    if (is.null(stat_names) || any(!nzchar(stat_names))) {
      stop("compile_logpmf_cs_bip(): change-stat matrix has no valid column names.", call. = FALSE)
    }

    if (is.null(stat_names_global)) {
      stat_names_global <- stat_names
    } else if (!identical(stat_names, stat_names_global)) {
      stop("compile_logpmf_cs_bip(): change-stat column names changed across events.", call. = FALSE)
    }

    # Ages for old perps (heads)
    times <- .cs_get_vertex_times(net2)
    names(times) <- network::network.vertex.names(net2)
    age <- t_k - times[heads]
    if (any(!is.finite(age))) stop("compile_logpmf_cs_bip(): non-finite ages.", call. = FALSE)
    if (any(age < 0)) stop("compile_logpmf_cs_bip(): negative ages encountered.", call. = FALSE)

    # Observed edges among candidates: (new_event, old_perp)
    edge_keys_all <- .edge_keys_from_df_undirected(edges_k)
    cand_keys <- paste(tails, heads, sep = "|")
    obs_keys <- paste(
      rep.int(new_event, length(heads)),
      heads,
      sep = "|"
    )
    # Mark an edge observed if undirected key exists in edges_k
    is_obs <- .edge_key_undirected(rep.int(new_event, length(heads)), heads) %in% edge_keys_all

    cache[[k]] <- list(
      type = "cs_bip",
      use_poisson = use_poisson,
      K_new = K_new,
      C = C,
      age = as.numeric(age),
      is_obs = is_obs
    )

    net <- net_step_fun(net, new_nodes, new_edges, t_k)
  }

  force(cache); force(stat_names_global); force(eps)

  function(params) {
    if (!is.list(params)) return(-Inf)

    beta_edges <- params[["beta_edges"]]
    lambda_new <- params[["lambda_new"]]

    if (!(is.numeric(beta_edges) && length(beta_edges) == 1L && is.finite(beta_edges) && beta_edges >= 0)) return(-Inf)
    if (!(is.numeric(lambda_new) && length(lambda_new) == 1L && is.finite(lambda_new) && lambda_new >= 0)) return(-Inf)

    # Extract CS_* thetas once
    if (is.null(stat_names_global) || length(stat_names_global) == 0L) {
      thetas <- numeric(0)
    } else {
      thetas <- .cs_extract_thetas(params, stat_names = stat_names_global)
    }

    ll <- 0.0

    for (k in seq_along(cache)) {
      c <- cache[[k]]

      if (identical(c$type, "invalid")) return(-Inf)

      # Poisson arrivals term for new perps
      if (isTRUE(c$use_poisson)) {
        ll <- ll + stats::dpois(c$K_new, lambda = lambda_new, log = TRUE)
      } else {
        if (c$K_new != 0L) return(-Inf)
      }

      if (!identical(c$type, "cs_bip")) next

      eta <- as.vector(c$C %*% unname(thetas))
      if (any(!is.finite(eta))) return(-Inf)

      p0 <- stats::plogis(eta)
      p <- p0 * exp(-beta_edges * c$age)
      p <- .clamp_prob(p, eps = eps)

      ll <- ll + .bern_ll(p, c$is_obs)
    }

    as.numeric(ll)
  }
}
