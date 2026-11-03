# R/pmf_cs.R
#
# CS (change-statistics) Bernoulli-product mark kernel (NETWORK-BASED STATE).
#
# Parameter convention (no nested vectors):
#   - beta_edges   : numeric scalar >= 0
#   - node_lambda  : numeric scalar > 0
#   - CS_*         : numeric scalars, one per ERNM change-stat column
#                   e.g. if change stats columns are c("edges","triangles","star.2"),
#                   then params must include c("CS_edges","CS_triangles","CS_star.2")
#
# Network requirements:
#   - net is a network::network
#   - vertex names stored in network::network.vertex.names(net)
#   - vertex attribute "time" stores node entrance times (numeric)

.cs_net_copy <- function(net) {
  if (exists("network.copy", where = asNamespace("network"), inherits = FALSE)) {
    return(network::network.copy(net))
  }
  unserialize(serialize(net, NULL))
}

# deterministic candidate set: each new node considers the `truncation` most
# recently arrived old nodes (no randomness in the candidate set).
# `old_nodes` must be in arrival order, as given by network.vertex.names().
.cs_candidates_new_old <- function(new_nodes, old_nodes, truncation = Inf) {
  stopifnot(is.character(new_nodes), is.character(old_nodes))
  if (length(new_nodes) == 0L || length(old_nodes) == 0L) {
    return(list(tail = character(0), head = character(0)))
  }

  if (is.finite(truncation)) {
    truncation <- as.integer(truncation)
    if (!(length(truncation) == 1L && is.finite(truncation) && truncation >= 1L)) {
      stop("CS truncation must be an integer >= 1 (or Inf).", call. = FALSE)
    }
    old_nodes <- old_nodes[seq.int(to = length(old_nodes), length.out = min(truncation, length(old_nodes)))]
  }

  heads_u <- sort(old_nodes)

  tail_out <- character(0)
  head_out <- character(0)

  for (u in new_nodes) {
    tail_out <- c(tail_out, rep.int(u, length(heads_u)))
    head_out <- c(head_out, heads_u)
  }

  list(tail = tail_out, head = head_out)
}

# Validate observed edges: must be new<->old only, no duplicates, no loops.
# Returns oriented endpoints (new_end, old_end) as character vectors (same length as nrow(edges_k)).
validate_edges_cs <- function(edges_k, new_nodes, old_nodes) {
  stopifnot(is.data.frame(edges_k), all(c("i", "j") %in% names(edges_k)))

  new_nodes <- as.character(new_nodes)
  old_nodes <- as.character(old_nodes)

  if (nrow(edges_k) == 0L) return(list(new_end = character(0), old_end = character(0)))

  i <- as.character(edges_k$i)
  j <- as.character(edges_k$j)

  if (anyNA(i) || anyNA(j) || any(!nzchar(i)) || any(!nzchar(j))) {
    stop("CS edges_k contains NA/empty node IDs.", call. = FALSE)
  }
  if (any(i == j)) stop("CS kernel does not allow self-loops.", call. = FALSE)

  key <- paste(pmin(i, j), pmax(i, j), sep = "|")
  if (anyDuplicated(key)) stop("CS kernel observed duplicate edges in one event.", call. = FALSE)

  is_i_new <- i %in% new_nodes
  is_j_new <- j %in% new_nodes
  is_i_old <- i %in% old_nodes
  is_j_old <- j %in% old_nodes

  ok <- (is_i_new & is_j_old) | (is_j_new & is_i_old)
  if (!all(ok)) {
    bad <- which(!ok)[1]
    stop(
      "CS kernel requires each observed edge to be between a new node and an old node.\n",
      "Bad edge: (", i[[bad]], ", ", j[[bad]], ").",
      call. = FALSE
    )
  }

  new_end <- ifelse(is_i_new, i, j)
  old_end <- ifelse(is_i_new, j, i)

  list(new_end = as.character(new_end), old_end = as.character(old_end))
}

# Old-code-compatible: ensure at least 4 vertices + remove 'na' attribute before C++
.cs_prepare_net_for_cpp <- function(net) {
  stopifnot(inherits(net, "network"))

  net2 <- net

  nv <- network::network.size(net2)
  if (nv < 4L) {
    network::add.vertices(net2, nv = 4L - nv)
  }

  # match old code: delete "na" attribute if present
  if ("na" %in% network::list.vertex.attributes(net2)) {
    network::delete.vertex.attribute(net2, "na")
  }

  net2
}

# Get/create (and cache) the ERNM cpp model, but ALWAYS setNetwork(as.BinaryNet(net)).
# model_cache is a list(rhs = <chr>|NULL, model = <cpp model>|NULL)
.cs_get_model <- function(net, formula_rhs, model_cache = NULL) {
  stopifnot(inherits(net, "network"))
  stopifnot(is.character(formula_rhs), length(formula_rhs) == 1L, nzchar(formula_rhs))

  if (is.null(model_cache)) model_cache <- list(rhs = NULL, model = NULL)

  # If RHS matches, reuse model object. Otherwise rebuild.
  if (!is.null(model_cache$model) && isTRUE(identical(model_cache$rhs, formula_rhs))) {
    mdl <- model_cache$model
  } else {
    # Build ERNM C++ model from formula; note: LHS dummy net will be replaced by setNetwork
    f <- stats::as.formula(paste0("net ~ ", formula_rhs))
    mdl <- ernm::createCppModel(f)
    model_cache$rhs <- formula_rhs
    model_cache$model <- mdl
  }

  # MUST set network each time (net changes as we augment with new nodes)
  net_cpp <- .cs_prepare_net_for_cpp(net)
  mdl$setNetwork(ernm::as.BinaryNet(net_cpp))

  list(mdl = mdl, model_cache = model_cache)
}

# Extract theta vector aligned to ERNM statistic names, using params$CS_<statname>
.cs_extract_thetas <- function(params, stat_names) {
  stopifnot(is.list(params))
  stopifnot(is.character(stat_names), length(stat_names) >= 1L)

  out <- numeric(length(stat_names))
  names(out) <- stat_names

  for (k in seq_along(stat_names)) {
    sn <- stat_names[[k]]
    nm <- paste0("CS_", sn)
    if (!(nm %in% names(params))) {
      stop("Missing CS parameter: params$", nm, call. = FALSE)
    }
    val <- params[[nm]]
    if (!(is.numeric(val) && length(val) == 1L && is.finite(val))) {
      stop("Invalid CS parameter: params$", nm, " must be finite numeric scalar.", call. = FALSE)
    }
    out[[k]] <- val
  }

  out
}

# Internal: ensure a net has a "time" vertex attribute with finite numeric values
.cs_get_vertex_times <- function(net) {
  tt <- network::get.vertex.attribute(net, "time")
  if (is.null(tt)) stop("CS kernel requires vertex attribute 'time' on net.", call. = FALSE)
  if (!is.numeric(tt)) stop("CS kernel: vertex attribute 'time' must be numeric.", call. = FALSE)
  if (any(!is.finite(tt))) stop("CS kernel: vertex attribute 'time' has non-finite values.", call. = FALSE)
  tt
}

# Internal: add new nodes to a *copy* of net (no edges), set names + time
.cs_augment_net_with_new_nodes <- function(net, new_nodes, t_k) {
  stopifnot(inherits(net, "network"))
  stopifnot(is.numeric(t_k), length(t_k) == 1L, is.finite(t_k))

  new_nodes <- unique(as.character(new_nodes))
  new_nodes <- new_nodes[!is.na(new_nodes) & nzchar(new_nodes)]
  if (length(new_nodes) == 0L) return(net)

  # Ensure existing nodes have time
  .cs_get_vertex_times(net)

  net2 <- .cs_net_copy(net)

  nv_old <- network::network.size(net2)
  network::add.vertices(net2, nv = length(new_nodes))
  nv_new <- network::network.size(net2)
  new_vids <- seq.int(nv_old + 1L, nv_new)

  # Set vertex names (network.vertex.names uses "vertex.names")
  network::set.vertex.attribute(net2, "vertex.names", value = new_nodes, v = new_vids)

  # Set time for new nodes
  network::set.vertex.attribute(net2, "time", value = rep(t_k, length(new_nodes)), v = new_vids)

  net2
}

# Internal: mapping node-id -> vertex index for a net
.cs_make_idx <- function(net) {
  nms <- network::network.vertex.names(net)
  if (is.null(nms)) nms <- character(0)
  if (anyNA(nms) || any(!nzchar(nms))) stop("CS kernel requires non-empty vertex names.", call. = FALSE)
  setNames(seq_along(nms), nms)
}

# Internal: compute CS attachment probabilities for candidate new-old edges.
# Returns named numeric vector with names "new|old".
edge_probs_cs <- function(net,
                          new_nodes,
                          old_nodes,
                          t_k,
                          params,
                          formula_rhs,
                          truncation = Inf,
                          model_cache = NULL,
                          beta_edges = NULL,
                          eps = 1e-12) {

  stopifnot(inherits(net, "network"))
  stopifnot(is.numeric(t_k), length(t_k) == 1L, is.finite(t_k))
  stopifnot(is.list(params))
  stopifnot(is.character(formula_rhs), length(formula_rhs) == 1L, nzchar(formula_rhs))
  stopifnot(is.numeric(eps), length(eps) == 1L, is.finite(eps), eps > 0, eps < 0.5)

  # allow passing beta_edges explicitly (useful for tests/back-compat), otherwise take from params
  if (is.null(beta_edges)) {
    if (!("beta_edges" %in% names(params))) stop("edge_probs_cs(): missing params$beta_edges", call. = FALSE)
    beta_edges <- params$beta_edges
  }
  stopifnot(is.numeric(beta_edges), length(beta_edges) == 1L, is.finite(beta_edges), beta_edges >= 0)

  new_nodes <- unique(as.character(new_nodes))
  old_nodes <- unique(as.character(old_nodes))

  if (length(new_nodes) == 0L || length(old_nodes) == 0L) {
    return(list(p = structure(numeric(0), names = character(0)), model_cache = model_cache))
  }

  cand <- .cs_candidates_new_old(new_nodes, old_nodes, truncation = truncation)
  tails <- cand$tail
  heads <- cand$head
  if (length(tails) == 0L) return(list(p = structure(numeric(0), names = character(0)), model_cache = model_cache))

  idx <- .cs_make_idx(net)

  tail_idx <- idx[tails]
  head_idx <- idx[heads]
  if (any(is.na(tail_idx)) || any(is.na(head_idx))) {
    stop("edge_probs_cs(): missing idx mapping for some candidate nodes. Did you augment net with arrivals first?",
         call. = FALSE)
  }

  res_mdl <- .cs_get_model(net, formula_rhs = formula_rhs, model_cache = model_cache)
  mdl <- res_mdl$mdl
  model_cache <- res_mdl$model_cache

  C <- mdl$computeChangeStats(as.integer(tail_idx), as.integer(head_idx))
  if (!is.matrix(C)) C <- as.matrix(C)

  cs_stats <- mdl$statistics()
  stat_names <- names(cs_stats)
  if (is.null(stat_names) || any(!nzchar(stat_names))) {
    stop("edge_probs_cs(): change-stat matrix has no valid column names; cannot map to CS_* params.",
         call. = FALSE)
  }

  thetas <- .cs_extract_thetas(params, stat_names = stat_names)

  eta <- as.vector(C %*% unname(thetas))
  if (any(!is.finite(eta))) stop("edge_probs_cs(): non-finite linear predictors.", call. = FALSE)

  p0 <- stats::plogis(eta)

  # decay uses OLD node entrance times (node arrival times), stored as vertex attr "time"
  times <- .cs_get_vertex_times(net)
  names(times) <- network::network.vertex.names(net)
  age <- t_k - times[heads]
  if (any(!is.finite(age))) stop("edge_probs_cs(): non-finite ages from time attribute.", call. = FALSE)
  if (any(age < 0)) stop("edge_probs_cs(): negative ages encountered (time after t_k).", call. = FALSE)

  p <- p0 * exp(-beta_edges * age)
  p <- pmin(pmax(p, eps), 1 - eps)

  names(p) <- paste(tails, heads, sep = "|")
  list(p = p, model_cache = model_cache)
}
