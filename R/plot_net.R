# Drawing the event-perp network.

.role_cols <- c(event = "#D55E00", perp = "#0072B2")

# Internal: weak components of a network, ranked by size (largest first; ties
# by the time of their first node). Returns membership as ranks.
.ranked_components <- function(net) {
  memb <- sna::component.dist(net, connected = "weak")$membership
  born <- network::get.vertex.attribute(net, "time")
  if (is.null(born) || all(is.na(born))) born <- seq_along(memb)
  size <- tabulate(memb)
  first <- vapply(seq_along(size), function(k) min(born[memb == k]), numeric(1))
  rank <- order(-size, first)
  match(memb, rank)
}

#' Components of an event-perp network
#'
#' Lists the connected components (groups of events and perps linked through
#' shared perps) of the network built from an events object, largest first.
#'
#' @param x An events object, or a network from [create_net_from_events()].
#' @return A data frame with one row per component: `component` (rank, 1 =
#'   largest), `n_events`, `n_perps`, `first` and `last` (the times of its
#'   first and last events), and attribute `"membership"` (each node's
#'   component).
#' @export
#'
#' @examples
#' set.seed(1)
#' p <- list(mu = 2, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
#' ev <- sim_hawkesNet(p, T_end = 20, mark_type = "ba_bip")$ev
#' head(net_components(ev))
net_components <- function(x) {
  net <- if (inherits(x, "events")) create_net_from_events(x) else x
  if (!inherits(net, "network")) stop("`x` must be an events object or a network.", call. = FALSE)
  memb <- .ranked_components(net)
  role <- network::get.vertex.attribute(net, "role")
  if (is.null(role) || all(is.na(role))) role <- rep("event", length(memb))
  born <- network::get.vertex.attribute(net, "time")
  k <- seq_len(max(memb))
  is_ev <- role == "event"
  out <- data.frame(
    component = k,
    n_events = tabulate(memb[is_ev], nbins = length(k)),
    n_perps = tabulate(memb[!is_ev], nbins = length(k)),
    first = vapply(k, function(j) min(born[memb == j & is_ev], Inf), numeric(1)),
    last = vapply(k, function(j) max(born[memb == j & is_ev], -Inf), numeric(1))
  )
  attr(out, "membership") <- memb
  out
}

#' Plot an event-perp network
#'
#' Draws the network built by [create_net_from_events()] (or returned by
#' [sim_hawkesNet()]) with events and perps told apart: events as orange
#' squares, perps as blue circles, each sized by its degree (perps per event,
#' or events per perp). Other arguments go to [network::plot.network()], so
#' any of its settings (layout `mode`, labels, colours) can be changed.
#'
#' Whole networks of thousands of nodes are slow to lay out and hard to read;
#' `component` draws one or more connected components instead, ranked by
#' size (see [net_components()]).
#'
#' @param x A network of class `"ppnet_net"`.
#' @param component `NULL` (default) for the whole network, `"largest"`, or
#'   component ranks (1 = largest), e.g. `1:4`, each drawn in its own panel.
#' @param size_by_degree Logical; scale node size by degree.
#' @param legend Logical; add a legend.
#' @param ... Passed to [network::plot.network()].
#'
#' @return Invisibly, the node coordinates (a list of them for several
#'   components).
#' @export
#'
#' @examples
#' set.seed(1)
#' p <- list(mu = 2, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
#' ev <- sim_hawkesNet(p, T_end = 20, mark_type = "ba_bip")$ev
#' net <- create_net_from_events(ev)
#' plot(net)
#' plot(net, component = 1:4)
plot.ppnet_net <- function(x, component = NULL, size_by_degree = TRUE, legend = TRUE, ...) {
  if (is.null(component)) {
    if (network::network.size(x) > 2000) {
      message("Laying out ", network::network.size(x), " nodes; `component = 1:4` draws the ",
              "largest components instead.")
    }
    return(invisible(.plot_bip_net(x, size_by_degree, legend, main = NULL, ...)))
  }
  memb <- .ranked_components(x)
  if (identical(component, "largest")) component <- 1L
  component <- as.integer(component)
  bad <- component[component < 1L | component > max(memb)]
  if (length(bad)) stop("There are ", max(memb), " components; no component ", bad[1], ".", call. = FALSE)

  if (length(component) > 1L) {
    old <- graphics::par(mfrow = grDevices::n2mfrow(length(component)))
    on.exit(graphics::par(old), add = TRUE)
  }
  coords <- lapply(component, function(k) {
    sub <- network::get.inducedSubgraph(x, which(memb == k))
    role <- network::get.vertex.attribute(sub, "role")
    main <- paste0("Component ", k, ": ", sum(role == "event"), " events, ",
                   sum(role == "perp"), " perps")
    .plot_bip_net(sub, size_by_degree, legend && k == component[1], main = main,
                  scale = sqrt(length(component)), ...)
  })
  invisible(if (length(coords) == 1L) coords[[1L]] else coords)
}

# Internal: draw one network with event/perp colours, shapes and sizes.
.plot_bip_net <- function(net, size_by_degree, legend, main, scale = 1, ...) {
  role <- network::get.vertex.attribute(net, "role")
  deg <- sna::degree(net, gmode = "graph")
  # Node sizes are relative to the panel, so smaller panels get larger nodes.
  cex <- scale * (if (size_by_degree) 0.7 + 0.4 * sqrt(deg) else rep(1, length(deg)))
  args <- utils::modifyList(list(
    x = net, vertex.col = unname(.role_cols[role]), vertex.sides = ifelse(role == "event", 4L, 50L),
    vertex.rot = ifelse(role == "event", 45, 0), vertex.cex = cex, vertex.border = NA,
    edge.col = "grey55", main = main
  ), list(...))
  coord <- do.call(network::plot.network, args)
  if (legend) {
    graphics::legend("bottomleft", legend = c("event", "perp"), pch = c(15, 16),
                     col = .role_cols, bty = "n", cex = 0.9)
  }
  invisible(coord)
}
