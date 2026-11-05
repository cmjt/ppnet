#' Build the network from an events object
#'
#' Replays the events in order, adding each event's nodes and edges.
#'
#' @param events An events object created by [make_events()].
#' @param net_init_fun Function returning the starting network (used when
#'   `net0` is `NULL`).
#' @param net0 Optional starting network.
#' @param net_step_fun Function adding one event's nodes and edges.
#' @param debug Logical; plot the network after each event and wait for Enter.
#'
#' @return The final `network` object (or `0` if there are no events).
#' @export
create_net_from_events <- function(events,
                                   net_init_fun = net_init,
                                   net0 = NULL,
                                   net_step_fun = net_add_event,
                                   debug = FALSE) {
  #
  arrivals <- nodes_by_event(events)
  edges <- edges_by_event(events)

  event_times <- events$times$t
  if (length(event_times) == 0L) return(0)

  net <- if (is.null(net0)) net_init_fun() else net0

  for (idx in seq_along(event_times)) {
    t_k <- event_times[idx]
    new_nodes <- arrivals[[idx]]
    new_edges <- edges[[idx]]

    net <- net_step_fun(net, new_nodes, new_edges, t_k)

    if (debug == TRUE) {
      plot(net)
      readline("Press Enter / Return to continue to next event...")
    }
  }

  return(net)
}

#' Convert a network with time attributes to an events object
#'
#' Preserves network attributes. Defaults to expecting "time" as edge / vertex
#' attribute for event / birth time.
#'
#' @param net A `network` object with node birth times and edge times.
#' @param node_time_attr Name of the vertex attribute holding birth times.
#' @param edge_time_attr Name of the edge attribute holding edge times.
#' @param ... Passed to [make_events()].
#'
#' @return An events object, as from [make_events()].
#' @export
create_events_from_net <- function(net,
                                   node_time_attr = "time",
                                   edge_time_attr = "time",
                                   ...) {
  if (!inherits(net, "network")) stop("`net` must be a `network` object.")

  # Create nodes data frame
  node_time_attr <- match.arg(node_time_attr)
  ids <- as.character(network::network.vertex.names(net))
  n <- length(ids)

  node_time <- network::get.vertex.attribute(net, node_time_attr)

  if (is.null(node_time) && n > 0L) {
    stop("`net` must have a vertex attribute 'time' or 'born' giving node birth times.")
  }

  nodes <- if (n == 0L) {
    data.frame(id = character(0), time = numeric(0))
  } else {
    data.frame(id = ids, time = as.numeric(node_time), stringsAsFactors = FALSE)
  }

  v_attrs <- setdiff(network::list.vertex.attributes(net), node_time_attr)
  for (a in v_attrs) nodes[[a]] <- network::get.vertex.attribute(net, a)

  # Create edges data frame
  edge_time_attr <- match.arg(edge_time_attr)

  edf <- network::as.data.frame.network(net, "edges")
  # edf typically has columns: tail, head, and any edge attrs
  if (nrow(edf) == 0L) {
    edges <- data.frame(i = character(0), j = character(0), time = numeric(0))
  } else {
    edges <- data.frame(
      i = as.character(edf$`.tail`),
      j = as.character(edf$`.head`),
      time = as.numeric(edf[[edge_time_attr]]),
      stringsAsFactors = FALSE
    )

    # keep extra edge attributes too (optional)
    extra <- setdiff(names(edf), c("tail", "head", edge_time_attr))
    for (a in extra) edges[[a]] <- edf[[a]]
  }

  make_events(edges = edges, nodes = nodes, ...)
}


# tiny helper (so we don't need rlang)
`%||%` <- function(x, y) if (is.null(x)) y else x

# Internal: transform parameters to working scale
to_working <- function(p_list, transform) {
  stopifnot(is.list(p_list))

  x <- numeric(length(p_list))
  names(x) <- names(p_list)

  for (nm in names(p_list)) {
    val <- p_list[[nm]]
    tf  <- transform[[nm]] %||% "log"

    if (tf == "log") {
      if (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val <= 0) {
        stop("Parameter '", nm, "' must be > 0 for log-transform.")
      }
      x[[nm]] <- log(val)
    } else if (tf == "none") {
      if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
        stop("Parameter '", nm, "' must be a finite scalar.")
      }
      x[[nm]] <- val
    } else {
      stop("Unknown transform '", tf, "' for parameter '", nm, "'.")
    }
  }

  x
}

# Internal: transform parameters back to original scale
from_working <- function(x, transform) {
  stopifnot(is.numeric(x))

  p <- as.list(x)
  names(p) <- names(x)

  for (nm in names(x)) {
    tf <- transform[[nm]] %||% "log"
    if (tf == "log") {
      p[[nm]] <- exp(x[[nm]])
    } else if (tf == "none") {
      p[[nm]] <- x[[nm]]
    } else {
      stop("Unknown transform '", tf, "' for parameter '", nm, "'.")
    }
  }

  p
}


