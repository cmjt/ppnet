#' Add one event to the network (nodes, then edges)
#'
#' @inheritParams new-event-params
#'
#' @return updated network object
#' @export
net_add_event <- function(net, new_nodes, new_edges, t_k) {
  if (!is.null(new_nodes) && nrow(new_nodes) > 0L) {
    net <- net_add_nodes(net, new_nodes, t_k)
  }
  if (!is.null(new_edges) && nrow(new_edges) > 0L) {
    net <- net_add_edges(net, new_edges, t_k)
  }
  net
}

#' Initialize the network-state object
#'
#' @return An empty undirected network.
#' @export
#'
#' @examples
#' net_init()
net_init <- function() {
  net <- network::network.initialize(0, directed = FALSE)

  return(net)
}

