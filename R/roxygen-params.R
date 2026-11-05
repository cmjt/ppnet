#' Arguments for adding one event to a network
#'
#' @param net network object, the current network at time t_k
#' @param new_nodes data.frame with at least column "id". Any additional
#'                  columns will be treated as vertex attributes.
#' @param new_edges data.frame with at least columns "i" and "j". Any additional
#'                  columns will be treated as edge attributes.
#' @param t_k numeric scalar, the event timestamp
#' @name new-event-params
#' @keywords internal
NULL

#' Common likelihood arguments
#'
#' @param events An events object created by `make_events()`.
#' @param T_end Optional end time for the likelihood window. Defaults to the
#'   last event time.
#' @param T0 Start time for the likelihood window.
#' @param mark_type One of `"ba"`, `"cs"`, `"ba_bip"`.
#' @param debug Logical; enable debug plotting (interactive; will pause each event).
#' @param ... Additional arguments forwarded to the mark PMF / edge-probability
#'   functions (e.g. `delta`, `model_cache`, etc.).
#'
#' @name hawkesnet-lik-args
NULL
