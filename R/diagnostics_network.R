# diagnostics_network.R

# Internal helper: degree vector for a `network` object
.degree_vec <- function(net) {
  # Prefer sna if available (fast)
  if (requireNamespace("sna", quietly = TRUE)) {
    return(as.numeric(sna::degree(net, gmode = "graph")))
  }

  # Fallback: adjacency matrix (OK for moderate size)
  A <- try(network::as.matrix.network(net, matrix.type = "adjacency"), silent = TRUE)
  if (inherits(A, "try-error")) {
    stop("degree fallback failed: install 'sna' or ensure network can be coerced to adjacency.",
         call. = FALSE)
  }
  as.numeric(rowSums(A != 0) + colSums(A != 0)) / 2  # undirected-ish fallback
}

# Internal helper: CCDF of a nonnegative vector
.ccdf <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(data.frame(k = numeric(0), ccdf = numeric(0)))
  tab <- table(x)
  k <- as.numeric(names(tab))
  p <- as.numeric(tab) / sum(tab)
  ord <- order(k)
  k <- k[ord]; p <- p[ord]
  ccdf <- rev(cumsum(rev(p)))
  data.frame(k = k, ccdf = ccdf)
}

#' Network growth diagnostics over events
#'
#' Replays the event sequence and tracks network size and simple network
#' summaries over time. This is kernel-agnostic as long as `net_step_fun`
#' updates the network using `(net, new_nodes, new_edges, t_k)`.
#'
#' Bipartite event-perp histories (node roles `"event"` and `"perp"`) are
#' summarised separately for events and perps, from the event walk used by
#' [event_summary()] (fast, no network objects): group sizes (perps per
#' event), returning and new perps per event, and perps' attendance (events
#' per perp), since event-node and perp degrees measure different things.
#'
#' @param events An events object created by `make_events()`.
#' @param net_init_fun Function that returns an initial `network` object.
#' @param net0 Optional initial network object (overrides `net_init_fun()`).
#' @param net_step_fun Function that updates the network after an event. Default
#'   is `net_add_event`.
#'
#' @return A list with:
#' \describe{
#'   \item{series}{Data frame with one row per event: time, nodes_n, edges_n,
#'   mean_degree, max_degree.}
#'   \item{deg_final}{Degree vector at the final network state.}
#'   \item{ccdf_final}{CCDF data frame for final degrees.}
#' }
#' For a bipartite history instead: `bipartite = TRUE`, `series` (one row per
#' event, after it: `t`, `events_n`, `perps_n`, `edges_n`, `group_size`,
#' `returning`, `new`), `group_size` and `attendance` (per event and per
#' perp) and `ccdf_attendance`.
#' @export
diag_network <- function(events,
                         net_init_fun = net_init,
                         net0 = NULL,
                         net_step_fun = net_add_event) {
  # validate_events(events)

  if (is.null(events$times) || !"t" %in% names(events$times)) {
    stop("diag_network(): `events$times$t` is required.", call. = FALSE)
  }

  t <- as.numeric(events$times$t)
  n_events <- length(t)

  if (.is_bipartite(events) && is.null(net0) && identical(net_init_fun, net_init) &&
      identical(net_step_fun, net_add_event)) {
    s <- event_summary(events)
    hp <- .history_perps(events)
    attendance <- as.vector(table(unlist(Map(c, hp$old, hp$new), use.names = FALSE)))
    series <- data.frame(t = t, events_n = seq_len(n_events), perps_n = cumsum(s$n_new),
                         edges_n = cumsum(s$n_perps), group_size = s$n_perps,
                         returning = s$n_old, new = s$n_new)
    return(list(bipartite = TRUE, series = series, group_size = s$n_perps,
                attendance = attendance, ccdf_attendance = .ccdf(attendance)))
  }

  arrivals <- nodes_by_event(events)
  edges <- edges_by_event(events)
  if (length(arrivals) != n_events || length(edges) != n_events) {
    stop("diag_network(): nodes_by_event/edges_by_event must return one element per event.",
         call. = FALSE)
  }

  net <- if (is.null(net0)) net_init_fun() else net0
  if (!inherits(net, "network")) stop("diag_network(): net must be a 'network' object.", call. = FALSE)

  nodes_n <- integer(n_events)
  edges_n <- integer(n_events)
  mean_deg <- numeric(n_events)
  max_deg <- numeric(n_events)

  for (k in seq_len(n_events)) {
    # Pre-update summaries (about current state)
    nodes_n[k] <- network::network.size(net)
    edges_n[k] <- network::network.edgecount(net)

    dv <- .degree_vec(net)
    mean_deg[k] <- if (length(dv)) mean(dv) else 0
    max_deg[k] <- if (length(dv)) max(dv) else 0

    # Update network with event k
    net <- net_step_fun(net, arrivals[[k]], edges[[k]], t[k])
  }

  # Final degrees (after last update)
  dv_final <- .degree_vec(net)
  ccdf_final <- .ccdf(dv_final)

  series <- data.frame(
    t = t,
    nodes_n = nodes_n,
    edges_n = edges_n,
    mean_degree = mean_deg,
    max_degree = max_deg
  )

  list(series = series, deg_final = dv_final, ccdf_final = ccdf_final)
}


#' Plot network diagnostics from `diag_network()`
#'
#' @param diag Output of `diag_network()`.
#' @param which Character vector of plots to show. Any of:
#'   `"size"`, `"edges_vs_nodes"`, `"degree_hist"`, `"degree_ccdf"`. For a
#'   bipartite history: `"size"` (events and perps over time),
#'   `"group_size"` (perps per event), `"returning"` (returning perps per event
#'   over time), `"attendance"` (events per perp) and `"attendance_ccdf"`;
#'   `"degree_hist"` and `"degree_ccdf"` then mean the attendance plots.
#' @param ask Logical; if `TRUE`, prompt between plots.
#' @param loglog Logical; for `"degree_ccdf"`, plot on log-log axes.
#' @export
plot_network_diag <- function(diag,
                              which = c("size", "edges_vs_nodes", "degree_hist", "degree_ccdf"),
                              ask = FALSE,
                              loglog = TRUE) {
  if (is.null(diag) || is.null(diag$series)) stop("plot_network_diag(): invalid `diag`.", call. = FALSE)
  if (isTRUE(diag$bipartite)) {
    if (missing(which)) which <- c("size", "group_size", "returning", "attendance_ccdf")
    return(.plot_bip_diag(diag, which, ask, loglog))
  }
  bip_only <- setdiff(which, c("size", "edges_vs_nodes", "degree_hist", "degree_ccdf"))
  if (length(bip_only)) {
    stop("Only for bipartite event-perp data: ", paste(bip_only, collapse = ", "), call. = FALSE)
  }
  ser <- diag$series
  which <- unique(which)

  old_ask <- par(ask = ask)
  on.exit(par(old_ask), add = TRUE)

  if ("size" %in% which) {
    plot(ser$t, ser$nodes_n, type = "l", xlab = "time", ylab = "count",
         main = "Network size over time")
    lines(ser$t, ser$edges_n, lty = 2)
    legend("topleft", legend = c("# nodes", "# edges"), lty = c(1, 2), bty = "n")
  }

  if ("edges_vs_nodes" %in% which) {
    plot(ser$nodes_n, ser$edges_n, xlab = "# nodes", ylab = "# edges",
         main = "Edges vs nodes", pch = 16)
  }

  if ("degree_hist" %in% which) {
    dv <- diag$deg_final
    if (!length(dv)) {
      plot.new(); text(0.5, 0.5, "No degrees to plot.")
    } else {
      hist(dv, breaks = "FD", xlab = "degree", main = "Final degree distribution",
           col = "grey90", border = "grey40")
    }
  }

  if ("degree_ccdf" %in% which) {
    cc <- diag$ccdf_final
    if (!nrow(cc)) {
      plot.new(); text(0.5, 0.5, "No CCDF to plot.")
    } else {
      if (isTRUE(loglog)) {
        # Avoid log(0)
        cc2 <- cc[cc$k > 0 & cc$ccdf > 0, , drop = FALSE]
        plot(cc2$k, cc2$ccdf, log = "xy", type = "p", pch = 16,
             xlab = "k", ylab = "P(Degree >= k)",
             main = "Degree CCDF (log-log)")
      } else {
        plot(cc$k, cc$ccdf, type = "p", pch = 16,
             xlab = "k", ylab = "P(Degree >= k)",
             main = "Degree CCDF")
      }
    }
  }

  invisible(NULL)
}

# Internal: plots of a bipartite diag_network().
.plot_bip_diag <- function(diag, which, ask, loglog) {
  which <- unique(sub("^degree_hist$", "attendance", sub("^degree_ccdf$", "attendance_ccdf", which)))
  bad <- setdiff(which, c("size", "group_size", "returning", "attendance", "attendance_ccdf"))
  if (length(bad)) stop("Not a bipartite network plot: ", paste(bad, collapse = ", "), call. = FALSE)
  ser <- diag$series
  old_ask <- par(ask = ask)
  on.exit(par(old_ask), add = TRUE)
  bars <- function(x, main, xlab) {
    tab <- table(factor(x, levels = seq.int(min(x), max(x))))
    graphics::barplot(tab / sum(tab), main = main, xlab = xlab, ylab = "proportion",
                      col = "grey80", border = "grey40")
  }
  if ("size" %in% which) {
    plot(ser$t, ser$perps_n, type = "l", xlab = "time", ylab = "count",
         ylim = c(0, max(ser$perps_n, ser$events_n)), main = "Events and perps over time")
    lines(ser$t, ser$events_n, lty = 2)
    legend("topleft", legend = c("# perps", "# events"), lty = c(1, 2), bty = "n")
  }
  if ("group_size" %in% which) bars(diag$group_size, "Group size", "perps per event")
  if ("returning" %in% which) {
    plot(ser$t, ser$returning, pch = 16, cex = 0.4, col = "grey50", xlab = "time",
         ylab = "returning perps", main = "Returning perps per event")
    k <- max(1L, round(nrow(ser) / 20))
    lines(ser$t, stats::filter(ser$returning, rep(1 / (2 * k + 1), 2 * k + 1)), lwd = 2)
  }
  if ("attendance" %in% which) bars(diag$attendance, "Attendance", "events per perp")
  if ("attendance_ccdf" %in% which) {
    cc <- diag$ccdf_attendance
    plot(cc$k, cc$ccdf, log = if (isTRUE(loglog)) "xy" else "", pch = 16, xlab = "k",
         ylab = "P(attendance >= k)", main = "Attendance CCDF")
  }
  invisible(NULL)
}
