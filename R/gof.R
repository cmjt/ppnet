# Goodness of fit by simulation, and the bipartite summaries it compares.

# Internal: TRUE for a bipartite event-perp history.
.is_bipartite <- function(ev) {
  !is.null(ev$nodes$role) && all(c("event", "perp") %in% ev$nodes$role)
}

# Internal: proportions of x in the categories 0/1..top-1 and "top+".
.binned_props <- function(x, from, top) {
  x <- pmin(x, top)
  levels <- from:top
  out <- tabulate(match(x, levels), nbins = length(levels)) / max(length(x), 1L)
  names(out) <- c(as.character(levels[-length(levels)]), paste0(top, "+"))
  out
}

# Internal: per-event perps (old and new) and summaries of a bipartite
# history, counting only events at or after t_from (the perps' earlier events
# still make them "old").
.bip_stats <- function(ev, t_from = -Inf, T_end = NULL, n_grid = 50L) {
  s <- event_summary(ev)
  hp <- .history_perps(ev)
  keep <- s$t >= t_from
  perps <- Map(c, hp$old, hp$new)[keep]
  t <- s$t[keep]

  # Attendance within the window, and how often co-offending pairs repeat.
  all_perps <- unlist(perps, use.names = FALSE)
  attendance <- as.vector(table(all_perps))
  pairs <- unlist(lapply(perps, function(p) {
    if (length(p) < 2L) return(character(0))
    p <- sort(p)
    cmb <- utils::combn(p, 2L)
    paste(cmb[1L, ], cmb[2L, ], sep = "\r")
  }), use.names = FALSE)
  pair_mult <- as.vector(table(pairs))

  lo <- if (is.finite(t_from)) t_from else min(s$t)
  hi <- T_end %||% max(s$t)
  grid <- seq(lo, hi, length.out = n_grid)
  gaps <- diff(t)
  gaps <- gaps[gaps > 0]

  list(
    count = stats::setNames(findInterval(grid, t), signif(grid, 4)),
    gaps = stats::setNames(if (length(gaps) > 1L) {
      stats::quantile(log10(gaps), c(0.05, 0.25, 0.5, 0.75, 0.95), names = FALSE)
    } else rep(NA_real_, 5L), c("5%", "25%", "50%", "75%", "95%")),
    size = .binned_props(s$n_perps[keep], 1L, 10L),
    returning = .binned_props(s$n_old[keep], 0L, 5L),
    new = .binned_props(s$n_new[keep], 0L, 5L),
    attendance = .binned_props(attendance, 1L, 10L),
    repeat_pairs = .binned_props(pair_mult, 1L, 4L)[-1L],
    scalars = c(events = length(t), perps = length(attendance),
                mean_size = mean(s$n_perps[keep]), mean_returning = mean(s$n_old[keep]),
                mean_new = mean(s$n_new[keep]),
                repeat_pair_share = if (length(pair_mult)) mean(pair_mult >= 2) else NA_real_)
  )
}

# Internal: time-only summaries, for unipartite histories.
.time_stats <- function(ev, t_from = -Inf, T_end = NULL, n_grid = 50L) {
  t <- ev$times$t
  t <- t[t >= t_from]
  lo <- if (is.finite(t_from)) t_from else min(t)
  grid <- seq(lo, T_end %||% max(t), length.out = n_grid)
  gaps <- diff(t)
  gaps <- gaps[gaps > 0]
  list(
    count = stats::setNames(findInterval(grid, t), signif(grid, 4)),
    gaps = stats::setNames(stats::quantile(log10(gaps), c(0.05, 0.25, 0.5, 0.75, 0.95),
                                           names = FALSE), c("5%", "25%", "50%", "75%", "95%")),
    scalars = c(events = length(t))
  )
}

#' Goodness of fit by simulation
#'
#' Simulates from a fitted model and compares summaries of the observed data
#' with their spread over the simulations, in the spirit of `ergm::gof()`.
#' For bipartite event-perp data the summaries cover both parts of the model:
#'
#' * **when events happen**: the cumulative count \eqn{N(t)} on a grid of
#'   times (`count`) and quantiles of the log10 gaps between events (`gaps`);
#' * **who takes part**: the distributions of group size (perps per event,
#'   `size`), returning perps per event (`returning`), new perps per event
#'   (`new`), perps' attendance (events per perp, `attendance`) and how often
#'   co-offending pairs recur (the share of pairs of perps who shared an event
#'   that shared 2, 3 or 4+, `repeat_pairs`), and their means.
#'
#' `print()` gives, for each scalar summary and for each distribution as a
#' whole, a two-sided Monte Carlo p-value. For a distribution the statistic
#' is its distance from the simulated mean, \eqn{\sum_j (x_j - \bar x_j)^2 /
#' s_j^2} over its categories, computed for the data and for each simulation.
#'
#' Unipartite data get the time summaries only.
#'
#' By default each simulation covers the fit's window from scratch, as the
#' fitted model does (perps seen first in the window are new). With
#' `burn_in`, each simulation starts that long before the window and only the
#' events in the window are summarised.
#'
#' @param object A fit from [fit_hawkesNet()]. A fit of one part needs the
#'   other part's parameters in `params`.
#' @param nsim Number of simulations.
#' @param params Parameters to add to or replace the fitted ones, e.g.
#'   `coef(time_fit)` for a fit with `part = "marks"`.
#' @param seed Optional seed, as for [simulate()].
#' @param burn_in Optional burn-in before the window, as in [sim_hawkesNet()].
#' @param n_grid Number of times at which \eqn{N(t)} is compared.
#' @param max_events,max_seconds Limits for each simulation (see
#'   [sim_hawkesNet()]); simulations that reach them are dropped and counted.
#' @param ... Passed to [sim_hawkesNet()].
#'
#' @return An object of class `"ppnet_gof"`: for each summary, the observed
#'   values (`obs`) and a matrix of simulated ones (`sim`, one row per
#'   simulation), plus `n_exploded`. `print()` shows the scalar summaries with
#'   simulation intervals and Monte Carlo p-values; `plot()` shows each
#'   summary's simulated spread (boxplots, or a band for \eqn{N(t)}) with the
#'   observed values.
#' @export
#'
#' @examples
#' set.seed(1)
#' p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
#' ev <- sim_hawkesNet(p, T_end = 30, mark_type = "ba_bip")$ev
#' fit <- fit_hawkesNet(ev, mark_type = "ba_bip")
#' g <- gof(fit, nsim = 20, seed = 1)
#' g
#' plot(g)
gof <- function(object, nsim = 50, params = NULL, seed = NULL, burn_in = 0, n_grid = 50L,
                max_events = 1e6, max_seconds = Inf, ...) {
  if (!inherits(object, "ppnet_fit")) stop("`object` must be a fit from fit_hawkesNet().", call. = FALSE)
  ev <- object$events
  bip <- .is_bipartite(ev)
  T0 <- object$model$T0
  T_end <- object$model$T_end
  stats_of <- function(e) {
    if (bip) .bip_stats(e, T0, T_end, n_grid) else .time_stats(e, T0, T_end, n_grid)
  }

  obs <- stats_of(ev)
  sims <- suppressWarnings(simulate(object, nsim = nsim, seed = seed, params = params,
                                    burn_in = burn_in, max_events = max_events,
                                    max_seconds = max_seconds, ...))
  exploded <- attr(sims, "exploded")
  sims <- sims[!exploded]
  if (!length(sims)) stop("Every simulation exploded; the fitted model may be explosive.", call. = FALSE)
  if (any(exploded)) {
    warning(sum(exploded), " of ", nsim, " simulations exploded and were dropped.", call. = FALSE)
  }
  sim_stats <- lapply(sims, stats_of)

  out <- lapply(stats::setNames(names(obs), names(obs)), function(nm) {
    list(obs = obs[[nm]], sim = do.call(rbind, lapply(sim_stats, `[[`, nm)))
  })
  structure(list(stats = out, nsim = length(sims), n_exploded = sum(exploded), bipartite = bip,
                 call = match.call()),
            class = "ppnet_gof")
}

# Internal: two-sided Monte Carlo p-value of obs among sims.
.mc_p <- function(obs, sims) {
  sims <- sims[is.finite(sims)]
  if (!is.finite(obs) || !length(sims)) return(NA_real_)
  min(1, 2 * min((1 + sum(sims >= obs)) / (length(sims) + 1),
                 (1 + sum(sims <= obs)) / (length(sims) + 1)))
}

# Internal: one-sided Monte Carlo p-value for a distribution: the distance of
# the observed and of each simulated vector from the simulated mean.
.mc_p_dist <- function(obs, sim) {
  m <- colMeans(sim, na.rm = TRUE)
  v <- apply(sim, 2L, stats::var, na.rm = TRUE)
  use <- is.finite(v) & v > 0
  if (!any(use)) return(NA_real_)
  d <- function(x) sum((x[use] - m[use])^2 / v[use])
  d_obs <- d(obs)
  d_sim <- apply(sim, 1L, d)
  (1 + sum(d_sim >= d_obs)) / (length(d_sim) + 1)
}

#' @export
print.ppnet_gof <- function(x, digits = 3L, ...) {
  sc <- x$stats$scalars
  sim <- sc$sim
  tab <- data.frame(observed = sc$obs, sim_mean = colMeans(sim, na.rm = TRUE),
                    sim_2.5 = apply(sim, 2L, stats::quantile, 0.025, na.rm = TRUE),
                    sim_97.5 = apply(sim, 2L, stats::quantile, 0.975, na.rm = TRUE),
                    p_value = vapply(seq_along(sc$obs), function(j) .mc_p(sc$obs[j], sim[, j]),
                                     numeric(1)))
  dists <- setdiff(names(x$stats), "scalars")
  dist_p <- vapply(dists, function(nm) .mc_p_dist(x$stats[[nm]]$obs, x$stats[[nm]]$sim),
                   numeric(1))
  cat("Goodness of fit:", x$nsim, "simulations")
  if (x$n_exploded) cat(" (", x$n_exploded, " exploded and dropped)", sep = "")
  cat("\n\nSummaries:\n")
  print(noquote(apply(as.matrix(tab), 2L, formatC, digits = digits, format = "g")))
  cat("\nDistributions (Monte Carlo p-value of the distance from the simulations):\n")
  print(data.frame(p_value = dist_p, row.names = dists), digits = digits)
  cat("\nplot() shows each distribution against its simulations.\n")
  invisible(x)
}

#' @param x A `"ppnet_gof"` object.
#' @param which Summaries to plot (default all but the scalars).
#' @rdname gof
#' @export
plot.ppnet_gof <- function(x, which = setdiff(names(x$stats), "scalars"), ...) {
  titles <- c(count = "Cumulative events N(t)", gaps = "Gaps between events (log10 quantiles)",
              size = "Group size (perps per event)", returning = "Returning perps per event",
              new = "New perps per event", attendance = "Attendance (events per perp)",
              repeat_pairs = "Co-offending pairs sharing k events")
  which <- intersect(which, names(x$stats))
  n <- length(which)
  old <- graphics::par(mfrow = grDevices::n2mfrow(n), mar = c(4, 4, 2.5, 1))
  on.exit(graphics::par(old), add = TRUE)
  for (nm in which) {
    st <- x$stats[[nm]]
    if (nm == "count") {
      tt <- as.numeric(names(st$obs))
      band <- apply(st$sim, 2L, stats::quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
      plot(tt, st$obs, type = "n", ylim = range(band, st$obs), xlab = "time", ylab = "N(t)",
           main = titles[[nm]])
      graphics::polygon(c(tt, rev(tt)), c(band[1, ], rev(band[3, ])), col = "grey85", border = NA)
      graphics::lines(tt, band[2, ], lty = 2)
      graphics::lines(tt, st$obs, lwd = 2)
    } else {
      graphics::boxplot(st$sim, names = names(st$obs), outline = FALSE, col = "grey90",
                        border = "grey40", main = titles[[nm]],
                        ylim = range(st$sim, st$obs, na.rm = TRUE),
                        xlab = if (nm == "gaps") "quantile" else "",
                        ylab = if (nm == "gaps") "log10 gap" else "proportion")
      graphics::lines(seq_along(st$obs), st$obs, type = "b", pch = 16, lwd = 2)
    }
  }
  invisible(x)
}
