# sim_perp_null.R

#' Simulate perp-driven events with no network effect
#'
#' A null model for the event-level connectedness effect `gamma` (see
#' [excitation_spec()]). Events here come from perps offending, so prolific
#' perps make events that are both well connected and followed by more
#' events, but no perp's rate depends on the network. Under this model the
#' true `gamma` is zero, and a fitted `gamma > 0` comes from confounding by
#' perp propensity.
#'
#' Perps enter at rate `entry_rate` from time `-burn_in` and stay active for
#' an exponential career with mean `career_mean`. While active, perp \eqn{v}
#' starts offences at rate
#' \deqn{r_v(t) = \text{rate}\, u_v + K_{self}\, \beta_{self}
#'   \sum_{s} e^{-\beta_{self}(t - s)},}
#' with a gamma random effect \eqn{u_v \sim \text{Gamma}(\text{shape}, \text{shape})}
#' (mean 1, variance 1 / shape) and the sum over the earlier events \eqn{s}
#' that \eqn{v} started (`excite_on = "own"`) or attended (`"all"`). With
#' `"all"`, being recruited into a group raises a perp's later rate, so group
#' events lead to more events through the perps at them. That is still an
#' individual effect, not a connectedness one, but it is no longer clean of
#' event size.
#'
#' Each offence is one event. Its initiator brings
#' \eqn{\text{Poisson}(\text{extra\_mean}\, u_v^{\text{size\_link}})}
#' co-offenders (`size_link > 0` links the random effect to group size). Each slot is
#' filled, with probability `p_past`, by one of the initiator's earlier
#' co-offenders who is still active and not already in the group, and
#' otherwise by an active perp chosen with probability proportional to
#' \eqn{u^{\text{recruit\_link}}}. A slot with no one left to fill it is
#' dropped.
#'
#' Events before time 0 build up the perps' histories and excitation but are
#' not returned, as in data observed from some date on. The returned
#' network therefore starts empty at time 0, as the fitted models assume.
#'
#' @param T_end End of the observation window `[0, T_end]`.
#' @param entry_rate New perps per unit time.
#' @param career_mean Mean career length; `Inf` for perps who never stop.
#' @param rate Base rate at which a perp with `u = 1` starts offences.
#' @param random_effect_shape Shape of the gamma random effect; `Inf` for none.
#' @param K_self,beta_self Per-perp self-excitation: expected number of
#'   offences each triggering event adds, and the decay rate.
#' @param excite_on `"own"` or `"all"`: which of a perp's events excite it.
#' @param extra_mean Mean number of co-offenders for a perp with `u = 1`.
#' @param size_link Power of `u` scaling the initiator's co-offender mean.
#' @param p_past Probability a co-offender slot goes to a past co-offender.
#' @param recruit_link Power of `u` weighting who is recruited from the
#'   active perps (0 for uniform).
#' @param burn_in Length of the unobserved run-in before time 0 (needs a
#'   value when `career_mean = Inf`).
#'
#' @return A list with `ev`, an events object for `mark_type = "ba_bip"`,
#'   and `perps`, a data frame with each perp's `id` (`NA` if never seen in
#'   the window), `u`, `start`, `end` and `n_events` in the window.
#' @export
#'
#' @examples
#' set.seed(1)
#' s <- sim_perp_events(T_end = 52, random_effect_shape = 1)
#' head(event_summary(s$ev))
sim_perp_events <- function(T_end,
                            entry_rate = 1,
                            career_mean = 104,
                            rate = 0.05,
                            random_effect_shape = Inf,
                            K_self = 0,
                            beta_self = 0.5,
                            excite_on = c("own", "all"),
                            extra_mean = 0.5,
                            size_link = 0,
                            p_past = 0.5,
                            recruit_link = 1,
                            burn_in = 3 * career_mean) {
  excite_on <- match.arg(excite_on)
  stopifnot(T_end > 0, entry_rate > 0, career_mean > 0, rate > 0, random_effect_shape > 0,
            K_self >= 0, K_self < 1, beta_self > 0, extra_mean >= 0,
            p_past >= 0, p_past <= 1, is.finite(burn_in), burn_in >= 0)

  # Perps: entry times, careers and random effects, all drawn up front.
  T0 <- -burn_in
  n_perp <- stats::rpois(1, entry_rate * (T_end - T0))
  start <- sort(stats::runif(n_perp, T0, T_end))
  end <- start + if (is.finite(career_mean)) stats::rexp(n_perp, 1 / career_mean) else Inf
  u <- if (is.finite(random_effect_shape)) {
    stats::rgamma(n_perp, random_effect_shape, random_effect_shape)
  } else {
    rep(1, n_perp)
  }
  w_size <- u^size_link
  w_recruit <- u^recruit_link

  # Times at which the active set changes; the background rate is constant
  # between them and the excitation only decays, so the current total rate
  # bounds the rate up to the next one.
  changes <- sort(c(start, end[end < T_end], T_end))
  next_change <- 1L

  active <- logical(n_perp)
  S <- numeric(n_perp)            # self-excitation sums at time t
  partners <- vector("list", n_perp)
  node_id <- rep(NA_character_, n_perp)
  n_seen <- 0L
  n_window <- integer(n_perp)

  ev_t <- numeric(0)
  rows_nodes <- list()
  rows_edges <- list()
  n_ev <- 0L

  t <- T0
  kb <- K_self * beta_self
  repeat {
    B <- rate * sum(u[active])
    M <- B + kb * sum(S[active])
    t_next <- if (M > 0) t + stats::rexp(1, M) else Inf
    if (t_next >= changes[next_change]) {
      # Move to the change point without an event.
      tc <- changes[next_change]
      S <- S * exp(-beta_self * (tc - t))
      t <- tc
      if (t >= T_end) break
      active <- start <= t & end > t
      next_change <- next_change + 1L
      next
    }
    S <- S * exp(-beta_self * (t_next - t))
    t <- t_next
    r_v <- (rate * u + kb * S) * active
    total <- sum(r_v)
    if (stats::runif(1) * M > total) next

    # An offence: pick its initiator, then the co-offenders.
    v <- .sample1(which(active), r_v[active])
    group <- v
    n_extra <- stats::rpois(1, extra_mean * w_size[v])
    for (s in seq_len(n_extra)) {
      past <- partners[[v]]
      past <- past[active[past] & !(past %in% group)]
      if (length(past) > 0L && stats::runif(1) < p_past) {
        w <- .sample1(past)
      } else {
        pool <- which(active)
        pool <- pool[!(pool %in% group)]
        if (length(pool) == 0L) break
        w <- .sample1(pool, w_recruit[pool])
      }
      group <- c(group, w)
    }

    for (g in group) partners[[g]] <- union(partners[[g]], setdiff(group, g))
    if (excite_on == "own") S[v] <- S[v] + 1 else S[group] <- S[group] + 1

    if (t >= 0) {
      n_ev <- n_ev + 1L
      ev_t[n_ev] <- t
      event_id <- paste0("e", n_ev)
      new <- group[is.na(node_id[group])]
      node_id[new] <- paste0("p", n_seen + seq_along(new))
      n_seen <- n_seen + length(new)
      n_window[group] <- n_window[group] + 1L
      rows_nodes[[n_ev]] <- data.frame(id = c(event_id, node_id[new]),
                                       role = c("event", rep("perp", length(new))),
                                       time = t, stringsAsFactors = FALSE)
      rows_edges[[n_ev]] <- data.frame(i = event_id, j = node_id[group], time = t,
                                       stringsAsFactors = FALSE)
    }
  }

  if (n_ev == 0L) stop("No events in the observation window.", call. = FALSE)
  ev <- make_events(nodes = .bind_rows_fast(rows_nodes), edges = .bind_rows_fast(rows_edges),
                    allow_implicit_birth = FALSE)
  list(ev = ev,
       perps = data.frame(id = node_id, u = u, start = start, end = end,
                          n_events = n_window, stringsAsFactors = FALSE))
}

# Internal: one draw from x (safe when x has length 1), optionally weighted.
.sample1 <- function(x, prob = NULL) {
  if (length(x) == 1L) return(x)
  x[sample.int(length(x), 1L, prob = prob)]
}
