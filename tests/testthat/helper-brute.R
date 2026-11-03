# Every likelihood here is recomputed from the raw event tables, without using
# the cached evaluators, so that those evaluators can be changed safely.

brute_time_ll <- function(t, mu, K, beta, T_end) {
  lam <- vapply(seq_along(t), function(i) {
    mu + K * sum(exp(-beta * (t[i] - t[seq_len(i - 1L)])))
  }, numeric(1))
  sum(log(lam)) - mu * T_end - (K / beta) * sum(1 - exp(-beta * (T_end - t)))
}

# Walk the events, tracking arrival time, role and degree of every node.
# `event_ll(k, t_k, born, deg, role, nodes_k, edges_k)` gives one event's mark term.
brute_mark_ll <- function(ev, event_ll) {
  arrivals <- nodes_by_event(ev)
  edges <- edges_by_event(ev)
  t <- ev$times$t
  born <- numeric(0); deg <- numeric(0); role <- character(0)

  ll <- 0
  for (k in seq_along(t)) {
    nodes_k <- arrivals[[k]]
    edges_k <- edges[[k]]
    ll <- ll + event_ll(k, t[k], born, deg, role, nodes_k, edges_k)

    if (!is.null(nodes_k) && nrow(nodes_k) > 0L) {
      ids <- as.character(nodes_k$id)
      born[ids] <- t[k]
      deg[ids] <- 0
      role[ids] <- if ("role" %in% names(nodes_k)) as.character(nodes_k$role) else "node"
    }
    if (!is.null(edges_k) && nrow(edges_k) > 0L) {
      ends <- table(c(as.character(edges_k$i), as.character(edges_k$j)))
      deg[names(ends)] <- deg[names(ends)] + as.numeric(ends)
    }
  }
  ll
}

neighbours_of <- function(id, edges_k) {
  if (is.null(edges_k) || nrow(edges_k) == 0L) return(character(0))
  i <- as.character(edges_k$i)
  j <- as.character(edges_k$j)
  c(j[i == id], i[j == id])
}

bern_ll <- function(p, obs) sum(ifelse(obs, log(p), log1p(-p)))

# sna::degree() counts each undirected edge twice, and the kernels use it as is.
ba_weights <- function(deg, age, beta_edges, delta = 0.001) {
  w <- (2 * deg + delta) * exp(-beta_edges * age)
  w / sum(w)
}


# BA-bip marks can have no perps, which size excitation cannot use. Redraw
# them until they have one; the time likelihood treats marks as observed
# history, so this is still a valid history for testing the time part.
local_min_one_perp <- function(env = parent.frame()) {
  orig <- sim_mark_ba_bip
  local_mocked_bindings(
    sim_mark_ba_bip = function(...) {
      repeat {
        mk <- orig(...)
        if (nrow(mk$edges) > 0L) return(mk)
      }
    },
    .env = env
  )
}
