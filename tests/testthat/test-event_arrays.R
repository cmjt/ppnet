# The array walks must reproduce the network-based evaluators, which are
# reached by passing an explicit starting network.

test_that("array and network BA mark likelihoods agree", {
  set.seed(31)
  ev <- sim_hawkesNet(list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3),
                      T_end = 15, mark_type = "ba")$ev
  arr <- compile_logpmf_ba_array(ev)
  net <- compile_logpmf_ba(ev, net0 = net_init())
  for (be in c(0, 0.3, 2)) {
    expect_equal(arr(list(beta_edges = be)), net(list(beta_edges = be)), tolerance = 1e-10)
  }
})

test_that("array and network BA-bip mark likelihoods agree", {
  set.seed(32)
  ev <- sim_hawkesNet(list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7),
                      T_end = 15, mark_type = "ba_bip")$ev
  arr <- compile_logpmf_ba_bip_array(ev)
  net <- compile_logpmf_ba_bip(ev, net0 = net_init())
  for (be in c(0, 0.3, 2)) {
    q <- list(beta_edges = be, lambda_new = 1.1)
    expect_equal(arr(q), net(q), tolerance = 1e-10)
  }
  expect_equal(arr(list(beta_edges = -1, lambda_new = 1)), -Inf)
})

test_that("loglik gives the same value on either path", {
  set.seed(33)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
  ev <- sim_hawkesNet(p, T_end = 12, mark_type = "ba_bip")$ev
  expect_equal(loglik(ev, p, mark_type = "ba_bip", T_end = 12),
               loglik(ev, p, mark_type = "ba_bip", T_end = 12, net0 = net_init()),
               tolerance = 1e-10)
})

# Event 1: e1-{a, b}; event 2: e2-{b, c}; event 3: e3-{a, c, d};
# event 4: e4-{d}; event 5: e5-{f}, which shares no perp.
bip_events <- function(drop_edge = FALSE) {
  edges <- data.frame(
    i = c("e1", "e1", "e2", "e2", "e3", "e3", "e3", "e4", "e5"),
    j = c("a", "b", "b", "c", "a", "c", "d", "d", "f"),
    time = c(1, 1, 2, 2, 3, 3, 3, 4, 5)
  )
  if (drop_edge) edges <- edges[-4, ]
  nodes <- data.frame(
    id = c("e1", "a", "b", "e2", "c", "e3", "d", "e4", "e5", "f"),
    time = c(1, 1, 1, 2, 2, 3, 3, 4, 5, 5),
    role = c("event", "perp", "perp", "event", "perp", "event", "perp", "event", "event", "perp")
  )
  make_events(edges = edges, nodes = nodes, allow_implicit_birth = FALSE)
}

test_that("event_summary counts perps and linked past events", {
  s <- event_summary(bip_events())
  expect_equal(s$n_perps, c(2, 2, 3, 1, 1))
  expect_equal(s$n_new, c(2, 1, 1, 0, 1))
  expect_equal(s$n_old, c(0, 1, 2, 1, 0))
  # e3 links to e1 (via a) and e2 (via c), and e1 twice is counted once.
  expect_equal(s$connectedness, c(0, 1, 2, 1, 0))
})

test_that("event_summary matches a brute-force count on simulated data", {
  set.seed(34)
  ev <- sim_hawkesNet(list(mu = 2, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7),
                      T_end = 15, mark_type = "ba_bip")$ev
  s <- event_summary(ev)

  arrivals <- nodes_by_event(ev)
  edges <- edges_by_event(ev)
  perps_of <- lapply(seq_along(edges), function(k) {
    focal <- as.character(arrivals[[k]]$id[arrivals[[k]]$role == "event"])
    neighbours_of(focal, edges[[k]])
  })
  brute <- vapply(seq_along(perps_of), function(k) {
    sum(vapply(seq_len(k - 1L), function(m) any(perps_of[[m]] %in% perps_of[[k]]), logical(1)))
  }, numeric(1))

  expect_equal(s$n_perps, lengths(perps_of))
  expect_equal(s$connectedness, brute)
})

test_that("a new perp without an edge to its event makes the likelihood -Inf", {
  ev <- bip_events(drop_edge = TRUE)
  expect_equal(compile_logpmf_ba_bip_array(ev)(list(beta_edges = 0.3, lambda_new = 1)), -Inf)
  expect_equal(compile_logpmf_ba_bip(ev, net0 = net_init())(list(beta_edges = 0.3, lambda_new = 1)), -Inf)
})
