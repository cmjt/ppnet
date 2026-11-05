test_that("BA loglik matches a brute-force calculation", {
  set.seed(11)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3)
  ev <- sim_hawkesNet(p, T_end = 12, mark_type = "ba")$ev
  skip_if(nrow(ev$times) < 5L)

  mark <- brute_mark_ll(ev, function(k, t_k, born, deg, role, nodes_k, edges_k) {
    if (length(born) == 0L) return(0)
    pr <- ba_weights(deg, t_k - born, p$beta_edges)
    bern_ll(pr, names(born) %in% neighbours_of(as.character(nodes_k$id), edges_k))
  })

  expect_equal(
    loglik(ev, p, mark_type = "ba", T_end = 12),
    brute_time_ll(ev$times$t, p$mu, p$K, p$beta, 12) + mark,
    tolerance = 1e-8
  )
})

test_that("BA-bip loglik matches a brute-force calculation", {
  set.seed(12)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
  ev <- sim_hawkesNet(p, T_end = 12, mark_type = "ba_bip")$ev
  skip_if(nrow(ev$times) < 5L)

  mark <- brute_mark_ll(ev, function(k, t_k, born, deg, role, nodes_k, edges_k) {
    event_id <- as.character(nodes_k$id[nodes_k$role == "event"])
    ll <- dpois(sum(nodes_k$role == "perp"), p$lambda_new, log = TRUE)
    old <- names(born)[role == "perp"]
    if (length(old) == 0L) return(ll)
    pr <- ba_weights(deg[old], t_k - born[old], p$beta_edges)
    ll + bern_ll(pr, old %in% neighbours_of(event_id, edges_k))
  })

  expect_equal(
    loglik(ev, p, mark_type = "ba_bip", T_end = 12),
    brute_time_ll(ev$times$t, p$mu, p$K, p$beta, 12) + mark,
    tolerance = 1e-8
  )
})

cs_params <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3,
                  node_lambda = 1.2, CS_edges = 0.4)

brute_cs_ll <- function(ev, p, T_end, truncation = Inf) {
  mark <- brute_mark_ll(ev, function(k, t_k, born, deg, role, nodes_k, edges_k) {
    new <- if (is.null(nodes_k)) character(0) else as.character(nodes_k$id)
    ll <- dpois(length(new), p$node_lambda, log = TRUE)
    if (length(born) == 0L || length(new) == 0L) return(ll)
    old <- utils::tail(names(born), truncation)
    pr <- plogis(p$CS_edges) * exp(-p$beta_edges * (t_k - born[old]))
    for (u in new) ll <- ll + bern_ll(pr, old %in% neighbours_of(u, edges_k))
    ll
  })
  brute_time_ll(ev$times$t, p$mu, p$K, p$beta, T_end) + mark
}

test_that("CS loglik matches a brute-force calculation", {
  set.seed(13)
  ev <- sim_hawkesNet(cs_params, T_end = 10, mark_type = "cs", formula_rhs = "edges")$ev
  skip_if(nrow(ev$times) < 5L)

  expect_equal(
    loglik(ev, cs_params, mark_type = "cs", T_end = 10, formula_rhs = "edges"),
    brute_cs_ll(ev, cs_params, 10),
    tolerance = 1e-8
  )
})

test_that("CS-bip loglik matches a brute-force calculation", {
  set.seed(14)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7, CS_edges = 0.4)
  ev <- sim_hawkesNet(p, T_end = 10, mark_type = "cs_bip", formula_rhs = "edges")$ev
  skip_if(nrow(ev$times) < 5L)

  mark <- brute_mark_ll(ev, function(k, t_k, born, deg, role, nodes_k, edges_k) {
    event_id <- as.character(nodes_k$id[nodes_k$role == "event"])
    ll <- dpois(sum(nodes_k$role == "perp"), p$lambda_new, log = TRUE)
    old <- names(born)[role == "perp"]
    if (length(old) == 0L) return(ll)
    pr <- plogis(p$CS_edges) * exp(-p$beta_edges * (t_k - born[old]))
    ll + bern_ll(pr, old %in% neighbours_of(event_id, edges_k))
  })

  expect_equal(
    loglik(ev, p, mark_type = "cs_bip", T_end = 10, formula_rhs = "edges"),
    brute_time_ll(ev$times$t, p$mu, p$K, p$beta, 10) + mark,
    tolerance = 1e-8
  )
})

test_that("the likelihood separates into time and mark parts", {
  set.seed(15)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
  ev <- sim_hawkesNet(p, T_end = 12, mark_type = "ba_bip")$ev
  ll <- function(q) loglik(ev, q, mark_type = "ba_bip", T_end = 12)

  # Changing a time parameter moves the total by the same amount for any mark
  # parameters, because no parameter appears in both parts.
  time2 <- modifyList(p, list(K = 0.2))
  mark2 <- modifyList(p, list(beta_edges = 1, lambda_new = 2))
  expect_equal(ll(time2) - ll(p),
               ll(modifyList(mark2, list(K = 0.2))) - ll(mark2),
               tolerance = 1e-8)
})

test_that("CS candidates are the most recently arrived old nodes", {
  # Arrival order differs from string order: "9" sorts after "10"
  old <- c("9", "10", "11", "12")
  cand <- .cs_candidates_new_old("13", old, truncation = 2)
  expect_setequal(cand$head, c("11", "12"))
  expect_equal(cand$tail, c("13", "13"))

  expect_setequal(.cs_candidates_new_old("13", old, truncation = 10)$head, old)
  expect_setequal(.cs_candidates_new_old("13", old)$head, old)
  expect_error(.cs_candidates_new_old("13", old, truncation = 0), "integer >= 1")
})

test_that("CS fitting uses the same truncation as the simulator", {
  set.seed(16)
  ev <- sim_hawkesNet(cs_params, T_end = 10, mark_type = "cs",
                      formula_rhs = "edges", truncation = 3)$ev
  skip_if(nrow(ev$times) < 5L)

  expect_equal(
    loglik(ev, cs_params, mark_type = "cs", T_end = 10, formula_rhs = "edges", truncation = 3),
    brute_cs_ll(ev, cs_params, 10, truncation = 3),
    tolerance = 1e-8
  )
})

test_that("CS fitting errors when an edge falls outside the truncation", {
  ev <- make_events(
    edges = data.frame(i = c("b", "c", "d"), j = c("a", "b", "a"), time = c(2, 3, 4)),
    nodes = data.frame(id = c("a", "b", "c", "d"), time = c(1, 2, 3, 4))
  )
  expect_error(
    loglik(ev, cs_params, mark_type = "cs", formula_rhs = "edges", truncation = 1),
    "increase `truncation`"
  )
  expect_true(is.finite(
    loglik(ev, cs_params, mark_type = "cs", formula_rhs = "edges", truncation = 3)
  ))
})

test_that("fit_hawkesNet improves on its starting likelihood", {
  set.seed(17)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
  ev <- sim_hawkesNet(p, T_end = 20, mark_type = "ba_bip")$ev
  init <- list(mu = 0.5, K = 0.2, beta = 2, beta_edges = 1, lambda_new = 2)

  fit <- suppressWarnings(fit_hawkesNet(ev, init, mark_type = "ba_bip", T_end = 20,
                                        control = list(maxit = 50)))
  expect_gt(fit$loglik, loglik(ev, init, mark_type = "ba_bip", T_end = 20))
})

test_that("fit_hawkesNet moves free parameters when others are fixed", {
  set.seed(18)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
  ev <- sim_hawkesNet(p, T_end = 20, mark_type = "ba_bip")$ev
  init <- modifyList(p, list(mu = 0.5, K = 0.2))

  fit <- suppressWarnings(fit_hawkesNet(ev, init, fixed_params = c("beta", "beta_edges", "lambda_new"),
                                        mark_type = "ba_bip", T_end = 20))
  expect_gt(fit$loglik, loglik(ev, init, mark_type = "ba_bip", T_end = 20))
  expect_equal(fit$par[c("beta", "beta_edges", "lambda_new")], init[c("beta", "beta_edges", "lambda_new")])
})

test_that("the cached mark part gives the same likelihood and gradient as a fresh one", {
  set.seed(19)
  p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
  ev <- sim_hawkesNet(p, T_end = 15, mark_type = "ba_bip")$ev
  p_time <- modifyList(p, list(mu = 1.3, K = 0.4))     # time parameters only
  p_mark <- modifyList(p, list(beta_edges = 0.6))      # a mark parameter

  fresh <- function(q, gradient = FALSE) {
    .build_loglik_fn(ev, T_end = 15, mark_type = "ba_bip")(q, gradient = gradient)
  }
  f <- .build_loglik_fn(ev, T_end = 15, mark_type = "ba_bip")
  # Alternate calls so each one either reuses or replaces the cached mark part,
  # with and without the gradient.
  for (q in list(p, p_time, p_mark, p_time, p)) {
    expect_equal(f(q), fresh(q))
    expect_equal(f(q, gradient = TRUE), fresh(q, gradient = TRUE))
  }
  expect_equal(f(p_mark, gradient = TRUE), fresh(p_mark, gradient = TRUE))
})
