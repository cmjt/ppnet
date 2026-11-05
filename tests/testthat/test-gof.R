sp_g <- attachment_spec(coffenders = TRUE, coffender_window = 20, new_perps = "negbin", min_perps = 1)
p_g <- list(mu = 3, K = 0.6, beta = 1.5, m = 0.5, beta_edges = 0.1, delta = 0.5, theta = 0.3,
            lambda_new = 0.6, lambda_new0 = 0.4, size_new = 2)
set.seed(401)
ev_g <- sim_hawkesNet(p_g, T_end = 60, mark_type = "ba_bip", attachment = sp_g, return_net = FALSE)$ev

test_that("bipartite network diagnostics keep events and perps apart", {
  d <- diag_network(ev_g)
  s <- event_summary(ev_g)
  expect_true(d$bipartite)
  expect_equal(d$group_size, s$n_perps)
  expect_equal(d$series$returning, s$n_old)
  expect_equal(sum(d$attendance), sum(s$n_perps))          # perp-event links
  expect_equal(length(d$attendance), sum(s$n_new))           # perps
  expect_equal(tail(d$series$perps_n, 1), sum(s$n_new))
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot_network_diag(d))
  for (type in c("group_size", "returning", "attendance", "attendance_ccdf", "degree_hist")) {
    expect_silent(plot(ev_g, type = type))
  }
  set.seed(402)
  ev_u <- sim_hawkesNet(list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3), T_end = 10)$ev
  expect_error(plot(ev_u, type = "group_size"), "bipartite")
})

test_that("summaries of a history are what they say", {
  st <- .bip_stats(ev_g)
  s <- event_summary(ev_g)
  expect_equal(unname(st$scalars["events"]), nrow(s))
  expect_equal(unname(st$scalars["mean_size"]), mean(s$n_perps))
  expect_equal(unname(st$size[1]), mean(s$n_perps == 1))
  expect_equal(unname(st$returning[1]), mean(s$n_old == 0))
  expect_equal(sum(st$size), 1)
  expect_equal(unname(tail(st$count, 1)), nrow(s))
  # Only events in the window count.
  st2 <- .bip_stats(ev_g, t_from = 30)
  expect_equal(unname(st2$scalars["events"]), sum(s$t >= 30))
})

test_that("gof() passes the true model and flags a wrong one", {
  fit <- fit_hawkesNet(ev_g, mark_type = "ba_bip", attachment = sp_g)
  g <- gof(fit, nsim = 30, seed = 1)
  expect_s3_class(g, "ppnet_gof")
  expect_equal(g$nsim, 30L)
  expect_equal(dim(g$stats$size$sim), c(30L, 10L))
  p_true <- vapply(setdiff(names(g$stats), "scalars"),
                   function(nm) .mc_p_dist(g$stats[[nm]]$obs, g$stats[[nm]]$sim), numeric(1))
  expect_true(all(p_true > 0.05))
  # Poisson new perps with events allowed to have none: the group sizes and
  # new perps per event do not match.
  fit0 <- fit_hawkesNet(ev_g, mark_type = "ba_bip", attachment = attachment_spec(degree = FALSE))
  g0 <- gof(fit0, nsim = 30, seed = 1)
  p_wrong <- vapply(c("size", "new"),
                    function(nm) .mc_p_dist(g0$stats[[nm]]$obs, g0$stats[[nm]]$sim), numeric(1))
  expect_true(all(p_wrong < 0.05))
  out <- capture.output(print(g))
  expect_true(any(grepl("Distributions", out)))
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(g))
})

test_that("gof() handles one-part fits, unipartite data and exploding simulations", {
  fm <- fit_hawkesNet(ev_g, mark_type = "ba_bip", attachment = sp_g, part = "marks")
  expect_error(gof(fm, nsim = 2), "give the time parameters")
  ft <- fit_hawkesNet(ev_g, mark_type = "ba_bip", part = "time")
  g <- gof(fm, nsim = 3, seed = 1, params = coef(ft))
  expect_equal(g$nsim, 3L)
  set.seed(403)
  ev_u <- sim_hawkesNet(list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3), T_end = 20)$ev
  gu <- gof(fit_hawkesNet(ev_u, mark_type = "ba"), nsim = 3, seed = 1)
  expect_false(gu$bipartite)
  expect_setequal(names(gu$stats), c("count", "gaps", "scalars"))
  expect_error(gof(fit_hawkesNet(ev_g, mark_type = "ba_bip", attachment = sp_g), nsim = 3,
                   seed = 1, max_events = 20), "Every simulation exploded")
})
