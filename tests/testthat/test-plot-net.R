set.seed(501)
ev_n <- sim_hawkesNet(list(mu = 2, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7),
                      T_end = 20, mark_type = "ba_bip")$ev

test_that("bipartite networks get the ppnet_net class and stay networks", {
  net <- create_net_from_events(ev_n)
  expect_s3_class(net, "ppnet_net")
  expect_s3_class(net, "network")
  expect_equal(network::network.size(net), nrow(ev_n$nodes))
  expect_s3_class(network::network.copy(net), "network")
  expect_equal(nrow(create_events_from_net(net)$times), nrow(ev_n$times))
  # Unipartite networks are left as they are.
  set.seed(502)
  ev_u <- sim_hawkesNet(list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3), T_end = 10)$ev
  expect_false(inherits(create_net_from_events(ev_u), "ppnet_net"))
})

test_that("net_components() ranks the components and counts their nodes", {
  comp <- net_components(ev_n)
  expect_equal(sum(comp$n_events), nrow(ev_n$times))
  expect_equal(sum(comp$n_perps), sum(ev_n$nodes$role == "perp"))
  size <- comp$n_events + comp$n_perps
  expect_false(is.unsorted(rev(size)))                        # largest first
  memb <- attr(comp, "membership")
  expect_equal(tabulate(memb), size)
  # Each component is connected and no edge crosses components.
  net <- create_net_from_events(ev_n)
  el <- network::as.edgelist(net)
  expect_true(all(memb[el[, 1]] == memb[el[, 2]]))
})

test_that("plot() draws the whole network or chosen components", {
  net <- create_net_from_events(ev_n)
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(net))
  expect_silent(plot(net, component = "largest"))
  co <- plot(net, component = 1:3)
  expect_length(co, 3L)
  expect_silent(plot(net, component = 1, size_by_degree = FALSE, legend = FALSE,
                     displaylabels = TRUE))
  expect_error(plot(net, component = 10000), "components")
})
