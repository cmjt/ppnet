#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(devtools)
  library(network)
  library(ernm)
  devtools::load_all()
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: run_cs_rep.R <rep_id> <out_dir>")

rep_id <- as.integer(args[[1]])
out_dir <- args[[2]]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- reproducible per-rep seed ----
seed0 <- 100000 + rep_id
set.seed(seed0)

# ---- CS sim-study settings (as provided) ----
T_end <- 15
formula_rhs <- "edges + star(c(2,3))"

params_true <- list(
  mu = 2.5,
  K = 0.5,
  beta = 1,
  beta_edges = 1,
  node_lambda = 4,
  CS_edges = -2.5,
  CS_star.2 = -1,
  CS_star.3 = -3
)

params_init <- list(
  mu = 1,
  K = 1,
  beta = 2,
  beta_edges = 2,
  node_lambda = 5,
  CS_edges = -1,
  CS_star.2 = 0,
  CS_star.3 = 0
)

# Explicitly set transforms to avoid warnings:
transform <- list(
  mu = "log",
  K = "log",
  beta = "log",
  beta_edges = "log",
  node_lambda = "log",
  CS_edges = "none",
  CS_star.2 = "none",
  CS_star.3 = "none"
)

# ---- small utilities ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_time <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  out <- eval.parent(substitute(expr))
  t1 <- proc.time()[["elapsed"]]
  list(value = out, elapsed = unname(t1 - t0))
}

# Compact network stats (avoid heavy ops if network is huge)
net_stats <- function(net) {
  n_nodes <- network::network.size(net)
  n_edges <- network::network.edgecount(net)

  deg <- tryCatch({
    if (requireNamespace("sna", quietly = TRUE)) {
      as.numeric(sna::degree(net, gmode = "graph"))
    } else {
      as.numeric(network::degree(net, cmode = "freeman"))
    }
  }, error = function(e) rep(NA_real_, n_nodes))

  list(
    n_nodes = n_nodes,
    n_edges = n_edges,
    mean_degree = if (all(is.na(deg))) NA_real_ else mean(deg),
    max_degree  = if (all(is.na(deg))) NA_real_ else max(deg),
    density = {
      d <- tryCatch(network::gden(net), error = function(e) NA_real_)
      if (!is.na(d)) d else {
        n <- n_nodes
        m <- n_edges
        directed <- isTRUE(network::is.directed(net))
        denom <- if (directed) n * (n - 1) else n * (n - 1) / 2
        if (denom > 0) m / denom else NA_real_
      }
    }
  )
}

# ---- result container ----
res <- list(
  rep_id = rep_id,
  seed = seed0,
  mark_type = "cs",
  T_end = T_end,
  formula_rhs = formula_rhs,
  params_true = params_true,
  params_init = params_init,
  status = "ok",
  error = NA_character_,
  trace = NA_character_
)

# ---- simulate ----
sim_out <- tryCatch({
  safe_time(
    sim_hawkesNet(
      params = params_true,
      T_end = T_end,
      mark_type = "cs",
      formula_rhs = formula_rhs
    )
  )
}, error = function(e) {
  msg <- paste0("[sim error] ", paste(class(e), collapse = "/"), ": ", conditionMessage(e))
  message(msg)

  trace <- paste(vapply(sys.calls(), deparse, ""), collapse = "\n")
  res$status <- "error_sim"
  res$error  <- msg
  res$trace  <- trace

  saveRDS(res, file = file.path(out_dir, sprintf("cs_rep_%03d.rds", rep_id)))
  quit(status = 1)
})

sim <- sim_out$value
res$sim_elapsed <- sim_out$elapsed

ev  <- sim$ev %||% NULL
net <- sim$net %||% NULL
# print(ev)
res$n_events <- tryCatch({
  if (!is.null(ev) && !is.null(ev$times) && "t" %in% names(ev$times)) nrow(ev$times) else NA_integer_
}, error = function(e) NA_integer_)

if (!is.null(net)) {
  st <- net_stats(net)
  res <- c(res, st)
} else {
  res$n_nodes <- NA_integer_
  res$n_edges <- NA_integer_
  res$mean_degree <- NA_real_
  res$max_degree <- NA_real_
  res$density <- NA_real_
}
# ---- fit ----
fit_out <- tryCatch({
  safe_time(
    fit_hawkesNet(
      ev = ev,
      params_init = params_init,
      mark_type = "cs",
      transform = transform,
      formula_rhs = formula_rhs
    )
  )
}, warning = function(w) {
  # record warning but keep going
  res$fit_warning <- paste0("[fit warning] ", conditionMessage(w))
  invokeRestart("muffleWarning")
}, error = function(e) {
  msg <- paste0("[fit error] ", paste(class(e), collapse = "/"), ": ", conditionMessage(e))
  message(msg)

  trace <- paste(vapply(sys.calls(), deparse, ""), collapse = "\n")
  res$status <- "error_fit"
  res$error  <- msg
  res$trace  <- trace

  saveRDS(res, file = file.path(out_dir, sprintf("cs_rep_%03d.rds", rep_id)))
  quit(status = 1)
})

fit <- fit_out$value
res$fit_elapsed <- fit_out$elapsed

# ---- convergence / loglik ----
res$convergence <- fit$convergence %||% NA_integer_
res$converged <- isTRUE(!is.na(res$convergence) && res$convergence == 0L)
res$loglik <- fit$loglik %||% NA_real_

res$fit_message <- fit$message %||% NA_character_
res$fit_counts  <- fit$counts %||% NULL

# ---- parameter estimates ----
p_hat <- fit$par %||% list()
res$mu_hat <- p_hat$mu %||% NA_real_
res$K_hat  <- p_hat$K %||% NA_real_
res$beta_hat <- p_hat$beta %||% NA_real_
res$beta_edges_hat <- p_hat$beta_edges %||% NA_real_
res$node_lambda_hat <- p_hat$node_lambda %||% NA_real_
res$CS_edges_hat <- p_hat$CS_edges %||% NA_real_
res$CS_star.2_hat <- p_hat$CS_star.2 %||% NA_real_
res$CS_star.3_hat <- p_hat$CS_star.3 %||% NA_real_

# ---- misc ----
res$has_net <- !is.null(net)
res$has_events <- !is.null(ev)

saveRDS(res, file = file.path(out_dir, sprintf("cs_rep_%03d.rds", rep_id)))
