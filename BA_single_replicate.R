#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(devtools)
  library(network)
  library(ernm)
  devtools::load_all()
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: run_ba_rep.R <rep_id> <out_dir>")

rep_id <- as.integer(args[[1]])
out_dir <- args[[2]]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- reproducible per-rep seed ----
seed0 <- 100000 + rep_id
set.seed(seed0)

# ---- paper params (as provided) ----
T_end <- 5
params_true <- list(
  mu = 10,
  K = 0.5,
  beta = 1,
  beta_edges = 1
)

# ---- empirical "blow-up" caps (for flags only) ----
cap_events <- 200000L
cap_edges  <- 200000L

# ---- small utilities ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_time <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- expr
  t1 <- proc.time()[["elapsed"]]
  list(value = val, seconds = t1 - t0)
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

# ---- result object (always saved) ----
res <- list(
  rep_id = rep_id,
  seed = seed0,
  T_end = T_end,
  params_true = params_true,
  status = "ok"
)

# ---- load package (edit if needed) ----
# If running inside the ppnet package repo, this is often enough:
# devtools::load_all()
# Otherwise assume it's installed:

# ---- simulate ----
sim_out <- tryCatch({
  safe_time(
    sim_hawkesNet(
      params = params_true,
      T_end = T_end,
      mark_type = "ba"
      # , max_events = cap_events, max_edges = cap_edges
    )
  )
}, error = function(e) {
  msg <- paste0("[sim error] ", paste(class(e), collapse = "/"), ": ", conditionMessage(e))
  message(msg)

  # capture a useful trace (works inside tryCatch)
  trace <- paste(vapply(sys.calls(), deparse, ""), collapse = "\n")

  res$status <- "error_sim"
  res$error  <- msg
  res$trace  <- trace

  # always write the rds so you can diagnose the failed reps
  saveRDS(res, file = file.path(out_dir, sprintf("ba_rep_%03d.rds", rep_id)))

  NULL
})

if (is.null(sim_out)) {
  print("HUHHSSF")
  saveRDS(res, file = file.path(out_dir, sprintf("ba_rep_%03d.rds", rep_id)))
  quit(status = 0)
}

sim <- sim_out$value
res$sim_seconds <- sim_out$seconds

# Extract events + net (defensive)
ev  <- sim$ev %||% sim$events %||% NULL
net <- sim$net %||% NULL

# ---- realised size ----
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

# ---- empirical blow-up flags ----
res$hit_cap_events <- isTRUE(!is.na(res$n_events) && res$n_events >= cap_events)
res$hit_cap_edges  <- isTRUE(!is.na(res$n_edges)  && res$n_edges  >= cap_edges)

# ---- temporal diagnostics summaries (cheap) ----
# Uses your internal diag function (exported as diag_temporal in diagnostics_temporal.R)
res$dt_mean <- NA_real_
res$dt_median <- NA_real_
res$lambda_mean <- NA_real_
res$lambda_max <- NA_real_

if (!is.null(ev)) {
  td <- tryCatch(
    diag_temporal(events = ev, params = params_true, T_end = T_end, T0 = 0),
    error = function(e) NULL
  )
  if (!is.null(td)) {
    res$dt_mean   <- mean(td$dt, na.rm = TRUE)
    res$dt_median <- stats::median(td$dt, na.rm = TRUE)
    res$lambda_mean <- mean(td$lambda, na.rm = TRUE)
    res$lambda_max  <- max(td$lambda, na.rm = TRUE)
  }
}

# ---- fit ----
# Choose reasonable starting values (slightly perturbed)
params_init <- list(
  mu = 8,
  K = 0.7,
  beta = 1.5,
  beta_edges = 0.7
)

fit_out <- tryCatch({
  safe_time(
    fit_hawkesNet(
      events = ev,
      params_init = params_init,
      mark_type = "ba",
      control = list(maxit = 2000)
    )
  )
}, error = function(e) {
  res$status <- "error_fit"
  res$error <- conditionMessage(e)
  NULL
})

if (is.null(fit_out)) {
  saveRDS(res, file = file.path(out_dir, sprintf("ba_rep_%03d.rds", rep_id)))
  quit(status = 0)
}

fit <- fit_out$value
res$fit_seconds <- fit_out$seconds

# ---- fit diagnostics ----
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

# ---- optional: save a compact fingerprint of the run ----
# (useful to re-link to raw sim object if you later decide to save them)
res$has_net <- !is.null(net)
res$has_events <- !is.null(ev)

saveRDS(res, file = file.path(out_dir, sprintf("ba_rep_%03d.rds", rep_id)))
