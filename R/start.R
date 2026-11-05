# Default starting values and the split of parameters into the time part and
# the mark part.

# Internal: names of the time-part parameters for a model.
.time_param_names <- function(kernel = "exp", baseline = NULL, excitation = NULL) {
  c("mu", "K", "beta", if (identical(kernel, "exp2")) c("K_fast", "beta_fast"),
    baseline_coef_names(baseline), excitation_coef_names(excitation))
}

# Internal: starting values worked out from the data, for every parameter of
# the model part(s) being fitted, except the CS_<statistic> coefficients of the
# ERGM-term mark types (those are left to the user).
#   - Time part: half the events from the baseline (mu = rate / 2), a total
#     branching ratio of 1/2, and decay times of 10 mean gaps (exp), or 50
#     and 0.1 mean gaps (exp2). Baseline and excitation coefficients 0.
#   - Mark part: m and the new-perp rates from the mean numbers of old and
#     new perps per event, beta_edges = 10 / (window length), delta = 1,
#     theta terms 0, size_new = 1, random_effect_var = 0.5.
.default_start <- function(events, part = "both", mark_type = "ba", kernel = "exp",
                           baseline = NULL, excitation = NULL, attachment = NULL,
                           T0 = 0, T_end = NULL) {
  t <- events$times$t
  len <- (T_end %||% max(t)) - T0
  rate <- length(t) / len
  out <- list()

  if (part != "marks") {
    out <- if (identical(kernel, "exp2")) {
      list(mu = rate / 2, K = 0.3 * rate / 50, beta = rate / 50,
           K_fast = 0.2 * 10 * rate, beta_fast = 10 * rate)
    } else {
      list(mu = rate / 2, K = 0.5 * rate / 10, beta = rate / 10)
    }
    zero <- c(baseline_coef_names(baseline), excitation_coef_names(excitation))
    out[zero] <- 0
  }

  if (part != "time") {
    pos <- function(x) if (is.finite(x) && x > 0.05) x else 0.05
    beta_edges <- 10 / len
    out$beta_edges <- beta_edges
    if (mark_type %in% c("ba_bip", "cs_bip")) {
      s <- event_summary(events)
      if (is.null(attachment)) {
        out$lambda_new <- pos(mean(s$n_new))
      } else {
        shifted <- attachment$min_perps == 1L & s$n_old == 0L
        out$beta_edges <- NULL          # re-added below in the spec's order
        defaults <- list(
          m = pos(mean(s$n_old)), beta_edges = beta_edges, delta = 1,
          lambda_new = pos(mean(s$n_new[!shifted])),
          lambda_new0 = pos(mean(s$n_new[shifted] - 1)), size_new = 1,
          random_effect_var = 0.5
        )
        for (nm in attachment_param_names(attachment)) {
          out[[nm]] <- defaults[[nm]] %||% 0      # theta terms start at 0
        }
      }
    } else if (mark_type == "cs") {
      new <- tabulate(match(events$nodes$event_id, events$times$event_id),
                      nbins = length(t))
      out$node_lambda <- pos(mean(new))
    }
  }
  out
}
