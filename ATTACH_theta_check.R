# Simulation check: does unmodelled perp frailty make the co-offender effect
# `theta` look positive when there is none?
#
# Each replicate simulates with theta = 0 and perp frailty u ~ Gamma(shape,
# shape) (variance 1 / shape), fits the attachment model (which has no
# frailty term), and records theta's estimate and Wald test. The time part is
# held at its true values: it shares no parameters with the marks, so fitting
# it would not change the theta estimates.
#
# Compare `attendance_dispersion` with the same quantity in your data,
#   a <- as.numeric(table(perp ids over event-perp edges)); var(a) / mean(a),
# to see which frailty level is realistic.
#
# Run from inside the ppnet repo, with the package installed.

library(ppnet)
library(parallel)

n_rep <- 20
frailty_shapes <- c(Inf, 2, 1, 0.5)   # frailty variance 0, 0.5, 1, 2
T_end <- 520                          # weeks
cores <- max(1L, detectCores() - 1L)

spec <- attachment_spec(degree = TRUE, coffenders = TRUE, coffender_window = 52,
                        new_perps = "negbin", min_perps = 1)
# Generic values; replace with your fitted estimates (keeping theta = 0) to
# match your data. Don't commit those.
truth <- list(mu = 10, K = 0.2, beta = 0.5, K_fast = 5, beta_fast = 50,
              m = 0.5, beta_edges = 0.05, delta = 0.5, theta = 0,
              lambda_new = 0.6, lambda_new0 = 0.5, size_new = 2)
time_names <- c("mu", "K", "beta", "K_fast", "beta_fast")

one_rep <- function(shape, rep) {
  set.seed(1000 * rep + round(10 * min(shape, 99)))
  ev <- suppressMessages(sim_hawkesNet(truth, T_end = T_end, mark_type = "ba_bip",
                                       attachment = spec, kernel = "exp2",
                                       frailty_shape = shape, return_net = FALSE))$ev
  att <- as.numeric(table(ev$edges$j))
  fit <- suppressWarnings(fit_hawkesNet(
    ev, modifyList(truth, list(m = 0.3, beta_edges = 0.1, delta = 1)),
    fixed_params = time_names, mark_type = "ba_bip", T_end = T_end,
    attachment = spec, kernel = "exp2", hessian = TRUE
  ))
  se <- tryCatch(sqrt(diag(solve(fit$hessian)))[["theta"]], error = function(e) NA_real_)
  data.frame(frailty_shape = shape, rep = rep, n_events = nrow(ev$times),
             attendance_dispersion = var(att) / mean(att),
             theta = fit$par$theta, se = se, converged = fit$convergence == 0)
}

grid <- expand.grid(shape = frailty_shapes, rep = seq_len(n_rep))
res <- do.call(rbind, mclapply(seq_len(nrow(grid)), function(i) {
  capture.output(r <- one_rep(grid$shape[i], grid$rep[i]))
  r
}, mc.cores = cores))

summary_tab <- do.call(rbind, lapply(split(res, res$frailty_shape), function(d) {
  data.frame(frailty_variance = 1 / d$frailty_shape[1],
             attendance_dispersion = mean(d$attendance_dispersion),
             mean_theta = mean(d$theta), sd_theta = sd(d$theta), mean_se = mean(d$se, na.rm = TRUE),
             reject_5pct = mean(abs(d$theta / d$se) > 1.96, na.rm = TRUE),
             converged = mean(d$converged))
}))
print(summary_tab[order(summary_tab$frailty_variance), ], digits = 3, row.names = FALSE)
saveRDS(res, "ATTACH_theta_check.rds")
