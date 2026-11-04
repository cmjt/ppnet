# Simulation check: does an unmodelled perp random effect make the co-offender effect
# `theta` look positive when there is none?
#
# Each replicate simulates with theta = 0 and perp random effect u ~ Gamma(shape,
# shape) (variance 1 / shape), then fits the attachment model without a
# random effect and with one (`random_effect = TRUE`, M4), and records theta's
# estimate and Wald test under each. The time part is held at its true
# values: it shares no parameters with the marks, so fitting it would not
# change the theta estimates.
#
# Compare `attendance_dispersion` with the same quantity in your data,
#   a <- as.numeric(table(perp ids over event-perp edges)); var(a) / mean(a),
# to see which random-effect variance is realistic.
#
# Run from inside the ppnet repo, with the package installed.

library(ppnet)
library(parallel)
`%||%` <- function(x, y) if (is.null(x)) y else x

n_rep <- 20
random_effect_shapes <- c(Inf, 2, 1, 0.5)   # random effect variance 0, 0.5, 1, 2
T_end <- 520                          # weeks
cores <- max(1L, detectCores() - 1L)

spec <- attachment_spec(degree = TRUE, coffenders = TRUE, coffender_window = 52,
                        new_perps = "negbin", min_perps = 1)
spec_re <- attachment_spec(degree = TRUE, coffenders = TRUE, coffender_window = 52,
                           new_perps = "negbin", min_perps = 1, random_effect = TRUE)
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
                                       random_effect_shape = shape, return_net = FALSE))$ev
  att <- as.numeric(table(ev$edges$j))
  init <- modifyList(truth, list(m = 0.3, beta_edges = 0.1, delta = 1))
  fit_one <- function(s, init) {
    fit <- suppressWarnings(fit_hawkesNet(
      ev, init, fixed_params = time_names, mark_type = "ba_bip", T_end = T_end,
      attachment = s, kernel = "exp2", hessian = TRUE
    ))
    se <- tryCatch(sqrt(diag(solve(fit$hessian)))[["theta"]], error = function(e) NA_real_)
    data.frame(theta = fit$par$theta, se = se, random_effect_var = fit$par$random_effect_var %||% 0,
               loglik = fit$loglik, converged = fit$convergence == 0)
  }
  r0 <- fit_one(spec, init)
  r1 <- fit_one(spec_re, c(init, random_effect_var = 0.5))
  data.frame(random_effect_shape = shape, rep = rep, n_events = nrow(ev$times),
             attendance_dispersion = var(att) / mean(att),
             model = c("no_re", "re"), rbind(r0, r1))
}

grid <- expand.grid(shape = random_effect_shapes, rep = seq_len(n_rep))
res <- do.call(rbind, mclapply(seq_len(nrow(grid)), function(i) {
  capture.output(r <- one_rep(grid$shape[i], grid$rep[i]))
  r
}, mc.cores = cores))

summary_tab <- do.call(rbind, lapply(split(res, list(res$model, res$random_effect_shape)), function(d) {
  data.frame(model = d$model[1], random_effect_var = 1 / d$random_effect_shape[1],
             attendance_dispersion = mean(d$attendance_dispersion),
             mean_theta = mean(d$theta), sd_theta = sd(d$theta), mean_se = mean(d$se, na.rm = TRUE),
             reject_5pct = mean(abs(d$theta / d$se) > 1.96, na.rm = TRUE),
             mean_random_effect_var = mean(d$random_effect_var), converged = mean(d$converged))
}))
print(summary_tab[order(summary_tab$model, summary_tab$random_effect_var), ], digits = 3,
      row.names = FALSE)
saveRDS(res, "ATTACH_theta_check.rds")
