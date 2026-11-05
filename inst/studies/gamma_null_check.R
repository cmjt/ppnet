# Simulation check: does perp propensity make the event-level connectedness
# effect `gamma` look positive when the network has no effect?
#
# Events come from perps offending (sim_perp_events()): each perp has its own
# rate, with a gamma random effect and optional self-excitation, over a finite career,
# and brings co-offenders along, partly past ones. No rate depends on the
# network, so the true gamma is 0, but prolific perps make events that are
# both well connected and followed by more events.
#
# Each replicate fits the event-time model with size and connectedness
# excitation and records gamma's estimate and Wald test. The mark part shares
# no parameters with the time part, so it is held fixed.
#
# Compare `attendance_dispersion` with the same quantity in your data,
#   a <- as.numeric(table(perp ids over event-perp edges)); var(a) / mean(a),
# to see which scenario is realistic.
#
# Run with the package installed, e.g.
#   source(system.file("studies", "gamma_null_check.R", package = "ppnet"))

library(ppnet)
library(parallel)

# Where to save the results (set `out_file` before sourcing to change it).
if (!exists("out_file")) out_file <- "gamma_null_check.rds"

n_rep <- 20
T_end <- 520                          # weeks
window <- 52
cores <- max(1L, detectCores() - 1L)

# Generic values; replace with ones matching your data. Don't commit those.
base <- list(entry_rate = 1, career_mean = 104, rate = 0.05, random_effect_shape = Inf,
             K_self = 0, beta_self = 0.5, excite_on = "own", extra_mean = 0.5,
             size_link = 0, p_past = 0.5, recruit_link = 1)
scenarios <- list(
  no_re = list(),
  re_var_1 = list(random_effect_shape = 1),
  re_var_2 = list(random_effect_shape = 0.5),
  re_var_1_self = list(random_effect_shape = 1, K_self = 0.3),
  re_var_1_self_all = list(random_effect_shape = 1, K_self = 0.3, excite_on = "all"),
  re_var_1_size = list(random_effect_shape = 1, size_link = 1)
)

spec <- attachment_spec(degree = TRUE, new_perps = "negbin", min_perps = 1)
exc <- excitation_spec(size = TRUE, connectedness = TRUE, connectedness_window = window)
marks <- list(m = 0.5, beta_edges = 0.05, delta = 0.5, lambda_new = 0.6,
              lambda_new0 = 0.5, size_new = 2)

one_rep <- function(scen, rep) {
  set.seed(1000 * rep + match(scen, names(scenarios)))
  args <- modifyList(base, scenarios[[scen]])
  ev <- do.call(sim_perp_events, c(list(T_end = T_end), args))$ev
  att <- as.numeric(table(ev$edges$j))
  n <- nrow(ev$times)
  init <- c(list(mu = n / T_end / 2, K = 0.3, beta = 0.5, alpha = 0, gamma = 0), marks)
  fit <- suppressWarnings(fit_hawkesNet(
    ev, init, fixed_params = names(marks), mark_type = "ba_bip", T_end = T_end,
    attachment = spec, excitation = exc, hessian = TRUE
  ))
  se <- tryCatch(sqrt(diag(solve(fit$hessian)))[c("alpha", "gamma")],
                 error = function(e) c(alpha = NA_real_, gamma = NA_real_))
  data.frame(scenario = scen, rep = rep, n_events = n,
             attendance_dispersion = var(att) / mean(att),
             K = fit$par$K, beta = fit$par$beta,
             alpha = fit$par$alpha, se_alpha = se[["alpha"]],
             gamma = fit$par$gamma, se_gamma = se[["gamma"]],
             converged = fit$convergence == 0)
}

grid <- expand.grid(scen = names(scenarios), rep = seq_len(n_rep), stringsAsFactors = FALSE)
res <- do.call(rbind, mclapply(seq_len(nrow(grid)), function(i) {
  capture.output(r <- one_rep(grid$scen[i], grid$rep[i]))
  r
}, mc.cores = cores))

summary_tab <- do.call(rbind, lapply(split(res, res$scenario), function(d) {
  z <- d$gamma / d$se_gamma
  data.frame(scenario = d$scenario[1], n_events = mean(d$n_events),
             attendance_dispersion = mean(d$attendance_dispersion),
             mean_gamma = mean(d$gamma), sd_gamma = sd(d$gamma),
             mean_se = mean(d$se_gamma, na.rm = TRUE),
             pos_signif = mean(z > 1.96, na.rm = TRUE),
             neg_signif = mean(z < -1.96, na.rm = TRUE),
             mean_alpha = mean(d$alpha), converged = mean(d$converged))
}))
print(summary_tab[match(names(scenarios), summary_tab$scenario), ], digits = 3,
      row.names = FALSE)
saveRDS(res, out_file)
