p_m <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)

set.seed(201)
ev_m <- sim_hawkesNet(p_m, T_end = 40, mark_type = "ba_bip")$ev
fit_m <- fit_hawkesNet(ev_m, p_m, mark_type = "ba_bip", hessian = TRUE)

test_that("fits and simulations are quiet by default", {
  expect_silent(fit_hawkesNet(ev_m, p_m, mark_type = "ba_bip"))
  expect_silent(sim_hawkesNet(p_m, T_end = 5, mark_type = "ba_bip"))
  expect_message(fit_hawkesNet(ev_m, p_m, mark_type = "ba_bip", verbose = TRUE), "Fitting took")
})

test_that("coef, vcov and confint are on the natural scale", {
  expect_s3_class(fit_m, "ppnet_fit")
  expect_equal(coef(fit_m), unlist(fit_m$par))

  Vw <- vcov(fit_m, scale = "working")
  expect_equal(unname(Vw), unname(solve(fit_m$hessian)), ignore_attr = TRUE)
  J <- unlist(fit_m$par)                          # all log-transformed here
  expect_equal(unname(vcov(fit_m)), unname(Vw * outer(J, J)), ignore_attr = TRUE)
  expect_identical(attr(vcov(fit_m), "type"), "hessian")

  ci <- confint(fit_m)
  expect_equal(rownames(ci), names(coef(fit_m)))
  expect_true(all(ci[, 1] > 0 & ci[, 1] < coef(fit_m) & coef(fit_m) < ci[, 2]))
  # Log scale, back-transformed: the interval is symmetric in log.
  se_w <- sqrt(diag(Vw))
  expect_equal(unname(log(ci[, 2]) - log(coef(fit_m))), unname(qnorm(0.975) * se_w))
  expect_equal(confint(fit_m, "mu"), ci["mu", , drop = FALSE])
})

test_that("the Hessian is computed when the fit has none", {
  fit0 <- fit_hawkesNet(ev_m, p_m, mark_type = "ba_bip")
  expect_null(fit0$hessian)
  expect_equal(vcov(fit0), vcov(fit_m), tolerance = 1e-4)
})

test_that("logLik, AIC, BIC, nobs and anova agree with compare_fits", {
  ll <- logLik(fit_m)
  expect_equal(as.numeric(ll), fit_m$loglik)
  expect_equal(attr(ll, "df"), 5L)
  expect_equal(nobs(fit_m), nrow(ev_m$times))
  fit_small <- fit_hawkesNet(ev_m, p_m, fixed_params = "beta", mark_type = "ba_bip")
  tab <- compare_fits(small = fit_small, full = fit_m)
  expect_equal(AIC(fit_small, fit_m)$AIC, tab$AIC)
  expect_equal(BIC(fit_m), tab$BIC[2])
  a <- anova(fit_small, fit_m)
  expect_equal(a$model, c("fit_small", "fit_m"))
  expect_equal(a$LR, tab$LR)
})

test_that("summary tests only the parameters that can be 0 and lists fixed ones", {
  p_c <- c(p_m, list(gamma = 0))
  set.seed(202)
  fit_c <- fit_hawkesNet(ev_m, p_c, fixed_params = "beta_edges", mark_type = "ba_bip",
                         excitation = excitation_spec(size = FALSE, connectedness = TRUE),
                         hessian = TRUE)
  s <- summary(fit_c)
  expect_s3_class(s, "summary.ppnet_fit")
  tab <- s$coefficients
  expect_equal(rownames(tab), c("mu", "K", "beta", "gamma", "lambda_new"))
  expect_true(all(is.na(tab[c("mu", "K", "beta", "lambda_new"), "z value"])))
  expect_equal(tab["gamma", "z value"], tab["gamma", "Estimate"] / tab["gamma", "Std. Error"])
  expect_equal(names(s$fixed), "beta_edges")
  out <- capture.output(print(s))
  expect_true(any(grepl("Fixed: beta_edges", out)))
  expect_true(any(grepl("AIC", out)))
  expect_true(any(grepl("Coefficients", capture.output(print(fit_c)))))
})

test_that("pseudo-likelihood fits use sandwich SEs and refuse likelihood comparisons", {
  p <- list(mu = 3, K = 0.4, beta = 1.5, m = 1.2, beta_edges = 0.3, delta = 1,
            lambda_new = 1, theta_cycles = 0.6)
  sp <- attachment_spec(degree = TRUE, cycles = TRUE, coffender_window = 10)
  set.seed(203)
  ev <- sim_hawkesNet(p, T_end = 100, mark_type = "ba_bip", attachment = sp, return_net = FALSE)$ev
  fit <- suppressWarnings(fit_hawkesNet(ev, p, fixed_params = c("mu", "K", "beta"),
                                        mark_type = "ba_bip", attachment = sp))
  V <- vcov(fit)
  expect_identical(attr(V, "type"), "sandwich")
  expect_equal(sqrt(diag(V)), attr(attachment_vcov(fit), "se_natural"))
  expect_error(logLik(fit), "pseudo-likelihood")
  expect_error(AIC(fit), "pseudo-likelihood")
  expect_error(compare_fits(a = fit, b = fit), "Pseudo-likelihood")
  expect_true(any(grepl("sandwich", capture.output(summary(fit)))))
})

test_that("simulate() uses the fitted model and handles the seed like stats::simulate", {
  s1 <- simulate(fit_m, nsim = 2, seed = 5)
  s2 <- simulate(fit_m, nsim = 2, seed = 5)
  expect_length(s1, 2L)
  expect_true(all(vapply(s1, inherits, logical(1), "events")))
  expect_equal(s1[[1]]$times, s2[[1]]$times)
  expect_false(identical(s1[[1]]$times, s1[[2]]$times))
  # The caller's random stream is left as it was.
  set.seed(9); a <- runif(1)
  set.seed(9); invisible(simulate(fit_m, seed = 1)); b <- runif(1)
  expect_equal(a, b)
  short <- simulate(fit_m, T_end = 5, seed = 1)[[1]]
  expect_lte(max(short$times$t), 5)
})

test_that("the time and mark parts add up to the full log-likelihood", {
  ll <- function(part) loglik(ev_m, p_m, mark_type = "ba_bip", part = part)
  expect_equal(ll("time") + ll("marks"), ll("both"))
  # Each part needs only its own parameters.
  expect_equal(loglik(ev_m, p_m[c("mu", "K", "beta")], mark_type = "ba_bip", part = "time"),
               ll("time"))
  expect_equal(loglik(ev_m, p_m[c("beta_edges", "lambda_new")], mark_type = "ba_bip",
                      part = "marks"), ll("marks"))
})

test_that("one-part fits match the joint fit and count only their own parameters", {
  ft <- fit_hawkesNet(ev_m, mark_type = "ba_bip", part = "time", hessian = TRUE)
  fm <- fit_hawkesNet(ev_m, mark_type = "ba_bip", part = "marks", hessian = TRUE)
  expect_equal(names(coef(ft)), c("mu", "K", "beta"))
  expect_equal(names(coef(fm)), c("beta_edges", "lambda_new"))
  expect_equal(c(coef(ft), coef(fm))[names(coef(fit_m))], coef(fit_m), tolerance = 1e-3)
  expect_equal(ft$loglik + fm$loglik, fit_m$loglik, tolerance = 1e-6)
  expect_equal(attr(logLik(ft), "df"), 3L)
  expect_equal(sqrt(diag(vcov(ft))), sqrt(diag(vcov(fit_m)))[1:3], tolerance = 1e-2)
  expect_true(any(grepl("time part only", capture.output(summary(ft)))))
  # Entries for the other part are ignored.
  ft2 <- fit_hawkesNet(ev_m, c(p_m, list(mu = 2)), mark_type = "ba_bip", part = "time")
  expect_equal(names(coef(ft2)), c("mu", "K", "beta"))
  expect_error(compare_fits(time = ft, marks = fm), "different parts")
  # simulate() needs the other part's parameters.
  expect_error(simulate(fm, seed = 1), "give the time parameters")
  sims <- simulate(fm, seed = 1, params = coef(ft))
  expect_s3_class(sims[[1]], "events")
})

test_that("default starting values reach the same fit as given ones", {
  f0 <- fit_hawkesNet(ev_m, mark_type = "ba_bip")
  expect_equal(coef(f0), coef(fit_m), tolerance = 1e-3)
  expect_equal(f0$loglik, fit_m$loglik, tolerance = 1e-6)
})

test_that("defaults cover the baseline, excitation, exp2 and attachment terms", {
  sp <- attachment_spec(coffenders = TRUE, paths = TRUE, coffender_window = 10, episode_gap = 0.2,
                        recent = TRUE, new_perps = "negbin", min_perps = 1, random_effect = TRUE)
  exc <- excitation_spec(size = FALSE, connectedness = TRUE, connectedness_window = 10)
  bl <- baseline_spec(periods = c(year = 20))
  st <- .default_start(ev_m, "both", "ba_bip", "exp2", bl, exc, sp)
  expect_setequal(names(st), c(.time_param_names("exp2", bl, exc), attachment_param_names(sp)))
  expect_true(st$beta_fast > st$beta)
  expect_true(all(unlist(st[c("theta", "theta_paths", "theta_recent", "gamma",
                              "b_year_cos1", "b_year_sin1")]) == 0))
  pos <- setdiff(names(st), c("theta", "theta_paths", "theta_recent", "gamma",
                              "b_year_cos1", "b_year_sin1"))
  expect_true(all(unlist(st[pos]) > 0))
  # And a time fit with them runs.
  f <- fit_hawkesNet(ev_m, mark_type = "ba_bip", part = "time", kernel = "exp2",
                     baseline = bl, excitation = exc)
  expect_equal(f$convergence, 0L)
})
