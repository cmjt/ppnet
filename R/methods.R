# Standard model methods for fits from fit_hawkesNet() (class "ppnet_fit").

# Internal: TRUE if the fit maximised a pseudo-likelihood (4-cycles).
.is_pseudo <- function(fit) isTRUE(fit$model$attachment$cycles)

.free_names <- function(fit) names(fit$fit$par)

# Internal: Hessian of the negative log-likelihood on the working scale, over
# the free parameters; the stored one, or computed from the fit's model.
.fit_hessian <- function(fit) {
  if (!is.null(fit$hessian)) return(fit$hessian)
  m <- fit$model
  ll_fn <- do.call(.build_loglik_fn, c(
    list(events = fit$events, T_end = m$T_end, T0 = m$T0, mark_type = m$mark_type,
         formula_rhs = m$formula_rhs, baseline = m$baseline, excitation = m$excitation,
         attachment = m$attachment, kernel = m$kernel),
    m$dots
  ))
  free <- .free_names(fit)
  x0_all <- to_working(fit$par, fit$transform)
  obj <- .make_objective(ll_fn, x0_all, free, fit$transform)
  has_grad <- !is.null(attr(tryCatch(ll_fn(fit$par, gradient = TRUE), error = function(e) NULL),
                            "gradient"))
  H <- stats::optimHess(x0_all[free], obj$fn, if (has_grad) obj$gr else NULL)
  dimnames(H) <- list(free, free)
  H
}

# Internal: d(natural)/d(working) for each free parameter.
.jacobian <- function(fit) {
  free <- .free_names(fit)
  est <- unlist(fit$par[free])
  stats::setNames(ifelse(fit$transform[free] == "log", est, 1), free)
}

#' Standard methods for fitted HawkesNet models
#'
#' Methods for fits from [fit_hawkesNet()], in the style of [stats::glm()]:
#' `coef()`, `vcov()`, `confint()`, `summary()`, `logLik()` (so [AIC()] and
#' [BIC()] work), `nobs()`, `anova()` and `simulate()`.
#'
#' Estimates and standard errors are on the natural scale of the parameters.
#' Positive parameters are fitted on the log scale, so their standard errors
#' come from the delta method and their confidence intervals are computed on
#' the log scale and back-transformed. Fixed parameters are not estimates, so
#' they are left out.
#'
#' Fits with 4-cycles (`attachment_spec(cycles = TRUE)`) maximise a
#' pseudo-likelihood: `vcov()` then uses the sandwich covariance from
#' [attachment_vcov()] for the attachment parameters, and `logLik()`,
#' `AIC()` and `anova()` give an error, as they are not valid.
#'
#' If the fit was made without `hessian = TRUE`, the Hessian is computed when
#' it is first needed (once per call).
#'
#' @param object,x A fit from [fit_hawkesNet()].
#' @param type For `vcov()`: `"auto"` (sandwich for pseudo-likelihood fits,
#'   otherwise the Hessian), `"hessian"` or `"sandwich"` (attachment
#'   parameters from [attachment_vcov()], the others from the Hessian; not
#'   available with the random effect).
#' @param scale For `vcov()`: `"natural"` (the parameters themselves) or
#'   `"working"` (the scale the fit was optimised on).
#' @param parm,level For `confint()`: parameters (names or indices; default
#'   all free ones) and confidence level.
#' @param nsim,seed For `simulate()`: number of simulations and an optional
#'   seed, as for [stats::simulate()].
#' @param T_end For `simulate()`: end of the simulation window (default: the
#'   fit's).
#' @param digits Number of significant digits to print.
#' @param ... For `anova()`: further fits to compare. For `simulate()`:
#'   passed to [sim_hawkesNet()]. Otherwise unused.
#'
#' @return `coef()`: named vector of estimates. `vcov()`: covariance matrix,
#'   with attribute `"type"`. `confint()`: matrix of lower and upper limits.
#'   `summary()`: an object of class `"summary.ppnet_fit"` with a coefficient
#'   table. `logLik()`: a `"logLik"` object. `anova()`: the table from
#'   [compare_fits()]. `simulate()`: a list of `nsim` events objects, with
#'   attribute `"seed"`.
#'
#' @examples
#' set.seed(1)
#' p <- list(mu = 1, K = 0.5, beta = 1, beta_edges = 0.3, lambda_new = 0.7)
#' ev <- sim_hawkesNet(p, T_end = 30, mark_type = "ba_bip")$ev
#' fit <- fit_hawkesNet(ev, p, mark_type = "ba_bip")
#' fit
#' summary(fit)
#' confint(fit)
#' AIC(fit)
#' sims <- simulate(fit, nsim = 2, seed = 2)
#' @name ppnet-methods
NULL

#' @rdname ppnet-methods
#' @export
coef.ppnet_fit <- function(object, ...) {
  unlist(object$par[.free_names(object)])
}

#' @rdname ppnet-methods
#' @export
vcov.ppnet_fit <- function(object, type = c("auto", "hessian", "sandwich"),
                           scale = c("natural", "working"), ...) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  free <- .free_names(object)
  if (type == "auto") type <- if (.is_pseudo(object)) "sandwich" else "hessian"

  H <- .fit_hessian(object)[free, free, drop = FALSE]
  V <- tryCatch(solve(H), error = function(e) {
    warning("The Hessian is singular; variances are not available.", call. = FALSE)
    matrix(NA_real_, length(free), length(free))
  })
  dimnames(V) <- list(free, free)

  if (type == "sandwich") {
    if (is.null(object$hessian)) object$hessian <- H
    Va <- attachment_vcov(object)
    a <- rownames(Va)
    V[a, ] <- 0
    V[, a] <- 0                    # time and mark parts share no parameters
    V[a, a] <- Va
  }

  if (scale == "natural") {
    J <- .jacobian(object)
    V <- V * outer(J, J)
  }
  attr(V, "type") <- type
  V
}

#' @rdname ppnet-methods
#' @export
confint.ppnet_fit <- function(object, parm, level = 0.95, ...) {
  free <- .free_names(object)
  if (missing(parm)) parm <- free
  if (is.numeric(parm)) parm <- free[parm]
  Vw <- vcov(object, scale = "working")
  x <- to_working(object$par[free], object$transform)
  z <- stats::qnorm(1 - (1 - level) / 2)
  se <- sqrt(diag(Vw))
  lo <- x - z * se
  hi <- x + z * se
  is_log <- object$transform[free] == "log"
  lo[is_log] <- exp(lo[is_log])
  hi[is_log] <- exp(hi[is_log])
  pct <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2), trim = TRUE), " %")
  out <- cbind(lo, hi)[parm, , drop = FALSE]
  dimnames(out) <- list(parm, pct)
  out
}

#' @rdname ppnet-methods
#' @export
logLik.ppnet_fit <- function(object, ...) {
  if (.is_pseudo(object)) {
    stop("This fit maximises a pseudo-likelihood (4-cycles), so logLik(), AIC() and ",
         "likelihood-ratio tests are not valid.", call. = FALSE)
  }
  structure(object$loglik, df = length(.free_names(object)), nobs = nobs(object),
            class = "logLik")
}

#' @rdname ppnet-methods
#' @export
nobs.ppnet_fit <- function(object, ...) nrow(object$events$times)

#' @rdname ppnet-methods
#' @export
summary.ppnet_fit <- function(object, ...) {
  free <- .free_names(object)
  est <- coef(object)
  V <- vcov(object)
  se <- sqrt(diag(V))
  # Wald tests against 0 only for parameters that can be 0 or negative.
  testable <- object$transform[free] == "none"
  z <- ifelse(testable, est / se, NA_real_)
  p <- ifelse(testable, 2 * stats::pnorm(-abs(z)), NA_real_)
  tab <- cbind(Estimate = est, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = p)
  rownames(tab) <- free
  fixed <- object$fixed %||% setdiff(names(object$par), free)
  structure(list(
    call = object$call, coefficients = tab, fixed = unlist(object$par[fixed]),
    pseudo = .is_pseudo(object), se_type = attr(V, "type"),
    loglik = object$loglik, n_par = length(free), nobs = nobs(object),
    convergence = object$convergence
  ), class = "summary.ppnet_fit")
}

#' @rdname ppnet-methods
#' @export
print.summary.ppnet_fit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients", if (x$pseudo) " (sandwich standard errors)", ":\n", sep = "")
  stats::printCoefmat(x$coefficients, digits = digits, na.print = "", has.Pvalue = TRUE,
                      P.values = TRUE)
  cat("Positive parameters are fitted on the log scale: their standard errors are",
      "from the delta method and they are not tested against 0.\n")
  if (length(x$fixed)) {
    cat("\nFixed: ", paste(names(x$fixed), format(x$fixed, digits = digits), sep = " = ",
                         collapse = ", "), "\n", sep = "")
  }
  cat("\n", x$nobs, " events, ", x$n_par, " free parameters. ", sep = "")
  if (x$pseudo) {
    cat("Log pseudo-likelihood: ", format(x$loglik, digits = digits + 2),
        " (AIC and likelihood-ratio tests are not valid).\n", sep = "")
  } else {
    cat("Log-likelihood: ", format(x$loglik, digits = digits + 2),
        ", AIC: ", format(-2 * x$loglik + 2 * x$n_par, digits = digits + 2), "\n", sep = "")
  }
  if (!identical(x$convergence, 0L) && !identical(x$convergence, 0)) {
    cat("Warning: optim() did not converge (code ", x$convergence, ").\n", sep = "")
  }
  cat("\n")
  invisible(x)
}

#' @rdname ppnet-methods
#' @export
print.ppnet_fit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n", nobs(x), " events. Log-", if (.is_pseudo(x)) "pseudo-" , "likelihood: ",
      format(x$loglik, digits = digits + 2), "\n\n", sep = "")
  invisible(x)
}

#' @rdname ppnet-methods
#' @export
anova.ppnet_fit <- function(object, ...) {
  fits <- c(list(object), list(...))
  nms <- vapply(as.list(substitute(list(object, ...)))[-1L],
                function(e) paste(deparse(e), collapse = ""), character(1))
  names(fits) <- nms
  do.call(compare_fits, fits)
}

#' @rdname ppnet-methods
#' @export
simulate.ppnet_fit <- function(object, nsim = 1, seed = NULL, T_end = object$model$T_end, ...) {
  # Seed handling as in stats:::simulate.lm.
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) stats::runif(1)
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  m <- object$model
  args <- c(list(params = object$par, T_end = T_end, T0 = m$T0, mark_type = m$mark_type,
                 baseline = m$baseline, excitation = m$excitation,
                 attachment = m$attachment, kernel = m$kernel, return_net = FALSE),
            if (!is.null(m$formula_rhs)) list(formula_rhs = m$formula_rhs),
            list(...))
  out <- lapply(seq_len(nsim), function(i) do.call(sim_hawkesNet, args)$ev)
  attr(out, "seed") <- RNGstate
  out
}
