# attachment_vcov.R

#' Sandwich covariance for the attachment parameters
#'
#' Robust (sandwich) covariance matrix \eqn{H^{-1} J H^{-1}} for the free
#' attachment parameters of a fit, on the working scale the fit used (log
#' for positive parameters, untransformed for `theta`, `theta_paths` and
#' `theta_cycles`). \eqn{H} is the Hessian of the fit and \eqn{J} the sum of
#' outer products of each event's score. Given the past, each event's score
#' has mean zero, so \eqn{J} estimates the variance of the total score.
#'
#' With `cycles = TRUE` the fit maximises a pseudo-likelihood and the
#' Hessian alone gives wrong standard errors, so use these. For a true
#' likelihood the two agree in large samples. The mark part shares no
#' parameters with the time part, so the time parameters are left out.
#'
#' @param fit A fit from [fit_hawkesNet()] with an `attachment` spec (the
#'   Hessian is computed if the fit has none). Not available with
#'   `random_effect = TRUE`. [vcov()] on the fit gives the same standard
#'   errors for pseudo-likelihood fits, together with the other parameters.
#'
#' @return The covariance matrix (working scale), with attributes `"se"`: the
#'   standard errors on the working scale, and `"se_natural"`: the standard
#'   errors of the parameters themselves (delta method, `estimate * se` for
#'   log-transformed parameters), comparable with the estimates in `fit$par`.
#' @export
attachment_vcov <- function(fit) {
  if (!inherits(fit, "hawkesNet_fit")) stop("`fit` must come from fit_hawkesNet().", call. = FALSE)
  spec <- fit$model$attachment
  if (is.null(spec)) stop("The fit has no attachment spec.", call. = FALSE)
  if (isTRUE(spec$random_effect)) {
    stop("attachment_vcov() is not available with the random effect.", call. = FALSE)
  }

  free <- intersect(names(fit$fit$par), attachment_param_names(spec))
  if (length(free) == 0L) stop("No free attachment parameters.", call. = FALSE)
  mark <- compile_logpmf_attach(fit$events, spec)
  ll <- mark(fit$par, scores = TRUE)
  S <- attr(ll, "scores")[, free, drop = FALSE]
  # Natural-scale scores to the working scale: d/d log p = p d/dp.
  tf <- fit$transform %||% stats::setNames(rep("log", length(free)), free)
  is_log <- vapply(free, function(nm) (tf[[nm]] %||% "log") == "log", logical(1))
  S[, is_log] <- sweep(S[, is_log, drop = FALSE], 2L, unlist(fit$par[free[is_log]]), `*`)

  Hinv <- solve(.fit_hessian(fit)[free, free, drop = FALSE])
  V <- Hinv %*% crossprod(S) %*% Hinv
  dimnames(V) <- list(free, free)
  se <- sqrt(diag(V))
  attr(V, "se") <- se
  attr(V, "se_natural") <- ifelse(is_log, unlist(fit$par[free]) * se, se)
  V
}
