# compare_fits.R

#' Compare fitted models
#'
#' Tabulates fits from [fit_hawkesNet()] side by side and, for models given in
#' order from smallest to largest, tests each against the one before it by a
#' likelihood ratio. The tests assume each model is nested in the next and
#' that all were fitted to the same events with the same `T0`, `T_end` and
#' fixed parameters; only the counts of free parameters are compared.
#'
#' @param ... Named `hawkesNet_fit` objects, smallest model first, e.g.
#'   `compare_fits(base = f0, size = f1, both = f2)`.
#' @param test Logical; add likelihood-ratio tests between successive models.
#'
#' @return A data.frame with one row per fit: `n_par` (free parameters),
#'   `loglik`, `AIC`, `BIC` (with the number of events as sample size),
#'   `converged`, `hessian_ok` (`NA` if fitted without `hessian = TRUE`), and,
#'   if `test`, `LR`, `df` and `p` against the previous row.
#' @export
compare_fits <- function(..., test = TRUE) {
  fits <- list(...)
  if (length(fits) == 0L) stop("compare_fits() needs at least one fit.", call. = FALSE)
  if (is.null(names(fits)) || any(!nzchar(names(fits)))) {
    names(fits) <- paste0("model", seq_along(fits))
  }
  if (!all(vapply(fits, inherits, logical(1), "hawkesNet_fit"))) {
    stop("All arguments must be fits from fit_hawkesNet().", call. = FALSE)
  }
  parts <- vapply(fits, function(f) f$model$part %||% "both", character(1))
  if (length(unique(parts)) > 1L) {
    stop("The fits are of different parts (", paste(unique(parts), collapse = ", "),
         "): their log-likelihoods cannot be compared.", call. = FALSE)
  }
  pseudo <- vapply(fits, .is_pseudo, logical(1))
  if (any(pseudo)) {
    stop("Pseudo-likelihood (4-cycle) fits cannot be compared by likelihood or AIC: ",
         paste(names(fits)[pseudo], collapse = ", "), ".", call. = FALSE)
  }
  n_ev <- vapply(fits, function(f) nrow(f$events$times), numeric(1))
  if (length(unique(n_ev)) > 1L) stop("The fits use different events.", call. = FALSE)

  n_par <- vapply(fits, function(f) if (is.null(f$fit)) 0 else length(f$fit$par), numeric(1))
  ll <- vapply(fits, function(f) f$loglik, numeric(1))
  hess_ok <- vapply(fits, function(f) {
    if (is.null(f$hessian)) return(NA)
    ev <- tryCatch(eigen(f$hessian, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NA_real_)
    all(is.finite(ev)) && all(ev > 0)
  }, logical(1))

  out <- data.frame(
    model = names(fits), n_par = n_par, loglik = ll,
    AIC = -2 * ll + 2 * n_par, BIC = -2 * ll + log(n_ev[1L]) * n_par,
    converged = vapply(fits, function(f) identical(f$convergence, 0L) || identical(f$convergence, 0),
                       logical(1)),
    hessian_ok = hess_ok, row.names = NULL
  )

  if (test && length(fits) > 1L) {
    df <- c(NA, diff(n_par))
    lr <- c(NA, 2 * diff(ll))
    if (any(df[-1L] <= 0)) {
      warning("Models are not in increasing order of free parameters; ",
              "tests against the previous row may be meaningless.", call. = FALSE)
    }
    out$LR <- lr
    out$df <- df
    out$p <- ifelse(df > 0, stats::pchisq(pmax(lr, 0), df, lower.tail = FALSE), NA)
  }
  out
}
