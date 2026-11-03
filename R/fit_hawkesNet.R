#' Fit a HawkesNet model by maximum likelihood
#'
#' @param events An events object created by `make_events()`.
#' @param params_init Named list of starting values for all model parameters.
#'   Must include the parameters required by the chosen `mark_type` (and any
#'   baseline/intensity parameters used by `loglik()`).
#' @param fixed_params Character vector of parameter names to hold fixed at
#'   their `params_init` values.
#' @param transform Named character vector giving the working-scale transform
#'   for each parameter in `params_init` (e.g. `"log"`). Missing entries default
#'   to `"log"`. If `NULL`, all parameters default to `"log"`.
#' @param hessian Logical; if `TRUE`, return the Hessian from `optim()`.
#' @param control List of control arguments passed to [stats::optim()].
#' @param mark_type One of `"ba"`, `"cs"`, `"ba_bip"`, `"cs_bip"`.
#' @param debug Logical; enable debug mode (passed through to likelihood/mark
#'   components where supported).
#' @param T_end Optional end time for the likelihood window.
#' @param T0 Start time for the likelihood window.
#' @param ... Additional arguments forwarded to the mark PMF / edge-probability
#'   functions used inside `loglik()`.
#' @param formula_rhs String, the CS kernel formula.
#' @param baseline Optional [baseline_spec()] for a time-varying baseline. Its
#'   coefficients (see [baseline_coef_names()]) must be in `params_init`, and
#'   default to the `"none"` transform.
#' @param excitation Optional [excitation_spec()]; its coefficients (`alpha`,
#'   `gamma`) must be in `params_init` and default to the `"none"` transform.
#' @param attachment Optional [attachment_spec()]; `theta` defaults to the
#'   `"none"` transform, the other attachment parameters to `"log"`.
#' @param kernel `"exp"` or `"exp2"`; see [loglik()]. With `"exp2"`,
#'   `params_init` must contain `K_fast` and `beta_fast`, with `beta_fast > beta`.
#'
#' @return A list with components:
#' \describe{
#'   \item{fit}{The `optim()` result (or `NULL` if all parameters are fixed).}
#'   \item{par}{Named list of parameter estimates on the natural scale.}
#'   \item{loglik}{Log-likelihood at `par`.}
#'   \item{convergence}{`optim()` convergence code.}
#'   \item{message}{`optim()` message (if any).}
#'   \item{hessian}{Hessian matrix (only if `hessian = TRUE`).}
#' }
#'
#' @examples
#' set.seed(1)
#'
#' params_true <- list(
#'   mu = 0.5,
#'   K = 0.5,
#'   beta = 0.5,
#'   beta_edges = 2
#' )
#'
#' sim <- sim_hawkesNet(
#'   params = params_true,
#'   T_end = 10,
#'   mark_type = "ba"
#' )
#'
#' params_init <- list(
#'   mu = 1,
#'   K = 1,
#'   beta = 1,
#'   beta_edges = 1
#' )
#'
#' fit <- fit_hawkesNet(
#'   events = sim$ev,
#'   params_init = params_init,
#'   mark_type = "ba"
#' )
#'
#' fit$par
#'
#' @export
fit_hawkesNet <- function(events,
                          params_init,
                          fixed_params = NULL,
                          transform = NULL,
                          hessian = FALSE,
                          control = list(),
                          mark_type = "ba",
                          debug = FALSE,
                          T_end = NULL,
                          T0 = 0,
                          formula_rhs = NULL,
                          baseline = NULL,
                          excitation = NULL,
                          attachment = NULL,
                          kernel = c("exp", "exp2"),
                          ...
                          ) {
  # 0) Input validation

  # Make sure starting values are provided for the optimizer
  stopifnot(is.list(params_init), length(params_init) > 0)

  start_time <- proc.time()

  # Check up front: inside the objective a missing coefficient would be caught
  # and look like an infeasible point, so the fit would fail silently.
  if (!is.null(baseline)) .baseline_coefs(baseline, params_init)
  if (!is.null(excitation)) .excitation_coefs(excitation, params_init)
  kernel <- match.arg(kernel)
  if (kernel == "exp2") {
    if (!all(c("K_fast", "beta_fast") %in% names(params_init))) {
      stop("kernel = 'exp2' needs K_fast and beta_fast in params_init.", call. = FALSE)
    }
    if (!(params_init$beta_fast > params_init$beta)) {
      stop("kernel = 'exp2' needs beta_fast > beta in params_init.", call. = FALSE)
    }
  }
  missing_attach <- setdiff(attachment_param_names(attachment), names(params_init))
  if (length(missing_attach) > 0L) {
    stop("Missing attachment parameters: ", paste(missing_attach, collapse = ", "), call. = FALSE)
  }

  # 1) Deal with the parameter transformations

  # If no transformation information was provided, assume we are fitting on the
  # log scale
  # Baseline and excitation coefficients can be negative, so they are fitted
  # untransformed unless the user says otherwise.
  b_names <- intersect(c(baseline_coef_names(baseline), excitation_coef_names(excitation),
                         if (!is.null(attachment) && attachment$coffenders) "theta"),
                       names(params_init))
  if (is.null(transform)) {
    transform <- rep("log", length(params_init))
    names(transform) <- names(params_init)
    transform[b_names] <- "none"
  }
  transform[setdiff(b_names, names(transform))] <- "none"
  # Make sure that each parameter has a valid transform, default to "log"
  # and warn user in the case any are missing
  missing_tf <- setdiff(names(params_init), names(transform))
  if (length(missing_tf) > 0) {
    warning(paste0('Missing `transform` entries for: ',
                  paste(missing_tf, collapse = ", "), '. Defaulting to "log"'))

    transform[missing_tf] <- "log"
  }

  # 2) Deal with any possible fixed parameters

  # Coerce fixed_params to character(0) if NULL
  fixed_params <- fixed_params %||% character(0)
  # Make sure to warn the user if fixed_params contains names not present in
  # params_init
  unknown_fixed <- setdiff(fixed_params, names(params_init))
  if (length(unknown_fixed) > 0) {
    warning(
      "The following fixed_params are not present in params_init and will be ignored: ",
      paste(unknown_fixed, collapse = ", "))
  }
  # Derive which parameters are free using fixed_params
  fixed_params <- intersect(fixed_params, names(params_init))
  free_params <- setdiff(names(params_init), fixed_params)


  # Build the compiled likelihood function ONCE (optim() will call the objective
  # many times; we don't want to rebuild mark caches each time).
  ll_fn <- .build_loglik_fn(
    events = events,
    T_end = T_end,
    T0 = T0,
    mark_type = mark_type,
    debug = debug,
    formula_rhs = formula_rhs,
    baseline = baseline,
    excitation = excitation,
    attachment = attachment,
    kernel = kernel,
    ...
  )

  # And in the case there are no free parameters, just return the likelihood
  # at the param_init values
  if (length(free_params) == 0L) {
    ll <- ll_fn(params_init)
    return(list(
      fit = NULL,
      par = params_init,
      loglik = ll,
      convergence = 0L,
      message = "No free parameters (all fixed)."
    ))
  }

  # 4) Optimizer setup

  # Transform the initial parameter values to working scale
  x0_all <- to_working(params_init, transform)

  # Define the objective function
  fn <- function(x) {
    # optim() only passes the free parameters, so fill in the fixed ones
    # before transforming back from working scale
    x_all <- x0_all
    x_all[free_params] <- x
    p <- from_working(x_all, transform)
    ll <- tryCatch(
      ll_fn(p),
      error = function(e) -Inf
    )
    if (!is.finite(ll)) return(1e100)
    -ll
  }

  # Make sure to only optimize the free parameters
  x0 <- x0_all[free_params]

  # Analytic gradient on the working scale, when the model provides one.
  # d/dx = p * d/dp for log-transformed parameters. Parameters the likelihood
  # does not use get gradient zero.
  gr <- function(x) {
    x_all <- x0_all
    x_all[free_params] <- x
    p <- from_working(x_all, transform)
    ll <- tryCatch(ll_fn(p, gradient = TRUE), error = function(e) -Inf)
    g <- attr(ll, "gradient")
    if (!is.finite(ll) || is.null(g)) return(rep(0, length(x)))
    g_free <- vapply(free_params, function(nm) if (nm %in% names(g)) g[[nm]] else 0, numeric(1))
    is_log <- vapply(free_params, function(nm) (transform[[nm]] %||% "log") == "log", logical(1))
    g_free[is_log] <- g_free[is_log] * unlist(p[free_params[is_log]])
    -g_free
  }
  p_init <- from_working(x0_all, transform)
  has_grad <- !is.null(attr(tryCatch(ll_fn(p_init, gradient = TRUE), error = function(e) NULL),
                            "gradient"))

  # 5) Fit the model
  fit <- stats::optim(
    par = x0,
    fn = fn,
    gr = if (has_grad) gr else NULL,
    method = "BFGS",
    control = control,
    hessian = hessian
  )

  # 6) Structure output
  x_all <- x0_all
  x_all[free_params] <- fit$par

  p_hat <- from_working(x_all, transform)
  ll_hat <- ll_fn(p_hat)

  # Prepare output object
  out <- list(
    fit = fit,
    par = p_hat,
    loglik = ll_hat,
    convergence = fit$convergence,
    message = fit$message,
    events = events
  )
  out$call = match.call()
  if (hessian) out$hessian <- fit$hessian
  class(out) <- c("hawkesNet_fit", class(out))

  print(paste0("Fitting took ", round((proc.time()-start_time)[3],2), " seconds"))

  out
}

