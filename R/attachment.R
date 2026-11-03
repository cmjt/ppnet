# attachment.R

#' Attachment specification for returning perps
#'
#' Replaces the BA-bip rule for which old perps join an event. At event
#' \eqn{k} each earlier perp \eqn{v} has weight
#' \deqn{w_v = (d_v + \delta)\, e^{\theta x_v}\, e^{-\beta_e (t_k - \ell_v)}}
#' where \eqn{d_v} is the number of events \eqn{v} has attended, \eqn{\ell_v}
#' the time of the latest, and \eqn{x_v} its co-offender feature, \eqn{h} of
#' the number of distinct perps it has shared an event with. It joins
#' independently with probability
#' \deqn{p_v = 1 - \exp(-m\, w_v / W), \qquad W = \textstyle\sum_u w_u,}
#' so `m` is close to the expected number of returning perps per event when
#' the \eqn{p_v} are small.
#'
#' The number of new perps is modelled given how many old perps joined, so
#' the mark likelihood factorises as \eqn{p(\text{old})\, p(\text{new} \mid
#' \text{old})}. With `min_perps = 0` it is `new_perps`(`lambda_new`) at every
#' event. With `min_perps = 1` every event has at least one perp: events with
#' an old perp have `new_perps`(`lambda_new`) new ones, and events with none
#' have \eqn{1 +} `new_perps`(`lambda_new0`). `"negbin"` uses mean
#' `lambda_new` and size `size_new`, so the variance is
#' \eqn{\lambda + \lambda^2 / \text{size}}.
#'
#' Parameters: `m`, `beta_edges` (\eqn{\beta_e}), `delta` (if `degree`) and
#' `theta` (if `coffenders`), then `lambda_new`, `lambda_new0` (if
#' `min_perps = 1`) and `size_new` (if `"negbin"`); see
#' [attachment_param_names()]. With `degree = FALSE` and `coffenders = FALSE`
#' the weights depend on recency only.
#'
#' Since \eqn{\sum_v m w_v / W = m}, each evaluation needs only the running
#' sum \eqn{W}, so its cost grows with the number of events and perp
#' appearances, not with events times perps.
#'
#' @param degree Logical; include \eqn{(d_v + \delta)}.
#' @param coffenders Logical; include \eqn{e^{\theta x_v}}.
#' @param coffender_scale `"log1p"` (default) for \eqn{h(n) = \log(1 + n)} or
#'   `"linear"` for \eqn{h(n) = n}.
#' @param coffender_window Only co-offenders sharing an event at most this long
#'   before count. `Inf` (the default) counts all.
#' @param new_perps Distribution of the new-perp count: `"poisson"` (default)
#'   or `"negbin"`.
#' @param min_perps `0` (default) allows events with no perps; `1` requires at
#'   least one, as in observed data, and adds `lambda_new0`.
#'
#' @return An object of class `attachment_spec`.
#' @export
attachment_spec <- function(degree = TRUE, coffenders = FALSE,
                            coffender_scale = c("log1p", "linear"),
                            coffender_window = Inf,
                            new_perps = c("poisson", "negbin"), min_perps = 0) {
  coffender_scale <- match.arg(coffender_scale)
  new_perps <- match.arg(new_perps)
  if (!(length(min_perps) == 1L && min_perps %in% c(0, 1))) {
    stop("`min_perps` must be 0 or 1.", call. = FALSE)
  }
  stopifnot(is.numeric(coffender_window), length(coffender_window) == 1L, coffender_window >= 0)
  structure(
    list(degree = isTRUE(degree), coffenders = isTRUE(coffenders),
         coffender_scale = coffender_scale, coffender_window = coffender_window,
         new_perps = new_perps, min_perps = as.integer(min_perps)),
    class = "attachment_spec"
  )
}

#' Names of the attachment parameters
#'
#' @param spec An `attachment_spec`, or `NULL`.
#' @return Character vector of parameter names.
#' @export
attachment_param_names <- function(spec) {
  if (is.null(spec)) return(character(0))
  c("m", "beta_edges", if (spec$degree) "delta", if (spec$coffenders) "theta", "lambda_new",
    if (spec$min_perps == 1L) "lambda_new0", if (spec$new_perps == "negbin") "size_new")
}

# Internal: log-likelihood of the new-perp counts given the old-perp counts,
# with its gradient (attribute "gradient") when asked.
.new_perp_count_ll <- function(spec, n_new, n_old, params, gradient = FALSE) {
  shifted <- spec$min_perps == 1L & n_old == 0L
  y <- n_new - shifted
  if (any(y < 0L)) return(-Inf)  # an event with no perps at all

  groups <- list(lambda_new = !shifted)
  if (spec$min_perps == 1L) groups$lambda_new0 <- shifted
  size <- if (spec$new_perps == "negbin") params[["size_new"]] else Inf

  ll <- 0
  g <- list()
  for (nm in names(groups)) {
    yy <- y[groups[[nm]]]
    mu <- params[[nm]]
    if (spec$new_perps == "poisson") {
      ll <- ll + sum(stats::dpois(yy, mu, log = TRUE))
      g[[nm]] <- sum(yy) / mu - length(yy)
    } else {
      ll <- ll + sum(stats::dnbinom(yy, size = size, mu = mu, log = TRUE))
      g[[nm]] <- sum(yy / mu - (yy + size) / (size + mu))
      g[["size_new"]] <- (g[["size_new"]] %||% 0) +
        sum(digamma(yy + size) - digamma(size) + log(size / (size + mu)) +
              1 - (yy + size) / (size + mu))
    }
  }
  if (gradient && is.finite(ll)) attr(ll, "gradient") <- unlist(g)
  ll
}

.coffender_feature <- function(spec, n) {
  if (spec$coffender_scale == "log1p") log1p(n) else as.numeric(n)
}

# Attachment mark evaluator for a bipartite history.
compile_logpmf_attach <- function(events, spec) {
  stopifnot(inherits(spec, "attachment_spec"))
  a <- .event_arrays(events)
  if (all(is.na(a$role))) stop("Attachment needs node roles 'event' and 'perp'.", call. = FALSE)
  f <- .bip_focal(a)
  invalid <- !all(f$valid)
  h <- if (invalid) NULL else attach_history_cpp(a$t, a$node_ptr, a$edge_ptr, a$edge_i, a$edge_j,
                                                 f$focal, length(a$ids), spec$coffender_window)
  feat <- if (spec$coffenders && !invalid) .coffender_feature(spec, h$upd_x) else numeric(length(h$upd_x))
  n_old <- if (invalid) integer(0) else diff(h$obs_ptr)

  # With `gradient = TRUE` the value carries attribute "gradient".
  function(params, gradient = FALSE) {
    get <- function(nm, default) if (nm %in% names(params)) params[[nm]] else default
    m <- get("m", NA_real_)
    beta_edges <- get("beta_edges", NA_real_)
    delta <- if (spec$degree) get("delta", NA_real_) else 0
    theta <- if (spec$coffenders) get("theta", NA_real_) else 0
    lambda_new <- get("lambda_new", NA_real_)
    lambda_new0 <- if (spec$min_perps == 1L) get("lambda_new0", NA_real_) else 1
    size_new <- if (spec$new_perps == "negbin") get("size_new", NA_real_) else 1

    ok <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
    if (!ok(m) || m <= 0 || !ok(beta_edges) || beta_edges < 0 || !ok(lambda_new) || lambda_new < 0 ||
        !ok(delta) || delta <= 0 && spec$degree || !ok(theta) ||
        !ok(lambda_new0) || lambda_new0 < 0 || !ok(size_new) || size_new <= 0) return(-Inf)
    if (invalid) return(-Inf)

    r <- attach_ll_cpp(a$t, h$upd_ptr, h$upd_node, h$upd_deg, h$upd_last, feat,
                       h$obs_ptr, h$obs_node, length(a$ids), m, beta_edges, delta, theta,
                       spec$degree)
    counts <- .new_perp_count_ll(spec, f$n_new, n_old, params, gradient = gradient)
    ll <- as.numeric(counts) + r[1L]
    if (gradient && is.finite(ll)) {
      g <- c(m = r[2L], beta_edges = r[3L], delta = r[4L], theta = r[5L],
             attr(counts, "gradient"))
      attr(ll, "gradient") <- g[intersect(names(g), attachment_param_names(spec))]
    }
    ll
  }
}

# Internal: a stateful mark simulator for sim_hawkesNet(). It is only called
# for accepted events, so it updates its own perp state on every call.
# With finite `frailty_shape` each perp gets a fixed propensity
# u ~ Gamma(shape, shape) (mean 1) multiplying its weight; the likelihood
# ignores it, so this simulates unmodelled differences between perps.
.make_sim_mark_attach <- function(spec, frailty_shape = Inf) {
  stopifnot(is.numeric(frailty_shape), length(frailty_shape) == 1L, frailty_shape > 0)
  st <- new.env()
  st$ids <- character(0)       # perps, in arrival order
  st$deg <- numeric(0)
  st$last <- numeric(0)
  st$u <- numeric(0)
  st$ties <- list()            # perp -> named vector: co-offender -> latest shared time
  st$n_nodes <- 0L

  function(net, t_k, params, ...) {
    n_old <- length(st$ids)
    joined <- character(0)
    if (n_old > 0L) {
      x <- vapply(st$ids, function(v) {
        tv <- st$ties[[v]]
        sum(t_k - tv <= spec$coffender_window)
      }, numeric(1))
      w <- st$u * exp(-params$beta_edges * (t_k - st$last))
      if (spec$degree) w <- w * (st$deg + params$delta)
      if (spec$coffenders) w <- w * exp(params$theta * .coffender_feature(spec, x))
      p <- -expm1(-params$m * w / sum(w))
      joined <- st$ids[stats::runif(n_old) < p]
    }

    shifted <- spec$min_perps == 1L && length(joined) == 0L
    mu_new <- if (shifted) params$lambda_new0 else params$lambda_new
    n_new <- shifted + if (spec$new_perps == "poisson") {
      stats::rpois(1L, mu_new)
    } else {
      stats::rnbinom(1L, size = params$size_new, mu = mu_new)
    }
    new_ids <- as.character(st$n_nodes + seq_len(1L + n_new))
    st$n_nodes <- st$n_nodes + 1L + n_new
    event_id <- new_ids[1L]
    new_perps <- new_ids[-1L]

    perps <- c(joined, new_perps)
    st$ids <- c(st$ids, new_perps)
    st$deg <- c(st$deg, numeric(n_new))
    st$last <- c(st$last, numeric(n_new))
    st$u <- c(st$u, if (is.finite(frailty_shape)) {
      stats::rgamma(n_new, shape = frailty_shape, rate = frailty_shape)
    } else {
      rep(1, n_new)
    })
    names(st$deg) <- names(st$last) <- st$ids
    st$deg[perps] <- st$deg[perps] + 1
    st$last[perps] <- t_k
    for (v in perps) {
      others <- setdiff(perps, v)
      tv <- st$ties[[v]] %||% numeric(0)
      tv[others] <- t_k
      st$ties[[v]] <- tv
    }

    list(
      arrivals = data.frame(id = new_ids, role = c("event", rep("perp", n_new)),
                            stringsAsFactors = FALSE),
      edges = data.frame(i = rep(event_id, length(perps)), j = perps, stringsAsFactors = FALSE)
    )
  }
}
