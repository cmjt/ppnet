# attachment.R

#' Attachment specification for returning perps
#'
#' Replaces the BA-bip rule for which old perps join an event. At event
#' \eqn{k} each earlier perp \eqn{v} has weight
#' \deqn{w_v = (d_v + \delta)\, e^{\theta x_v + \theta_p y_v}\,
#'   e^{-\beta_e (t_k - \ell_v)}}
#' where \eqn{d_v} is the number of events \eqn{v} has attended, \eqn{\ell_v}
#' the time of the latest, and \eqn{x_v} and \eqn{y_v} are change statistics
#' for adding the edge (event \eqn{k}, \eqn{v}), from \eqn{v}'s own history:
#' \itemize{
#'   \item \eqn{x_v = h(\cdot)} of the number of distinct perps \eqn{v} has
#'     shared an event with (`coffenders`);
#'   \item \eqn{y_v = h(\cdot)} of the number of bipartite 3-paths
#'     \eqn{k - v - j - u} the edge completes, one for each earlier event
#'     \eqn{j} of \eqn{v} and other perp \eqn{u} at it: co-offenders counted
#'     with multiplicity, \eqn{\sum_j (n_j - 1)} (`paths`).
#' }
#' Both count only events at most `coffender_window` before \eqn{t_k}.
#'
#' Records of one episode (a spree, or one incident logged as several
#' events) can sit minutes apart. An episode gap \eqn{h} (`episode_gap`)
#' separates continuing an episode from coming back on a separate occasion:
#' `recent = TRUE` adds \eqn{e^{\theta_r r_v}} with \eqn{r_v = h(\cdot)} of
#' the number of \eqn{v}'s own events at most \eqn{h} before \eqn{t_k}
#' (`theta_recent`), and with `cycles` the 4-cycles below are split into those
#' through a shared event at most \eqn{h} before (`theta_cycles_short`) and
#' those more than \eqn{h} before (`theta_cycles`).
#' These depend only on \eqn{v}'s own history, so joins stay independent
#' given the past. It joins
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
#' With `random_effect = TRUE` each perp also has a gamma random effect (a
#' frailty, in survival-analysis terms)
#' \eqn{u_v \sim \text{Gamma}(1/\sigma^2, 1/\sigma^2)} (mean 1, variance
#' \eqn{\sigma^2} = `random_effect_var`) that multiplies its rate,
#' \eqn{p_v = 1 - \exp(-u_v m w_v / W)}, with \eqn{W} still the sum of the
#' weights without the random effect. This separates perps who return often for their
#' own reasons from the effects of attendance and co-offending. \eqn{u_v} is
#' integrated out of the likelihood perp by perp: exactly for perps that
#' never return, and otherwise by quadrature in \eqn{\log u_v} over the range
#' that holds the integrand, located from its mode (Gauss-Legendre panels of
#' `random_effect_nodes` nodes, widening away from the mode). The cost stays
#' linear in events and joins. `random_effect_var = 0` gives the model
#' without the random effect exactly.
#'
#' Parameters: `m`, `beta_edges` (\eqn{\beta_e}), `delta` (if `degree`),
#' `theta` (if `coffenders`) and `theta_paths` (\eqn{\theta_p}, if `paths`), then `lambda_new`, `lambda_new0` (if
#' `min_perps = 1`), `size_new` (if `"negbin"`) and `random_effect_var` (if
#' `random_effect`); see [attachment_param_names()]. With `degree`,
#' `coffenders` and `paths` all `FALSE` the weights depend on recency only.
#'
#' With `cycles = TRUE` joins at the same event depend on each other through
#' the 4-cycle change statistic \eqn{z_{vk} = h(\cdot)} of the number of
#' 4-cycles \eqn{k - v - j - u - k} that \eqn{v} joining closes: one for each
#' other old perp \eqn{u} at event \eqn{k} and earlier event \eqn{j} (within
#' `coffender_window`) that \eqn{u} and \eqn{v} both attended, i.e. past
#' co-offenders turning up together. Each old perp is scored by its
#' conditional given the others at the event,
#' \deqn{\mathrm{logit}\, P(v \text{ joins} \mid \text{rest}) =
#'   \log(e^{q_v} - 1) + \theta_c z_{vk} + \theta_{cs} z^{s}_{vk},
#'   \qquad q_v = m w_v / W,}
#' where \eqn{z^s_{vk}} counts the 4-cycles through shared events within the
#' episode gap (only with `episode_gap > 0`; \eqn{z_{vk}} then counts the
#' rest),
#' the full conditionals of an autologistic model for each event, so the
#' fit maximises a pseudo-likelihood. With the 4-cycle coefficients at 0
#' this is the model above exactly. A pseudo-likelihood is not a likelihood:
#' standard errors from its Hessian are wrong (typically too small), so use
#' [attachment_vcov()], and likelihood-ratio tests and AIC do not apply to
#' models with `cycles`. It cannot be combined with `random_effect`. Given
#' each perp's propensity joins are independent, so a fixed propensity does
#' not by itself make past co-offenders turn up together.
#'
#' Since \eqn{\sum_v m w_v / W = m}, each evaluation needs only the running
#' sum \eqn{W}, so its cost grows with the number of events and perp
#' appearances, not with events times perps.
#'
#' @param degree Logical; include \eqn{(d_v + \delta)}. Old perps have
#'   \eqn{d_v \ge 1}, so `delta = 0` (joins proportional to past attendance)
#'   is allowed; fix it there with `fixed_params` if the fit runs to 0.
#' @param coffenders Logical; include \eqn{e^{\theta x_v}} (distinct
#'   co-offenders).
#' @param coffender_scale `"log1p"` (default) for \eqn{h(n) = \log(1 + n)} or
#'   `"linear"` for \eqn{h(n) = n}, for both `coffenders` and `paths`.
#' @param coffender_window Only events at most this long before count towards
#'   `coffenders` and `paths`. `Inf` (the default) counts all.
#' @param paths Logical; include \eqn{e^{\theta_p y_v}} (3-paths, i.e.
#'   co-offenders with multiplicity).
#' @param cycles Logical; let joins depend on past co-offenders at the same
#'   event through 4-cycles (`theta_cycles`, and `theta_cycles_short` with an
#'   `episode_gap`), fitted by pseudo-likelihood.
#' @param episode_gap The episode gap \eqn{h}: shared events at most this long
#'   before count as the same episode. `0` (the default) for none. Must be
#'   less than `coffender_window`.
#' @param recent Logical; include \eqn{e^{\theta_r r_v}}, the perp's own
#'   events within `episode_gap` (`theta_recent`). Needs `episode_gap > 0`.
#' @param new_perps Distribution of the new-perp count: `"poisson"` (default)
#'   or `"negbin"`.
#' @param min_perps `0` (default) allows events with no perps; `1` requires at
#'   least one, as in observed data, and adds `lambda_new0`.
#' @param random_effect Logical; give each perp a gamma random effect with
#'   variance `random_effect_var`.
#' @param random_effect_nodes Gauss-Legendre nodes per panel for the
#'   random-effect integrals.
#'
#' @return An object of class `attachment_spec`.
#' @export
attachment_spec <- function(degree = TRUE, coffenders = FALSE,
                            coffender_scale = c("log1p", "linear"),
                            coffender_window = Inf, paths = FALSE, cycles = FALSE,
                            episode_gap = 0, recent = FALSE,
                            new_perps = c("poisson", "negbin"), min_perps = 0,
                            random_effect = FALSE, random_effect_nodes = 8L) {
  coffender_scale <- match.arg(coffender_scale)
  new_perps <- match.arg(new_perps)
  if (!(length(min_perps) == 1L && min_perps %in% c(0, 1))) {
    stop("`min_perps` must be 0 or 1.", call. = FALSE)
  }
  stopifnot(is.numeric(coffender_window), length(coffender_window) == 1L, coffender_window >= 0)
  stopifnot(is.numeric(episode_gap), length(episode_gap) == 1L, episode_gap >= 0)
  if (episode_gap > 0 && episode_gap >= coffender_window) {
    stop("`episode_gap` must be less than `coffender_window`.", call. = FALSE)
  }
  if (isTRUE(recent) && !(episode_gap > 0)) {
    stop("`recent` needs `episode_gap > 0`.", call. = FALSE)
  }
  if (isTRUE(cycles) && isTRUE(random_effect)) {
    stop("`cycles` cannot be combined with `random_effect`.", call. = FALSE)
  }
  structure(
    list(degree = isTRUE(degree), coffenders = isTRUE(coffenders),
         coffender_scale = coffender_scale, coffender_window = coffender_window,
         paths = isTRUE(paths), cycles = isTRUE(cycles), episode_gap = episode_gap,
         recent = isTRUE(recent), new_perps = new_perps, min_perps = as.integer(min_perps),
         random_effect = isTRUE(random_effect), random_effect_nodes = as.integer(random_effect_nodes)),
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
  c("m", "beta_edges", if (spec$degree) "delta", .attach_theta_names(spec), "lambda_new",
    if (spec$min_perps == 1L) "lambda_new0", if (spec$new_perps == "negbin") "size_new",
    if (isTRUE(spec$random_effect)) "random_effect_var", .attach_cycle_names(spec))
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

# Internal: each event's contribution to the gradient of .new_perp_count_ll(),
# one row per event.
.new_perp_count_scores <- function(spec, n_new, n_old, params) {
  shifted <- spec$min_perps == 1L & n_old == 0L
  y <- n_new - shifted
  nb <- spec$new_perps == "negbin"
  nms <- c("lambda_new", if (spec$min_perps == 1L) "lambda_new0", if (nb) "size_new")
  S <- matrix(0, length(y), length(nms), dimnames = list(NULL, nms))
  mu <- ifelse(shifted, params[["lambda_new0"]] %||% NA_real_, params[["lambda_new"]])
  if (!nb) {
    d_mu <- y / mu - 1
  } else {
    size <- params[["size_new"]]
    d_mu <- y / mu - (y + size) / (size + mu)
    S[, "size_new"] <- digamma(y + size) - digamma(size) + log(size / (size + mu)) +
      1 - (y + size) / (size + mu)
  }
  S[!shifted, "lambda_new"] <- d_mu[!shifted]
  if (spec$min_perps == 1L) S[shifted, "lambda_new0"] <- d_mu[shifted]
  S
}

.coffender_feature <- function(spec, n) {
  if (spec$coffender_scale == "log1p") log1p(n) else as.numeric(n)
}

# Internal: names of the own-history change-statistic coefficients, in
# column order.
.attach_theta_names <- function(spec) {
  c(if (spec$coffenders) "theta", if (isTRUE(spec$paths)) "theta_paths",
    if (isTRUE(spec$recent)) "theta_recent")
}

# Internal: names of the 4-cycle coefficients, in column order.
.attach_cycle_names <- function(spec) {
  if (!isTRUE(spec$cycles)) return(character(0))
  c("theta_cycles", if ((spec$episode_gap %||% 0) > 0) "theta_cycles_short")
}

# Internal: the change statistics for each state snapshot, one column per
# coefficient in .attach_theta_names().
.attach_features <- function(spec, h) {
  cols <- list(theta = h$upd_x, theta_paths = h$upd_y,
               theta_recent = h$upd_r)[.attach_theta_names(spec)]
  X <- matrix(vapply(cols, function(n) .coffender_feature(spec, n), numeric(length(h$upd_x))),
              nrow = length(h$upd_x))
  colnames(X) <- names(cols)
  X
}

# Attachment mark evaluator for a bipartite history.
compile_logpmf_attach <- function(events, spec) {
  stopifnot(inherits(spec, "attachment_spec"))
  a <- .event_arrays(events)
  if (all(is.na(a$role))) stop("Attachment needs node roles 'event' and 'perp'.", call. = FALSE)
  f <- .bip_focal(a)
  invalid <- !all(f$valid)
  cycles <- isTRUE(spec$cycles)
  h <- if (invalid) NULL else attach_history_cpp(a$t, a$node_ptr, a$edge_ptr, a$edge_i, a$edge_j,
                                                 f$focal, length(a$ids), spec$coffender_window,
                                                 cycles, spec$episode_gap %||% 0)
  cyc_names <- .attach_cycle_names(spec)
  cyc_z <- if (cycles && !invalid) {
    cbind(theta_cycles = .coffender_feature(spec, h$cyc_z),
          theta_cycles_short = .coffender_feature(spec, h$cyc_zs))[, cyc_names, drop = FALSE]
  } else {
    matrix(0, 0, 0)
  }
  no_cyc <- integer(0)
  th_names <- .attach_theta_names(spec)
  feat <- if (invalid) NULL else .attach_features(spec, h)
  n_old <- if (invalid) integer(0) else diff(h$obs_ptr)
  gl <- .gauss_legendre(spec$random_effect_nodes %||% 8L)

  # With `gradient = TRUE` the value carries attribute "gradient", and with
  # `scores = TRUE` attribute "scores": each event's contribution to it, one
  # row per event (not available with the random effect).
  function(params, gradient = FALSE, scores = FALSE) {
    get <- function(nm, default) if (nm %in% names(params)) params[[nm]] else default
    m <- get("m", NA_real_)
    beta_edges <- get("beta_edges", NA_real_)
    delta <- if (spec$degree) get("delta", NA_real_) else 0
    theta <- vapply(th_names, get, numeric(1), default = NA_real_)
    lambda_new <- get("lambda_new", NA_real_)
    lambda_new0 <- if (spec$min_perps == 1L) get("lambda_new0", NA_real_) else 1
    size_new <- if (spec$new_perps == "negbin") get("size_new", NA_real_) else 1
    random_effect_var <- if (isTRUE(spec$random_effect)) get("random_effect_var", NA_real_) else 0
    theta_c <- vapply(cyc_names, get, numeric(1), default = NA_real_)

    ok <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
    if (!ok(m) || m <= 0 || !ok(beta_edges) || beta_edges < 0 || !ok(lambda_new) || lambda_new < 0 ||
        !ok(delta) || delta < 0 || !all(is.finite(theta)) ||
        !ok(lambda_new0) || lambda_new0 < 0 || !ok(size_new) || size_new <= 0 ||
        !ok(random_effect_var) || random_effect_var < 0 || !all(is.finite(theta_c))) return(-Inf)
    if (invalid) return(-Inf)

    out <- attach_ll_cpp(a$t, h$upd_ptr, h$upd_node, h$upd_deg, h$upd_last, feat,
                         h$obs_ptr, h$obs_node, length(a$ids), m, beta_edges, delta, theta,
                         spec$degree, random_effect_var, gl$x, gl$w,
                         if (cycles) h$cyc_ptr else no_cyc, if (cycles) h$cyc_node else no_cyc,
                         cyc_z, if (cycles) h$cyc_y else no_cyc, theta_c, scores)
    r <- out$value
    counts <- .new_perp_count_ll(spec, f$n_new, n_old, params, gradient = gradient)
    ll <- as.numeric(counts) + r[1L]
    keep <- attachment_param_names(spec)
    P <- length(th_names)
    walk_names <- c("m", "beta_edges", "delta", th_names, "random_effect_var", cyc_names)
    if (gradient && is.finite(ll)) {
      g <- c(stats::setNames(r[-1L], walk_names), attr(counts, "gradient"))
      attr(ll, "gradient") <- g[intersect(keep, names(g))]
    }
    if (scores && is.finite(ll)) {
      S <- out$scores
      colnames(S) <- walk_names
      S <- cbind(S, .new_perp_count_scores(spec, f$n_new, n_old, params))
      attr(ll, "scores") <- S[, keep, drop = FALSE]
    }
    ll
  }
}

# Internal: a stateful mark simulator for sim_hawkesNet(), backed by the C++
# simulator in src/attach_sim.cpp. It is only called for accepted events, so
# each call advances the perp state. Each perp gets a fixed propensity
# u ~ Gamma(shape, shape) (mean 1) multiplying its rate, with shape
# 1 / random_effect_var if the spec has a random effect, else
# `random_effect_shape` (a random effect the fitted model ignores, for checking what other effects pick up).
.make_sim_mark_attach <- function(spec, random_effect_shape = Inf) {
  stopifnot(is.numeric(random_effect_shape), length(random_effect_shape) == 1L, random_effect_shape > 0)
  if (isTRUE(spec$random_effect) && is.finite(random_effect_shape)) {
    stop("Give `random_effect_shape` only when the attachment spec has no random effect.", call. = FALSE)
  }
  sim <- NULL
  perp_ids <- character(0)   # perp id by C++ index
  n_ids <- 0L                # node ids handed out so far (events and perps)
  id_prefix <- ""            # "s" after a history, so new ids cannot clash

  start <- function(params, t0) {
    if (isTRUE(spec$random_effect)) {
      fv <- params[["random_effect_var"]]
      random_effect_shape <- if (fv > 0) 1 / fv else Inf
    }
    sim <<- attach_sim_new(params, spec$degree, spec$coffenders, isTRUE(spec$paths),
                           isTRUE(spec$cycles), isTRUE(spec$recent), spec$episode_gap %||% 0,
                           spec$coffender_scale == "log1p", spec$coffender_window,
                           spec$new_perps == "negbin", spec$min_perps == 1L,
                           random_effect_shape, t0)
  }

  # With `history` (an events object), load the observed events into the
  # state instead of simulating: their perps' attendance, recency,
  # co-offenders and shared events. Random effects of history perps are
  # drawn from their distribution, not conditioned on their history.
  function(net, t_k, params, ..., history = NULL) {
    if (!is.null(history)) {
      per_event <- .history_perps(history)
      start(params, if (length(per_event$t)) per_event$t[1] else t_k)
      index <- integer(0)
      for (k in seq_along(per_event$t)) {
        old <- per_event$old[[k]]
        new <- per_event$new[[k]]
        attach_sim_observe(sim, per_event$t[k], unname(index[old]) - 1L, length(new))
        index[new] <- length(perp_ids) + seq_along(new)
        perp_ids <<- c(perp_ids, new)
      }
      id_prefix <<- "s"
      return(invisible(NULL))
    }
    if (is.null(sim)) start(params, t_k)
    r <- attach_sim_step(sim, t_k)

    new_ids <- paste0(id_prefix, n_ids + seq_len(1L + r$n_new))
    n_ids <<- n_ids + 1L + r$n_new
    event_id <- new_ids[1L]
    new_perps <- new_ids[-1L]
    joined <- perp_ids[r$joined + 1L]
    perp_ids <<- c(perp_ids, new_perps)

    list(
      arrivals = data.frame(id = new_ids, role = c("event", rep("perp", r$n_new)),
                            stringsAsFactors = FALSE),
      edges = data.frame(i = rep(event_id, length(joined) + r$n_new), j = c(joined, new_perps),
                         stringsAsFactors = FALSE)
    )
  }
}

# Internal: for each event of a bipartite history, in time order, its time,
# its old perps (seen at an earlier event) and its new perps.
.history_perps <- function(ev) {
  roles <- stats::setNames(ev$nodes$role, as.character(ev$nodes$id))
  if (is.null(ev$nodes$role) || !all(c("event", "perp") %in% roles)) {
    stop("A history needs node roles \"event\" and \"perp\".", call. = FALSE)
  }
  born <- stats::setNames(ev$nodes$event_id, as.character(ev$nodes$id))
  e <- ev$edges
  ei <- as.character(e$i)
  ej <- as.character(e$j)
  perp <- ifelse(roles[ei] == "perp", ei, ej)
  by_event <- split(perp, factor(e$event_id, levels = ev$times$event_id))
  old <- new <- vector("list", nrow(ev$times))
  for (k in seq_len(nrow(ev$times))) {
    p <- unique(by_event[[k]])
    is_new <- born[p] == ev$times$event_id[k]
    old[[k]] <- p[!is_new]
    new[[k]] <- p[is_new]
  }
  list(t = ev$times$t, old = old, new = new)
}
