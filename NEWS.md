# ppnet 0.1.0

## Models

* Attachment model for bipartite event-perp data (`attachment_spec()`):
  degree and recency, co-offenders, 3-paths, 4-cycles (pseudo-likelihood,
  with sandwich standard errors from `attachment_vcov()`), an episode gap
  with recent own events, negative-binomial new perps, and a perp random
  effect.
* Time part: seasonal and trend baselines (`baseline_spec()`), a two-scale
  kernel (`kernel = "exp2"`), and size, connectedness and mark-path dependent
  excitation (`excitation_spec()`).

## Fitting and inference

* Fits have class `"ppnet_fit"` with `print()`, `summary()`, `coef()`,
  `vcov()`, `confint()`, `logLik()`, `nobs()`, `anova()` and `simulate()`
  methods. Estimates and standard errors are on the natural scale; pseudo-
  likelihood fits use sandwich standard errors and refuse likelihood-based
  comparisons.
* `params_init` is optional (starting values default to values worked out
  from the data), and `part = "time"` or `"marks"` fits one part of the model
  on its own.
* The mark part of the likelihood is cached while its parameters are
  unchanged, which makes time-part fits much faster.
* `fit_hawkesNet()` and `sim_hawkesNet()` report timings only with
  `verbose = TRUE`.

## Simulation and goodness of fit

* `sim_hawkesNet()` can continue an observed history (`history =`), use a
  burn-in (`burn_in =`) and stop runaway simulations (`max_events`,
  `max_seconds`); `simulate(fit, from = "end")` continues the data.
* `gof()` compares observed and simulated summaries of event times, group
  sizes, returning and new perps, attendance and repeat co-offending.

## Plots and diagnostics

* `diag_network()` and its plots summarise bipartite data with events and
  perps kept apart (group size, returning perps, attendance).
* `plot()` on networks from bipartite data colours events and perps, and
  `component =` draws the largest connected components; `net_components()`
  lists them.

## Documentation and structure

* Vignettes: "Getting started with ppnet" and "The likelihoods in ppnet and
  hawkesNet" (each likelihood set against the HawkesNet paper and package).
* Simulation studies and null checks moved to `inst/studies/`.
* The undocumented `net` data set was removed.
