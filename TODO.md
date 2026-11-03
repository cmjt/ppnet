# TODO

## Model checks (step 6)

- [ ] **Simulation check of the co-offender effect `theta`.** Simulate with
  `attachment_spec(degree = TRUE, coffenders = TRUE)`, `theta = 0` and perp
  frailty (`sim_hawkesNet(frailty_shape = ...)`), fit the model, and record how
  often `theta` comes out positive and significant. Repeat over a few frailty
  variances, fitted to match the observed spread in attendance counts.
- [ ] **Perp-driven null for the event-level `gamma`.** In the current simulator
  event times are generated first and perps are assigned to events afterwards,
  so event timing cannot depend on who offends, and `gamma` is unbiased under it
  by construction. The real-data confound is the reverse: prolific perps both
  make events connected and cause further events. Write a simulator where each
  perp has its own offending process (rate with frailty, optional
  self-excitation) and events are formed when perps offend, alone or in groups,
  with no network effect. Fit the event model and check how often `gamma > 0`.
- [ ] **Move the attachment simulator to C++.** `.make_sim_mark_attach()` loops
  over every perp at every event in R (co-offender counts especially), which is
  too slow for simulation checks at the size of real data.
- [ ] **Window sensitivity.** Refit with `connectedness_window` and
  `coffender_window` at half and double the chosen value, and at `Inf`.

## Model extensions

- [ ] **Individual effects (M4).** Add a perp-level effect to the attachment
  weights, e.g. a gamma frailty integrated out per perp. Integrating it couples
  a perp's joins across events, so the running-sum evaluator needs rethinking:
  `W` would depend on each perp's posterior frailty given its history.
- [ ] **Comparison tools.** `compare_fits()` for nested likelihood-ratio tests and
  AIC/BIC tables, replacing the hand-written comparisons.
- [ ] **Temporal diagnostics.** `R/diagnostics_temporal.R` still assumes a
  constant baseline and a single unweighted kernel; update it for
  `baseline_spec()`, `excitation_spec()` and `kernel = "exp2"` (compensator
  residuals, time-rescaling QQ plot).

## Towards a change-statistic (CS) version of the attachment model

The current CS-bip kernel uses `ernm` change statistics on `network` objects, a
logistic link with decay since a perp's first event, and stores a matrix per
(event, old perp) pair. To bring it in line with the attachment model:

- [ ] **Change-statistic weights.** Generalise the attachment weight to
  `log w_v = theta' s_v(t) - beta_e (t - last_v)`, where `s_v(t)` are bipartite
  change statistics for adding the edge (event, `v`): attendance count,
  co-offender count (optionally windowed), perp covariates. The current model
  is the case `s_v = (log(d_v + delta), h(x_v))` with the first coefficient
  fixed at 1, so it stays nested.
- [ ] **Compute the statistics in the C++ walk**, not through `ernm`, keeping
  `ernm` only as a reference in tests (statistics match `computeChangeStats()`
  on small networks).
- [ ] **Decide how to handle event-dependent statistics.** Statistics that depend
  only on `v`'s own history keep joins independent given the past, so the
  `p_v = 1 - exp(-m w_v / W)` link and the running sum `W` still give linear
  cost. Statistics that depend on who else joins the same event (e.g. "`v` has
  co-offended with another perp at this event", shared partners) make joins
  within an event dependent: they need an ordering of joins within the event,
  or an ERGM-style within-event model with its own normalising constant (or a
  pseudo-likelihood), and the linear-cost trick no longer applies.
- [ ] **Link function.** Keep the normalised `m w_v / W` link (linear cost,
  `m` interpretable as returning perps per event) or the current unnormalised
  logistic link (absolute probabilities, but every evaluation touches every
  (event, old perp) pair). Compare both on simulated data.
- [ ] **Tests as for the attachment model:** brute-force likelihood, gradients
  against central differences, recovery from simulation.
