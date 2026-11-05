# TODO

## Model checks (step 6)

- [x] **Simulation check of the co-offender effect `theta`** (`inst/studies/theta_null_check.R`). Simulate with
  `attachment_spec(degree = TRUE, coffenders = TRUE)`, `theta = 0` and perp
  random effect (`sim_hawkesNet(random_effect_shape = ...)`), fit the model, and record how
  often `theta` comes out positive and significant. Repeat over a few random-effect
  variances, fitted to match the observed spread in attendance counts.
  Result: with any random effect, `theta` comes out clearly positive and
  significant in every replicate, so a positive `theta` is not evidence of a
  co-offender effect without individual effects (M4).
- [x] **Perp-driven null for the event-level `gamma`** (`sim_perp_events()`,
  `inst/studies/gamma_null_check.R`). Perps enter over time with exponential careers and
  offend at rate `rate * u_v` (gamma random effect), optionally self-exciting on
  their own events or on every event they attend. Each offence is an event
  whose initiator brings `Poisson(extra_mean * u^size_link)` co-offenders,
  partly past co-offenders (`p_past`), otherwise active perps weighted by
  `u^recruit_link`. No rate depends on the network; a burn-in builds history
  before the window, which starts with an empty network.
  Result (20 replicates per scenario, generic values, `exc` window 52, `exp`
  kernel): `gamma` is not pushed positive in any scenario (random-effect variance
  0, 1, 2; with self-excitation on own or all events; random effect linked to
  group size). Rejections for `gamma > 0` were 0-16% and for `gamma < 0`
  5-20%. Without self-excitation the fitted `K` is near 0 and `gamma` is
  barely identified (huge standard errors). With self-excitation it is well
  identified and close to 0, slightly negative with own-event excitation.
  Still to do: rerun at values matched to the real data (dispersion of
  attendance counts, events per week), and with the `exp2` kernel and the
  seasonal baseline used in the real fits (`exp2` collapsed to
  `beta_fast = beta` on these data, so it was left out).
- [x] **Move the attachment simulator to C++** (`src/attach_sim.cpp`, Fenwick-tree sampling). `.make_sim_mark_attach()` loops
  over every perp at every event in R (co-offender counts especially), which is
  too slow for simulation checks at the size of real data.
- [ ] **Window sensitivity.** Refit with `connectedness_window` and
  `coffender_window` at half and double the chosen value, and at `Inf`.

## Model extensions

- [x] **Individual effects (M4).** `attachment_spec(random_effect = TRUE)` adds
  `random_effect_var`: each perp has `u_v ~ Gamma(1/sigma2, 1/sigma2)` and joins
  with `1 - exp(-u_v q_vk)`, `q_vk = m w_vk / W_k`, with `W_k` from the
  weights without the random effect. `u_v` is integrated out per perp:
  - `Q_v = sum_k q_vk` over the events where `v` is a candidate is linear-time
    from suffix sums of `exp(-beta t_k) / W_k` within reference-time epochs
    (prefix sums lost late stretches to cancellation).
  - Perps that never return are closed form. For the others the integral in
    `log u` uses the exact join probabilities: Newton for the mode (the log
    integrand is concave), then Gauss-Legendre panels widening away from it.
  - The Poisson-hit closed form first proposed was dropped: it cannot give a
    join probability near 1, so for prolific perps it made `theta` worse
    (1.55 against 0.22 without the random effect, truth 0).
  - `random_effect_var = 0` (or below 1e-8) is the model without the random effect, so it can
    be tested by a likelihood ratio (on the boundary: halve the chi-squared(1)
    p-value) with `fixed_params = "random_effect_var"`.
  - The simulator uses the same model (`random_effect_var` in a spec with `random_effect = TRUE`, or
    `random_effect_shape` for a random effect the fitted model ignores).
- [x] **Comparison tools.** `compare_fits()` for nested likelihood-ratio tests and
  AIC/BIC tables, replacing the hand-written comparisons.
- [x] **Temporal diagnostics.** `R/diagnostics_temporal.R` still assumes a
  constant baseline and a single unweighted kernel; update it for
  `baseline_spec()`, `excitation_spec()` and `kernel = "exp2"` (compensator
  residuals, time-rescaling QQ plot).

## Towards a change-statistic (CS) version of the attachment model

Change statistics are the change in counts of bipartite configurations when
the edge (event k, perp v) is added. Decisions so far: own-history
configurations first, normalised link `p_v = 1 - exp(-m w_v / W)` (linear
cost, nests the current model, works with the random effect); within-event
configurations, when added, by pseudo-likelihood.

- [x] **Vector of change statistics.** `log w_v = log(d_v + delta) +
  sum_j theta_j s_vj - beta_e (t - last_v)`; the C++ walk, gradients and
  random-effect integrals take a matrix of statistics, so new own-history
  statistics only need a column.
- [x] **3-paths** (`attachment_spec(paths = TRUE)`, `theta_paths`): the
  3-paths k - v - j - u the edge completes, i.e. v's co-offenders counted with
  multiplicity, `sum_j (n_j - 1)` over v's events in `coffender_window`, on
  the `coffender_scale`. Computed in the C++ history pass and the simulator;
  brute-force, gradient, simulator-calibration and recovery tests.
  `theta_paths` is strongly correlated with `delta` (3-paths are roughly
  attendance times group size - 1), so it is only well determined when group
  sizes vary.
- [ ] **Fit to the real data**: `paths` alone, with `coffenders` (distinct
  ties vs repeat co-offending with the same partners), and with the random
  effect.
- [x] **4-cycles by pseudo-likelihood** (`attachment_spec(cycles = TRUE)`,
  `theta_cycles`). `z_vk` = 4-cycles k - v - j - u - k closed by v joining k
  (shared past events, within `coffender_window`, with the other old perps at
  k). Each old perp is scored by its conditional given the rest of the event,
  `logit P = log(expm1(q_v)) + theta_cycles h(z_vk)`: the full conditionals
  of a per-event autologistic model, nesting the current model at 0. Only
  perps with z > 0 change from the independent terms, so the cost stays
  linear. Sandwich SEs from per-event scores (`attachment_vcov()`); Hessian
  SEs, LR tests and AIC are not valid for these fits. The simulator draws
  joins by Gibbs sampling (10 sweeps from the independent draw); tests check
  its joins against the conditionals.
  Like ERGMs with closure terms the model degenerates for large positive
  `theta_cycles` (past co-offenders pull each other in, groups grow, more
  shared events): with m = 1.5 it runs away by `theta_cycles` = 0.9, so
  simulations from a fitted value need checking.
- [x] **Episode gap** for records minutes apart (sprees, or incidents
  logged as several events), keeping continuous times. `episode_gap = h` in
  `attachment_spec()`: `recent = TRUE` adds `theta_recent` (own events within
  h), and with `cycles` the 4-cycles split into `theta_cycles_short` (shared
  event within h) and `theta_cycles` (beyond h, up to the window).
  `excitation_spec(connectedness_min_gap = h)` leaves same-episode events
  out of connectedness, so the fast kernel carries within-episode excitation.
- [x] **Mark-path dependent excitation** (`excitation_spec(path_dependent =
  TRUE)`, the hawkesNet paper's Definition 4.1). `c_i(t)` also counts later
  events sharing a perp with event i (same window and minimum gap, after
  `t_i`), so an event's excitation changes after it happens. Each change acts
  as an extra jump at the later event, discounted by `exp(-beta * lag)`, so
  the likelihood stays an exact O(events + links) recursion; the jumps now
  depend on `beta`. Brute-force, gradient, simulator-calibration and
  recovery tests. Much less stable than the fixed-at-birth form: with
  `gamma = 0.8`, window 5, `lambda_new = 0.5` the simulator ran away, and
  the paper's contraction condition fails for it. `gamma` here is partly
  defined by later offending, so it is not the connectedness effect.
- [ ] **Refit the real data with the episode gap** (h = 1 day, and 1 hour /
  1 week as checks): does `theta_cycles` beyond h stay out of the
  degenerate range, and do `gamma`, `theta` and the fast kernel change?
- [ ] **4-cycles with the random effect.** Disallowed for now: the per-perp
  integral would be over a product of conditionals (a composite marginal
  likelihood). A fixed propensity does not by itself make co-offenders turn
  up together (joins are independent given propensities), so this matters
  less than for `theta`.
- [ ] **Event stars** (perps already at event k) by the same route, if wanted.
- [ ] **Perp covariates** as nodal change statistics, if wanted.
- [ ] **Check the old `cs_bip` kernel**: does it score each candidate edge on
  the pre-event network (so within-event configurations never register)?
