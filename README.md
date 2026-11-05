# ppnet

Self-exciting point process models for growing networks, including bipartite networks.

## Usage

```r
library(ppnet)

# Simulate a bipartite event-perp history
truth <- list(mu = 2, K = 0.6, beta = 1.5, beta_edges = 0.3, lambda_new = 0.7)
ev <- sim_hawkesNet(truth, T_end = 100, mark_type = "ba_bip")$ev

# Fit, then use the usual model methods
fit <- fit_hawkesNet(ev, truth, mark_type = "ba_bip")
summary(fit)
confint(fit)
AIC(fit)
sims <- simulate(fit, nsim = 10, seed = 1)
```

See `vignette("ppnet")` for building up the time part and the mark part of
a model, and `?attachment_spec`, `?excitation_spec` and `?baseline_spec` for
the model terms.

## Acknowledgements

`ppnet` derives from `hawkesNet`. It continues
[Joseph Reps's work](https://github.com/JosephReps/hawkesNet)
([documentation](https://josephreps.github.io/hawkesNet/)), and all of the
development up to the move to this repository is his. `hawkesNet` is by
[Duncan Clark](https://github.com/duncan-clark/hawkesNet).
