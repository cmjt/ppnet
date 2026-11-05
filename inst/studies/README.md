# Simulation studies

Scripts that check ppnet's estimators by simulation. They are installed with
the package; find them with

```r
system.file("studies", package = "ppnet")
```

They need ppnet installed (they call `library(ppnet)`).

## Simulation studies (one replicate per run, for a cluster)

Each simulates one data set from a fixed model, fits it, and saves the result
to `<out_dir>`:

```sh
Rscript ba_replicate.R <rep_id> <out_dir>       # BA, unipartite
Rscript ba_bip_replicate.R <rep_id> <out_dir>   # BA, bipartite
Rscript cs_replicate.R <rep_id> <out_dir>       # change statistics (ERGM terms)
```

`rep_id` sets the seed, so replicates are reproducible. Results are
`<out_dir>/<model>_rep_<rep_id>.rds`.

## Null checks

Each simulates data with no network effect and checks how often the fitted
model finds one anyway. Run them with `source()`; they save their results to
`out_file` (set it before sourcing to change the default file name in the
working directory).

* `theta_null_check.R`: does an unmodelled perp random effect make the
  co-offender effect `theta` look positive? Fits the attachment model with
  and without the random effect.
* `gamma_null_check.R`: does perps' own offending (`sim_perp_events()`) make
  the connectedness effect `gamma` look positive?
