# Diagnostic plots for drlate fits

Diagnostic plots for drlate fits

## Usage

``` r
# S3 method for class 'drlate'
plot(
  x,
  type = c("overlap", "balance", "balance_density", "weights"),
  bins = 30,
  geom = c("histogram", "density"),
  var = NULL,
  ...
)
```

## Arguments

- x:

  A fitted
  [`drlate()`](https://kvenkita.github.io/drlate/reference/drlate.md)
  object (with `keep_data = TRUE`).

- type:

  One of:

  - `"overlap"`: histograms (or kernel densities, see `geom`) of the
    estimated instrument propensity score by instrument arm, with the
    `pstolerance` bounds marked. Mass piling up near 0 or 1 signals
    overlap problems.

  - `"balance"`: a love plot of standardized mean differences from
    [`balance()`](https://kvenkita.github.io/drlate/reference/balance.md),
    unweighted vs IPW-weighted, with the conventional \|SMD\| = 0.1
    reference lines.

  - `"balance_density"`: kernel densities of the covariates by
    instrument arm, raw versus IPW-weighted (the Stata
    `latebalance density` display). Weighting that balances a covariate
    brings the two arm densities together in the weighted panel.

  - `"weights"`: distributions of the implied IPW weights by arm; a long
    right tail means a few observations dominate the estimate.

- bins:

  Number of histogram bins for `"overlap"` and `"weights"`.

- geom:

  For `type = "overlap"`, either `"histogram"` (default) or `"density"`
  (a kernel-density overlap matching Stata `lateoverlap`).

- var:

  For `type = "balance_density"`, an optional character vector selecting
  covariates to plot; defaults to all model covariates.

- ...:

  Currently unused.

## Value

A `ggplot` object.
