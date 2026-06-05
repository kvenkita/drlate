# Diagnostic plots for drlate fits

Diagnostic plots for drlate fits

## Usage

``` r
# S3 method for class 'drlate'
plot(x, type = c("overlap", "balance", "weights"), bins = 30, ...)
```

## Arguments

- x:

  A fitted
  [`drlate()`](https://kvenkita.github.io/drlate/reference/drlate.md)
  object (with `keep_data = TRUE`).

- type:

  One of:

  - `"overlap"`: histograms of the estimated instrument propensity score
    by instrument arm, with the `pstolerance` bounds marked. Mass piling
    up near 0 or 1 signals overlap problems.

  - `"balance"`: a love plot of standardized mean differences from
    [`balance()`](https://kvenkita.github.io/drlate/reference/balance.md),
    unweighted vs IPW-weighted, with the conventional \|SMD\| = 0.1
    reference lines.

  - `"weights"`: distributions of the implied IPW weights by arm; a long
    right tail means a few observations dominate the estimate.

- bins:

  Number of histogram bins for `"overlap"` and `"weights"`.

- ...:

  Currently unused.

## Value

A `ggplot` object.
