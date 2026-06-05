# Covariate balance across instrument arms

Computes standardized mean differences (SMDs) of the model covariates
between the two instrument arms, before and after weighting by the
inverse of the estimated instrument propensity score. Well-balanced
weighted covariates (conventionally, absolute SMD below 0.1) indicate
that the propensity score model is doing its job.

## Usage

``` r
balance(object, ...)

# S3 method for class 'drlate'
balance(object, ...)
```

## Arguments

- object:

  A fitted
  [`drlate()`](https://kvenkita.github.io/drlate/reference/drlate.md)
  object (with `keep_data = TRUE`).

- ...:

  Currently unused.

## Value

A data frame with one row per covariate and columns `variable`,
`smd_unweighted`, and `smd_weighted`.

## Details

The covariate set is the union of the columns of the instrument,
outcome, and treatment model matrices (the intercept is dropped). The
SMD denominator is the unweighted pooled standard deviation
\\\sqrt{(s_1^2 + s_0^2)/2}\\ in both columns, so the two columns are
directly comparable. Weighted arm means are Hájek means using the
inverse-propensity weights implied by the fit (for `estimand = "latt"`,
the Z=0 arm uses the ATT odds weights \\p/(1-p)\\, matching the
estimator).

## See also

[`plot.drlate()`](https://kvenkita.github.io/drlate/reference/plot.drlate.md)
with `type = "balance"` for the love plot.
